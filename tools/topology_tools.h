// created by Sebastian Reiter
// s.b.reiter@gmail.com
// Jun 21, 2013 (d,m,y)

#ifndef __H__UG__topology_tools__
#define __H__UG__topology_tools__

#include <vector>
#include "../mesh.h"
#include "lib_grid/algorithms/remeshing/resolve_intersections.h"
#include "lib_grid/algorithms/geom_obj_util/vertex_util.h"
#include "lib_grid/algorithms/geom_obj_util/edge_util.h"
#include "lib_grid/algorithms/remove_duplicates_util.h"
#include "lib_grid/file_io/file_io.h"

namespace ug{
namespace promesh{

inline void EraseSelectedElements(Mesh* obj, bool eraseUnusedVrts,
						   bool eraseUnusedEdges, bool eraseUnusedFaces)
{
//	adjust selection
	Selector& sel = obj->selector();
	SelectAssociatedEdges(sel, sel.begin<Vertex>(), sel.end<Vertex>());
	SelectAssociatedFaces(sel, sel.begin<Edge>(), sel.end<Edge>());
	SelectAssociatedVolumes(sel, sel.begin<Face>(), sel.end<Face>());

	if(eraseUnusedFaces)
		SelectInnerSelectionFaces(sel);

	if(eraseUnusedEdges)
		SelectInnerSelectionEdges(sel);

	if(eraseUnusedVrts)
		SelectInnerSelectionVertices(sel);

//	erase selected elements
	Grid& grid = obj->grid();
	grid.erase(sel.begin<Volume>(), sel.end<Volume>());
	grid.erase(sel.begin<Face>(), sel.end<Face>());
	grid.erase(sel.begin<Edge>(), sel.end<Edge>());
	grid.erase(sel.begin<Vertex>(), sel.end<Vertex>());
}

///	returns the number of removed vertices
inline size_t RemoveDoubles(Mesh* obj, number threshold)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();

	size_t numVrts = grid.num<Vertex>();
	ug::RemoveDoubles<3>(grid, sel.begin<Vertex>(), sel.end<Vertex>(),
					 	 obj->position_attachment(), threshold);
	return numVrts - grid.num<Vertex>();
}

inline size_t RemoveDoubleEdges(Mesh* obj)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();

	size_t numEdges= grid.num<Edge>();
	RemoveDuplicates(grid, sel.begin<Edge>(), sel.end<Edge>());
	return numEdges - grid.num<Edge>();
}

inline size_t RemoveDoubleFaces(Mesh* obj)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();

	size_t numFaces= grid.num<Face>();
	RemoveDuplicates(grid, sel.begin<Face>(), sel.end<Face>());
	return numFaces - grid.num<Face>();
}

inline void MergeAtFirst(Mesh* obj)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	Mesh::position_accessor_t& aaPos = obj->position_accessor();

	SelectAssociatedGridObjects(sel);

	vector3 first = aaPos[*sel.vertices_begin()];
	Vertex* vrt = MergeMultipleVertices(grid, sel.vertices_begin(), sel.vertices_end());
	if(vrt)
		aaPos[vrt] = first;
}

inline void MergeAtCenter(Mesh* obj)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	Mesh::position_accessor_t& aaPos = obj->position_accessor();

	SelectAssociatedGridObjects(sel);

	vector3 center;
	CalculateCenter(center, sel, aaPos);
	Vertex* vrt = MergeMultipleVertices(grid, sel.vertices_begin(), sel.vertices_end());
	if(vrt)
		aaPos[vrt] = center;
}

inline void MergeAtLast(Mesh* obj)
{
	Mesh::position_accessor_t& aaPos = obj->position_accessor();
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	SelectAssociatedGridObjects(sel);

//	todo: This iteration shouldn't be necessary!
	VertexIterator vrtIter = sel.begin<Vertex>();
	Vertex* lastVrt = *vrtIter;
	for(; vrtIter != sel.end<Vertex>(); ++vrtIter)
		lastVrt = *vrtIter;

	vector3 last = aaPos[lastVrt];
	Vertex* vrt = MergeMultipleVertices(grid, sel.vertices_begin(), sel.vertices_end());
	if(vrt)
		aaPos[vrt] = last;
}

inline void CollapseEdge(Mesh* obj)
{
	using namespace std;
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	Mesh::position_accessor_t& aaPos = obj->position_accessor();

	vector<Face*> vFaces;
	vector<Edge*> vEdges;
	while(sel.num<Edge>() > 0){
		Edge* e = *sel.begin<Edge>();
	//	to make sure that all selected edges are collapsed,
	//	we have to check the adjacent triangles
		CollectFaces(vFaces, grid, e);
		for(size_t i = 0; i < vFaces.size(); ++i){
			Face* f = vFaces[i];
			if(f->num_edges() == 3){
				CollectEdges(vEdges, grid, f);
				int counter = 0;
				for(size_t j = 0; j < 3; ++j){
					if(sel.is_selected(vEdges[j]))
						++counter;
				}
			//	if two edges are selected, we have
			//	to mark the unselected edge, too (since we
			//	don't know which will be removed).
				if(counter == 2){
					for(size_t j = 0; j < 3; ++j)
						sel.select(vEdges[j]);
				}
			}
		}

	//	calculate the center
		VecAdd(aaPos[e->vertex(0)], aaPos[e->vertex(0)], aaPos[e->vertex(1)]);
		VecScale(aaPos[e->vertex(0)], aaPos[e->vertex(0)], 0.5);

	//	perform collapse
		CollapseEdge(grid, e, e->vertex(0));
	}
}

inline void SplitEdge(Mesh* obj)
{
	using namespace std;
//	collect all edges that shall be splitted in a vector
//	since new edges will be automatically selected again.
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	Mesh::position_accessor_t& aaPos = obj->position_accessor();

	vector<Edge*> vEdges;
	for(EdgeIterator iter = sel.begin<Edge>();
		iter != sel.end<Edge>(); ++iter)
	{
		vEdges.push_back(*iter);
	}

//	iterate through all edges in the vector and split them
	for(size_t i = 0; i < vEdges.size(); ++i){
		vector3 center = CalculateCenter(vEdges[i], aaPos);
		RegularVertex* vrt = ug::SplitEdge<RegularVertex>(grid, vEdges[i]);
		aaPos[vrt] = center;
	}
}

inline void SwapEdge(Mesh* obj)
{
	using namespace std;
//	collect all edges that shall be swapped in a vector
//	since new edges will be automatically selected again.
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	vector<Edge*> vEdges;
	for(EdgeIterator iter = sel.begin<Edge>();
		iter != sel.end<Edge>(); ++iter)
	{
		vEdges.push_back(*iter);
	}

//	iterate through all edges in the vector and swap them
//	if they are adjacent to two triangles
	Face* faces[2];
	for(size_t i = 0; i < vEdges.size(); ++i){
		int numFaces = GetAssociatedFaces(faces, grid, vEdges[i], 2);
		if(numFaces == 2){
			if(faces[0]->num_vertices() == 3 && faces[1]->num_vertices() == 3){
				SwapEdge(grid, vEdges[i]);
			}
		}
	}
}

inline void PlaneCut(Mesh* obj, const vector3& p, const vector3& n)
{
	Selector& sel = obj->selector();
	CutEdgesWithPlane(sel, p, n);
}

inline void AdjustEdgeOrientation(Mesh* obj)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();

	AdjustEdgeOrientationToFaceOrientation(grid, sel.begin<Edge>(),
												sel.end<Edge>());
}

inline void FixFaceOrientation(Mesh* obj)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();

	FixFaceOrientation(grid, sel.begin<Face>(),
							sel.end<Face>());
}

inline void FixFaceSubsetOrientations(Mesh* obj)
{
	Grid& grid = obj->grid();
	SubsetHandler& sh = obj->subset_handler();

	for(int i = 0; i < sh.num_subsets(); ++i){
		FixFaceOrientation(grid, sh.begin<Face>(i), sh.end<Face>(i));
	}
}

inline int FixVolumeOrientation(Mesh* obj)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	Mesh::position_accessor_t& aaPos = obj->position_accessor();

	return FixOrientation(grid, sel.begin<Volume>(), sel.end<Volume>(), aaPos);
}


inline void InvertFaceOrientation(Mesh* obj)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();

	InvertOrientation(grid, sel.begin<Face>(), sel.end<Face>());
}

inline void ResolveEdgeIntersection(Mesh* obj, number snapThreshold)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	Mesh::position_accessor_t& aaPos = obj->position_accessor();

//	the grid may now contain some degenerated triangles. We'll try to
//	remove most of them by projecting vertices onto close edges
	SelectAssociatedVertices(sel, sel.begin<Edge>(), sel.end<Edge>());


	ProjectVerticesToCloseEdges(grid, sel.get_grid_objects(),
								obj->position_attachment(), snapThreshold);

	IntersectCloseEdges(grid, sel, aaPos, snapThreshold);

//	remove doubles now
	ug::RemoveDoubles<3>(grid, sel.begin<Vertex>(), sel.end<Vertex>(),
					 obj->position_attachment(), snapThreshold);
}

inline void ResolveTriangleIntersections(Mesh* obj, number snapThreshold)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	// Mesh::position_accessor_t& aaPos = obj->position_accessor();

//	first we'll resolve triangle-triangle intersections
	ug::ResolveTriangleIntersections(grid, sel.begin<Triangle>(),
							 	 sel.end<Triangle>(), snapThreshold,
							 	 obj->position_attachment());

//	the grid may now contain some degenerated triangles. We'll try to
//	remove most of them by projecting vertices onto close edges
	SelectAssociatedVertices(sel, sel.begin<Triangle>(), sel.end<Triangle>());
	SelectAssociatedEdges(sel, sel.begin<Triangle>(), sel.end<Triangle>());
	ProjectVerticesToCloseEdges(grid, sel.get_grid_objects(),
								obj->position_attachment(), snapThreshold);

//	remove doubles now
	ug::RemoveDoubles<3>(grid, sel.begin<Vertex>(), sel.end<Vertex>(),
					 obj->position_attachment(), snapThreshold);
}

inline void ProjectVerticesToCloseEdges(Mesh* obj, number snapThreshold)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	// Mesh::position_accessor_t& aaPos = obj->position_accessor();

	ProjectVerticesToCloseEdges(grid, sel.get_grid_objects(),
								obj->position_attachment(), snapThreshold);

//	remove doubles now
	ug::RemoveDoubles<3>(grid, sel.begin<Vertex>(), sel.end<Vertex>(),
					 obj->position_attachment(), snapThreshold);
}

inline void ProjectVerticesToCloseFaces(Mesh* obj, number snapThreshold)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	//Mesh::position_accessor_t& aaPos = obj->position_accessor();

	ProjectVerticesToCloseFaces(grid, sel,
								obj->position_attachment(), snapThreshold);

//	remove doubles now
	ug::RemoveDoubles<3>(grid, sel.begin<Vertex>(), sel.end<Vertex>(),
					 obj->position_attachment(), snapThreshold);
}

inline void IntersectCloseEdges(Mesh* obj, number snapThreshold)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	Mesh::position_accessor_t& aaPos = obj->position_accessor();

	IntersectCloseEdges(grid, sel, aaPos, snapThreshold);
}


inline void ResolveSelfIntersections(Mesh* obj, number snapThreshold)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	Mesh::position_accessor_t& aaPos = obj->position_accessor();

	Triangulate(grid, sel.begin<Quadrilateral>(), sel.end<Quadrilateral>(), &aaPos);
	
	bool intersectFaces = sel.num<Face>() > 0;

	SelectAssociatedEdges(sel, sel.begin<Face>(), sel.end<Face>());
	SelectAssociatedVertices(sel, sel.begin<Face>(), sel.end<Face>());

	//	remove doubles again
	ug::RemoveDoubles<3>(grid, sel.begin<Vertex>(), sel.end<Vertex>(),
					 obj->position_attachment(), snapThreshold);

	if(intersectFaces)
		ProjectVerticesToCloseFaces(grid, sel,
									obj->position_attachment(), snapThreshold);

	ProjectVerticesToCloseEdges(grid, sel.get_grid_objects(),
								obj->position_attachment(), snapThreshold);

//	resolve face intersections
	if(intersectFaces){
		ug::ResolveTriangleIntersections(grid, sel.begin<Triangle>(),
								 	 sel.end<Triangle>(), snapThreshold,
								 	 obj->position_attachment());

	//	the grid may now contain some degenerated Faces. We'll try to
	//	remove most of them by projecting vertices onto close edges
		SelectAssociatedEdges(sel, sel.begin<Face>(), sel.end<Face>());
		SelectAssociatedVertices(sel, sel.begin<Face>(), sel.end<Face>());
		ProjectVerticesToCloseEdges(grid, sel.get_grid_objects(),
									obj->position_attachment(), snapThreshold);
	}
	else{
		Mesh::position_accessor_t& aaPos = obj->position_accessor();
		IntersectCloseEdges(grid, sel, aaPos, snapThreshold);
	}


//	remove doubles again
	ug::RemoveDoubles<3>(grid, sel.begin<Vertex>(), sel.end<Vertex>(),
					 obj->position_attachment(), snapThreshold);
}

}}// end of namespace

#endif
