// created by Sebastian Reiter
// s.b.reiter@gmail.com
// Jun 21, 2013 (d,m,y)

#ifndef __H__UG__topology_tools__
#define __H__UG__topology_tools__

#include <vector>
#include "../mesh_object.h"
#include "lib_grid/algorithms/remeshing/resolve_intersections.h"
#include "lib_grid/algorithms/geom_obj_util/vertex_util.h"
#include "lib_grid/algorithms/geom_obj_util/edge_util.h"

namespace ug{
namespace promesh{

void EraseSelectedElements(MeshObject* obj, bool eraseUnusedVrts,
						   bool eraseUnusedEdges, bool eraseUnusedFaces)
{
//	adjust selection
	Selector& sel = obj->get_selector();
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
	Grid& grid = obj->get_grid();
	grid.erase(sel.begin<Volume>(), sel.end<Volume>());
	grid.erase(sel.begin<Face>(), sel.end<Face>());
	grid.erase(sel.begin<Edge>(), sel.end<Edge>());
	grid.erase(sel.begin<Vertex>(), sel.end<Vertex>());
}

///	returns the number of removed vertices
size_t RemoveDoubles(MeshObject* obj, number threshold)
{
	Grid& grid = obj->get_grid();
	Selector& sel = obj->get_selector();

	size_t numVrts = grid.num<Vertex>();
	ug::RemoveDoubles<3>(grid, sel.begin<Vertex>(), sel.end<Vertex>(),
					 	 obj->position_attachment(), threshold);
	return numVrts - grid.num<Vertex>();
}

size_t RemoveDoubleEdges(MeshObject* obj)
{
	Grid& grid = obj->get_grid();
	Selector& sel = obj->get_selector();

	size_t numEdges= grid.num<Edge>();
	RemoveDoubleEdges(grid, sel.begin<Edge>(), sel.end<Edge>());
	return numEdges - grid.num<Edge>();
}

void MergeAtFirst(MeshObject* obj)
{
	Grid& grid = obj->get_grid();
	Selector& sel = obj->get_selector();
	MeshObject::position_accessor_t& aaPos = obj->position_accessor();

	SelectAssociatedGridObjects(sel);

	vector3 first = aaPos[*sel.vertices_begin()];
	Vertex* vrt = MergeMultipleVertices(grid, sel.vertices_begin(), sel.vertices_end());
	if(vrt)
		aaPos[vrt] = first;
}

void MergeAtCenter(MeshObject* obj)
{
	Grid& grid = obj->get_grid();
	Selector& sel = obj->get_selector();
	MeshObject::position_accessor_t& aaPos = obj->position_accessor();

	SelectAssociatedGridObjects(sel);

	vector3 center;
	CalculateCenter(center, sel, aaPos);
	Vertex* vrt = MergeMultipleVertices(grid, sel.vertices_begin(), sel.vertices_end());
	if(vrt)
		aaPos[vrt] = center;
}

void MergeAtLast(MeshObject* obj)
{
	MeshObject::position_accessor_t& aaPos = obj->position_accessor();
	Grid& grid = obj->get_grid();
	Selector& sel = obj->get_selector();
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

void CollapseEdge(MeshObject* obj)
{
	using namespace std;
	Grid& grid = obj->get_grid();
	Selector& sel = obj->get_selector();
	MeshObject::position_accessor_t& aaPos = obj->position_accessor();

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

void SplitEdge(MeshObject* obj)
{
	using namespace std;
//	collect all edges that shall be splitted in a vector
//	since new edges will be automatically selected again.
	Grid& grid = obj->get_grid();
	Selector& sel = obj->get_selector();
	MeshObject::position_accessor_t& aaPos = obj->position_accessor();

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

void SwapEdge(MeshObject* obj)
{
	using namespace std;
//	collect all edges that shall be swapped in a vector
//	since new edges will be automatically selected again.
	Grid& grid = obj->get_grid();
	Selector& sel = obj->get_selector();
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

void PlaneCut(MeshObject* obj, const vector3& p, const vector3& n)
{
	Selector& sel = obj->get_selector();
	CutEdgesWithPlane(sel, p, n);
}

void AdjustEdgeOrientation(MeshObject* obj)
{
	Grid& grid = obj->get_grid();
	Selector& sel = obj->get_selector();

	AdjustEdgeOrientationToFaceOrientation(grid, sel.begin<Edge>(),
												sel.end<Edge>());
}

void FixFaceOrientation(MeshObject* obj)
{
	Grid& grid = obj->get_grid();
	Selector& sel = obj->get_selector();

	FixFaceOrientation(grid, sel.begin<Face>(),
							sel.end<Face>());
}

void FixFaceSubsetOrientations(MeshObject* obj)
{
	Grid& grid = obj->get_grid();
	SubsetHandler& sh = obj->get_subset_handler();

	for(int i = 0; i < sh.num_subsets(); ++i){
		FixFaceOrientation(grid, sh.begin<Face>(i), sh.end<Face>(i));
	}
}

int FixVolumeOrientation(MeshObject* obj)
{
	Grid& grid = obj->get_grid();
	Selector& sel = obj->get_selector();
	MeshObject::position_accessor_t& aaPos = obj->position_accessor();

	return FixOrientation(grid, sel.begin<Volume>(), sel.end<Volume>(), aaPos);
}


void InvertFaceOrientation(MeshObject* obj)
{
	Grid& grid = obj->get_grid();
	Selector& sel = obj->get_selector();

	InvertOrientation(grid, sel.begin<Face>(), sel.end<Face>());
}

void ResolveEdgeIntersection(MeshObject* obj, number snapThreshold)
{
	Grid& grid = obj->get_grid();
	Selector& sel = obj->get_selector();
	MeshObject::position_accessor_t& aaPos = obj->position_accessor();

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

void ResolveTriangleIntersections(MeshObject* obj, number snapThreshold)
{
	Grid& grid = obj->get_grid();
	Selector& sel = obj->get_selector();
	// MeshObject::position_accessor_t& aaPos = obj->position_accessor();

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

void ProjectVerticesToCloseEdges(MeshObject* obj, number snapThreshold)
{
	Grid& grid = obj->get_grid();
	Selector& sel = obj->get_selector();
	// MeshObject::position_accessor_t& aaPos = obj->position_accessor();

	ProjectVerticesToCloseEdges(grid, sel.get_grid_objects(),
								obj->position_attachment(), snapThreshold);

//	remove doubles now
	ug::RemoveDoubles<3>(grid, sel.begin<Vertex>(), sel.end<Vertex>(),
					 obj->position_attachment(), snapThreshold);
}

void ProjectVerticesToCloseFaces(MeshObject* obj, number snapThreshold)
{
	Grid& grid = obj->get_grid();
	Selector& sel = obj->get_selector();
	//MeshObject::position_accessor_t& aaPos = obj->position_accessor();

	ProjectVerticesToCloseFaces(grid, sel,
								obj->position_attachment(), snapThreshold);

//	remove doubles now
	ug::RemoveDoubles<3>(grid, sel.begin<Vertex>(), sel.end<Vertex>(),
					 obj->position_attachment(), snapThreshold);
}

void IntersectCloseEdges(MeshObject* obj, number snapThreshold)
{
	Grid& grid = obj->get_grid();
	Selector& sel = obj->get_selector();
	MeshObject::position_accessor_t& aaPos = obj->position_accessor();

	IntersectCloseEdges(grid, sel, aaPos, snapThreshold);
}


void ResolveSelfIntersections(MeshObject* obj, number snapThreshold)
{
	Grid& grid = obj->get_grid();
	Selector& sel = obj->get_selector();
	MeshObject::position_accessor_t& aaPos = obj->position_accessor();

	Triangulate(grid, sel.begin<Quadrilateral>(), sel.end<Quadrilateral>(), &aaPos);
	
	bool intersectFaces = sel.num<Face>() > 0;

	SelectAssociatedEdges(sel, sel.begin<Face>(), sel.end<Face>());
	SelectAssociatedVertices(sel, sel.begin<Face>(), sel.end<Face>());

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
		MeshObject::position_accessor_t& aaPos = obj->position_accessor();
		IntersectCloseEdges(grid, sel, aaPos, snapThreshold);
	}


//	remove doubles again
	ug::RemoveDoubles<3>(grid, sel.begin<Vertex>(), sel.end<Vertex>(),
					 obj->position_attachment(), snapThreshold);
}

}}// end of namespace

#endif
