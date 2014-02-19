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
	SelectAssociatedEdges(sel, sel.begin<VertexBase>(), sel.end<VertexBase>());
	SelectAssociatedFaces(sel, sel.begin<EdgeBase>(), sel.end<EdgeBase>());
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
	grid.erase(sel.begin<EdgeBase>(), sel.end<EdgeBase>());
	grid.erase(sel.begin<VertexBase>(), sel.end<VertexBase>());
}

///	returns the number of removed vertices
size_t RemoveDoubles(MeshObject* obj, number threshold)
{
	Grid& grid = obj->get_grid();
	Selector& sel = obj->get_selector();

	size_t numVrts = grid.num<VertexBase>();
	ug::RemoveDoubles<3>(grid, sel.begin<VertexBase>(), sel.end<VertexBase>(),
					 	 obj->position_attachment(), threshold);
	return numVrts - grid.num<VertexBase>();
}

size_t RemoveDoubleEdges(MeshObject* obj)
{
	Grid& grid = obj->get_grid();
	Selector& sel = obj->get_selector();

	size_t numEdges= grid.num<EdgeBase>();
	RemoveDoubleEdges(grid, sel.begin<EdgeBase>(), sel.end<EdgeBase>());
	return numEdges - grid.num<EdgeBase>();
}

void MergeAtFirst(MeshObject* obj)
{
	Grid& grid = obj->get_grid();
	Selector& sel = obj->get_selector();
	MeshObject::position_accessor_t& aaPos = obj->position_accessor();

	SelectAssociatedGridObjects(sel);

	vector3 first = aaPos[*sel.vertices_begin()];
	VertexBase* vrt = MergeMultipleVertices(grid, sel.vertices_begin(), sel.vertices_end());
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
	VertexBase* vrt = MergeMultipleVertices(grid, sel.vertices_begin(), sel.vertices_end());
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
	VertexBaseIterator vrtIter = sel.begin<VertexBase>();
	VertexBase* lastVrt = *vrtIter;
	for(; vrtIter != sel.end<VertexBase>(); ++vrtIter)
		lastVrt = *vrtIter;

	vector3 last = aaPos[lastVrt];
	VertexBase* vrt = MergeMultipleVertices(grid, sel.vertices_begin(), sel.vertices_end());
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
	vector<EdgeBase*> vEdges;
	while(sel.num<EdgeBase>() > 0){
		EdgeBase* e = *sel.begin<EdgeBase>();
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

	vector<EdgeBase*> vEdges;
	for(EdgeBaseIterator iter = sel.begin<EdgeBase>();
		iter != sel.end<EdgeBase>(); ++iter)
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
	vector<EdgeBase*> vEdges;
	for(EdgeBaseIterator iter = sel.begin<EdgeBase>();
		iter != sel.end<EdgeBase>(); ++iter)
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

	AdjustEdgeOrientationToFaceOrientation(grid, sel.begin<EdgeBase>(),
												sel.end<EdgeBase>());
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
	SelectAssociatedVertices(sel, sel.begin<EdgeBase>(), sel.end<EdgeBase>());


	ProjectVerticesToCloseEdges(grid, sel.get_grid_objects(),
								aaPos, snapThreshold);

	IntersectCloseEdges(grid, sel.get_grid_objects(),
						aaPos, snapThreshold);
//	remove doubles now
	ug::RemoveDoubles<3>(grid, sel.begin<VertexBase>(), sel.end<VertexBase>(),
					 obj->position_attachment(), snapThreshold);
}

void ResolveTriangleIntersections(MeshObject* obj, number snapThreshold)
{
	Grid& grid = obj->get_grid();
	Selector& sel = obj->get_selector();
	MeshObject::position_accessor_t& aaPos = obj->position_accessor();

//	first we'll resolve triangle-triangle intersections
	ug::ResolveGridIntersections(grid, sel.begin<Triangle>(),
							 	 sel.end<Triangle>(), snapThreshold, aaPos);

//	the grid may now contain some degenerated triangles. We'll try to
//	remove most of them by projecting vertices onto close edges
	SelectAssociatedVertices(sel, sel.begin<Triangle>(), sel.end<Triangle>());
	SelectAssociatedEdges(sel, sel.begin<Triangle>(), sel.end<Triangle>());
	ProjectVerticesToCloseEdges(grid, sel.get_grid_objects(),
								aaPos, snapThreshold);

//	remove doubles now
	ug::RemoveDoubles<3>(grid, sel.begin<VertexBase>(), sel.end<VertexBase>(),
					 obj->position_attachment(), snapThreshold);
}

void ProjectVerticesToCloseEdges(MeshObject* obj, number snapThreshold)
{
	Grid& grid = obj->get_grid();
	Selector& sel = obj->get_selector();
	MeshObject::position_accessor_t& aaPos = obj->position_accessor();

	ProjectVerticesToCloseEdges(grid, sel.get_grid_objects(),
								aaPos, snapThreshold);

//	remove doubles now
	ug::RemoveDoubles<3>(grid, sel.begin<VertexBase>(), sel.end<VertexBase>(),
					 obj->position_attachment(), snapThreshold);
}

void ProjectVerticesToCloseFaces(MeshObject* obj, number snapThreshold)
{
	Grid& grid = obj->get_grid();
	Selector& sel = obj->get_selector();
	MeshObject::position_accessor_t& aaPos = obj->position_accessor();

	ProjectVerticesToCloseFaces(grid, sel.get_grid_objects(),
								aaPos, snapThreshold);

//	remove doubles now
	ug::RemoveDoubles<3>(grid, sel.begin<VertexBase>(), sel.end<VertexBase>(),
					 obj->position_attachment(), snapThreshold);
}

void IntersectCloseEdges(MeshObject* obj, number snapThreshold)
{
	Grid& grid = obj->get_grid();
	Selector& sel = obj->get_selector();
	MeshObject::position_accessor_t& aaPos = obj->position_accessor();

	IntersectCloseEdges(grid, sel.get_grid_objects(),
						aaPos, snapThreshold);
}

}}// end of namespace

#endif
