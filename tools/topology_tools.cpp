/*
 * Copyright (c) 2013-2017:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#include "topology_tools.h"

namespace ug{
namespace promesh{

void EraseSelectedElements(
			Mesh* obj,
			bool eraseUnusedVrts,
			bool eraseUnusedEdges,
			bool eraseUnusedFaces)
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
size_t RemoveDoubles(Mesh* obj, number threshold)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();

	size_t numVrts = grid.num<Vertex>();
	ug::RemoveDoubles<3>(grid, sel.begin<Vertex>(), sel.end<Vertex>(),
					 	 obj->position_attachment(), threshold);
	return numVrts - grid.num<Vertex>();
}

size_t RemoveDoubleEdges(Mesh* obj)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();

	size_t numEdges= grid.num<Edge>();
	RemoveDuplicates(grid, sel.begin<Edge>(), sel.end<Edge>());
	return numEdges - grid.num<Edge>();
}

size_t RemoveDoubleFaces(Mesh* obj)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();

	size_t numFaces= grid.num<Face>();
	RemoveDuplicates(grid, sel.begin<Face>(), sel.end<Face>());
	return numFaces - grid.num<Face>();
}

void MergeAtFirst(Mesh* obj)
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

void MergeAtCenter(Mesh* obj)
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

void MergeAtLast(Mesh* obj)
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

void CollapseEdge(Mesh* obj)
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

void SplitEdge(Mesh* obj)
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

void SwapEdge(Mesh* obj)
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

void PlaneCut(Mesh* obj, const vector3& p, const vector3& n)
{
	Selector& sel = obj->selector();
	CutEdgesWithPlane(sel, p, n);
}

void AdjustEdgeOrientation(Mesh* obj)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();

	AdjustEdgeOrientationToFaceOrientation(grid, sel.begin<Edge>(),
												sel.end<Edge>());
}

void FixFaceOrientation(Mesh* obj)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();

	FixFaceOrientation(grid, sel.begin<Face>(),
							sel.end<Face>());
}

void FixFaceSubsetOrientations(Mesh* obj)
{
	Grid& grid = obj->grid();
	SubsetHandler& sh = obj->subset_handler();

	for(int i = 0; i < sh.num_subsets(); ++i){
		FixFaceOrientation(grid, sh.begin<Face>(i), sh.end<Face>(i));
	}
}

int FixVolumeOrientation(Mesh* obj)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	Mesh::position_accessor_t& aaPos = obj->position_accessor();

	return FixOrientation(grid, sel.begin<Volume>(), sel.end<Volume>(), aaPos);
}


void InvertFaceOrientation(Mesh* obj)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();

	InvertOrientation(grid, sel.begin<Face>(), sel.end<Face>());
}

void ResolveEdgeIntersection(Mesh* obj, number snapThreshold)
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

void ResolveTriangleIntersections(Mesh* obj, number snapThreshold)
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

void ProjectVerticesToCloseEdges(Mesh* obj, number snapThreshold)
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

void ProjectVerticesToCloseFaces(Mesh* obj, number snapThreshold)
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

void IntersectCloseEdges(Mesh* obj, number snapThreshold)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	Mesh::position_accessor_t& aaPos = obj->position_accessor();

	IntersectCloseEdges(grid, sel, aaPos, snapThreshold);
}


void ResolveSelfIntersections(Mesh* obj, number snapThreshold)
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

}}//	end of namespace
