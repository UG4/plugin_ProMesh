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

#include <vector>
#include "selection_tools.h"
#include "common/util/index_list_util.h"
#include "lib_grid/algorithms/orientation_util.h"
#include "lib_grid/algorithms/selection_util.h"

using namespace std;

namespace ug{
namespace promesh{

void ClearSelection(Mesh* obj)
{
	obj->selector().clear();
}

void SelectAll(Mesh* obj)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	sel.select(grid.vertices_begin(), grid.vertices_end());
	sel.select(grid.edges_begin(), grid.edges_end());
	sel.select(grid.faces_begin(), grid.faces_end());
	sel.select(grid.volumes_begin(), grid.volumes_end());
}

void SelectElementsByIndexRange (Mesh* obj,
                            const char* vrtRanges,
                            const char* edgeRanges,
                            const char* faceRanges,
                            const char* volRanges,
                            bool clearSelection)
{
	Selector& sel = obj->selector();
	if(clearSelection)
		sel.clear();

	vector<size_t> inds;
	RangeStringToIndexList (inds, vrtRanges);
	ug::SelectElementsByIndex<Vertex>(sel, inds);
	RangeStringToIndexList (inds, edgeRanges);
	ug::SelectElementsByIndex<Edge>(sel, inds);
	RangeStringToIndexList (inds, faceRanges);
	ug::SelectElementsByIndex<Face>(sel, inds);
	RangeStringToIndexList (inds, volRanges);
	ug::SelectElementsByIndex<Volume>(sel, inds);
}

void ExtendSelection(Mesh* obj, int neighborhoodSize)
{
	ExtendSelection(obj->selector(), (size_t)neighborhoodSize);
}

void ExtendSelectionInDirection(
		Mesh* obj,
		int neighborhoodSize,
		const vector3& dir,
		number minAngle,
		number maxAngle)
{
	ExtendSelectionInDirection(
			obj->selector(),
			(size_t)neighborhoodSize,
			dir,
			minAngle,
			maxAngle,
			obj->position_accessor());
}


void SelectSubset(
			Mesh* obj,
			int si,
			bool selVrts,
			bool selEdges,
			bool selFaces,
			bool selVols)
{
	SubsetHandler& sh = obj->subset_handler();
	Selector& sel = obj->selector();

	if(si >= 0){
		if(selVrts)
			sel.select(sh.begin<Vertex>(si), sh.end<Vertex>(si));
		if(selEdges)
			sel.select(sh.begin<Edge>(si), sh.end<Edge>(si));
		if(selFaces)
			sel.select(sh.begin<Face>(si), sh.end<Face>(si));
		if(selVols)
			sel.select(sh.begin<Volume>(si), sh.end<Volume>(si));
	}
	else{
		Grid& grid = obj->grid();
	//	subset -1 has to be selected. Those are not directly accessible.
		if(selVrts){
			for(VertexIterator iter = grid.vertices_begin();
				iter != grid.vertices_end(); ++iter)
			{
				if(sh.get_subset_index(*iter) == -1)
					sel.select(*iter);
			}
		}
		if(selEdges){
			for(EdgeIterator iter = grid.edges_begin();
				iter != grid.edges_end(); ++iter)
			{
				if(sh.get_subset_index(*iter) == -1)
					sel.select(*iter);
			}
		}
		if(selFaces){
			for(FaceIterator iter = grid.faces_begin();
				iter != grid.faces_end(); ++iter)
			{
				if(sh.get_subset_index(*iter) == -1)
					sel.select(*iter);
			}
		}
		if(selVols){
			for(VolumeIterator iter = grid.volumes_begin();
				iter != grid.volumes_end(); ++iter)
			{
				if(sh.get_subset_index(*iter) == -1)
					sel.select(*iter);
			}
		}
	}
}

void SelectSubsetBoundary(
			Mesh* obj,
			int si,
			bool edgeBnds,
			bool faceBnds,
			bool volBnds)
{
	SubsetHandler& sh = obj->subset_handler();
	Selector& sel = obj->selector();

	if(edgeBnds)
		SelectAreaBoundary(sel, sh.begin<Edge>(si), sh.end<Edge>(si));
	if(faceBnds)
		SelectAreaBoundary(sel, sh.begin<Face>(si), sh.end<Face>(si));
	if(volBnds)
		SelectAreaBoundary(sel, sh.begin<Volume>(si), sh.end<Volume>(si));
}

template <class TGeomObj>
static void SelectUnassignedElementsHelper(
					Grid& grid,
					SubsetHandler& sh,
					Selector& sel)
{
	typedef typename geometry_traits<TGeomObj>::iterator	iterator;
	for(iterator iter = grid.begin<TGeomObj>(); iter != grid.end<TGeomObj>(); ++iter)
	{
		if(sh.get_subset_index(*iter) == -1){
			sel.select(*iter);
		}
	}
}

template <class TGeomObj>
static void DeselectUnassignedElementsHelper(
					Grid& grid,
					SubsetHandler& sh,
					Selector& sel)
{
	typedef typename geometry_traits<TGeomObj>::iterator	iterator;
	for(iterator iter = grid.begin<TGeomObj>(); iter != grid.end<TGeomObj>(); ++iter)
	{
		if(sh.get_subset_index(*iter) == -1){
			sel.deselect(*iter);
		}
	}
}

void SelectUnassignedElements(
			Mesh* obj,
			bool selVrts,
			bool selEdges,
			bool selFaces,
			bool selVols)
{
	Grid& grid = obj->grid();
	SubsetHandler& sh = obj->subset_handler();
	Selector& sel = obj->selector();

	if(selVrts)
		SelectUnassignedElementsHelper<Vertex>(grid, sh, sel);
	if(selEdges)
		SelectUnassignedElementsHelper<Edge>(grid, sh, sel);
	if(selFaces)
		SelectUnassignedElementsHelper<Face>(grid, sh, sel);
	if(selVols)
		SelectUnassignedElementsHelper<Volume>(grid, sh, sel);
}

void InvertSelection(
			Mesh* obj,
			bool invVrts,
			bool invEdges,
			bool invFaces,
			bool invVols)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();

	if(invVrts)
		InvertSelection(sel, grid.begin<Vertex>(),
							grid.end<Vertex>());

	if(invEdges)
		InvertSelection(sel, grid.begin<Edge>(),
							grid.end<Edge>());

	if(invFaces)
		InvertSelection(sel, grid.begin<Face>(),
							grid.end<Face>());

	if(invVols)
		InvertSelection(sel, grid.begin<Volume>(),
							grid.end<Volume>());
}

void SelectSelectionBoundary(Mesh* obj)
{
	Selector& sel = obj->selector();

	if(sel.num<Volume>() > 0)
		SelectAreaBoundary(sel, sel.begin<Volume>(), sel.end<Volume>());
	if(sel.num<Face>() > 0)
		SelectAreaBoundary(sel, sel.begin<Face>(), sel.end<Face>());
}

void CloseSelection(Mesh* obj)
{
	CloseSelection (obj->selector());
}

template <class TElem>
void RestrictSelectionToSubset (Selector& sel, const SubsetHandler& sh, int si)
{
	typedef typename geometry_traits<TElem>::iterator	iterator;

	for (iterator iter = sel.begin<TElem>(); iter != sel.end<TElem>();){
		TElem* elem = *iter;
		++iter;
		if (sh.get_subset_index(elem) != si)
			sel.deselect(elem);
	}
}

void RestrictSelectionToSubset(Mesh* obj, int si)
{
	SubsetHandler& sh = obj->subset_handler();
	Selector& sel = obj->selector();

	RestrictSelectionToSubset<Vertex>(sel, sh, si);
	RestrictSelectionToSubset<Edge>(sel, sh, si);
	RestrictSelectionToSubset<Face>(sel, sh, si);
	RestrictSelectionToSubset<Volume>(sel, sh, si);
}

////////////////////////////////////////////////////////////////////////////////
//	VERTICES
void SelectBoundaryVertices(Mesh* obj)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	ug::SelectBoundaryElements(sel, grid.begin<Vertex>(), grid.end<Vertex>());
}

void SelectInnerVertices(Mesh* obj)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	SelectInnerElements(sel, grid.begin<Vertex>(), grid.end<Vertex>());
}

void SelectAssociatedVertices(Mesh* obj)
{
	Selector& sel = obj->selector();
	SelectAssociatedVertices(sel, sel.begin<Edge>(), sel.end<Edge>());
	SelectAssociatedVertices(sel, sel.begin<Face>(), sel.end<Face>());
	SelectAssociatedVertices(sel, sel.begin<Volume>(), sel.end<Volume>());
}

void SelectAllVertices(Mesh* obj)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	sel.select(grid.vertices_begin(), grid.vertices_end());
}

void DeselectAllVertices(Mesh* obj)
{
	Selector& sel = obj->selector();
	sel.deselect(sel.vertices_begin(), sel.vertices_end());
}

void SelectMarkedVertices(Mesh* obj)
{
	obj->selector().select(
		obj->crease_handler().begin<Vertex>(REM_FIXED),
		obj->crease_handler().end<Vertex>(REM_FIXED));
}

bool SelectVertexByIndex(Mesh* obj, int index)
{
	Grid& grid = obj->grid();
	int counter = 0;

	for(VertexIterator iter = grid.begin<Vertex>();
		iter != grid.end<Vertex>(); ++iter, ++counter)
	{
		if(counter == index){
			obj->selector().select(*iter);
			return true;
		}
	}
	return false;
}

template <class TElem>
static size_t
SelectUnconnectedVerticesHelper(Grid& grid, Selector& sel)
{
	using namespace ug;
	typename Grid::traits<TElem>::secure_container	elems;

	size_t numUnconnected = 0;
	for(VertexIterator iter = grid.vertices_begin();
		iter != grid.vertices_end(); ++iter)
	{
		if(!sel.is_selected(*iter)){
			grid.associated_elements(elems, *iter);
			if(elems.size() == 0){
				sel.select(*iter);
				numUnconnected++;
			}
		}
	}
	return numUnconnected;
}

size_t SelectUnconnectedVertices(
				Mesh* obj,
				bool edgeCons,
				bool faceCons,
				bool volCons)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();

	size_t numUnconnected = 0;
	if(edgeCons)
		numUnconnected += SelectUnconnectedVerticesHelper<Edge>(grid, sel);
	if(faceCons)
		numUnconnected += SelectUnconnectedVerticesHelper<Face>(grid, sel);
	if(volCons)
		numUnconnected += SelectUnconnectedVerticesHelper<Volume>(grid, sel);
	return numUnconnected;
}

size_t SelectSelectionKinkVertices(
				Mesh* obj,
				number thresholdAngle,
				bool selectDarts)
{
	Selector& sel = obj->selector();
	std::vector<Vertex*> candidates;
	CollectVerticesTouchingSelection(candidates, sel);
	size_t numSel = sel.num<Vertex>();
	
	SelectKinkVertices(sel, candidates.begin(), candidates.end(),
						 thresholdAngle, selectDarts, obj->position_accessor(),
						 IsSelected(sel));

	return sel.num<Vertex>() - numSel;
}

size_t SelectSubsetKinkVertices(
				Mesh* obj,
				int subsetIndex,
				number thresholdAngle,
				bool selectDarts)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	SubsetHandler& sh = obj->subset_handler();

	std::vector<Vertex*> candidates;
	grid.begin_marking();
	for(EdgeIterator eiter = sh.begin<Edge>(subsetIndex);
		eiter != sh.end<Edge>(subsetIndex); ++eiter)
	{
		Edge* e = *eiter;
		for(size_t i = 0; i < 2; ++i){
			if(!grid.is_marked(e->vertex(i))){
				candidates.push_back(e->vertex(i));
				grid.mark(e->vertex(i));
			}
		}
	}
	grid.end_marking();

	size_t numSel = sel.num<Vertex>();
	
	SelectKinkVertices(sel, candidates.begin(), candidates.end(),
						 thresholdAngle, selectDarts, obj->position_accessor(),
						 IsInSubset(obj->subset_handler(), subsetIndex));

	return sel.num<Vertex>() - numSel;
}

////////////////////////////////////////////////////////////////////////////////
//	EDGES
void SelectBoundaryEdges(Mesh* obj)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	SelectBoundaryElements(sel, grid.begin<Edge>(), grid.end<Edge>());
}

void SelectInnerEdges(Mesh* obj)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	SelectInnerElements(sel, grid.begin<Edge>(), grid.end<Edge>());
}

void SelectNonManifoldEdges(Mesh* obj)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();

//	iterate over all edges and check how many adjacent faces each has.
//	if there are more than 2, the edge will be selected
	for(EdgeIterator iter = grid.edges_begin();
		iter != grid.edges_end(); ++iter)
	{
		if(NumAssociatedFaces(grid, *iter) != 2)
			sel.select(*iter);
	}
}

void SelectSmoothEdgePath(
			Mesh* obj,
			number maxDeviation,
			number normalWeight,
			bool stopAtSelectedVrts)
{
	SelectSmoothEdgePath(obj->selector(), maxDeviation, normalWeight, stopAtSelectedVrts);
}

void SelectShortEdges(Mesh* obj, number maxLength)
{
	number maxLengthSq = maxLength * maxLength;
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	Mesh::position_accessor_t& aaPos = obj->position_accessor();

	for(EdgeIterator iter = grid.begin<Edge>();
		iter != grid.end<Edge>(); ++iter)
	{
		Edge* e = *iter;
		if(VecDistanceSq(aaPos[e->vertex(0)], aaPos[e->vertex(1)]) < maxLengthSq)
			sel.select(e);
	}
}

void SelectLongEdges(Mesh* obj, number minLength)
{
	number minLengthSq = minLength * minLength;
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	Mesh::position_accessor_t& aaPos = obj->position_accessor();

	for(EdgeIterator iter = grid.begin<Edge>();
		iter != grid.end<Edge>(); ++iter)
	{
		Edge* e = *iter;
		if(VecDistanceSq(aaPos[e->vertex(0)], aaPos[e->vertex(1)]) > minLengthSq)
			sel.select(e);
	}
}

void SelectCreaseEdges(Mesh* obj, number minAngle)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	SelectCreaseEdges(sel, grid.begin<Edge>(), grid.end<Edge>(),
					  minAngle, obj->position_attachment());
}

void SelectLinkedBoundaryEdges(Mesh* obj, bool stopAtSelectedVrts)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();

	if(stopAtSelectedVrts)
		SelectLinkedElements<Edge>(sel, IsOnBoundary(grid), IsNotSelected(sel));
	else
		SelectLinkedElements<Edge>(sel, IsOnBoundary(grid));
}

void SelectAssociatedEdges(Mesh* obj)
{
	Selector& sel = obj->selector();
	SelectAssociatedEdges(sel,sel.begin<Face>(), sel.end<Face>());
	SelectAssociatedEdges(sel, sel.begin<Volume>(), sel.end<Volume>());
}

void SelectAllEdges(Mesh* obj)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	sel.select(grid.edges_begin(), grid.edges_end());
}

void DeselectAllEdges(Mesh* obj)
{
	Selector& sel = obj->selector();
	sel.deselect(sel.edges_begin(), sel.edges_end());
}

void SelectMarkedEdges(Mesh* obj)
{
	obj->selector().select(
		obj->crease_handler().begin<Edge>(REM_CREASE),
		obj->crease_handler().end<Edge>(REM_CREASE));
}

bool SelectEdgeByIndex(Mesh* obj, int index)
{
	Grid& grid = obj->grid();
	int counter = 0;

	for(EdgeIterator iter = grid.begin<Edge>();
		iter != grid.end<Edge>(); ++iter, ++counter)
	{
		if(counter == index){
			obj->selector().select(*iter);
			return true;
		}
	}
	return false;
}

void EdgeSelectionFill(Mesh* obj)
{
	Selector& sel = obj->selector();
	SelectionFill<Edge>(sel, IsSelected(sel));
}

void SelectShortPolychains(Mesh* m, number maxChainLength, bool closedChainsOnly)
{
	SelectShortPolychains(m->selector(), maxChainLength,
						  closedChainsOnly, m->position_accessor());
}

void SelectEdgesByDirection(
				Mesh* m,
				const vector3& dir,
				number minDeviationAngle,
				number maxDeviationAngle,
				bool selectFlipped)
{
	Selector& sel = m->selector();
	Mesh::position_accessor_t aaPos = m->position_accessor();
	SelectEdgesByDirection(sel, aaPos, dir, minDeviationAngle,
	                       maxDeviationAngle, selectFlipped);
}

void SelectSubsetEdgesByDirection(
				Mesh* m,
				int subsetIndex,
				const vector3& dir,
				number minDeviationAngle,
				number maxDeviationAngle,
				bool selectFlipped)
{
	Selector& sel = m->selector();
	SubsetHandler& sh = m->subset_handler();
	Mesh::position_accessor_t aaPos = m->position_accessor();
	SelectSubsetEdgesByDirection(sel, sh, subsetIndex, aaPos, dir, minDeviationAngle,
	                       maxDeviationAngle, selectFlipped);
}

////////////////////////////////////////////////////////////////////////////////
//	FACES
void SelectBoundaryFaces(Mesh* obj)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	SelectBoundaryElements(sel, grid.begin<Face>(), grid.end<Face>());
}

void SelectInnerFaces(Mesh* obj)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	SelectInnerElements(sel, grid.begin<Face>(), grid.end<Face>());
}

void SelectMarkedFaces(Mesh* obj)
{
	obj->selector().select(
		obj->crease_handler().begin<Face>(REM_CREASE),
		obj->crease_handler().end<Face>(REM_CREASE));
}

void SelectLinkedManifoldFaces(Mesh* obj)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
//	select associated faces of selected elements
	SelectAssociatedFaces(sel, sel.vertices_begin(), sel.vertices_end());
	SelectAssociatedFaces(sel, sel.edges_begin(), sel.edges_end());
	SelectAssociatedFaces(sel, sel.volumes_begin(), sel.volumes_end());

//	push all selected faces to a stack
	std::stack<Face*> stk;
	for(FaceIterator iter = sel.faces_begin(); iter != sel.faces_end(); ++iter){
		stk.push(*iter);
	}

//	while there are faces in the stack, get their associated edges.
//	if those edges are adjacent to exactly two faces, then select the
//	neighboured face and push it to the stack (if it was unselected previously)
	std::vector<Edge*> edges;
	while(!stk.empty()){
		Face* f = stk.top();
		stk.pop();
		CollectEdges(edges, grid, f);

		for(size_t i = 0; i < edges.size(); ++i){
			Face* faces[2];
			if(GetAssociatedFaces(faces, grid, edges[i], 2) == 2){
				for(size_t j = 0; j < 2; ++j){
					if(!sel.is_selected(faces[j])){
						sel.select(faces[j]);
						stk.push(faces[j]);
					}
				}
			}
		}
	}
}


void SelectLinkedBoundaryFaces(Mesh* obj, bool stopAtSelectedEdges)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();

	if(stopAtSelectedEdges)
		SelectLinkedElements<Face>(sel, IsOnBoundary(grid), IsNotSelected(sel));
	else
		SelectLinkedElements<Face>(sel, IsOnBoundary(grid));
}

void SelectDegenerateFaces(Mesh* obj, number maxHeight)
{
	number maxHeightSq = maxHeight * maxHeight;
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	Mesh::position_accessor_t& aaPos = obj->position_accessor();

	for(FaceIterator iter = grid.begin<Face>();
		iter != grid.end<Face>(); ++iter)
	{
		Face* f = *iter;
	//	iterate over the edges and check which is the longest.
		number maxLenSq = -1;
		size_t bestInd = -1;
		EdgeDescriptor ed;
		for(size_t i = 0; i < f->num_edges(); ++i){
			f->edge_desc(i, ed);
			number lenSq = VecDistanceSq(aaPos[ed.vertex(0)], aaPos[ed.vertex(1)]);
			if(lenSq > maxLenSq){
				maxLenSq = lenSq;
				bestInd = i;
			}
		}

		if(maxLenSq < maxHeightSq)
			sel.select(f);
		else{
		//	project the other vertices to the line and check the height
		//todo: this is not enough for quadrilaterals
			vector3 p;
			vector3& v = aaPos[f->vertex((bestInd + 2) % f->num_vertices())];
			f->edge_desc(bestInd, ed);
			number t;
			if(DistancePointToLine(t, v, aaPos[ed.vertex(0)], aaPos[ed.vertex(1)]) < maxHeight)
				sel.select(f);
		}
	}
}

void SelectLinkedFlatFaces(
			Mesh* obj,
			number maxDeviationAngle,
			bool ignoreOrientation,
			bool traverseDegeneratedFaces,
			bool stopAtSelectedEdges)
{
	Selector& sel = obj->selector();
	if(traverseDegeneratedFaces)
		ug::SelectLinkedFlatAndDegeneratedFaces(sel, maxDeviationAngle,
												ignoreOrientation,
												stopAtSelectedEdges);
	else
		ug::SelectLinkedFlatFaces(sel, maxDeviationAngle, ignoreOrientation,
								  stopAtSelectedEdges);
}

void SelectIntersectingTriangles(Mesh* obj)
{
	using namespace ug::node_tree;

	Grid& grid = obj->grid();
	Selector& sel = obj->selector();

//	make sure that triangles are present.
	if(grid.num<Triangle>() == 0){
		UG_LOG("Given grid does not contain any triangles. Aborting 'Select Intersecting Triangles'.\n");
		return;
	}

	if(grid.num<Quadrilateral>() > 0){
		UG_LOG("WARNING: Quadrilateral intersections are ignored during 'Select Intersecting Triangles'.\n");
	}

//	create an octree
//	sort the triangles of grid into an octree to speed-up projection performance
	SPOctree octree;
	octree = CreateOctree(grid, grid.begin<Triangle>(),
								grid.end<Triangle>(),
								10, 30, false, obj->position_attachment());

	if(!octree.valid()){
		UG_LOG("  Octree creation failed in 'Select Intersecting Triangles'. Aborting.\n");
		return;
	}

//	access the position attachment of the grid
	Mesh::position_accessor_t& aaPos = obj->position_accessor();

//	we'll use a traverser to find the intersections.
	Traverser_IntersectFaces intersectionTraverser;

//	for each face the face itself and direct neighbors shall be ignored
	std::vector<Face*> ignoreList;

	size_t triCount = 0;
//	now iterate over all triangles of the grid and find intersections
	for(TriangleIterator iter = grid.begin<Triangle>();
		iter != grid.end<Triangle>(); ++iter, ++triCount)
	{
		Triangle* t = *iter;

	//	add neighbors and the face itself to the ignore list
		CollectNeighbors(ignoreList, t, grid, NHT_VERTEX_NEIGHBORS);
		intersectionTraverser.clear_ignore_list();
		for(size_t i = 0; i < ignoreList.size(); ++i){
			intersectionTraverser.ignore_element(ignoreList[i]);
		}
		intersectionTraverser.ignore_element(t);

		if(intersectionTraverser.intersect_tri(aaPos[t->vertex(0)], aaPos[t->vertex(1)],
												aaPos[t->vertex(2)], octree))
		{
		//	check the intersecting face:
			/*
			const std::vector<CollisionElementID>& faces =
				intersectionTraverser.get_intersected_element_ids();
			*/
		//	an intersection occurred. Log the index and select the triangle.
			sel.select(t);
		}
	}

}

void SelectAssociatedFaces(Mesh* obj)
{
	Selector& sel = obj->selector();
	SelectAssociatedFaces(sel, sel.begin<Volume>(), sel.end<Volume>());
}

void SelectAllFaces(Mesh* obj)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	sel.select(grid.faces_begin(), grid.faces_end());
}

void DeselectAllFaces(Mesh* obj)
{
	Selector& sel = obj->selector();
	sel.deselect(sel.faces_begin(), sel.faces_end());
}

bool SelectFaceByIndex(Mesh* obj, int index)
{
	Grid& grid = obj->grid();
	int counter = 0;

	for(FaceIterator iter = grid.begin<Face>();
		iter != grid.end<Face>(); ++iter, ++counter)
	{
		if(counter == index){
			obj->selector().select(*iter);
			return true;
		}
	}
	return false;
}

void SelectFacesByNormal(
			Mesh* obj,
			const vector3& refNormal,
			number maxDeviationAngle)
{
	number dotThreshold = cos(deg_to_rad(maxDeviationAngle));
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	Grid::AttachmentAccessor<Vertex, APosition> aaPos(grid, aPosition);

	vector3 n;
	VecNormalize(n, refNormal);

	for(FaceIterator iter = grid.faces_begin(); iter != grid.faces_end(); ++iter){
		vector3 fn;
		CalculateNormal(fn, *iter, aaPos);
		if(VecDot(fn, n) >= dotThreshold)
			sel.select(*iter);
	}
}


void SelectFacesByNormal(
			Mesh* obj,
			const vector3& refNormal,
			number minDeviationAngle,
			number maxDeviationAngle,
			bool noInnerFaces)
{
	number minDot = cos(deg_to_rad(maxDeviationAngle));
	number maxDot = cos(deg_to_rad(minDeviationAngle));

	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	Grid::AttachmentAccessor<Vertex, APosition> aaPos(grid, aPosition);

	vector3 n;
	VecNormalize(n, refNormal);

	for(FaceIterator iter = grid.faces_begin(); iter != grid.faces_end(); ++iter){
		vector3 fn;
		CalculateNormal(fn, *iter, aaPos);
		const number d = VecDot(fn, n);
		if((d > minDot - SMALL && d < maxDot + SMALL)
		   && (!noInnerFaces || (NumAssociatedVolumes(grid, *iter) < 2)))
		{
			sel.select(*iter);
		}
	}
}


void FaceSelectionFill(Mesh* obj)
{
	Selector& sel = obj->selector();
	SelectionFill<Face>(sel, IsSelected(sel));
}

size_t SelectBentQuadrilaterals(Mesh* obj, number dotThreshold)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	Grid::AttachmentAccessor<Vertex, APosition> aaPos(grid, aPosition);

//	iterate over all quadrilaterals and search for bent ones
	size_t selCount = 0;
	for(QuadrilateralIterator iter = grid.begin<Quadrilateral>();
		iter != grid.end<Quadrilateral>(); ++iter)
	{
		Quadrilateral* q = *iter;

	//	we'll compare the dot-product of the normals
		vector3 n1, n2;
		CalculateTriangleNormal(n1, aaPos[q->vertex(0)],
								aaPos[q->vertex(1)], aaPos[q->vertex(2)]);
		CalculateTriangleNormal(n2, aaPos[q->vertex(2)],
								aaPos[q->vertex(3)], aaPos[q->vertex(0)]);

		number d1 = VecDot(n1, n2);

		if(d1 < dotThreshold){
			++selCount;
			sel.select(q);
		}
	}

	return selCount;
}


////////////////////////////////////////////////////////////////////////////////
//	VOLUMES
void SelectAllVolumes(Mesh* obj)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	sel.select(grid.volumes_begin(), grid.volumes_end());
}

void DeselectAllVolumes(Mesh* obj)
{
	Selector& sel = obj->selector();
	sel.deselect(sel.volumes_begin(), sel.volumes_end());
}

int SelectUnorientableVolumes(Mesh* obj)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();

	Grid::AttachmentAccessor<Vertex, APosition> aaPos(grid, aPosition);

	int numUnorientable = 0;

//	iterate over all volumes. Check whether the orientation can be determined.
	for(VolumeIterator iter = grid.volumes_begin();
		iter != grid.volumes_end(); ++iter)
	{
		Volume* v = *iter;
	//	get orientation of original volume
		bool bOriented = CheckOrientation(v, aaPos);
	//	flip the volume
		grid.flip_orientation(v);

	//	if orientations match, then the volume can not be oriented.
		if(bOriented == CheckOrientation(v, aaPos)){
			sel.select(v);
			++numUnorientable;
		}

	//	reflip orientation
		grid.flip_orientation(v);
	}

	return numUnorientable;
}

int SelectSlivers(Mesh* obj, number thresholdRatio)
{
	std::vector<Tetrahedron*> slivers;
	Grid& g = obj->grid();
	FindSlivers(slivers, g.begin<Tetrahedron>(), g.end<Tetrahedron>(),
			    thresholdRatio, obj->position_accessor());
	obj->selector().select(slivers.begin(), slivers.end());
	return static_cast<int>(slivers.size());
}

bool SelectVolumeByIndex(Mesh* obj, int index)
{
	Grid& grid = obj->grid();
	int counter = 0;

	for(VolumeIterator iter = grid.begin<Volume>();
		iter != grid.end<Volume>(); ++iter, ++counter)
	{
		if(counter == index){
			obj->selector().select(*iter);
			return true;
		}
	}
	return false;
}

void SelectVolumesByType(
			Mesh* obj,
			bool selHexahedra,
			bool selOctahedra,
			bool selPrisms,
			bool selPyramids,
			bool selTetrahedra)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();

	if(selHexahedra)
		sel.select(grid.begin<Hexahedron>(), grid.end<Hexahedron>());
	if(selOctahedra)
		sel.select(grid.begin<Octahedron>(), grid.end<Octahedron>());
	if(selPrisms)
		sel.select(grid.begin<Prism>(), grid.end<Prism>());
	if(selPyramids)
		sel.select(grid.begin<Pyramid>(), grid.end<Pyramid>());
	if(selTetrahedra)
		sel.select(grid.begin<Tetrahedron>(), grid.end<Tetrahedron>());
}

void VolumeSelectionFill(Mesh* obj)
{
	Selector& sel = obj->selector();
	SelectionFill<Volume>(sel, IsSelected(sel));
}


void ClearMarks(Mesh* obj)
{
	obj->crease_handler().clear();
}

void MarkCornersOfMarkedEdges(Mesh* obj, number angle)
{
	Grid& grid = obj->grid();
	MarkCorners(obj->grid(), obj->crease_handler(),
				grid.begin<Vertex>(), grid.end<Vertex>(),
				IsInSubset(obj->crease_handler(), ug::REM_CREASE),
				ug::REM_FIXED, angle, obj->position_attachment());
}


void MarkSelection(Mesh* obj)
{
	obj->crease_handler().assign_subset(
		obj->selector().begin<Vertex>(),
		obj->selector().end<Vertex>(),
		ug::REM_FIXED);
	obj->crease_handler().assign_subset(
		obj->selector().begin<Edge>(),
		obj->selector().end<Edge>(),
		ug::REM_CREASE);
	obj->crease_handler().assign_subset(
		obj->selector().begin<Face>(),
		obj->selector().end<Face>(),
		ug::REM_CREASE);
}

void UnmarkSelection(Mesh* obj)
{
	obj->crease_handler().assign_subset(
		obj->selector().begin<Vertex>(),
		obj->selector().end<Vertex>(),
		ug::REM_NONE);

	obj->crease_handler().assign_subset(
		obj->selector().begin<Edge>(),
		obj->selector().end<Edge>(),
		ug::REM_NONE);
}

void MarkCreaseEdges(Mesh* obj, number minAngle, bool clearMarks)
{
	if(clearMarks)
		obj->crease_handler().clear();

	MarkCreaseEdges(obj->grid(),
					obj->crease_handler(),
					obj->grid().begin<Edge>(),
					obj->grid().end<Edge>(),
					REM_CREASE, minAngle);
	MarkFixedCreaseVertices(obj->grid(),
							obj->crease_handler(),
							REM_CREASE, REM_FIXED);
}


template <class elem_t, class vector_t, class AAPos>
static void
SelectElementsBySplitPlane_IMPL(
		Grid& g,
		Selector& sel,
		const vector_t& pivot,
		const vector_t& normal,
		AAPos aaPos)
{
	typedef typename Grid::traits<elem_t>::iterator	iter_t;

	for(iter_t iter = g.begin<elem_t>(); iter != g.end<elem_t>(); ++ iter) {
		elem_t* elem = *iter;
		vector_t p = CalculateCenter(elem, aaPos);
		vector_t dir;
		VecSubtract(dir, p, pivot);
		if(VecDot(dir, normal) >= 0)
			sel.select(elem);
	}
}

void SelectElementsBySplitPlane(
			Mesh* obj,
			bool selectVrts,
			bool selectEdges,
			bool selectFaces,
			bool selectVols,
			const vector3& pivot,
			const vector3& normal)
{
	Grid& g = obj->grid();
	Selector& sel = obj->selector();
	Mesh::position_accessor_t aaPos = obj->position_accessor();

	if(selectVrts)
		SelectElementsBySplitPlane_IMPL<Vertex>(g, sel, pivot, normal, aaPos);
	if(selectEdges)
		SelectElementsBySplitPlane_IMPL<Edge>(g, sel, pivot, normal, aaPos);
	if(selectFaces)
		SelectElementsBySplitPlane_IMPL<Face>(g, sel, pivot, normal, aaPos);
	if(selectVols)
		SelectElementsBySplitPlane_IMPL<Volume>(g, sel, pivot, normal, aaPos);
}

}}//	end of namespace
