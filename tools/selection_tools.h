/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__selection_tools__
#define __H__UG__selection_tools__

#include <vector>
#include <stack>
#include "../mesh.h"
#include "common/node_tree/node_tree.h"
#include "lib_grid/algorithms/crease_util.h"
#include "lib_grid/algorithms/mark_util.h"
#include "lib_grid/algorithms/problem_detection_util.h"
#include "lib_grid/algorithms/refinement_mark_util.h"
#include "lib_grid/algorithms/selection_util.h"
#include "lib_grid/algorithms/remeshing/edge_length_adjustment.h"
#include "lib_grid/algorithms/trees/octree.h"
#include "lib_grid/callbacks/callbacks.h"

//selection tools
#define	TOOLTIP_SELECT_LINKED_MANIFOLD_FACES "Selects faces linked with the selection, not crossing non-manifold edges." //?
#define	TOOLTIP_SELECT_NON_MANIFOLD_EDGES "Selects edges with more than 2 associated faces."
#define	TOOLTIP_CLEAR_SELECTION "Clears the selection"
#define	TOOLTIP_SELECT_SMOOTH_EDGE_PATH "Selects a smooth edge path."
#define	TOOLTIP_SELECT_BOUNDARY_VERTICES "Selects vertices that lie on the boundary of the geometry"
#define	TOOLTIP_SELECT_INNER_VERTICES "Selects vertices that do not lie on the boundary of the geometry"
#define	TOOLTIP_SELECT_BOUNDARY_EDGES "Selects edges that lie on the boundary of the geometry"
#define	TOOLTIP_SELECT_INNER_EDGES "Selects edges that do not lie on the boundary of the geometry"
#define	TOOLTIP_SELECT_BOUNDARY_FACES "Selects faces that lie on the boundary of the geometry"
#define	TOOLTIP_SELECT_INNER_FACES "Selects faces that do not lie on the boundary of the geometry"
#define	TOOLTIP_SELECT_SHORT_EDGES "Selects edges that are shorter than a given threshold."
#define	TOOLTIP_SELECT_LONG_EDGES "Selects edges that are longer than a given threshold."
#define	TOOLTIP_SELECT_CREASE_EDGES "Selects edges that at which triangles meet in a given angle."
#define	TOOLTIP_SELECT_DEGENERATE_FACES "Selects faces that have a height shorter than a given threshold."
#define	TOOLTIP_SELECT_LINKED_FLAT_FACES "Selects linked faces of selected faces that have a similar normal."
#define	TOOLTIP_SELECT_LINKED_BOUNDARY_EDGES "Selects linked boundary edges of selected edges."
#define	TOOLTIP_SELECT_LINKED_BOUNDARY_FACES "Selects linked boundary faces of selected faces."
#define	TOOLTIP_SELECT_INTERSECTING_TRIANGLES "Selects intersecting triangles. Neighbors are ignored."
#define	TOOLTIP_SELECT_ASSOCIATED_VERTICES "Selects vertices that belong to selected edges, faces and volumes."
#define	TOOLTIP_SELECT_ASSOCIATED_EDGES "Selects edges that belong to selected faces and volumes."
#define	TOOLTIP_SELECT_ASSOCIATED_FACES "Selects faces that belong to selected volumes."
#define	TOOLTIP_SELECT_ALL "Selects all vertices, edges ,faces and volumes of the current grid"
#define	TOOLTIP_SELECT_ALL_VERTICES "Selects all vertices of the current grid"
#define	TOOLTIP_DESELECT_ALL_VERTICES "Deselects all vertices of the current grid"
#define	TOOLTIP_SELECT_ALL_EDGES "Selects all edges of the current grid"
#define	TOOLTIP_DESELECT_ALL_EDGES "Deselects all edges of the current grid"
#define	TOOLTIP_SELECT_ALL_FACES "Selects all faces of the current grid"
#define	TOOLTIP_DESELECT_ALL_FACES "Deselects all faces of the current grid"
#define	TOOLTIP_SELECT_ALL_VOLUMES "Selects all volumes of the current grid"
#define	TOOLTIP_DESELECT_ALL_VOLUMES "Deselects all volumes of the current grid"
#define	TOOLTIP_SELECT_MARKED_VERTICES "Selects vertices which are marked."
#define	TOOLTIP_SELECT_MARKED_EDGES "Selects edges which are marked."
#define	TOOLTIP_SELECT_UNORIENTABLE_VOLUMES "Selects all volumes whose orientation can not be determined"
#define	TOOLTIP_EXTEND_SELECTION "Selects neighbors of selected elements."
#define	TOOLTIP_SELECT_VERTEX_BY_INDEX "Selects a vertex given its index."
#define	TOOLTIP_SELECT_EDGE_BY_INDEX "Selects a edge given its index."
#define	TOOLTIP_SELECT_FACE_BY_INDEX "Selects a face given its index."
#define	TOOLTIP_SELECT_FACES_BY_NORMAL "Selects faces given a normal and a maximum deviation angle."
#define	TOOLTIP_SELECT_VOLUME_BY_INDEX "Selects a volume given its index."
#define	TOOLTIP_SELECT_VERTEX_BY_COORDINATE "Selects a vertex given a coordinate."
#define	TOOLTIP_SELECT_EDGE_BY_COORDINATE  "Selects the edge whose center is closest to the specified coordinate."
#define	TOOLTIP_SELECT_FACE_BY_COORDINATE "Selects the face whose center is closest to the specified coordinate."
#define	TOOLTIP_SELECT_VOLUME_BY_COORDINATE "Selects the volume whose center is closest to the specified coordinate."
#define	TOOLTIP_SELECT_VERTEX_BY_CYL_COORDINATE "Selects a vertex given a cylindrical coordinate."
#define	TOOLTIP_SELECT_EDGE_BY_CYL_COORDINATE  "Selects the edge whose center is closest to the specified cylindrical coordinate."
#define	TOOLTIP_SELECT_FACE_BY_CYL_COORDINATE "Selects the face whose center is closest to the specified cylindrical coordinate."
#define	TOOLTIP_SELECT_VOLUME_BY_CYL_COORDINATE "Selects the volume whose center is closest to the specified cylindrical coordinate."
#define	TOOLTIP_SELECT_UNCONNECTED_VERTICES "Selects vertices which are not connected to the given element type."
#define	TOOLTIP_SELECT_SUBSET "Selects all elements of a subset."
#define	TOOLTIP_SELECT_SUBSET_BOUNDARY "Selects the boundary of a subset."
#define	TOOLTIP_SELECT_UNASSIGNED_ELEMENTS "Selects all elements not assigned to any subset."
#define	TOOLTIP_INVERT_SELECTION "Inverts current selection."
#define	TOOLTIP_EDGE_SELECTION_FILL "Selects neighbours of selected edges over non-selected vertices."
#define	TOOLTIP_FACE_SELECTION_FILL "Selects neighbours of selected faces over non-selected edges."
#define	TOOLTIP_VOLUME_SELECTION_FILL "Selects neighbours of selected volumes over non-selected faces."
#define	TOOLTIP_SELECT_SELECTION_BOUNDARY "Selects the boundary of the current selection."
#define	TOOLTIP_SELECT_BENT_QUADRILATERALS "Selects quadrilaterals which do not lie in a plane."
#define	TOOLTIP_CLOSE_SELECTION "Selects all associated elements of lower dimensions."
#define TOOLTIP_SELECT_SLIVERS "Selects flat tetrahedrons. Threshold-ratio specifies the minimal ratio between the distance of two opposing edges to the length of the longest edge."
#define TOOLTIP_SELECT_SELECTION_KINK_VERTICES "Selects kink vertices in selected paths"
#define TOOLTIP_SELECT_SUBSET_KINK_VERTICES "Selects kink vertices in subset-paths"
#define TOOLTIP_SELECT_LINKED_EDGES "Repeatedly selects all edges which are vertex-neighbors of selected edges."
#define TOOLTIP_SELECT_LINKED_FACES "Repeatedly selects all faces which are edge-neighbors of selected faces."
#define TOOLTIP_SELECT_LINKED_VOLUMES "Repeatedly selects all volumes which are face-neighbors of selected volumes."
#define TOOLTIP_SELECT_SHORT_POLYCHAINS "Selects polygonal chains which are shorter than the given threshold."
#define TOOLTIP_SELECT_INTERFACE_ELEMENTS "Selects elements which are adjacent to higher dimensional elements of different subsets."
#define TOOLTIP_SELECT_ANISOTROPIC_ELEMENTS "Selects elements and associated long edges wich have a shortest-to-longest edge ratio smaller than the specified one."

namespace ug{
namespace promesh{

inline void ClearSelection(Mesh* obj)
{
	obj->selector().clear();
}

inline void SelectAll(Mesh* obj)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	sel.select(grid.vertices_begin(), grid.vertices_end());
	sel.select(grid.edges_begin(), grid.edges_end());
	sel.select(grid.faces_begin(), grid.faces_end());
	sel.select(grid.volumes_begin(), grid.volumes_end());
}

inline void ExtendSelection(Mesh* obj, int neighborhoodSize)
{
	Selector& sel = obj->selector();
	ExtendSelection(sel, (size_t)neighborhoodSize);
}

template <class TElem>
TElem* SelectElemByCoordinate(Mesh* obj, const vector3& coord)
{
	Grid& grid = obj->grid();
	Mesh::position_accessor_t& aaPos = obj->position_accessor();

	TElem* e = FindClosestByCoordinate<TElem>(coord,
											 grid.begin<TElem>(),
											 grid.end<TElem>(),
											 aaPos);

	if(e)
		obj->selector().select(e);

	return e;
}

template <class TElem>
TElem* SelectElemByCylindricalCoordinate(Mesh* obj, number rho, number phi, number z)
{
	number x = rho * cos(phi);
	number y = rho * sin(phi);

	vector3 coord = vector3(x,y,z);
	return SelectElemByCoordinate<TElem>(obj, coord);
}

inline void SelectSubset(Mesh* obj, int si, bool selVrts, bool selEdges,
				  		 bool selFaces, bool selVols)
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

inline void SelectSubsetBoundary(Mesh* obj, int si, bool edgeBnds,
						  bool faceBnds, bool volBnds)
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
static inline void SelectUnassignedElementsHelper(Grid& grid, SubsetHandler& sh, Selector& sel)
{
	typedef typename geometry_traits<TGeomObj>::iterator	iterator;
	for(iterator iter = grid.begin<TGeomObj>(); iter != grid.end<TGeomObj>(); ++iter)
	{
		if(sh.get_subset_index(*iter) == -1){
			sel.select(*iter);
		}
	}
}

inline void SelectUnassignedElements(Mesh* obj, bool selVrts, bool selEdges,
							  bool selFaces, bool selVols)
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

inline void InvertSelection(Mesh* obj, bool invVrts, bool invEdges,
					 bool invFaces, bool invVols)
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

inline void SelectSelectionBoundary(Mesh* obj)
{
	Selector& sel = obj->selector();

	if(sel.num<Volume>() > 0)
		SelectAreaBoundary(sel, sel.begin<Volume>(), sel.end<Volume>());
	if(sel.num<Face>() > 0)
		SelectAreaBoundary(sel, sel.begin<Face>(), sel.end<Face>());
}

inline void CloseSelection(Mesh* obj)
{
	Selector& sel = obj->selector();
	SelectAssociatedFaces(sel, sel.begin<Volume>(), sel.end<Volume>());
	SelectAssociatedEdges(sel, sel.begin<Face>(), sel.end<Face>());
	SelectAssociatedVertices(sel, sel.begin<Edge>(), sel.end<Edge>());
}


template <class TElem>
inline void SelectInterfaceElements(Mesh* obj, bool regardSelectedNbrsOnly)
{
	Grid& g = obj->grid();
	SelectInterfaceElements(obj->selector(), obj->subset_handler(),
							g.begin<TElem>(), g.end<TElem>(),
							regardSelectedNbrsOnly);
}


template <class TElem>
inline void SelectAnisotropicElements(Mesh* obj, number minEdgeRatio)
{
	Grid& g = obj->grid();
	MarkForAnisotropicRefinement(
		g,
		obj->selector(),
		minEdgeRatio,
		g.begin<TElem>(),
		g.end<TElem>(),
		obj->position_accessor());
}


////////////////////////////////////////////////////////////////////////////////
//	VERTICES
inline void SelectBoundaryVertices(Mesh* obj)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	ug::SelectBoundaryElements(sel, grid.begin<Vertex>(), grid.end<Vertex>());
}

inline void SelectInnerVertices(Mesh* obj)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	SelectInnerElements(sel, grid.begin<Vertex>(), grid.end<Vertex>());
}

inline void SelectAssociatedVertices(Mesh* obj)
{
	Selector& sel = obj->selector();
	SelectAssociatedVertices(sel, sel.begin<Edge>(), sel.end<Edge>());
	SelectAssociatedVertices(sel, sel.begin<Face>(), sel.end<Face>());
	SelectAssociatedVertices(sel, sel.begin<Volume>(), sel.end<Volume>());
}

inline void SelectAllVertices(Mesh* obj)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	sel.select(grid.vertices_begin(), grid.vertices_end());
}

inline void DeselectAllVertices(Mesh* obj)
{
	Selector& sel = obj->selector();
	sel.deselect(sel.vertices_begin(), sel.vertices_end());
}

inline void SelectMarkedVertices(Mesh* obj)
{
	obj->selector().select(
		obj->crease_handler().begin<Vertex>(REM_FIXED),
		obj->crease_handler().end<Vertex>(REM_FIXED));
}

inline bool SelectVertexByIndex(Mesh* obj, int index)
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
inline size_t SelectUnconnectedVerticesHelper(Grid& grid, Selector& sel)
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

inline size_t SelectUnconnectedVertices(Mesh* obj, bool edgeCons, bool faceCons,
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

inline size_t SelectSelectionKinkVertices(Mesh* obj, number thresholdAngle,
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

inline size_t SelectSubsetKinkVertices(Mesh* obj, int subsetIndex,
									   number thresholdAngle, bool selectDarts)
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
inline void SelectBoundaryEdges(Mesh* obj)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	SelectBoundaryElements(sel, grid.begin<Edge>(), grid.end<Edge>());
}

inline void SelectInnerEdges(Mesh* obj)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	SelectInnerElements(sel, grid.begin<Edge>(), grid.end<Edge>());
}

inline void SelectNonManifoldEdges(Mesh* obj)
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

inline void SelectSmoothEdgePath(Mesh* obj, number maxDeviation,
								 number normalWeight, bool stopAtSelectedVrts)
{
	SelectSmoothEdgePath(obj->selector(), maxDeviation, normalWeight, stopAtSelectedVrts);
}

inline void SelectShortEdges(Mesh* obj, number maxLength)
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

inline void SelectLongEdges(Mesh* obj, number minLength)
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

inline void SelectCreaseEdges(Mesh* obj, number minAngle)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	SelectCreaseEdges(sel, grid.begin<Edge>(), grid.end<Edge>(),
					  minAngle, obj->position_attachment());
}

inline void SelectLinkedBoundaryEdges(Mesh* obj, bool stopAtSelectedVrts)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();

	if(stopAtSelectedVrts)
		SelectLinkedElements<Edge>(sel, IsOnBoundary(grid), IsNotSelected(sel));
	else
		SelectLinkedElements<Edge>(sel, IsOnBoundary(grid));
}

inline void SelectAssociatedEdges(Mesh* obj)
{
	Selector& sel = obj->selector();
	SelectAssociatedEdges(sel,sel.begin<Face>(), sel.end<Face>());
	SelectAssociatedEdges(sel, sel.begin<Volume>(), sel.end<Volume>());
}

inline void SelectAllEdges(Mesh* obj)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	sel.select(grid.edges_begin(), grid.edges_end());
}

inline void DeselectAllEdges(Mesh* obj)
{
	Selector& sel = obj->selector();
	sel.deselect(sel.edges_begin(), sel.edges_end());
}

inline void SelectMarkedEdges(Mesh* obj)
{
	obj->selector().select(
		obj->crease_handler().begin<Edge>(REM_CREASE),
		obj->crease_handler().end<Edge>(REM_CREASE));
}

inline bool SelectEdgeByIndex(Mesh* obj, int index)
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

inline void EdgeSelectionFill(Mesh* obj)
{
	Selector& sel = obj->selector();
	SelectionFill<Edge>(sel, IsSelected(sel));
}

inline void SelectShortPolychains(Mesh* m, number maxChainLength, bool closedChainsOnly)
{
	SelectShortPolychains(m->selector(), maxChainLength,
						  closedChainsOnly, m->position_accessor());
}

////////////////////////////////////////////////////////////////////////////////
//	FACES
inline void SelectBoundaryFaces(Mesh* obj)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	SelectBoundaryElements(sel, grid.begin<Face>(), grid.end<Face>());
}

inline void SelectInnerFaces(Mesh* obj)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	SelectInnerElements(sel, grid.begin<Face>(), grid.end<Face>());
}

inline void SelectLinkedManifoldFaces(Mesh* obj)
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

template <class TElem>
inline void SelectLinkedElements(Mesh* obj)
{
	SelectLinkedElements<TElem>(obj->selector());
}

inline void SelectLinkedBoundaryFaces(Mesh* obj, bool stopAtSelectedEdges)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();

	if(stopAtSelectedEdges)
		SelectLinkedElements<Face>(sel, IsOnBoundary(grid), IsNotSelected(sel));
	else
		SelectLinkedElements<Face>(sel, IsOnBoundary(grid));
}

inline void SelectDegenerateFaces(Mesh* obj, number maxHeight)
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

inline void SelectLinkedFlatFaces(Mesh* obj, number maxDeviationAngle, bool traverseFlipped,
						   bool traverseDegeneratedFaces, bool stopAtSelectedEdges)
{
	Selector& sel = obj->selector();
	if(traverseDegeneratedFaces)
		ug::SelectLinkedFlatAndDegeneratedFaces(sel, maxDeviationAngle,
												traverseFlipped,
												stopAtSelectedEdges);
	else
		ug::SelectLinkedFlatFaces(sel, maxDeviationAngle, traverseFlipped,
								  stopAtSelectedEdges);
}

inline void SelectIntersectingTriangles(Mesh* obj)
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

inline void SelectAssociatedFaces(Mesh* obj)
{
	Selector& sel = obj->selector();
	SelectAssociatedFaces(sel, sel.begin<Volume>(), sel.end<Volume>());
}

inline void SelectAllFaces(Mesh* obj)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	sel.select(grid.faces_begin(), grid.faces_end());
}

inline void DeselectAllFaces(Mesh* obj)
{
	Selector& sel = obj->selector();
	sel.deselect(sel.faces_begin(), sel.faces_end());
}

inline bool SelectFaceByIndex(Mesh* obj, int index)
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

inline void SelectFacesByNormal(Mesh* obj, const vector3& refNormal, number maxDeviationAngle)
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

inline void FaceSelectionFill(Mesh* obj)
{
	Selector& sel = obj->selector();
	SelectionFill<Face>(sel, IsSelected(sel));
}

inline size_t SelectBentQuadrilaterals(Mesh* obj, number dotThreshold)
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
inline void SelectAllVolumes(Mesh* obj)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	sel.select(grid.volumes_begin(), grid.volumes_end());
}

inline void DeselectAllVolumes(Mesh* obj)
{
	Selector& sel = obj->selector();
	sel.deselect(sel.volumes_begin(), sel.volumes_end());
}

inline int SelectUnorientableVolumes(Mesh* obj)
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

inline int SelectSlivers(Mesh* obj, number thresholdRatio)
{
	std::vector<Tetrahedron*> slivers;
	Grid& g = obj->grid();
	FindSlivers(slivers, g.begin<Tetrahedron>(), g.end<Tetrahedron>(),
			    thresholdRatio, obj->position_accessor());
	obj->selector().select(slivers.begin(), slivers.end());
	return static_cast<int>(slivers.size());
}

inline bool SelectVolumeByIndex(Mesh* obj, int index)
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

inline void VolumeSelectionFill(Mesh* obj)
{
	Selector& sel = obj->selector();
	SelectionFill<Volume>(sel, IsSelected(sel));
}


inline void ClearMarks(Mesh* obj)
{
	obj->crease_handler().clear();
}

inline void MarkCornersOfMarkedEdges(Mesh* obj, number angle)
{
	Selector& sel = obj->selector();
	MarkCorners(obj->grid(), obj->crease_handler(),
				sel.begin<Vertex>(), sel.end<Vertex>(),
				IsInSubset(obj->crease_handler(), ug::REM_CREASE),
				ug::REM_FIXED, angle, obj->position_attachment());
}


inline void MarkSelection(Mesh* obj)
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

inline void UnmarkSelection(Mesh* obj)
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

}}// end of namespace

#endif
