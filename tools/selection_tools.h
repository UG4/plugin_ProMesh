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
#define	TOOLTIP_CLEAR_SELECTION "Clears the selection."
#define TOOLTIP_SELECT_ELEMENTS_BY_INDEX_RANGE "Selects elements by the given index ranges (e.g.: \"0,1,2-5,9,11-23\")."
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
#define	TOOLTIP_SELECT_MARKED_FACES "Selects faces which are marked."
#define	TOOLTIP_SELECT_UNORIENTABLE_VOLUMES "Selects all volumes whose orientation can not be determined"
#define	TOOLTIP_EXTEND_SELECTION "Selects neighbors of selected elements."
#define	TOOLTIP_EXTEND_SELECTION_IN_DIRECTION "Selects neighbors of selected elements whose center can be reached in the given direction from the center of already selected elements."
#define	TOOLTIP_SELECT_VERTEX_BY_INDEX "Selects a vertex given its index."
#define	TOOLTIP_SELECT_EDGE_BY_INDEX "Selects a edge given its index."
#define	TOOLTIP_SELECT_FACE_BY_INDEX "Selects a face given its index."
#define	TOOLTIP_SELECT_FACES_BY_NORMAL "Selects faces given a normal and a maximum deviation angle."
#define	TOOLTIP_SELECT_VOLUME_BY_INDEX "Selects a volume given its index."
#define	TOOLTIP_SELECT_VOLUMES_BY_TYPE "Selects all volumes of a given type."
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
#define TOOLTIP_SELECT_ELEMENTS_BY_SPLIT_PLANE "Selects elements whose center lies in front of the specified plane."
#define TOOLTIP_SELECT_EDGES_BY_DIRECTION "Selects all edges which do not deviate further from the specified direction than the given angle. A minimal required deviation angle can also be specified."
#define TOOLTIP_SELECT_SUBSET_EDGES_BY_DIRECTION "Selects all subset edges which do not deviate further from the specified direction than the given angle. A minimal required deviation angle can also be specified."
#define TOOLTIP_SELECT_ELEMENTS_IN_COORDINATE_RANGE "Selects all elements whose center lies in the specified range."
#define TOOLTIP_DESELECT_ELEMENTS_IN_COORDINATE_RANGE "Deselects all elements whose center lies in the specified range."
#define TOOLTIP_SELECT_VERTEX_IN_BOX "Selects all vertices in the given box"
#define TOOLTIP_SELECT_EDGE_IN_BOX "Selects all edges in the given box"
#define TOOLTIP_SELECT_FACE_IN_BOX "Selects all faces in the given box"
#define TOOLTIP_SELECT_VOLUME_IN_BOX "Selects all volumes in the given box"
#define TOOLTIP_SELECT_VERTEX_IN_CYLINDER "Selects all vertices in the given cylinder"
#define TOOLTIP_SELECT_EDGE_IN_CYLINDER "Selects all edges in the given cylinder"
#define TOOLTIP_SELECT_FACE_IN_CYLINDER "Selects all faces in the given cylinder"
#define TOOLTIP_SELECT_VOLUME_IN_CYLINDER "Selects all volumes in the given cylinder"

#define	TOOLTIP_CLEAR_MARKS "Clears all marks"
#define	TOOLTIP_MARK_CREASE_EDGES "Marks edges whose associated faces have a certain angle as crease-edge."
#define	TOOLTIP_MARK_SELECTION "Marks selected vertices and edges."
#define	TOOLTIP_UNMARK_SELECTION "Unmarks selected elements."
#define TOOLTIP_MARK_CORNERS_OF_MARKED_EDGES "Marks selected vertices as fixed, if they lie at a sharp corner of a marked path or if they are at endpoints or at junctions of marked edges."


namespace ug{
namespace promesh{

/// \addtogroup promesh
/// \{

void ClearSelection(Mesh* obj);

void SelectAll(Mesh* obj);

void SelectElementsByIndexRange (Mesh* obj,
                            const char* vrtRanges,
                            const char* edgeRanges,
                            const char* faceRanges,
                            const char* volRanges,
                            bool clearSelection);

void ExtendSelection(Mesh* obj, int neighborhoodSize);

void ExtendSelectionInDirection(
		Mesh* obj,
		int neighborhoodSize,
		const vector3& dir,
		number minAngle,
		number maxAngle);

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

void SelectSubset(
			Mesh* obj,
			int si,
			bool selVrts,
			bool selEdges,
			bool selFaces,
			bool selVols);

void SelectSubsetBoundary(
			Mesh* obj,
			int si,
			bool edgeBnds,
			bool faceBnds,
			bool volBnds);

void SelectUnassignedElements(
			Mesh* obj,
			bool selVrts,
			bool selEdges,
			bool selFaces,
			bool selVols);

void InvertSelection(
			Mesh* obj,
			bool invVrts,
			bool invEdges,
			bool invFaces,
			bool invVols);

void SelectSelectionBoundary(Mesh* obj);

void CloseSelection(Mesh* obj);


template <class TElem>
void SelectInterfaceElements(Mesh* obj, bool regardSelectedNbrsOnly)
{
	Grid& g = obj->grid();
	SelectInterfaceElements(obj->selector(), obj->subset_handler(),
							g.begin<TElem>(), g.end<TElem>(),
							regardSelectedNbrsOnly);
}


template <class TElem>
void SelectAnisotropicElements(Mesh* obj, number minEdgeRatio)
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
void SelectBoundaryVertices(Mesh* obj);

void SelectInnerVertices(Mesh* obj);

void SelectAssociatedVertices(Mesh* obj);

void SelectAllVertices(Mesh* obj);

void DeselectAllVertices(Mesh* obj);

void SelectMarkedVertices(Mesh* obj);

bool SelectVertexByIndex(Mesh* obj, int index);

size_t SelectUnconnectedVertices(
				Mesh* obj,
				bool edgeCons,
				bool faceCons,
				bool volCons);

size_t SelectSelectionKinkVertices(
				Mesh* obj,
				number thresholdAngle,
				bool selectDarts);

size_t SelectSubsetKinkVertices(
				Mesh* obj,
				int subsetIndex,
				number thresholdAngle,
				bool selectDarts);

void SelectBoundaryEdges(Mesh* obj);

void SelectInnerEdges(Mesh* obj);

void SelectNonManifoldEdges(Mesh* obj);

void SelectSmoothEdgePath(
			Mesh* obj,
			number maxDeviation,
			number normalWeight,
			bool stopAtSelectedVrts);

void SelectShortEdges(Mesh* obj, number maxLength);

void SelectLongEdges(Mesh* obj, number minLength);

void SelectCreaseEdges(Mesh* obj, number minAngle);

void SelectLinkedBoundaryEdges(Mesh* obj, bool stopAtSelectedVrts);

void SelectAssociatedEdges(Mesh* obj);

void SelectAllEdges(Mesh* obj);

void DeselectAllEdges(Mesh* obj);

void SelectMarkedEdges(Mesh* obj);

bool SelectEdgeByIndex(Mesh* obj, int index);

void EdgeSelectionFill(Mesh* obj);

void SelectShortPolychains(Mesh* m, number maxChainLength, bool closedChainsOnly);

void SelectEdgesByDirection(
				Mesh* m,
				const vector3& dir,
				number minDeviationAngle,
				number maxDeviationAngle,
				bool selectFlipped);

void SelectSubsetEdgesByDirection(
				Mesh* m,
				int subsetIndex,
				const vector3& dir,
				number minDeviationAngle,
				number maxDeviationAngle,
				bool selectFlipped);

void SelectBoundaryFaces(Mesh* obj);

void SelectInnerFaces(Mesh* obj);

void SelectMarkedFaces(Mesh* obj);

void SelectLinkedManifoldFaces(Mesh* obj);

template <class TElem>
void SelectLinkedElements(Mesh* obj)
{
	SelectLinkedElements<TElem>(obj->selector());
}

void SelectLinkedBoundaryFaces(Mesh* obj, bool stopAtSelectedEdges);

void SelectDegenerateFaces(Mesh* obj, number maxHeight);

void SelectLinkedFlatFaces(
			Mesh* obj,
			number maxDeviationAngle,
			bool ignoreOrientation,
			bool traverseDegeneratedFaces,
			bool stopAtSelectedEdges);

void SelectIntersectingTriangles(Mesh* obj);

void SelectAssociatedFaces(Mesh* obj);

void SelectAllFaces(Mesh* obj);

void DeselectAllFaces(Mesh* obj);

bool SelectFaceByIndex(Mesh* obj, int index);

void SelectFacesByNormal(
			Mesh* obj,
			const vector3& refNormal,
			number maxDeviationAngle);

void SelectFacesByNormal(
			Mesh* obj,
			const vector3& refNormal,
			number minDeviationAngle,
			number maxDeviationAngle,
			bool noInnerFaces);

void FaceSelectionFill(Mesh* obj);

size_t SelectBentQuadrilaterals(Mesh* obj, number dotThreshold);

void SelectAllVolumes(Mesh* obj);

void DeselectAllVolumes(Mesh* obj);

int SelectUnorientableVolumes(Mesh* obj);

int SelectSlivers(Mesh* obj, number thresholdRatio);

bool SelectVolumeByIndex(Mesh* obj, int index);

void SelectVolumesByType(
			Mesh* obj,
			bool selHexahedra,
			bool selOctahedra,
			bool selPrisms,
			bool selPyramids,
			bool selTetrahedra);

void VolumeSelectionFill(Mesh* obj);

void ClearMarks(Mesh* obj);

void MarkCornersOfMarkedEdges(Mesh* obj, number angle);

void MarkSelection(Mesh* obj);

void UnmarkSelection(Mesh* obj);

void MarkCreaseEdges(Mesh* obj, number minAngle, bool clearMarks);

void SelectElementsBySplitPlane(
			Mesh* obj,
			bool selectVrts,
			bool selectEdges,
			bool selectFaces,
			bool selectVols,
			const vector3& pivot,
			const vector3& normal);


template <int icoord, int selFlag, class TElem>
void SelectElementsInCoordinateRange(
			Mesh* obj,
			number min,
			number max)
{
	Grid& g = obj->grid();
	Selector& sel = obj->selector();
	Mesh::position_accessor_t aaPos = obj->position_accessor();

	for(typename Grid::traits<TElem>::iterator iter = g.begin<TElem>();
		iter != g.end<TElem>(); ++iter)
	{
		vector3 c = CalculateCenter(*iter, aaPos);
		if(c[icoord] >= min && c[icoord] <= max)
			sel.select(*iter, selFlag);
	}
}

template <int icoord, int selFlag>
void SelectElementsInCoordinateRange(
			Mesh* mesh,
			number min,
			number max,
			bool vrts,
			bool edges,
			bool faces,
			bool vols)
{
	if(vrts)
		SelectElementsInCoordinateRange<icoord, selFlag, Vertex>(mesh, min, max);
	if(edges)
		SelectElementsInCoordinateRange<icoord, selFlag, Edge>(mesh, min, max);
	if(faces)
		SelectElementsInCoordinateRange<icoord, selFlag, Face>(mesh, min, max);
	if(vols)
		SelectElementsInCoordinateRange<icoord, selFlag, Volume>(mesh, min, max);
}


///	Selects elements whose center lie in a box
template <class TElem>
void SelectElementsInBox(Mesh* obj, const vector3& min, const vector3& max)
{
	Mesh::position_accessor_t& aaPos = obj->position_accessor();

	Grid& grid = obj->grid();
	Selector& sel = obj->selector();

	for(typename Grid::traits<TElem>::iterator iter = grid.begin<TElem>();
		iter != grid.end<TElem>(); ++iter)
	{
		if(BoxBoundProbe(CalculateCenter(*iter, aaPos), min, max))
			sel.select(*iter);
	}
}


///	Selects elements whose center lie in a cylinder
template <class TElem>
void SelectElementsInCylinder(
			Mesh* obj,
			const vector3& cylBase,
			const vector3& cylTop,
			number radius)
{
	Mesh::position_accessor_t& aaPos = obj->position_accessor();

	Grid& grid = obj->grid();
	Selector& sel = obj->selector();

	vector3 from = cylBase;
	vector3 dir;
	VecSubtract(dir, cylTop, cylBase);

	for(typename Grid::traits<TElem>::iterator iter = grid.begin<TElem>();
		iter != grid.end<TElem>(); ++iter)
	{
		vector3 c = CalculateCenter(*iter, aaPos);
		vector3 p;
		number s = ProjectPointToRay(p, c, from, dir);
		if((s > -SMALL) && (s < (1. + SMALL))){
			if(VecDistanceSq(p, c) <= sq(radius))
				sel.select(*iter);
		}
	}
}

/// \}

}}// end of namespace

#endif
