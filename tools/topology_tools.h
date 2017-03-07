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

#ifndef __H__UG__topology_tools__
#define __H__UG__topology_tools__

#include <vector>
#include "../mesh.h"
#include "lib_grid/algorithms/remeshing/resolve_intersections.h"
#include "lib_grid/algorithms/geom_obj_util/vertex_util.h"
#include "lib_grid/algorithms/geom_obj_util/edge_util.h"
#include "lib_grid/algorithms/remove_duplicates_util.h"
#include "lib_grid/file_io/file_io.h"

#define TOOLTIP_RESOLVE_SELF_INTERSECTIONS "Resolves self intersections of faces and edges."
#define	TOOLTIP_RESOLVE_EDGE_INTERSECTIONS "Makes sure that all edge intersections are represented by a vertex."
#define	TOOLTIP_RESOLVE_TRIANGLE_INTERSECTIONS "Makes sure that all triangle intersections are represented by an edge and vertices."
#define	TOOLTIP_PROJECT_VERTICES_TO_CLOSE_EDGES "Projects selected vertices to selected close edges."
#define	TOOLTIP_PROJECT_VERTICES_TO_CLOSE_FACES "Projects selected vertices to selected close faces."
#define	TOOLTIP_INTERSECT_CLOSE_EDGES "Performs intersections between selected close edges."
#define	TOOLTIP_ERASE_SELECTED_ELEMENTS "Erases selected elements and associated unreferenced geometry."
#define	TOOLTIP_REMOVE_DOUBLES "Removes selected vertices that are close to each other"
#define	TOOLTIP_REMOVE_DOUBLE_EDGES "Removes selected duplicates of selected edges."
#define	TOOLTIP_REMOVE_DOUBLE_FACES "Removes selected duplicates of selected faces."
#define	TOOLTIP_MERGE_AT_FIRST "Merges all selected objects into a single vertex at the first selected vertex."
#define	TOOLTIP_MERGE_AT_CENTER "Merges all selected objects into a single vertex at the center of the selection."
#define	TOOLTIP_MERGE_AT_LAST "Merges all selected objects into a single vertex at the last selected vertex."
#define	TOOLTIP_COLLAPSE_EDGE "Collapses the edge and removes adjacent triangles."
#define	TOOLTIP_SPLIT_EDGE "Splits the edge and inserts new triangles."
#define	TOOLTIP_SWAP_EDGE "Swaps selected edges that are adjacent to exactly two triangles."
#define	TOOLTIP_PLANE_CUT "Cuts selected edges along the given plane."
#define	TOOLTIP_ADJUST_EDGE_ORIENTATION "Adjusts the orientation of boundary edges to associated faces."
#define	TOOLTIP_FIX_FACE_ORIENTATION "Tries to change orientation of selected faces so that all neighbouring faces point into the same direction. Only works correctly for manifold selections."
#define	TOOLTIP_FIX_FACE_SUBSET_ORIENTATIONS "Iterates over all subset and tries to fix face orientation for each. Only works correctly for manifold subsets."
#define	TOOLTIP_FIX_VOLUME_ORIENTATION "Changes orientation of selected volumes, so that they are oriented correctly."
#define	TOOLTIP_INVERT_FACE_ORIENTATION "Inverts the orientation of all selected faces."


namespace ug{
namespace promesh{

/// \addtogroup promesh
/// \{

void EraseSelectedElements(
			Mesh* obj,
			bool eraseUnusedVrts,
			bool eraseUnusedEdges,
			bool eraseUnusedFaces);

size_t RemoveDoubles(Mesh* obj, number threshold);

size_t RemoveDoubleEdges(Mesh* obj);

size_t RemoveDoubleFaces(Mesh* obj);

void MergeAtFirst(Mesh* obj);

void MergeAtCenter(Mesh* obj);

void MergeAtLast(Mesh* obj);

void CollapseEdge(Mesh* obj);

void SplitEdge(Mesh* obj);

void SwapEdge(Mesh* obj);

void PlaneCut(Mesh* obj, const vector3& p, const vector3& n);

void AdjustEdgeOrientation(Mesh* obj);

void FixFaceOrientation(Mesh* obj);

void FixFaceSubsetOrientations(Mesh* obj);

int FixVolumeOrientation(Mesh* obj);

void InvertFaceOrientation(Mesh* obj);

void ResolveEdgeIntersection(Mesh* obj, number snapThreshold);

void ResolveTriangleIntersections(Mesh* obj, number snapThreshold);

void ProjectVerticesToCloseEdges(Mesh* obj, number snapThreshold);

void ProjectVerticesToCloseFaces(Mesh* obj, number snapThreshold);

void IntersectCloseEdges(Mesh* obj, number snapThreshold);

void ResolveSelfIntersections(Mesh* obj, number snapThreshold);

/// \}

}}// end of namespace

#endif
