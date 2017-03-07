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

#ifndef __H__UG_PROMESH__subset_tools__
#define __H__UG_PROMESH__subset_tools__

#include <vector>
#include <queue>
#include "../mesh.h"
#include "lib_grid/algorithms/subset_util.h"
#include "lib_grid/algorithms/grid_statistics.h"

#define TOOLTIP_SET_SUBSET_NAME ""
#define	TOOLTIP_ASSIGN_SUBSET "Assigns the selected elements to a subset."
#define TOOLTIP_ASSIGN_NEW_SUBSET "Assigns selected elements to a new subset."
#define TOOLTIP_UNASSIGN_SUBSET "Unassigns selected elements from any subset."
#define	TOOLTIP_ASSIGN_SUBSET_COLORS "assigns subset colors by a procedural scheme."
#define	TOOLTIP_SEPARATE_FACES_BY_EDGE_SUBSETS "Assigns faces that are surrounded by a set of edge-subsets to a common subset."
#define	TOOLTIP_SEPARATE_FACES_BY_SELECTED_EDGES "Assigns faces that are surrounded by a set of selected edges to a common subset."
#define	TOOLTIP_SEPARATE_VOLUMES_BY_FACE_SUBSETS "Assigns volumes that are surrounded by a set of face-subsets to a common subset."
#define	TOOLTIP_SEPARATE_VOLUMES_BY_SELECTED_FACES "Assigns volumes that are surrounded by a set of selected faces to a common subset."
#define	TOOLTIP_SEPARATE_IRREGULAR_MANIFOLD_SUBSETS "After this algorithm all face-subsets are regular manifolds."
#define	TOOLTIP_MOVE_SUBSET "Moves a subset to another index."
#define	TOOLTIP_SWAP_SUBSETS "Swaps two subsets"
#define	TOOLTIP_JOIN_SUBSETS "Joins two subsets"
#define	TOOLTIP_ERASE_SUBSET "Erases a subset, but not its associated geometry."
#define	TOOLTIP_ERASE_EMPTY_SUBSETS "Erases Subsets, which do not contain any elements at all."
#define	TOOLTIP_ADJUST_SUBSETS_FOR_UG3 "Assigns face and edge indices so that the geometry can be used with ug3."
#define	TOOLTIP_ADJUST_SUBSETS_FOR_UG4 "Adjusts subsets for simulation with ug4."
#define	TOOLTIP_SEPARATE_FACE_SUBSETS_BY_NORMAL "Collects faces of each subset that have a similar normal and assigns them to new subsets."
#define	TOOLTIP_SEPARATE_FACE_SUBSET_BY_NORMAL "Collects faces of a given subset that have a similar normal and assigns them to new subsets."
#define	TOOLTIP_ASSIGN_SUBSETS_BY_QUALITY "Assigns the selected to a subset depending on their quality."
#define	TOOLTIP_SEPARATE_DEGENERATED_BOUNDARY_FACE_SUBSETS "Separates degenerated boundary face subsets at sharp creases."
#define	TOOLTIP_COPY_SUBSET_INDICES_TO_SIDES "Copies subset indices of selected elements to sides of those elements."
#define	TOOLTIP_ASSIGN_SUBSETS_BY_ELEMENT_TYPE "Assigns elemets to subsets based on their concrete type."

namespace ug{
namespace promesh{

/// \addtogroup promesh
/// \{

void AssignSubset(Mesh* obj, int newIndex);

void AssignSubset(
			Mesh* obj,
			int newIndex,
			bool vertices,
			bool edges,
			bool faces,
			bool volumes);

void AssignNewSubset(
			Mesh* obj,
			const char* name,
			bool vertices,
			bool edges,
			bool faces,
			bool volumes);

void UnassignSubset(
			Mesh* obj,
			bool vertices,
			bool edges,
			bool faces,
			bool volumes);

void SetSubsetName(Mesh* obj, int si, const char* name);

void AssignSubsetColors(Mesh* obj);

void SeparateFacesByEdgeSubsets(Mesh* obj);

void SeparateFacesBySelectedEdges(Mesh* obj);

void SeparateVolumesByFaceSubsets(Mesh* obj);

void SeparateVolumesBySelectedFaces(Mesh* obj);

void SeparateIrregularManifoldSubsets(Mesh* obj);

void MoveSubset(Mesh* obj, int oldIndex, int newIndex);

void SwapSubsets(Mesh* obj, int oldIndex, int newIndex);

void JoinSubsets(Mesh* obj, int target, int si1, int si2, bool eraseUnused);

void EraseSubset(Mesh* obj, int si, bool eraseGeometry);

void EraseEmptySubsets(Mesh* obj);

void AdjustSubsetsForUG3(Mesh* obj, bool keepIntfSubs);

void AdjustSubsetsForUG4(Mesh* obj, bool preserveExistingSubsets);

void SeparateFaceSubsetsByNormal(Mesh* obj);

void SeparateFaceSubsetByNormal(Mesh* obj, int si);

void AssignSubsetsByQuality(Mesh* obj, int numSections);

void SeparateDegeneratedBoundaryFaceSubsets(Mesh* obj, number angle);

void AssignSubsetsByElementType(Mesh* obj);

void CopySubsetIndicesToSides(
			Mesh* obj,
			bool selectionOnly,
			bool toUnassignedOnly);

/// \}

}}// end of namespace

#endif
