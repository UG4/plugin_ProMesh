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

#include "registration_routines.h"
#include "registry/registry.h"
#include "bridge/util.h"
#include "tools/subset_tools.h"

using namespace std;
using namespace ug::bridge;

namespace ug{
namespace promesh{
void RegisterSubsetTools(ProMeshRegistry& reg, std::string baseGrp)
{
	baseGrp.append("/Subsets");
	string grp = baseGrp;
	reg.add_function("AssignSubset",
			static_cast<void (*)(Mesh*, int)>(&AssignSubset), grp, "",
			"mesh # subset index", TOOLTIP_ASSIGN_SUBSET, "", RT_NO_PROMESH)
		.add_function("AssignSubset",
			static_cast<void (*)(Mesh*, int, bool, bool, bool, bool)>(&AssignSubset), grp, "",
			"mesh #"
			"subset index #"
			"assign vertices || value=true #"
			"assign edges || value=true #"
			"assign faces || value=true #"
			"assign volumes || value=true ", TOOLTIP_ASSIGN_SUBSET)
		.add_function("AssignNewSubset", &AssignNewSubset, grp, "",
			"mesh #"
			"name || value=\"subset\" #"
			"assign vertices || value=true #"
			"assign edges || value=true #"
			"assign faces || value=true #"
			"assign volumes || value=true ", TOOLTIP_ASSIGN_NEW_SUBSET)
		.add_function("UnassignSubsets", &UnassignSubsets, grp, "",
			"mesh #"
			"unassign vertices || value=true #"
			"unassign edges || value=true #"
			"unassign faces || value=true #"
			"unassign volumes || value=true ", TOOLTIP_UNASSIGN_SUBSETS)
		.add_function("SetSubsetName", &SetSubsetName, grp, "",
			"mesh # subset index # name ", TOOLTIP_SET_SUBSET_NAME, "", RT_NO_PROMESH)
		.add_function("MoveSubset", &MoveSubset, grp, "",
			"mesh # old subset index # new subset index", TOOLTIP_MOVE_SUBSET)
		.add_function("SwapSubsets", &SwapSubsets, grp, "",
			"mesh # subset index 1 # subset index 2", TOOLTIP_SWAP_SUBSETS)
		.add_function("JoinSubsets", &JoinSubsets, grp, "",
			"mesh # target subset index # subset index 1 # subset index 2", TOOLTIP_JOIN_SUBSETS)
		.add_function("EraseSubset", &EraseSubset, grp, "",
			"mesh # subset index # erase geometry", TOOLTIP_ERASE_SUBSET)
		.add_function("AdjustSubsetsForUG3", &AdjustSubsetsForUG3, grp, "",
			"mesh # keep interface subsets", TOOLTIP_ADJUST_SUBSETS_FOR_UG3, "", RT_NO_PROMESH)
		.add_function("AdjustSubsetsForUG4", &AdjustSubsetsForUG4, grp, "",
			"mesh # preserve existing subsets", TOOLTIP_ADJUST_SUBSETS_FOR_UG4, "", RT_NO_PROMESH)
		.add_function("AssignSubsetsByQuality", &AssignSubsetsByQuality, grp, "",
			"mesh # num sections || value=10", TOOLTIP_ASSIGN_SUBSETS_BY_QUALITY)
		.add_function("CopySubsetIndicesToSides", &CopySubsetIndicesToSides, grp, "",
			"mesh #"
			"selection only #"
			"to unassigned elements only || value=true", TOOLTIP_COPY_SUBSET_INDICES_TO_SIDES);

	grp = baseGrp + "/Separate";
	reg.add_function("SeparateFacesByEdgeSubsets", &SeparateFacesByEdgeSubsets, grp, "",
			"mesh", TOOLTIP_SEPARATE_FACES_BY_EDGE_SUBSETS)
		.add_function("SeparateFacesBySelectedEdges", &SeparateFacesBySelectedEdges, grp, "",
			"mesh", TOOLTIP_SEPARATE_FACES_BY_SELECTED_EDGES)
		.add_function("SeparateVolumesByFaceSubsets", &SeparateVolumesByFaceSubsets, grp, "",
			"mesh", TOOLTIP_SEPARATE_VOLUMES_BY_FACE_SUBSETS)
		.add_function("SeparateVolumesBySelectedFaces", &SeparateVolumesBySelectedFaces, grp, "",
			"mesh", TOOLTIP_SEPARATE_VOLUMES_BY_SELECTED_FACES)
		.add_function("SeparateIrregularManifoldSubsets", &SeparateIrregularManifoldSubsets, grp, "",
			"mesh", TOOLTIP_SEPARATE_IRREGULAR_MANIFOLD_SUBSETS)
		.add_function("SeparateFaceSubsetsByNormal", &SeparateFaceSubsetsByNormal, grp, "",
			"mesh", TOOLTIP_SEPARATE_FACE_SUBSETS_BY_NORMAL)
		.add_function("SeparateFaceSubsetByNormal", &SeparateFaceSubsetByNormal, grp, "",
			"mesh # subset index", TOOLTIP_SEPARATE_FACE_SUBSET_BY_NORMAL)
		.add_function("SeparateDegeneratedBoundaryFaceSubsets", &SeparateDegeneratedBoundaryFaceSubsets, grp, "",
			"mesh # angle", TOOLTIP_SEPARATE_DEGENERATED_BOUNDARY_FACE_SUBSETS);

	grp = baseGrp;
	reg.add_function("EraseEmptySubsets", &EraseEmptySubsets, grp, "",
			"mesh", TOOLTIP_ERASE_EMPTY_SUBSETS)
		.add_function("AssignSubsetColors", &AssignSubsetColors, grp, "", "mesh", TOOLTIP_ASSIGN_SUBSET_COLORS)
		.add_function("AssignSubsetsByElementType", &AssignSubsetsByElementType, grp, "",
			"mesh", TOOLTIP_ASSIGN_SUBSETS_BY_ELEMENT_TYPE);
}

}}//	end of namespace
