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

#include "mesh.h"
#include "promesh_plugin.h"
#include "promesh_registry.h"
#include "registration_routines.h"
#include "tooltips.h"

#include "tools/file_io_tools.h"
#include "tools/new_tools.h"
#include "tools/subset_tools.h"
#include "tools/measure_tools.h"
#include "bridge/util.h"

using namespace std;
using namespace ug::bridge;

namespace ug{
namespace promesh{

static ProMeshRegistry* g_promeshRegistry = NULL;

ProMeshRegistry&
GetProMeshRegistry()
{
	return *g_promeshRegistry;
}

static void RegisterMisc(ProMeshRegistry& reg, string baseGrp)
{
	string grp;

//	file io
	grp = baseGrp + "/File IO";
	reg.add_function("LoadMesh", &LoadMesh, grp, "",
			"mesh # filename", TOOLTIP_LOAD_MESH, "", RT_NO_PROMESH)
		.add_function("SaveMesh", &SaveMesh, grp, "",
			"mesh # filename", TOOLTIP_SAVE_MESH, "", RT_NO_PROMESH)
		.add_function("ExportToUG3", &ExportToUG3, grp, "",
			"mesh # filenamePrefix # lgmName # problemName", TOOLTIP_EXPORT_TO_UG3, "", RT_NO_PROMESH);

//	subset tools
	grp = baseGrp + "/Subsets";
	reg.add_function("AssignSubset",
			static_cast<void (*)(Mesh*, int)>(&AssignSubset), grp, "",
			"mesh # subset index", TOOLTIP_ASSIGN_SUBSET)
		.add_function("AssignSubset",
			static_cast<void (*)(Mesh*, int, bool, bool, bool, bool)>(&AssignSubset), grp, "",
			"mesh # subset index # assign vertices # assign edges # assign faces # assign volumes", TOOLTIP_ASSIGN_SUBSET)
		.add_function("SetSubsetName", &SetSubsetName, grp, "",
			"mesh # subset index # name ", TOOLTIP_SET_SUBSET_NAME)
		.add_function("AssignSubsetColors", &AssignSubsetColors, grp, "", "mesh", TOOLTIP_ASSIGN_SUBSET_COLORS)
		.add_function("MoveSubset", &MoveSubset, grp, "",
			"mesh # old subset index # new subset index", TOOLTIP_MOVE_SUBSET)
		.add_function("SwapSubsets", &SwapSubsets, grp, "",
			"mesh # subset index 1 # subset index 2", TOOLTIP_SWAP_SUBSETS)
		.add_function("JoinSubsets", &JoinSubsets, grp, "",
			"mesh # target subset index # subset index 1 # subset index 2", TOOLTIP_JOIN_SUBSETS)
		.add_function("EraseSubset", &EraseSubset, grp, "",
			"mesh # subset index # erase geometry", TOOLTIP_ERASE_SUBSET)
		.add_function("EraseEmptySubsets", &EraseEmptySubsets, grp, "",
			"mesh", TOOLTIP_ERASE_EMPTY_SUBSETS)
		.add_function("AdjustSubsetsForUG3", &AdjustSubsetsForUG3, grp, "",
			"mesh # keep interface subsets", TOOLTIP_ADJUST_SUBSETS_FOR_UG3)
		.add_function("AdjustSubsetsForUG4", &AdjustSubsetsForUG4, grp, "",
			"mesh # preserve existing subsets", TOOLTIP_ADJUST_SUBSETS_FOR_UG4)
		.add_function("AssignSubsetsByQuality", &AssignSubsetsByQuality, grp, "",
			"mesh # num sections", TOOLTIP_ASSIGN_SUBSETS_BY_QUALITY)
		.add_function("CopySubsetIndicesToSides", &CopySubsetIndicesToSides, grp, "",
			"mesh # selection only # to unassigned elements only", TOOLTIP_COPY_SUBSET_INDICES_TO_SIDES)
		.add_function("AssignSubsetsByElementType", &AssignSubsetsByElementType, grp, "",
			"mesh", TOOLTIP_ASSIGN_SUBSETS_BY_ELEMENT_TYPE);

	grp = baseGrp + "/Subsets/Separate";
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

//	info tools
	grp = baseGrp + string("/Info/Measure length, area, volume");
	reg.add_function("MeasureGridLength", &MeasureGridLength, grp, "length",
			"mesh", TOOLTIP_MEASURE_GRID_LENGTH)
	   .add_function("MeasureGridArea", &MeasureGridArea, grp, "area",
	   		"mesh", TOOLTIP_MEASURE_GRID_AREA)
	   .add_function("MeasureGridVolume", &MeasureGridVolume, grp, "volume",
	   		"mesh", TOOLTIP_MEASURE_GRID_VOLUME)
	   .add_function("MeasureSubsetLength", &MeasureSubsetLength, grp, "length",
	   		"mesh#subset", TOOLTIP_MEASURE_SUBSET_LENGTH)
	   .add_function("MeasureSubsetArea", &MeasureSubsetArea, grp, "area",
	   		"mesh#subset", TOOLTIP_MEASURE_SUBSET_AREA)
	   .add_function("MeasureSubsetVolume", &MeasureSubsetVolume, grp, "volume",
	   		"mesh#subset", TOOLTIP_MEASURE_SUBSET_VOLUME)
	   .add_function("MeasureSelectionLength", &MeasureSelectionLength, grp, "length",
	   		"mesh", TOOLTIP_MEASURE_SELECTION_LENGTH)
	   .add_function("MeasureSelectionArea", &MeasureSelectionArea, grp, "area",
	   		"mesh", TOOLTIP_MEASURE_SELECTION_AREA)
	   .add_function("MeasureSelectionVolume", &MeasureSelectionVolume, grp, "volume",
	   		"mesh", TOOLTIP_MEASURE_SELECTION_VOLUME);

//	new tools
	grp = baseGrp + "/Util";
	reg.add_class_<Box>("Box", grp)
		.add_method("set_min", &Box::set_min)
		.add_method("set_max", &Box::set_max)
		.add_method("min", &Box::get_min)
		.add_method("max", &Box::get_max)
		.add_method("local_to_global", &Box::local_to_global)
		.add_method("global_to_local", &Box::global_to_local);

	reg.add_function("GetBoundingBox", &GetBoundingBox, grp, "", "",
					 TOOLTIP_GET_BOUNDING_BOX, "", RT_NO_PROMESH);


	// reg.add_projector_<RefinementProjector>("default",
	// 		"Places new vertices at the center of their parent element.")
	// 	.add_projector_<SphereProjectorNew>("sphere",
	// 		"Places new vertices on a sphere. The distance of the new vertex to the center "
	// 		"is thereby calculated as the average distance of the parent's corners "
	// 		"to the center. The radius property is thereby ignored. The radius property "
	// 		"is only used during reprojection of associated vertices.")
	// 	.add_projector_<CylinderProjectorNew>("cylinder",
	// 		"Places new vertices on a cylinder. The distance of the new vertex to the axis "
	// 		"is thereby calculated as the average distance of the parent's corners "
	// 		"to the axis. The radius property is thereby ignored. The radius property "
	// 		"is only used during reprojection of associated vertices.");
}

} // end namespace promesh


/**
 * This function is called when the plugin is loaded.
 */
extern "C" void
InitUGPlugin_ProMesh(Registry* reg, string grp)
{
	using namespace ug::promesh;
	grp.append("promesh");

	try{
		if(g_promeshRegistry == NULL)
			g_promeshRegistry = new ProMeshRegistry(reg);

		ProMeshRegistry& pmreg = GetProMeshRegistry();
		RegisterMesh(pmreg, grp);
		RegisterCoordinateTransformTools(pmreg, grp);
		RegisterSelectionTools(pmreg, grp);
		RegisterMeshingTools(pmreg, grp);
		RegisterMisc(pmreg, grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

}// end of namespace

