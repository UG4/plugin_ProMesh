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

#include "mesh.h"
#include "promesh_plugin.h"
#include "promesh_registry.h"
#include "registration_routines.h"
#include "tooltips.h"

#include "tools/file_io_tools.h"
#include "tools/new_tools.h"
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
		RegisterSubsetTools(pmreg, grp);
		RegisterMeshingTools(pmreg, grp);
		RegisterMisc(pmreg, grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

}// end of namespace

