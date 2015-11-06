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

#ifndef __H__UG_PROMESH__file_io_tools__
#define __H__UG_PROMESH__file_io_tools__

#include "../mesh.h"
#include "lib_grid/file_io/file_io.h"
#include "lib_grid/file_io/file_io_ug.h"
#include "lib_grid/file_io/file_io_ugx.h"

//	file io
#define TOOLTIP_LOAD_MESH "Loads a Mesh from File. The format is automatically recognized by the filname's suffix. Supported formats are: (ugx, vtu (ascii), obj, stl (ascii and binary), lgm, ng, ele, msh, asc, net, art)."
#define TOOLTIP_SAVE_MESH "Saves a Mesh to File. The format is automatically recognized by the filname's suffix. Supported formats are: (ugx, vtu (ascii), obj, stl (ascii), ncdf, smesh, ele, tikz/tex, asc, net, art)."
#define TOOLTIP_EXPORT_TO_UG3 "Writes a mesh to the UG3 legacy format (.lgm and .ng)"

namespace ug{
namespace promesh{

inline bool LoadMesh(Mesh* obj, const char* filename)
{
	const char* pSuffix = strrchr(filename, '.');
	if(!pSuffix)
		return false;

	if(strcmp(pSuffix, ".ugx") == 0)
	{
	//	load from ugx
		GridReaderUGX ugxReader;
		if(!ugxReader.parse_file(filename)){
			UG_LOG("ERROR in LoadGridFromUGX: File not found: " << filename << std::endl);
			return false;
		}
		else{
			if(ugxReader.num_grids() < 1){
				UG_LOG("ERROR in LoadGridFromUGX: File contains no grid.\n");
				return false;
			}
			else{
				ugxReader.grid(obj->grid(), 0, obj->position_attachment());

				if(ugxReader.num_subset_handlers(0) > 0)
					ugxReader.subset_handler(obj->subset_handler(), 0, 0);

				if(ugxReader.num_subset_handlers(0) > 1)
					ugxReader.subset_handler(obj->crease_handler(), 1, 0);

				if(ugxReader.num_selectors(0) > 0)
					ugxReader.selector(obj->selector(), 0, 0);

				return true;
			}
		}
	}
	else{
		return LoadGridFromFile(obj->grid(), obj->subset_handler(), filename,
								obj->position_attachment());
	}
	return false;
}

inline bool SaveMesh(Mesh* obj, const char* filename)
{
	const char* pSuffix = strrchr(filename, '.');
	if(!pSuffix)
		return false;

	if(strcmp(pSuffix, ".ugx") == 0){
		GridWriterUGX ugxWriter;
		ugxWriter.add_grid(obj->grid(), "defGrid", obj->position_attachment());
		ugxWriter.add_subset_handler(obj->subset_handler(), "defSH", 0);
		ugxWriter.add_subset_handler(obj->crease_handler(), "markSH", 0);
		ugxWriter.add_selector(obj->selector(), "defSel", 0);

		return ugxWriter.write_to_file(filename);
	}
	else
		return SaveGridToFile(obj->grid(), obj->subset_handler(), filename,
							  obj->position_attachment());
}

inline bool ExportToUG3(Mesh* obj, const char* filenamePrefix, const char* lgmName,
				 const char* problemName)
{
	bool saveOk = false;

	Grid& g = obj->grid();
	SubsetHandler& sh = obj->subset_handler();

	if(g.num_volumes() > 0){
	//	Create the subset-handlers
		SubsetHandler shFaces(g, SHE_FACE);
		SubsetHandler shVolumes(g, SHE_VOLUME);

		for(int i = 0; i < sh.num_subsets(); ++i){
			shFaces.assign_subset(sh.begin<Face>(i), sh.end<Face>(i), i);
			shVolumes.assign_subset(sh.begin<Volume>(i), sh.end<Volume>(i), i);
		}

		saveOk = ExportGridToUG(g, shFaces, shVolumes, filenamePrefix,
								lgmName, problemName, 0);
	}
	else{
		saveOk = ExportGridToUG_2D(g, filenamePrefix, lgmName, problemName, 0, &sh);
	}
	return saveOk;
}

}}// end of namespace

#endif
