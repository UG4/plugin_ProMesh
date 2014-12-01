// created by Sebastian Reiter
// s.b.reiter@gmail.com
// Jun 25, 2013 (d,m,y)

#ifndef __H__UG_PROMESH__file_io_tools__
#define __H__UG_PROMESH__file_io_tools__

#include "../mesh.h"
#include "lib_grid/file_io/file_io.h"
#include "lib_grid/file_io/file_io_ug.h"
#include "lib_grid/file_io/file_io_ugx.h"

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
