// created by Sebastian Reiter
// s.b.reiter@gmail.com
// Jun 25, 2013 (d,m,y)

#ifndef __H__UG_PROMESH__file_io_tools__
#define __H__UG_PROMESH__file_io_tools__

#include "../mesh_object.h"
#include "lib_grid/file_io/file_io.h"

namespace ug{
namespace promesh{

bool LoadMesh(MeshObject* obj, const char* filename)
{
	return LoadGridFromFile(obj->get_grid(), obj->get_subset_handler(), filename,
							obj->position_attachment());
}

bool SaveMesh(MeshObject* obj, const char* filename)
{
	return SaveGridToFile(obj->get_grid(), obj->get_subset_handler(), filename,
						  obj->position_attachment());
}

bool ExportToUG3(MeshObject* obj, const char* filenamePrefix, const char* lgmName,
				 const char* problemName)
{
	bool saveOk = false;

	Grid& g = obj->get_grid();
	SubsetHandler& sh = obj->get_subset_handler();

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
