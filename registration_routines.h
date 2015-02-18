// created by Sebastian Reiter
// s.b.reiter@gmail.com

#ifndef __H__UG_registration_routines
#define __H__UG_registration_routines

#include <string>
#include "registry/registry.h"

namespace ug{
namespace promesh{
	void RegisterMesh(bridge::Registry& reg, std::string grp);
	void RegisterCoordinateTransformTools(bridge::Registry& reg, std::string grp);
	void RegisterSelectionTools(bridge::Registry& reg, std::string grp);
	void RegisterMeshingTools(bridge::Registry& reg, std::string grp);
}//	end of namespace
}//	end of namespace

#endif	//__H__registration_routines
