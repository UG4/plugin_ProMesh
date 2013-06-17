// created by Sebastian Reiter
// s.b.reiter@gmail.com
// Jun 14, 2013 (d,m,y)

#include "mesh_object.h"
#include "tools/grid_generation_tools.h"
#include "tools/coordinate_transform_tools.h"
#include "bridge/util.h"
#include "bridge/util_domain_dependent.h"

using namespace std;
using namespace ug::bridge;

namespace ug{
namespace promesh{

/**
 *  \defgroup promesh_plugin ProMesh Plugin
 *  \ingroup plugins_experimental
 *  The promesh plugin gives access to many functions and tools which are contained
 *  in promesh.
 *  \{
 */

/**
 * Class exporting the functionality. All functionality that is to
 * be used in scripts or visualization must be registered here.
 */
struct Functionality
{

/**
 * Function called for the registration of Domain dependent parts.
 * All Functions and Classes depending on the Domain
 * are to be placed here when registering. The method is called for all
 * available Domain types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <typename TDomain>
static void Domain(Registry& reg, string grp)
{
	string suffix = GetDomainSuffix<TDomain>();
	string tag = GetDomainTag<TDomain>();

}

/**
 * Function called for the registration of Dimension dependent parts.
 * All Functions and Classes depending on the Dimension
 * are to be placed here when registering. The method is called for all
 * available Dimension types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <int dim>
static void Dimension(Registry& reg, string grp)
{
	string suffix = GetDimensionSuffix<dim>();
	string tag = GetDimensionTag<dim>();

}


/**
 * Function called for the registration of Domain and Algebra independent parts.
 * All Functions and Classes not depending on Domain and Algebra
 * are to be placed here when registering.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
static void Common(Registry& reg, string grp)
{
	reg.add_class_<MeshObject>("PM_MeshObject", grp)
		.add_constructor()
		.add_method("get_grid", &MeshObject::get_grid, grp)
		.add_method("get_subset_handler", &MeshObject::get_subset_handler, grp)
		.add_method("get_crease_handler", &MeshObject::get_crease_handler, grp)
		.add_method("get_selector", &MeshObject::get_selector, grp)
		.add_method("set_pivot", &MeshObject::set_pivot, grp)
		.add_method("get_pivot", &MeshObject::get_pivot, grp)
		.set_construct_as_smart_pointer(true);

//	grid generation tools
	reg.add_function("PM_CreateVertex", &CreateVertex, grp)
		.add_function("PM_CreateEdge", &CreateEdge, grp)
		.add_function("PM_CreateFace", &CreateFace, grp)
		.add_function("PM_CreateVolume", &CreateVolume, grp)
		.add_function("PM_CreatePlane", &CreatePlane, grp)
		.add_function("PM_CreateCircle", &CreateCircle, grp)
		.add_function("PM_CreateBox", &CreateBox, grp)
		.add_function("PM_CreateSphere", &CreateSphere, grp)
		.add_function("PM_CreateTetrahedron", &CreateTetrahedron, grp)
		.add_function("PM_CreatePyramid", &CreatePyramid, grp)
		.add_function("PM_CreatePrism", &CreatePrism, grp);

//	coordinate transform tools
	reg.add_function("PM_GetSelectionCenter", &GetSelectionCenter, grp)
		.add_function("PM_SetSelectionCenter", &SetSelectionCenter, grp)
		.add_function("PM_Move", &Move, grp)
		.add_function("PM_MoveAlongNormal", &MoveAlongNormal, grp)
		.add_function("PM_ScaleAroundCenter", &ScaleAroundCenter, grp)
		.add_function("PM_ScaleAroundPivot", &ScaleAroundPivot, grp)
		.add_function("PM_RotateAroundCenter", &RotateAroundCenter, grp)
		.add_function("PM_RotateAroundPivot", &RotateAroundPivot, grp)
		.add_function("PM_ConeTransform", &ConeTransform, grp)
		.add_function("PM_LaplacianSmooth", &LaplacianSmooth, grp)
		.add_function("PM_ProjectToLimitPLoop", &ProjectToLimitPLoop, grp)
		.add_function("PM_ProjectToLimitSmoothBoundary", &ProjectToLimitSmoothBoundary, grp)
		.add_function("PM_SetPivot", &SetPivot, grp)
		.add_function("PM_SetPivotToCenter", &SetPivotToCenter, grp)
		.add_function("PM_FlattenBentQuadrilaterals", &FlattenBentQuadrilaterals, grp);
}

}; // end Functionality

// end group promesh plugin
/// \}

} // end namespace promesh


/**
 * This function is called when the plugin is loaded.
 */
extern "C" void
InitUGPlugin_ProMesh(Registry* reg, string grp)
{
	grp.append("promesh");
	typedef promesh::Functionality Functionality;

	try{
		RegisterCommon<Functionality>(*reg,grp);
		RegisterDimensionDependent<Functionality>(*reg,grp);
		RegisterDomainDependent<Functionality>(*reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

}// end of namespace
