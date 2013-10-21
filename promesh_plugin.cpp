// created by Sebastian Reiter
// s.b.reiter@gmail.com
// Jun 14, 2013 (d,m,y)

#include "mesh_object.h"
#include "tools/coordinate_transform_tools.h"
#include "tools/file_io_tools.h"
#include "tools/grid_generation_tools.h"
#include "tools/new_tools.h"
#include "tools/refinement_tools.h"
#include "tools/remeshing_tools.h"
#include "tools/selection_tools.h"
#include "tools/subset_tools.h"
#include "tools/topology_tools.h"
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

//	file io
	reg.add_function("PM_LoadMesh", &LoadMesh, grp)
		.add_function("PM_SaveMesh", &SaveMesh, grp)
		.add_function("PM_ExportToUG3", &ExportToUG3, grp);

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

//	refinement
	reg.add_function("PM_Refine", &Refine, grp)
		.add_function("PM_HangingNodeRefine", &HangingNodeRefine, grp)
		.add_function("PM_RefineSmooth", &RefineSmooth, grp)
		.add_function("PM_RefineSmoothBoundary2D", &RefineSmoothBoundary2D, grp)
		.add_function("PM_CreateFractal", &CreateFractal, grp)
		.add_function("PM_InsertCenter", &InsertCenter, grp);

//	remeshing
	reg.add_function("PM_ConvertToTriangles", &ConvertToTriangles, grp)
		.add_function("PM_TriangleFill", &TriangleFill, grp)
		.add_function("PM_Retriangulate", &Retriangulate, grp)
		.add_function("PM_AdjustEdgeLength", &AdjustEdgeLength, grp)
		.add_function("PM_AdaptSurfaceToCylinder", &AdaptSurfaceToCylinder, grp)
		.add_function("PM_Tetrahedralize", &Tetrahedralize, grp)
		.add_function("PM_AssignVolumeConstraints", &AssignVolumeConstraints, grp)
		.add_function("PM_ClearVolumeConstraints", &ClearVolumeConstraints, grp)
		.add_function("PM_Retetrahedralize", &Retetrahedralize, grp)
		.add_function("PM_Duplicate", &Duplicate, grp)
		.add_function("PM_Extrude", &Extrude, grp)
		.add_function("PM_ExtrudeCylinders", &ExtrudeCylinders, grp);

//	selection tools
	reg.add_function("PM_ClearSelection", &ClearSelection, grp)
		.add_function("PM_SelectAll", &SelectAll, grp)
		.add_function("PM_ExtendSelection", &ExtendSelection, grp)
		.add_function("PM_SelectSubset", &SelectSubset, grp)
		.add_function("PM_SelectSubsetBoundary", &SelectSubsetBoundary, grp)
		.add_function("PM_SelectUnassignedElements", &SelectUnassignedElements, grp)
		.add_function("PM_InvertSelection", &InvertSelection, grp)
		.add_function("PM_SelectSelectionBoundary", &SelectSelectionBoundary, grp)
		.add_function("PM_CloseSelection", &CloseSelection, grp)

		.add_function("PM_SelectVertexByCoordinate", &SelectElemByCoordinate<VertexBase>, grp)
		.add_function("PM_SelectBoundaryVertices", &SelectBoundaryVertices, grp)
		.add_function("PM_SelectInnerVertices", &SelectInnerVertices, grp)
		.add_function("PM_SelectAssociatedVertices", &SelectAssociatedVertices, grp)
		.add_function("PM_SelectAllVertices", &SelectAllVertices, grp)
		.add_function("PM_DeselectAllVertices", &DeselectAllVertices, grp)
		.add_function("PM_SelectMarkedVertices", &SelectMarkedVertices, grp)
		.add_function("PM_SelectVertexByIndex", &SelectVertexByIndex, grp)
		.add_function("PM_SelectUnconnectedVertices", &SelectUnconnectedVertices, grp)

		.add_function("PM_SelectEdgeByCoordinate", &SelectElemByCoordinate<EdgeBase>, grp)
		.add_function("PM_SelectBoundaryEdges", &SelectBoundaryEdges, grp)
		.add_function("PM_SelectInnerEdges", &SelectInnerEdges, grp)
		.add_function("PM_SelectNonManifoldEdges", &SelectNonManifoldEdges, grp)
		.add_function("PM_SelectSmoothEdgePath", &SelectSmoothEdgePath, grp)
		.add_function("PM_SelectShortEdges", &SelectShortEdges, grp)
		.add_function("PM_SelectLongEdges", &SelectLongEdges, grp)
		.add_function("PM_SelectLinkedBoundaryEdges", &SelectLinkedBoundaryEdges, grp)
		.add_function("PM_SelectAssociatedEdges", &SelectAssociatedEdges, grp)
		.add_function("PM_SelectAllEdges", &SelectAllEdges, grp)
		.add_function("PM_DeselectAllEdges", &DeselectAllEdges, grp)
		.add_function("PM_SelectMarkedEdges", &SelectMarkedEdges, grp)
		.add_function("PM_SelectEdgeByIndex", &SelectEdgeByIndex, grp)
		.add_function("PM_EdgeSelectionFill", &EdgeSelectionFill, grp)

		.add_function("PM_SelectFaceByCoordinate", &SelectElemByCoordinate<Face>, grp)
		.add_function("PM_SelectBoundaryFaces", &SelectBoundaryFaces, grp)
		.add_function("PM_SelectInnerFaces", &SelectInnerFaces, grp)
		.add_function("PM_SelectLinkedManifoldFaces", &SelectLinkedManifoldFaces, grp)
		.add_function("PM_SelectLinkedBoundaryFaces", &SelectLinkedBoundaryFaces, grp)
		.add_function("PM_SelectDegenerateFaces", &SelectDegenerateFaces, grp)
		.add_function("PM_SelectLinkedFlatFaces", &SelectLinkedFlatFaces, grp)
		.add_function("PM_SelectIntersectingTriangles", &SelectIntersectingTriangles, grp)
		.add_function("PM_SelectAssociatedFaces", &SelectAssociatedFaces, grp)
		.add_function("PM_SelectAllFaces", &SelectAllFaces, grp)
		.add_function("PM_DeselectAllFaces", &DeselectAllFaces, grp)
		.add_function("PM_SelectFaceByIndex", &SelectFaceByIndex, grp)
		.add_function("PM_FaceSelectionFill", &FaceSelectionFill, grp)
		.add_function("PM_SelectBentQuadrilaterals", &SelectBentQuadrilaterals, grp)

		.add_function("PM_SelectVolumeByCoordinate", &SelectElemByCoordinate<Volume>, grp)
		.add_function("PM_SelectAllVolumes", &SelectAllVolumes, grp)
		.add_function("PM_DeselectAllVolumes", &DeselectAllVolumes, grp)
		.add_function("PM_SelectUnorientableVolumes", &SelectUnorientableVolumes, grp)
		.add_function("PM_SelectVolumeByIndex", &SelectVolumeByIndex, grp)
		.add_function("PM_VolumeSelectionFill", &VolumeSelectionFill, grp);

//	subset tools
	reg.add_function("PM_AssignSubset", &AssignSubset, grp)
		.add_function("PM_SetSubsetName", &SetSubsetName, grp)
		.add_function("PM_AssignSubsetColors", &AssignSubsetColors, grp)
		.add_function("PM_SeparateFacesByEdgeSubsets", &SeparateFacesByEdgeSubsets, grp)
		.add_function("PM_SeparateFacesBySelectedEdges", &SeparateFacesBySelectedEdges, grp)
		.add_function("PM_SeparateVolumesByFaceSubsets", &SeparateVolumesByFaceSubsets, grp)
		.add_function("PM_SeparateVolumesBySelectedFaces", &SeparateVolumesBySelectedFaces, grp)
		.add_function("PM_SeparateIrregularManifoldSubsets", &SeparateIrregularManifoldSubsets, grp)
		.add_function("PM_MoveSubset", &MoveSubset, grp)
		.add_function("PM_SwapSubsets", &SwapSubsets, grp)
		.add_function("PM_JoinSubsets", &JoinSubsets, grp)
		.add_function("PM_EraseSubset", &EraseSubset, grp)
		.add_function("PM_EraseEmptySubsets", &EraseEmptySubsets, grp)
		.add_function("PM_AdjustSubsetsForUG3", &AdjustSubsetsForUG3, grp)
		.add_function("PM_AdjustSubsetsForUG4", &AdjustSubsetsForUG4, grp)
		.add_function("PM_SeparateFaceSubsetsByNormal", &SeparateFaceSubsetsByNormal, grp)
		.add_function("PM_SeparateFaceSubsetByNormal", &SeparateFaceSubsetByNormal, grp)
		.add_function("PM_AssignSubsetsByQuality", &AssignSubsetsByQuality, grp)
		.add_function("PM_SeparateDegeneratedBoundaryFaceSubsets", &SeparateDegeneratedBoundaryFaceSubsets, grp)
		.add_function("PM_AssignSubsetsByElementType", &AssignSubsetsByElementType, grp);

//	topology
	reg.add_function("PM_EraseSelectedElements", &EraseSelectedElements, grp)
		.add_function("PM_RemoveDoubles", &RemoveDoubles, grp)
		.add_function("PM_RemoveDoubleEdges", &RemoveDoubleEdges, grp)
		.add_function("PM_MergeAtFirst", &MergeAtFirst, grp)
		.add_function("PM_MergeAtCenter", &MergeAtCenter, grp)
		.add_function("PM_MergeAtLast", &MergeAtLast, grp)
		.add_function("PM_CollapseEdge", &CollapseEdge, grp)
		.add_function("PM_SplitEdge", &SplitEdge, grp)
		.add_function("PM_SwapEdge", &SwapEdge, grp)
		.add_function("PM_PlaneCut", &PlaneCut, grp)
		.add_function("PM_AdjustEdgeOrientation", &AdjustEdgeOrientation, grp)
		.add_function("PM_FixFaceOrientation", &FixFaceOrientation, grp)
		.add_function("PM_FixFaceSubsetOrientations", &FixFaceSubsetOrientations, grp)
		.add_function("PM_FixVolumeOrientation", &FixVolumeOrientation, grp)
		.add_function("PM_InvertFaceOrientation", &InvertFaceOrientation, grp)
		.add_function("PM_ResolveEdgeIntersection", &ResolveEdgeIntersection, grp)
		.add_function("PM_ResolveTriangleIntersections", &ResolveTriangleIntersections, grp)
		.add_function("PM_ProjectVerticesToCloseEdges", &ProjectVerticesToCloseEdges, grp)
		.add_function("PM_ProjectVerticesToCloseFaces", &ProjectVerticesToCloseFaces, grp)
		.add_function("PM_IntersectCloseEdges", &IntersectCloseEdges, grp);

//	new tools
	reg.add_function("PM_SelectVerticesInBox", &SelectElementsInBox<VertexBase>, grp)
		.add_function("PM_SelectEdgesInBox", &SelectElementsInBox<EdgeBase>, grp)
		.add_function("PM_SelectFacesInBox", &SelectElementsInBox<Face>, grp)
		.add_function("PM_SelectVolumesInBox", &SelectElementsInBox<Volume>, grp);

	reg.add_function("PM_ScaleAroundPoint", &ScaleAroundPoint, grp);
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

