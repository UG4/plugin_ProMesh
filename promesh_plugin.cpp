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
#include "tools/measure_tools.h"
#include "bridge/util.h"
#include "bridge/util_domain_dependent.h"
#include "tooltips.h"

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
	
	reg.add_function("PM_GetSelectionCenter", &GetSelectionCenter, grp, "", "", TOOLTIP_GET_SELECTION_CENTER)
		.add_function("PM_SetSelectionCenter", &SetSelectionCenter, grp, "", "", TOOLTIP_SET_SELECTION_CENTER)
		.add_function("PM_Move", &Move, grp, "", "", TOOLTIP_MOVE)
		.add_function("PM_MoveMeshTo", &MoveMeshTo, grp, "", "", TOOLTIP_MOVE_MESH_TO) 
		.add_function("PM_MoveAlongNormal", &MoveAlongNormal, grp, "", "", TOOLTIP_MOVE_ALONG_NORMAL)  
		.add_function("PM_ScaleAroundCenter", &ScaleAroundCenter, grp, "", "", TOOLTIP_SCALE_AROUND_CENTER) 
		.add_function("PM_ScaleAroundPivot", &ScaleAroundPivot, grp, "", "", TOOLTIP_SCALE_AROUND_PIVOT)	
		.add_function("PM_RotateAroundCenter", &RotateAroundCenter, grp, "", "", TOOLTIP_ROTATE_AROUND_CENTER) 
		.add_function("PM_RotateAroundPivot", &RotateAroundPivot, grp, "", "", TOOLTIP_ROTATE_AROUND_PIVOT)  
		.add_function("PM_ConeTransform", &ConeTransform, grp, "", "", TOOLTIP_CONE_TRANSFORM)  
		.add_function("PM_LaplacianSmooth", &LaplacianSmooth, grp, "", "", TOOLTIP_LAPLACIAN_SMOOTH)  
		.add_function("PM_TangentialSmooth", &TangentialSmooth, grp, "", "", TOOLTIP_TANGENTIAL_SMOOTH) 
		.add_function("PM_ProjectToLimitPLoop", &ProjectToLimitPLoop, grp, "", "", TOOLTIP_PROJECT_TO_LIMIT_PLOOP)  
		.add_function("PM_ProjectToLimitSmoothBoundary", &ProjectToLimitSmoothBoundary, grp, "", "", TOOLTIP_PROJECT_TO_LIMIT_SMOOTH_BOUNDARY) 
		.add_function("PM_SetPivot", &SetPivot, grp, "", "", TOOLTIP_SET_PIVOT) 
		.add_function("PM_SetPivotToSelectionCenter", &SetPivotToSelectionCenter, grp, "", "", TOOLTIP_SET_PIVOT_TO_SELECTION_CENTER)
		.add_function("PM_SetPivotToMeshCenter", &SetPivotToMeshCenter, grp, "", "", TOOLTIP_SET_PIVOT_TO_MESH_CENTER)
		.add_function("PM_FlattenBentQuadrilaterals", &FlattenBentQuadrilaterals, grp, "", "", TOOLTIP_FLATTEN_BENT_QUADRILATERALS);  

//	file io
	reg.add_function("PM_LoadMesh", &LoadMesh, grp, "", "", TOOLTIP_LOAD_MESH)
		.add_function("PM_SaveMesh", &SaveMesh, grp,TOOLTIP_SAVE_MESH)
		.add_function("PM_ExportToUG3", &ExportToUG3, grp, "", "", TOOLTIP_EXPORT_TO_UG3);

//	grid generation tools
	reg.add_function("PM_CreateVertex", &CreateVertex, grp, "", "", TOOLTIP_CREATE_VERTEX)
		.add_function("PM_CreateEdge", &CreateEdge, grp, "", "", TOOLTIP_CREATE_EDGE)
		.add_function("PM_CreateFace", &CreateFace, grp, "", "", TOOLTIP_CREATE_FACE)
		.add_function("PM_CreateVolume", &CreateVolume, grp, "", "", TOOLTIP_CREATE_VOLUME)
		.add_function("PM_CreatePlane", &CreatePlane, grp, "", "", TOOLTIP_CREATE_PLANE)
		.add_function("PM_CreateCircle", &CreateCircle, grp, "", "", TOOLTIP_CREATE_CIRCLE)
		.add_function("PM_CreateBox", &CreateBox, grp, "", "", TOOLTIP_CREATE_BOX)
		.add_function("PM_CreateSphere", &CreateSphere, grp, "", "", TOOLTIP_CREATE_SPHERE)
		.add_function("PM_CreateTetrahedron", &CreateTetrahedron, grp, "", "", TOOLTIP_CREATE_TETRAHEDRON)
		.add_function("PM_CreatePyramid", &CreatePyramid, grp, "", "", TOOLTIP_CREATE_PYRAMID)
		.add_function("PM_CreatePrism", &CreatePrism, grp, "", "", TOOLTIP_CREATE_PRISM);

//	refinement
	reg.add_function("PM_Refine", &Refine, grp, "", "", TOOLTIP_REFINE)
		.add_function("PM_HangingNodeRefine", &HangingNodeRefine, grp, "", "", TOOLTIP_HANGING_NODE_REFINE)
		.add_function("PM_RefineSmooth", &RefineSmooth, grp, "", "", TOOLTIP_REFINE_SMOOTH)
		.add_function("PM_RefineSmoothBoundary2D", &RefineSmoothBoundary2D, grp, "", "", TOOLTIP_REFINE_SMOOTH_BOUNDARY_2D)
		.add_function("PM_CreateFractal", &CreateFractal, grp, "", "", TOOLTIP_CREATE_FRACTAL)
		.add_function("PM_InsertCenter", &InsertCenter, grp, "", "", TOOLTIP_INSERT_CENTER);

//	remeshing
	reg.add_function("PM_ConvertToTriangles", &ConvertToTriangles, grp, "", "", TOOLTIP_CONVERT_TO_TRIANGLES)
		.add_function("PM_TriangleFill", &TriangleFill, grp, "", "", TOOLTIP_TRIANGLE_FILL)
		.add_function("PM_Retriangulate", &Retriangulate, grp, "", "", TOOLTIP_RETRIANGULATE)
		.add_function("PM_AdjustEdgeLength", &AdjustEdgeLength, grp, "", "", TOOLTIP_ADJUST_EDGE_LENGTH)
		.add_function("PM_AdaptSurfaceToCylinder", &AdaptSurfaceToCylinder, grp, "", "", TOOLTIP_ADAPT_SURFACE_TO_CYLINDER)
		.add_function("PM_Tetrahedralize", &Tetrahedralize, grp, "", "", TOOLTIP_TETRAHEDRALIZE)
		.add_function("PM_AssignVolumeConstraints", &AssignVolumeConstraints, grp, "", "", TOOLTIP_ASSIGN_VOLUME_CONSTRAINTS)
		.add_function("PM_ClearVolumeConstraints", &ClearVolumeConstraints, grp, "", "", TOOLTIP_CLEAR_VOLUME_CONSTRAINTS)
		.add_function("PM_Retetrahedralize", &Retetrahedralize, grp, "", "", TOOLTIP_RETETRAHEDRALIZE)
		.add_function("PM_Duplicate", &Duplicate, grp, "", "", TOOLTIP_DUPLICATE)
		.add_function("PM_Extrude", &Extrude, grp, "", "", TOOLTIP_EXTRUDE)
		.add_function("PM_ExtrudeCylinders", &ExtrudeCylinders, grp, "", "", TOOLTIP_EXTRUDE_CYLINDERS)
		.add_function("PM_CreateShrinkGeometry", &CreateShrinkGeometry, grp, "", "", TOOLTIP_CREATE_SHRINK_GEOMETRY)
	    .add_function("PM_ExtrudeFacesWithTets", &ExtrudeFacesWithTets, grp, "", "", TOOLTIP_EXTRUDE_FACES_WITH_TETS);

//	selection tools
	reg.add_function("PM_ClearSelection", &ClearSelection, grp, "", "", TOOLTIP_CLEAR_SELECTION)
		.add_function("PM_SelectAll", &SelectAll, grp, "", "", TOOLTIP_SELECT_ALL) 
		.add_function("PM_ExtendSelection", &ExtendSelection, grp, "", "", TOOLTIP_EXTEND_SELECTION)
		.add_function("PM_SelectSubset", &SelectSubset, grp, "", "", TOOLTIP_SELECT_SUBSET)
		.add_function("PM_SelectSubsetBoundary", &SelectSubsetBoundary, grp, "", "", TOOLTIP_SELECT_SUBSET_BOUNDARY)
		.add_function("PM_SelectUnassignedElements", &SelectUnassignedElements, grp, "", "", TOOLTIP_SELECT_UNASSIGNED_ELEMENTS)
		.add_function("PM_InvertSelection", &InvertSelection, grp, "", "", TOOLTIP_INVERT_SELECTION)
		.add_function("PM_SelectSelectionBoundary", &SelectSelectionBoundary, grp, "", "", TOOLTIP_SELECT_SELECTION_BOUNDARY)
		.add_function("PM_CloseSelection", &CloseSelection, grp, "", "", TOOLTIP_CLOSE_SELECTION)

		.add_function("PM_SelectVertexByCoordinate", &SelectElemByCoordinate<Vertex>, grp,TOOLTIP_SELECT_VERTEX_BY_COORDINATE)
		.add_function("PM_SelectVertexByCylindricalCoordinate", &SelectElemByCylindricalCoordinate<Vertex>, grp,TOOLTIP_SELECT_VERTEX_BY_CYL_COORDINATE)
		.add_function("PM_SelectBoundaryVertices", &SelectBoundaryVertices, grp, "", "", TOOLTIP_SELECT_BOUNDARY_VERTICES)
		.add_function("PM_SelectInnerVertices", &SelectInnerVertices, grp,TOOLTIP_SELECT_INNER_VERTICES)
		.add_function("PM_SelectAssociatedVertices", &SelectAssociatedVertices, grp, "", "", TOOLTIP_SELECT_ASSOCIATED_VERTICES)
		.add_function("PM_SelectAllVertices", &SelectAllVertices, grp, "", "", TOOLTIP_SELECT_ALL_VERTICES)
		.add_function("PM_DeselectAllVertices", &DeselectAllVertices, grp, "", "", TOOLTIP_DESELECT_ALL_VERTICES)
		.add_function("PM_SelectMarkedVertices", &SelectMarkedVertices, grp, "", "", TOOLTIP_SELECT_MARKED_VERTICES)
		.add_function("PM_SelectVertexByIndex", &SelectVertexByIndex, grp, "", "", TOOLTIP_SELECT_VERTEX_BY_INDEX)
		.add_function("PM_SelectUnconnectedVertices", &SelectUnconnectedVertices, grp, "", "", TOOLTIP_SELECT_UNCONNECTED_VERTICES)

		.add_function("PM_SelectEdgeByCoordinate", &SelectElemByCoordinate<Edge>, grp, "", "", TOOLTIP_SELECT_EDGE_BY_COORDINATE)
		.add_function("PM_SelectEdgeByCylindricalCoordinate", &SelectElemByCylindricalCoordinate<Edge>, grp, "", "", TOOLTIP_SELECT_EDGE_BY_CYL_COORDINATE)
		.add_function("PM_SelectBoundaryEdges", &SelectBoundaryEdges, grp, "", "", TOOLTIP_SELECT_BOUNDARY_EDGES)
		.add_function("PM_SelectInnerEdges", &SelectInnerEdges, grp, "", "", TOOLTIP_SELECT_INNER_EDGES)
		.add_function("PM_SelectNonManifoldEdges", &SelectNonManifoldEdges, grp, "", "", TOOLTIP_SELECT_NON_MANIFOLD_EDGES)
		.add_function("PM_SelectSmoothEdgePath", &SelectSmoothEdgePath, grp, "", "", TOOLTIP_SELECT_SMOOTH_EDGE_PATH)
		.add_function("PM_SelectShortEdges", &SelectShortEdges, grp, "", "", TOOLTIP_SELECT_SHORT_EDGES)
		.add_function("PM_SelectLongEdges", &SelectLongEdges, grp, "", "", TOOLTIP_SELECT_LONG_EDGES)
		.add_function("PM_SelectCreaseEdges", &SelectCreaseEdges, grp, "", "", TOOLTIP_SELECT_CREASE_EDGES)
		.add_function("PM_SelectLinkedBoundaryEdges", &SelectLinkedBoundaryEdges, grp, "", "", TOOLTIP_SELECT_LINKED_BOUNDARY_EDGES)
		.add_function("PM_SelectAssociatedEdges", &SelectAssociatedEdges, grp, "", "", TOOLTIP_SELECT_ASSOCIATED_EDGES)
		.add_function("PM_SelectAllEdges", &SelectAllEdges, grp, "", "", TOOLTIP_SELECT_ALL_EDGES)
		.add_function("PM_DeselectAllEdges", &DeselectAllEdges, grp, "", "", TOOLTIP_DESELECT_ALL_EDGES)
		.add_function("PM_SelectMarkedEdges", &SelectMarkedEdges, grp, "", "", TOOLTIP_SELECT_MARKED_EDGES)
		.add_function("PM_SelectEdgeByIndex", &SelectEdgeByIndex, grp, "", "", TOOLTIP_SELECT_EDGE_BY_INDEX)
		.add_function("PM_EdgeSelectionFill", &EdgeSelectionFill, grp, "", "", TOOLTIP_EDGE_SELECTION_FILL)

		.add_function("PM_SelectFaceByCoordinate", &SelectElemByCoordinate<Face>, grp, "", "", TOOLTIP_SELECT_FACE_BY_COORDINATE)
		.add_function("PM_SelectFaceByCylindricalCoordinate", &SelectElemByCylindricalCoordinate<Face>, grp, "", "", TOOLTIP_SELECT_FACE_BY_CYL_COORDINATE)
		.add_function("PM_SelectBoundaryFaces", &SelectBoundaryFaces, grp, "", "", TOOLTIP_SELECT_BOUNDARY_FACES)
		.add_function("PM_SelectInnerFaces", &SelectInnerFaces, grp, "", "", TOOLTIP_SELECT_INNER_FACES)
		.add_function("PM_SelectLinkedManifoldFaces", &SelectLinkedManifoldFaces, grp, "", "", TOOLTIP_SELECT_LINKED_MANIFOLD_FACES)
		.add_function("PM_SelectLinkedBoundaryFaces", &SelectLinkedBoundaryFaces, grp, "", "", TOOLTIP_SELECT_LINKED_BOUNDARY_FACES)
		.add_function("PM_SelectDegenerateFaces", &SelectDegenerateFaces, grp, "", "", TOOLTIP_SELECT_DEGENERATE_FACES)
		.add_function("PM_SelectLinkedFlatFaces", &SelectLinkedFlatFaces, grp, "", "", TOOLTIP_SELECT_LINKED_FLAT_FACES)
		.add_function("PM_SelectIntersectingTriangles", &SelectIntersectingTriangles, grp, "", "", TOOLTIP_SELECT_INTERSECTING_TRIANGLES)
		.add_function("PM_SelectAssociatedFaces", &SelectAssociatedFaces, grp, "", "", TOOLTIP_SELECT_ASSOCIATED_FACES)
		.add_function("PM_SelectAllFaces", &SelectAllFaces, grp, "", "", TOOLTIP_SELECT_ALL_FACES)
		.add_function("PM_DeselectAllFaces", &DeselectAllFaces, grp, "", "", TOOLTIP_DESELECT_ALL_FACES)
		.add_function("PM_SelectFaceByIndex", &SelectFaceByIndex, grp, "", "", TOOLTIP_SELECT_FACE_BY_INDEX)
		.add_function("PM_SelectFacesByNormal", &SelectFacesByNormal, grp, "", "", TOOLTIP_SELECT_FACES_BY_NORMAL)
		.add_function("PM_FaceSelectionFill", &FaceSelectionFill, grp, "", "", TOOLTIP_FACE_SELECTION_FILL)
		.add_function("PM_SelectBentQuadrilaterals", &SelectBentQuadrilaterals, grp, "", "", TOOLTIP_SELECT_BENT_QUADRILATERALS)

		.add_function("PM_SelectVolumeByCoordinate", &SelectElemByCoordinate<Volume>, grp, "", "", TOOLTIP_SELECT_VOLUME_BY_COORDINATE)
		.add_function("PM_SelectVolumeByCylindricalCoordinate", &SelectElemByCylindricalCoordinate<Volume>, grp, "", "", TOOLTIP_SELECT_VOLUME_BY_CYL_COORDINATE)
		.add_function("PM_SelectAllVolumes", &SelectAllVolumes, grp, "", "", TOOLTIP_SELECT_ALL_VOLUMES)
		.add_function("PM_DeselectAllVolumes", &DeselectAllVolumes, grp, "", "", TOOLTIP_DESELECT_ALL_VOLUMES)
		.add_function("PM_SelectUnorientableVolumes", &SelectUnorientableVolumes, grp, "", "", TOOLTIP_SELECT_UNORIENTABLE_VOLUMES)
		.add_function("PM_SelectVolumeByIndex", &SelectVolumeByIndex, grp, "", "", TOOLTIP_SELECT_VOLUME_BY_INDEX)
		.add_function("PM_VolumeSelectionFill", &VolumeSelectionFill, grp, "", "", TOOLTIP_VOLUME_SELECTION_FILL)

		.add_function("PM_MarkSelection", &MarkSelection, grp, "", "", TOOLTIP_MARK_SELECTION)
		.add_function("PM_UnmarkSelection", &UnmarkSelection, grp, "", "", TOOLTIP_UNMARK_SELECTION);

//	subset tools
	reg.add_function("PM_AssignSubset", &AssignSubset, grp, "", "", TOOLTIP_ASSIGN_SUBSET)
		.add_function("PM_SetSubsetName", &SetSubsetName, grp, "", "", TOOLTIP_SET_SUBSET_NAME)
		.add_function("PM_AssignSubsetColors", &AssignSubsetColors, grp, "", "", TOOLTIP_ASSIGN_SUBSET_COLORS)
		.add_function("PM_SeparateFacesByEdgeSubsets", &SeparateFacesByEdgeSubsets, grp, "", "", TOOLTIP_SEPARATE_FACES_BY_EDGE_SUBSETS)
		.add_function("PM_SeparateFacesBySelectedEdges", &SeparateFacesBySelectedEdges, grp, "", "", TOOLTIP_SEPARATE_FACES_BY_SELECTED_EDGES)
		.add_function("PM_SeparateVolumesByFaceSubsets", &SeparateVolumesByFaceSubsets, grp, "", "", TOOLTIP_SEPARATE_VOLUMES_BY_FACE_SUBSETS)
		.add_function("PM_SeparateVolumesBySelectedFaces", &SeparateVolumesBySelectedFaces, grp, "", "", TOOLTIP_SEPARATE_VOLUMES_BY_SELECTED_FACES)
		.add_function("PM_SeparateIrregularManifoldSubsets", &SeparateIrregularManifoldSubsets, grp, "", "", TOOLTIP_SEPARATE_IRREGULAR_MANIFOLD_SUBSETS)
		.add_function("PM_MoveSubset", &MoveSubset, grp, "", "", TOOLTIP_MOVE_SUBSET)
		.add_function("PM_SwapSubsets", &SwapSubsets, grp, "", "", TOOLTIP_SWAP_SUBSETS)
		.add_function("PM_JoinSubsets", &JoinSubsets, grp, "", "", TOOLTIP_JOIN_SUBSETS)
		.add_function("PM_EraseSubset", &EraseSubset, grp, "", "", TOOLTIP_ERASE_SUBSET)
		.add_function("PM_EraseEmptySubsets", &EraseEmptySubsets, grp, "", "", TOOLTIP_ERASE_EMPTY_SUBSETS)
		.add_function("PM_AdjustSubsetsForUG3", &AdjustSubsetsForUG3, grp, "", "", TOOLTIP_ADJUST_SUBSETS_FOR_UG3)
		.add_function("PM_AdjustSubsetsForUG4", &AdjustSubsetsForUG4, grp, "", "", TOOLTIP_ADJUST_SUBSETS_FOR_UG4)
		.add_function("PM_SeparateFaceSubsetsByNormal", &SeparateFaceSubsetsByNormal, grp,TOOLTIP_SEPARATE_FACE_SUBSETS_BY_NORMAL )
		.add_function("PM_SeparateFaceSubsetByNormal", &SeparateFaceSubsetByNormal, grp, "", "", TOOLTIP_SEPARATE_FACE_SUBSET_BY_NORMAL)
		.add_function("PM_AssignSubsetsByQuality", &AssignSubsetsByQuality, grp, "", "", TOOLTIP_ASSIGN_SUBSETS_BY_QUALITY)
		.add_function("PM_SeparateDegeneratedBoundaryFaceSubsets", &SeparateDegeneratedBoundaryFaceSubsets, grp, "", "", TOOLTIP_SEPARATE_DEGENERATED_BOUNDARY_FACE_SUBSETS)
		.add_function("PM_AssignSubsetsByElementType", &AssignSubsetsByElementType, grp, "", "", TOOLTIP_ASSIGN_SUBSETS_BY_ELEMENT_TYPE);

//	topology
	reg.add_function("PM_EraseSelectedElements", &EraseSelectedElements, grp, "", "", TOOLTIP_ERASE_SELECTED_ELEMENTS)
		.add_function("PM_RemoveDoubles", &RemoveDoubles, grp, "", "", TOOLTIP_REMOVE_DOUBLES)
		.add_function("PM_RemoveDoubleEdges", &RemoveDoubleEdges, grp, "", "", TOOLTIP_REMOVE_DOUBLE_EDGES)
		.add_function("PM_MergeAtFirst", &MergeAtFirst, grp, "", "", TOOLTIP_MERGE_AT_FIRST)
		.add_function("PM_MergeAtCenter", &MergeAtCenter, grp, "", "", TOOLTIP_MERGE_AT_CENTER)
		.add_function("PM_MergeAtLast", &MergeAtLast, grp, "", "", TOOLTIP_MERGE_AT_LAST)
		.add_function("PM_CollapseEdge", &CollapseEdge, grp, "", "", TOOLTIP_COLLAPSE_EDGE)
		.add_function("PM_SplitEdge", &SplitEdge, grp, "", "", TOOLTIP_SPLIT_EDGE)
		.add_function("PM_SwapEdge", &SwapEdge, grp, "", "", TOOLTIP_SWAP_EDGE)
		.add_function("PM_PlaneCut", &PlaneCut, grp, "", "", TOOLTIP_PLANE_CUT)
		.add_function("PM_AdjustEdgeOrientation", &AdjustEdgeOrientation, grp, "", "", TOOLTIP_ADJUST_EDGE_ORIENTATION)
		.add_function("PM_FixFaceOrientation", &FixFaceOrientation, grp, "", "", TOOLTIP_FIX_FACE_ORIENTATION)
		.add_function("PM_FixFaceSubsetOrientations", &FixFaceSubsetOrientations, grp, "", "", TOOLTIP_FIX_FACE_SUBSET_ORIENTATIONS)
		.add_function("PM_FixVolumeOrientation", &FixVolumeOrientation, grp, "", "", TOOLTIP_FIX_VOLUME_ORIENTATION)
		.add_function("PM_InvertFaceOrientation", &InvertFaceOrientation, grp, "", "", TOOLTIP_INVERT_FACE_ORIENTATION)
		.add_function("PM_ResolveSelfIntersections", &ResolveSelfIntersections, grp, "", "", TOOLTIP_RESOLVE_SELF_INTERSECTIONS)
		.add_function("PM_ResolveEdgeIntersection", &ResolveEdgeIntersection, grp, "", "", TOOLTIP_RESOLVE_EDGE_INTERSECTIONS) 
		.add_function("PM_ResolveTriangleIntersections", &ResolveTriangleIntersections, grp, "", "", TOOLTIP_RESOLVE_TRIANGLE_INTERSECTIONS)
		.add_function("PM_ProjectVerticesToCloseEdges", &ProjectVerticesToCloseEdges, grp, "", "", TOOLTIP_PROJECT_VERTICES_TO_CLOSE_EDGES)
		.add_function("PM_ProjectVerticesToCloseFaces", &ProjectVerticesToCloseFaces, grp, "", "", TOOLTIP_PROJECT_VERTICES_TO_CLOSE_FACES)
		.add_function("PM_IntersectCloseEdges", &IntersectCloseEdges, grp, "", "", TOOLTIP_INTERSECT_CLOSE_EDGES);

//	info tools
		reg.add_function("PM_MeasureGridLength", &MeasureGridLength, grp, "", "", TOOLTIP_MEASURE_GRID_LENGTH)
		   .add_function("PM_MeasureGridArea", &MeasureGridArea, grp, "", "", TOOLTIP_MEASURE_GRID_AREA)
		   .add_function("PM_MeasureGridVolume", &MeasureGridVolume, grp, "", "", TOOLTIP_MEASURE_GRID_VOLUME)
		   .add_function("PM_MeasureSubsetLength", &MeasureSubsetLength, grp, "", "", TOOLTIP_MEASURE_SUBSET_LENGTH)
		   .add_function("PM_MeasureSubsetArea", &MeasureSubsetArea, grp, "", "", TOOLTIP_MEASURE_SUBSET_AREA)
		   .add_function("PM_MeasureSubsetVolume", &MeasureSubsetVolume, grp, "", "", TOOLTIP_MEASURE_SUBSET_VOLUME)
		   .add_function("PM_MeasureSelectionLength", &MeasureSelectionLength, grp, "", "", TOOLTIP_MEASURE_SELECTION_LENGTH)
		   .add_function("PM_MeasureSelectionArea", &MeasureSelectionArea, grp, "", "", TOOLTIP_MEASURE_SELECTION_AREA)
		   .add_function("PM_MeasureSelectionVolume", &MeasureSelectionVolume, grp, "", "", TOOLTIP_MEASURE_SELECTION_VOLUME);

//	new tools
	reg.add_class_<Box>("PM_Box", grp)
		.add_method("set_min", &Box::set_min)
		.add_method("set_max", &Box::set_max)
		.add_method("min", &Box::get_min)
		.add_method("max", &Box::get_max);

	reg.add_function("PM_GetBoundingBox", &GetBoundingBox, grp, "", "", TOOLTIP_GET_BOUNDING_BOX)
		.add_function("PM_SelectVerticesInBox", &SelectElementsInBox<Vertex>, grp, "", "", TOOLTIP_SELECT_VERTEX_IN_BOX)
		.add_function("PM_SelectEdgesInBox", &SelectElementsInBox<Edge>, grp, "", "", TOOLTIP_SELECT_EDGE_IN_BOX)  
		.add_function("PM_SelectFacesInBox", &SelectElementsInBox<Face>, grp, "", "", TOOLTIP_SELECT_FACE_IN_BOX)
		.add_function("PM_SelectVolumesInBox", &SelectElementsInBox<Volume>, grp, "", "", TOOLTIP_SELECT_VOLUME_IN_BOX)
		.add_function("PM_SelectVerticesInCylinder", &SelectElementsInCylinder<Vertex>, grp, "", "", TOOLTIP_SELECT_VERTEX_IN_CYLINDER)
		.add_function("PM_SelectEdgesInCylinder", &SelectElementsInCylinder<Edge>, grp, "", "", TOOLTIP_SELECT_EDGE_IN_CYLINDER)
		.add_function("PM_SelectFacesInCylinder", &SelectElementsInCylinder<Face>, grp, "", "", TOOLTIP_SELECT_FACE_IN_CYLINDER)
		.add_function("PM_SelectVolumesInCylinder", &SelectElementsInCylinder<Volume>, grp, "", "", TOOLTIP_SELECT_VOLUME_IN_CYLINDER);

	reg.add_function("PM_ScaleAroundPoint", &ScaleAroundPoint, grp, "", "", TOOLTIP_SCALE_AROUND_PIVOT);
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

