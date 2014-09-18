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
static void Common(Registry& reg, string baseGrp)
{
	string grp = baseGrp;

	reg.add_class_<MeshObject>("MeshObject", grp)
		.add_constructor()
		.add_method("get_grid", &MeshObject::get_grid, grp)
		.add_method("get_subset_handler", &MeshObject::get_subset_handler, grp)
		.add_method("get_crease_handler", &MeshObject::get_crease_handler, grp)
		.add_method("get_selector", &MeshObject::get_selector, grp)
		.add_method("set_pivot", &MeshObject::set_pivot, grp)
		.add_method("get_pivot", &MeshObject::get_pivot, grp)
		.set_construct_as_smart_pointer(true);

//	coordinate transform tools
	grp = baseGrp + string("/Coordinate Transform");
	reg.add_function("GetSelectionCenter", &GetSelectionCenter, grp, "", "", TOOLTIP_GET_SELECTION_CENTER)
		.add_function("SetSelectionCenter", &SetSelectionCenter, grp, "", "", TOOLTIP_SET_SELECTION_CENTER)
		.add_function("Move", &Move, grp, "", "", TOOLTIP_MOVE)
		.add_function("MoveMeshTo", &MoveMeshTo, grp, "", "", TOOLTIP_MOVE_MESH_TO) 
		.add_function("MoveAlongNormal", &MoveAlongNormal, grp, "", "", TOOLTIP_MOVE_ALONG_NORMAL)  
		.add_function("ScaleAroundCenter", &ScaleAroundCenter, grp, "", "", TOOLTIP_SCALE_AROUND_CENTER) 
		.add_function("ScaleAroundPivot", &ScaleAroundPivot, grp, "", "", TOOLTIP_SCALE_AROUND_PIVOT)	
		.add_function("RotateAroundCenter", &RotateAroundCenter, grp, "", "", TOOLTIP_ROTATE_AROUND_CENTER) 
		.add_function("RotateAroundPivot", &RotateAroundPivot, grp, "", "", TOOLTIP_ROTATE_AROUND_PIVOT)  
		.add_function("ConeTransform", &ConeTransform, grp, "", "", TOOLTIP_CONE_TRANSFORM)  
		.add_function("LaplacianSmooth", &LaplacianSmooth, grp, "", "", TOOLTIP_LAPLACIAN_SMOOTH)  
		.add_function("WeightedEdgeSmooth", &WeightedEdgeSmooth, grp, "",
					  "mesh#"
					  "alpha | default | min=0D; max=1D; value=0.25D#"
					  "num iterations | default | min=0; value=10",
					  TOOLTIP_WEIGHTED_EDGE_SMOOTH)  
		.add_function("WeightedFaceSmooth", &WeightedFaceSmooth, grp, "",
					  "mesh#"
					  "alpha | default | min=0D; max=1D; value=0.25D#"
					  "num iterations | default | min=0; value=10",
					  TOOLTIP_WEIGHTED_FACE_SMOOTH)  
		.add_function("TangentialSmooth", &TangentialSmooth, grp, "", "", TOOLTIP_TANGENTIAL_SMOOTH) 
		.add_function("ProjectToLimitPLoop", &ProjectToLimitPLoop, grp, "", "", TOOLTIP_PROJECT_TO_LIMIT_PLOOP)  
		.add_function("ProjectToLimitSmoothBoundary", &ProjectToLimitSmoothBoundary, grp, "", "", TOOLTIP_PROJECT_TO_LIMIT_SMOOTH_BOUNDARY) 
		.add_function("SetPivot", &SetPivot, grp, "", "", TOOLTIP_SET_PIVOT) 
		.add_function("SetPivotToSelectionCenter", &SetPivotToSelectionCenter, grp, "", "", TOOLTIP_SET_PIVOT_TO_SELECTION_CENTER)
		.add_function("SetPivotToMeshCenter", &SetPivotToMeshCenter, grp, "", "", TOOLTIP_SET_PIVOT_TO_MESH_CENTER)
		.add_function("FlattenBentQuadrilaterals", &FlattenBentQuadrilaterals, grp, "", "", TOOLTIP_FLATTEN_BENT_QUADRILATERALS);  

	
//	file io
	grp = baseGrp;
	reg.add_function("LoadMesh", &LoadMesh, grp, "", "", TOOLTIP_LOAD_MESH)
		.add_function("SaveMesh", &SaveMesh, grp,TOOLTIP_SAVE_MESH)
		.add_function("ExportToUG3", &ExportToUG3, grp, "", "", TOOLTIP_EXPORT_TO_UG3);

//	grid generation tools
	reg.add_function("CreateVertex", &CreateVertex, grp, "", "", TOOLTIP_CREATE_VERTEX)
		.add_function("CreateEdge", &CreateEdge, grp, "", "", TOOLTIP_CREATE_EDGE)
		.add_function("CreateFace", &CreateFace, grp, "", "", TOOLTIP_CREATE_FACE)
		.add_function("CreateVolume", &CreateVolume, grp, "", "", TOOLTIP_CREATE_VOLUME)
		.add_function("CreatePlane", &CreatePlane, grp, "", "", TOOLTIP_CREATE_PLANE)
		.add_function("CreateCircle", &CreateCircle, grp, "", "", TOOLTIP_CREATE_CIRCLE)
		.add_function("CreateBox", &CreateBox, grp, "", "", TOOLTIP_CREATE_BOX)
		.add_function("CreateSphere", &CreateSphere, grp, "", "", TOOLTIP_CREATE_SPHERE)
		.add_function("CreateTetrahedron", &CreateTetrahedron, grp, "", "", TOOLTIP_CREATE_TETRAHEDRON)
		.add_function("CreatePyramid", &CreatePyramid, grp, "", "", TOOLTIP_CREATE_PYRAMID)
		.add_function("CreatePrism", &CreatePrism, grp, "", "", TOOLTIP_CREATE_PRISM);

//	refinement
	reg.add_function("Refine", &Refine, grp, "", "", TOOLTIP_REFINE)
		.add_function("HangingNodeRefine", &HangingNodeRefine, grp, "", "", TOOLTIP_HANGING_NODE_REFINE)
		.add_function("RefineSmooth", &RefineSmooth, grp, "", "", TOOLTIP_REFINE_SMOOTH)
		.add_function("RefineSmoothBoundary2D", &RefineSmoothBoundary2D, grp, "", "", TOOLTIP_REFINE_SMOOTH_BOUNDARY_2D)
		.add_function("CreateFractal", &CreateFractal, grp, "", "", TOOLTIP_CREATE_FRACTAL)
		.add_function("InsertCenter", &InsertCenter, grp, "", "", TOOLTIP_INSERT_CENTER);

//	remeshing
	reg.add_function("ConvertToTriangles", &ConvertToTriangles, grp, "", "", TOOLTIP_CONVERT_TO_TRIANGLES)
		.add_function("TriangleFill", &TriangleFill, grp, "", "", TOOLTIP_TRIANGLE_FILL)
		.add_function("Retriangulate", &Retriangulate, grp, "", "", TOOLTIP_RETRIANGULATE)
		.add_function("AdjustEdgeLength", &AdjustEdgeLength, grp, "", "", TOOLTIP_ADJUST_EDGE_LENGTH)
		.add_function("AdaptSurfaceToCylinder", &AdaptSurfaceToCylinder, grp, "", "", TOOLTIP_ADAPT_SURFACE_TO_CYLINDER)
		.add_function("Tetrahedralize", &Tetrahedralize, grp, "", "", TOOLTIP_TETRAHEDRALIZE)
		.add_function("AssignVolumeConstraints", &AssignVolumeConstraints, grp, "", "", TOOLTIP_ASSIGN_VOLUME_CONSTRAINTS)
		.add_function("ClearVolumeConstraints", &ClearVolumeConstraints, grp, "", "", TOOLTIP_CLEAR_VOLUME_CONSTRAINTS)
		.add_function("Retetrahedralize", &Retetrahedralize, grp, "", "", TOOLTIP_RETETRAHEDRALIZE)
		.add_function("Duplicate", &Duplicate, grp, "", "", TOOLTIP_DUPLICATE)
		.add_function("Extrude", &Extrude, grp, "", "", TOOLTIP_EXTRUDE)
		.add_function("ExtrudeCylinders", &ExtrudeCylinders, grp, "", "", TOOLTIP_EXTRUDE_CYLINDERS)
		.add_function("CreateShrinkGeometry", &CreateShrinkGeometry, grp, "", "", TOOLTIP_CREATE_SHRINK_GEOMETRY)
	    .add_function("ExtrudeFacesWithTets", &ExtrudeFacesWithTets, grp, "", "", TOOLTIP_EXTRUDE_FACES_WITH_TETS);

//	selection tools
	reg.add_function("ClearSelection", &ClearSelection, grp, "", "", TOOLTIP_CLEAR_SELECTION)
		.add_function("SelectAll", &SelectAll, grp, "", "", TOOLTIP_SELECT_ALL) 
		.add_function("ExtendSelection", &ExtendSelection, grp, "", "", TOOLTIP_EXTEND_SELECTION)
		.add_function("SelectSubset", &SelectSubset, grp, "", "", TOOLTIP_SELECT_SUBSET)
		.add_function("SelectSubsetBoundary", &SelectSubsetBoundary, grp, "", "", TOOLTIP_SELECT_SUBSET_BOUNDARY)
		.add_function("SelectUnassignedElements", &SelectUnassignedElements, grp, "", "", TOOLTIP_SELECT_UNASSIGNED_ELEMENTS)
		.add_function("InvertSelection", &InvertSelection, grp, "", "", TOOLTIP_INVERT_SELECTION)
		.add_function("SelectSelectionBoundary", &SelectSelectionBoundary, grp, "", "", TOOLTIP_SELECT_SELECTION_BOUNDARY)
		.add_function("CloseSelection", &CloseSelection, grp, "", "", TOOLTIP_CLOSE_SELECTION)

		.add_function("SelectVertexByCoordinate", &SelectElemByCoordinate<Vertex>, grp,TOOLTIP_SELECT_VERTEX_BY_COORDINATE)
		.add_function("SelectVertexByCylindricalCoordinate", &SelectElemByCylindricalCoordinate<Vertex>, grp,TOOLTIP_SELECT_VERTEX_BY_CYL_COORDINATE)
		.add_function("SelectBoundaryVertices", &SelectBoundaryVertices, grp, "", "", TOOLTIP_SELECT_BOUNDARY_VERTICES)
		.add_function("SelectInnerVertices", &SelectInnerVertices, grp,TOOLTIP_SELECT_INNER_VERTICES)
		.add_function("SelectAssociatedVertices", &SelectAssociatedVertices, grp, "", "", TOOLTIP_SELECT_ASSOCIATED_VERTICES)
		.add_function("SelectAllVertices", &SelectAllVertices, grp, "", "", TOOLTIP_SELECT_ALL_VERTICES)
		.add_function("DeselectAllVertices", &DeselectAllVertices, grp, "", "", TOOLTIP_DESELECT_ALL_VERTICES)
		.add_function("SelectMarkedVertices", &SelectMarkedVertices, grp, "", "", TOOLTIP_SELECT_MARKED_VERTICES)
		.add_function("SelectVertexByIndex", &SelectVertexByIndex, grp, "", "", TOOLTIP_SELECT_VERTEX_BY_INDEX)
		.add_function("SelectUnconnectedVertices", &SelectUnconnectedVertices, grp, "", "", TOOLTIP_SELECT_UNCONNECTED_VERTICES)

		.add_function("SelectEdgeByCoordinate", &SelectElemByCoordinate<Edge>, grp, "", "", TOOLTIP_SELECT_EDGE_BY_COORDINATE)
		.add_function("SelectEdgeByCylindricalCoordinate", &SelectElemByCylindricalCoordinate<Edge>, grp, "", "", TOOLTIP_SELECT_EDGE_BY_CYL_COORDINATE)
		.add_function("SelectBoundaryEdges", &SelectBoundaryEdges, grp, "", "", TOOLTIP_SELECT_BOUNDARY_EDGES)
		.add_function("SelectInnerEdges", &SelectInnerEdges, grp, "", "", TOOLTIP_SELECT_INNER_EDGES)
		.add_function("SelectNonManifoldEdges", &SelectNonManifoldEdges, grp, "", "", TOOLTIP_SELECT_NON_MANIFOLD_EDGES)
		.add_function("SelectSmoothEdgePath", &SelectSmoothEdgePath, grp, "", "", TOOLTIP_SELECT_SMOOTH_EDGE_PATH)
		.add_function("SelectShortEdges", &SelectShortEdges, grp, "", "", TOOLTIP_SELECT_SHORT_EDGES)
		.add_function("SelectLongEdges", &SelectLongEdges, grp, "", "", TOOLTIP_SELECT_LONG_EDGES)
		.add_function("SelectCreaseEdges", &SelectCreaseEdges, grp, "", "", TOOLTIP_SELECT_CREASE_EDGES)
		.add_function("SelectLinkedBoundaryEdges", &SelectLinkedBoundaryEdges, grp, "", "", TOOLTIP_SELECT_LINKED_BOUNDARY_EDGES)
		.add_function("SelectAssociatedEdges", &SelectAssociatedEdges, grp, "", "", TOOLTIP_SELECT_ASSOCIATED_EDGES)
		.add_function("SelectAllEdges", &SelectAllEdges, grp, "", "", TOOLTIP_SELECT_ALL_EDGES)
		.add_function("DeselectAllEdges", &DeselectAllEdges, grp, "", "", TOOLTIP_DESELECT_ALL_EDGES)
		.add_function("SelectMarkedEdges", &SelectMarkedEdges, grp, "", "", TOOLTIP_SELECT_MARKED_EDGES)
		.add_function("SelectEdgeByIndex", &SelectEdgeByIndex, grp, "", "", TOOLTIP_SELECT_EDGE_BY_INDEX)
		.add_function("EdgeSelectionFill", &EdgeSelectionFill, grp, "", "", TOOLTIP_EDGE_SELECTION_FILL)

		.add_function("SelectFaceByCoordinate", &SelectElemByCoordinate<Face>, grp, "", "", TOOLTIP_SELECT_FACE_BY_COORDINATE)
		.add_function("SelectFaceByCylindricalCoordinate", &SelectElemByCylindricalCoordinate<Face>, grp, "", "", TOOLTIP_SELECT_FACE_BY_CYL_COORDINATE)
		.add_function("SelectBoundaryFaces", &SelectBoundaryFaces, grp, "", "", TOOLTIP_SELECT_BOUNDARY_FACES)
		.add_function("SelectInnerFaces", &SelectInnerFaces, grp, "", "", TOOLTIP_SELECT_INNER_FACES)
		.add_function("SelectLinkedManifoldFaces", &SelectLinkedManifoldFaces, grp, "", "", TOOLTIP_SELECT_LINKED_MANIFOLD_FACES)
		.add_function("SelectLinkedBoundaryFaces", &SelectLinkedBoundaryFaces, grp, "", "", TOOLTIP_SELECT_LINKED_BOUNDARY_FACES)
		.add_function("SelectDegenerateFaces", &SelectDegenerateFaces, grp, "", "", TOOLTIP_SELECT_DEGENERATE_FACES)
		.add_function("SelectLinkedFlatFaces", &SelectLinkedFlatFaces, grp, "", "", TOOLTIP_SELECT_LINKED_FLAT_FACES)
		.add_function("SelectIntersectingTriangles", &SelectIntersectingTriangles, grp, "", "", TOOLTIP_SELECT_INTERSECTING_TRIANGLES)
		.add_function("SelectAssociatedFaces", &SelectAssociatedFaces, grp, "", "", TOOLTIP_SELECT_ASSOCIATED_FACES)
		.add_function("SelectAllFaces", &SelectAllFaces, grp, "", "", TOOLTIP_SELECT_ALL_FACES)
		.add_function("DeselectAllFaces", &DeselectAllFaces, grp, "", "", TOOLTIP_DESELECT_ALL_FACES)
		.add_function("SelectFaceByIndex", &SelectFaceByIndex, grp, "", "", TOOLTIP_SELECT_FACE_BY_INDEX)
		.add_function("SelectFacesByNormal", &SelectFacesByNormal, grp, "", "", TOOLTIP_SELECT_FACES_BY_NORMAL)
		.add_function("FaceSelectionFill", &FaceSelectionFill, grp, "", "", TOOLTIP_FACE_SELECTION_FILL)
		.add_function("SelectBentQuadrilaterals", &SelectBentQuadrilaterals, grp, "", "", TOOLTIP_SELECT_BENT_QUADRILATERALS)

		.add_function("SelectVolumeByCoordinate", &SelectElemByCoordinate<Volume>, grp, "", "", TOOLTIP_SELECT_VOLUME_BY_COORDINATE)
		.add_function("SelectVolumeByCylindricalCoordinate", &SelectElemByCylindricalCoordinate<Volume>, grp, "", "", TOOLTIP_SELECT_VOLUME_BY_CYL_COORDINATE)
		.add_function("SelectAllVolumes", &SelectAllVolumes, grp, "", "", TOOLTIP_SELECT_ALL_VOLUMES)
		.add_function("DeselectAllVolumes", &DeselectAllVolumes, grp, "", "", TOOLTIP_DESELECT_ALL_VOLUMES)
		.add_function("SelectUnorientableVolumes", &SelectUnorientableVolumes, grp, "", "", TOOLTIP_SELECT_UNORIENTABLE_VOLUMES)
		.add_function("SelectVolumeByIndex", &SelectVolumeByIndex, grp, "", "", TOOLTIP_SELECT_VOLUME_BY_INDEX)
		.add_function("VolumeSelectionFill", &VolumeSelectionFill, grp, "", "", TOOLTIP_VOLUME_SELECTION_FILL)

		.add_function("MarkSelection", &MarkSelection, grp, "", "", TOOLTIP_MARK_SELECTION)
		.add_function("UnmarkSelection", &UnmarkSelection, grp, "", "", TOOLTIP_UNMARK_SELECTION);

//	subset tools
	reg.add_function("AssignSubset", &AssignSubset, grp, "", "", TOOLTIP_ASSIGN_SUBSET)
		.add_function("SetSubsetName", &SetSubsetName, grp, "", "", TOOLTIP_SET_SUBSET_NAME)
		.add_function("AssignSubsetColors", &AssignSubsetColors, grp, "", "", TOOLTIP_ASSIGN_SUBSET_COLORS)
		.add_function("SeparateFacesByEdgeSubsets", &SeparateFacesByEdgeSubsets, grp, "", "", TOOLTIP_SEPARATE_FACES_BY_EDGE_SUBSETS)
		.add_function("SeparateFacesBySelectedEdges", &SeparateFacesBySelectedEdges, grp, "", "", TOOLTIP_SEPARATE_FACES_BY_SELECTED_EDGES)
		.add_function("SeparateVolumesByFaceSubsets", &SeparateVolumesByFaceSubsets, grp, "", "", TOOLTIP_SEPARATE_VOLUMES_BY_FACE_SUBSETS)
		.add_function("SeparateVolumesBySelectedFaces", &SeparateVolumesBySelectedFaces, grp, "", "", TOOLTIP_SEPARATE_VOLUMES_BY_SELECTED_FACES)
		.add_function("SeparateIrregularManifoldSubsets", &SeparateIrregularManifoldSubsets, grp, "", "", TOOLTIP_SEPARATE_IRREGULAR_MANIFOLD_SUBSETS)
		.add_function("MoveSubset", &MoveSubset, grp, "", "", TOOLTIP_MOVE_SUBSET)
		.add_function("SwapSubsets", &SwapSubsets, grp, "", "", TOOLTIP_SWAP_SUBSETS)
		.add_function("JoinSubsets", &JoinSubsets, grp, "", "", TOOLTIP_JOIN_SUBSETS)
		.add_function("EraseSubset", &EraseSubset, grp, "", "", TOOLTIP_ERASE_SUBSET)
		.add_function("EraseEmptySubsets", &EraseEmptySubsets, grp, "", "", TOOLTIP_ERASE_EMPTY_SUBSETS)
		.add_function("AdjustSubsetsForUG3", &AdjustSubsetsForUG3, grp, "", "", TOOLTIP_ADJUST_SUBSETS_FOR_UG3)
		.add_function("AdjustSubsetsForUG4", &AdjustSubsetsForUG4, grp, "", "", TOOLTIP_ADJUST_SUBSETS_FOR_UG4)
		.add_function("SeparateFaceSubsetsByNormal", &SeparateFaceSubsetsByNormal, grp,TOOLTIP_SEPARATE_FACE_SUBSETS_BY_NORMAL )
		.add_function("SeparateFaceSubsetByNormal", &SeparateFaceSubsetByNormal, grp, "", "", TOOLTIP_SEPARATE_FACE_SUBSET_BY_NORMAL)
		.add_function("AssignSubsetsByQuality", &AssignSubsetsByQuality, grp, "", "", TOOLTIP_ASSIGN_SUBSETS_BY_QUALITY)
		.add_function("SeparateDegeneratedBoundaryFaceSubsets", &SeparateDegeneratedBoundaryFaceSubsets, grp, "", "", TOOLTIP_SEPARATE_DEGENERATED_BOUNDARY_FACE_SUBSETS)
		.add_function("AssignSubsetsByElementType", &AssignSubsetsByElementType, grp, "", "", TOOLTIP_ASSIGN_SUBSETS_BY_ELEMENT_TYPE);

//	topology
	reg.add_function("EraseSelectedElements", &EraseSelectedElements, grp, "", "", TOOLTIP_ERASE_SELECTED_ELEMENTS)
		.add_function("RemoveDoubles", &RemoveDoubles, grp, "", "", TOOLTIP_REMOVE_DOUBLES)
		.add_function("RemoveDoubleEdges", &RemoveDoubleEdges, grp, "", "", TOOLTIP_REMOVE_DOUBLE_EDGES)
		.add_function("MergeAtFirst", &MergeAtFirst, grp, "", "", TOOLTIP_MERGE_AT_FIRST)
		.add_function("MergeAtCenter", &MergeAtCenter, grp, "", "", TOOLTIP_MERGE_AT_CENTER)
		.add_function("MergeAtLast", &MergeAtLast, grp, "", "", TOOLTIP_MERGE_AT_LAST)
		.add_function("CollapseEdge", &CollapseEdge, grp, "", "", TOOLTIP_COLLAPSE_EDGE)
		.add_function("SplitEdge", &SplitEdge, grp, "", "", TOOLTIP_SPLIT_EDGE)
		.add_function("SwapEdge", &SwapEdge, grp, "", "", TOOLTIP_SWAP_EDGE)
		.add_function("PlaneCut", &PlaneCut, grp, "", "", TOOLTIP_PLANE_CUT)
		.add_function("AdjustEdgeOrientation", &AdjustEdgeOrientation, grp, "", "", TOOLTIP_ADJUST_EDGE_ORIENTATION)
		.add_function("FixFaceOrientation", &FixFaceOrientation, grp, "", "", TOOLTIP_FIX_FACE_ORIENTATION)
		.add_function("FixFaceSubsetOrientations", &FixFaceSubsetOrientations, grp, "", "", TOOLTIP_FIX_FACE_SUBSET_ORIENTATIONS)
		.add_function("FixVolumeOrientation", &FixVolumeOrientation, grp, "", "", TOOLTIP_FIX_VOLUME_ORIENTATION)
		.add_function("InvertFaceOrientation", &InvertFaceOrientation, grp, "", "", TOOLTIP_INVERT_FACE_ORIENTATION)
		.add_function("ResolveSelfIntersections", &ResolveSelfIntersections, grp, "", "", TOOLTIP_RESOLVE_SELF_INTERSECTIONS)
		.add_function("ResolveEdgeIntersection", &ResolveEdgeIntersection, grp, "", "", TOOLTIP_RESOLVE_EDGE_INTERSECTIONS) 
		.add_function("ResolveTriangleIntersections", &ResolveTriangleIntersections, grp, "", "", TOOLTIP_RESOLVE_TRIANGLE_INTERSECTIONS)
		.add_function("ProjectVerticesToCloseEdges", &ProjectVerticesToCloseEdges, grp, "", "", TOOLTIP_PROJECT_VERTICES_TO_CLOSE_EDGES)
		.add_function("ProjectVerticesToCloseFaces", &ProjectVerticesToCloseFaces, grp, "", "", TOOLTIP_PROJECT_VERTICES_TO_CLOSE_FACES)
		.add_function("IntersectCloseEdges", &IntersectCloseEdges, grp, "", "", TOOLTIP_INTERSECT_CLOSE_EDGES);

//	info tools
		grp = baseGrp + string("/Info/Measure length, area, volume");
		reg.add_function("MeasureGridLength", &MeasureGridLength, grp, "length", "", TOOLTIP_MEASURE_GRID_LENGTH)
		   .add_function("MeasureGridArea", &MeasureGridArea, grp, "area", "", TOOLTIP_MEASURE_GRID_AREA)
		   .add_function("MeasureGridVolume", &MeasureGridVolume, grp, "volume", "", TOOLTIP_MEASURE_GRID_VOLUME)
		   .add_function("MeasureSubsetLength", &MeasureSubsetLength, grp, "length", "mesh#subset", TOOLTIP_MEASURE_SUBSET_LENGTH)
		   .add_function("MeasureSubsetArea", &MeasureSubsetArea, grp, "area", "mesh#subset", TOOLTIP_MEASURE_SUBSET_AREA)
		   .add_function("MeasureSubsetVolume", &MeasureSubsetVolume, grp, "volume", "mesh#subset", TOOLTIP_MEASURE_SUBSET_VOLUME)
		   .add_function("MeasureSelectionLength", &MeasureSelectionLength, grp, "length", "", TOOLTIP_MEASURE_SELECTION_LENGTH)
		   .add_function("MeasureSelectionArea", &MeasureSelectionArea, grp, "area", "", TOOLTIP_MEASURE_SELECTION_AREA)
		   .add_function("MeasureSelectionVolume", &MeasureSelectionVolume, grp, "volume", "", TOOLTIP_MEASURE_SELECTION_VOLUME);

//	new tools
	reg.add_class_<Box>("Box", grp)
		.add_method("set_min", &Box::set_min)
		.add_method("set_max", &Box::set_max)
		.add_method("min", &Box::get_min)
		.add_method("max", &Box::get_max);

	reg.add_function("GetBoundingBox", &GetBoundingBox, grp, "", "", TOOLTIP_GET_BOUNDING_BOX)
		.add_function("SelectVerticesInBox", &SelectElementsInBox<Vertex>, grp, "", "", TOOLTIP_SELECT_VERTEX_IN_BOX)
		.add_function("SelectEdgesInBox", &SelectElementsInBox<Edge>, grp, "", "", TOOLTIP_SELECT_EDGE_IN_BOX)  
		.add_function("SelectFacesInBox", &SelectElementsInBox<Face>, grp, "", "", TOOLTIP_SELECT_FACE_IN_BOX)
		.add_function("SelectVolumesInBox", &SelectElementsInBox<Volume>, grp, "", "", TOOLTIP_SELECT_VOLUME_IN_BOX)
		.add_function("SelectVerticesInCylinder", &SelectElementsInCylinder<Vertex>, grp, "", "", TOOLTIP_SELECT_VERTEX_IN_CYLINDER)
		.add_function("SelectEdgesInCylinder", &SelectElementsInCylinder<Edge>, grp, "", "", TOOLTIP_SELECT_EDGE_IN_CYLINDER)
		.add_function("SelectFacesInCylinder", &SelectElementsInCylinder<Face>, grp, "", "", TOOLTIP_SELECT_FACE_IN_CYLINDER)
		.add_function("SelectVolumesInCylinder", &SelectElementsInCylinder<Volume>, grp, "", "", TOOLTIP_SELECT_VOLUME_IN_CYLINDER);

	reg.add_function("ScaleAroundPoint", &ScaleAroundPoint, grp, "", "", TOOLTIP_SCALE_AROUND_PIVOT);
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

