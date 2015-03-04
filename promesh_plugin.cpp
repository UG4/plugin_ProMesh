// created by Sebastian Reiter
// s.b.reiter@gmail.com
// Jun 14, 2013 (d,m,y)

#include "mesh.h"
#include "tools/file_io_tools.h"
#include "tools/new_tools.h"
#include "tools/subset_tools.h"
#include "tools/measure_tools.h"
#include "bridge/util.h"
#include "tooltips.h"
#include "registration_routines.h"

using namespace std;
using namespace ug::bridge;

namespace ug{
namespace promesh{

static void RegisterMisc(Registry& reg, string baseGrp)
{
	string grp = baseGrp;

//	file io
	grp = baseGrp;
	reg.add_function("LoadMesh", &LoadMesh, grp, "", "", TOOLTIP_LOAD_MESH)
		.add_function("SaveMesh", &SaveMesh, grp,TOOLTIP_SAVE_MESH)
		.add_function("ExportToUG3", &ExportToUG3, grp, "", "", TOOLTIP_EXPORT_TO_UG3);

//	subset tools
	grp = baseGrp;
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

	reg.add_function("ScaleAroundPoint", &ScaleAroundPoint, grp, "", "", TOOLTIP_SCALE_AROUND_POINT);
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
		RegisterMesh(*reg, grp);
		RegisterCoordinateTransformTools(*reg, grp);
		RegisterSelectionTools(*reg, grp);
		RegisterMeshingTools(*reg, grp);
		RegisterMisc(*reg, grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

}// end of namespace

