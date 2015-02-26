// created by Sebastian Reiter
// s.b.reiter@gmail.com

#include "lib_grid/algorithms/remove_duplicates_util.h"
#include "registration_routines.h"
#include "registry/registry.h"
#include "bridge/util.h"
#include "tooltips.h"
#include "tools/grid_generation_tools.h"
#include "tools/refinement_tools.h"
#include "tools/remeshing_tools.h"
#include "tools/topology_tools.h"

using namespace std;
using namespace ug::bridge;

namespace ug{
namespace promesh{

void RegisterMeshingTools(Registry& reg, string baseGrp)
{
	try{
		string grp = baseGrp;
		reg.add_function("CopySelection", &CopySelection, grp, "", "", TOOLTIP_COPY_SELECTION)
			.add_function("CreateVertex", &CreateVertex, grp, "", "", TOOLTIP_CREATE_VERTEX)
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

	//	layer meshing
		reg.add_function("MeshLayers", &MeshLayers, grp, "",
				"mesh # layers",
				TOOLTIP_MESH_LAYERS)
			.add_function("MeshLayerBoundaries", &MeshLayerBoundaries, grp, "",
				"mesh # layers",
				TOOLTIP_MESH_LAYER_BOUNDARIES);
		
	//	refinement
		reg.add_function("Refine", &Refine, grp, "", "", TOOLTIP_REFINE)
			.add_function("HangingNodeRefine", &HangingNodeRefine, grp, "", "", TOOLTIP_HANGING_NODE_REFINE)
			.add_function("RefineSmooth", &RefineSmooth, grp, "", "", TOOLTIP_REFINE_SMOOTH)
			.add_function("RefineSmoothBoundary2D", &RefineSmoothBoundary2D, grp, "", "", TOOLTIP_REFINE_SMOOTH_BOUNDARY_2D)
			.add_function("CreateFractal", &CreateFractal, grp, "", "", TOOLTIP_CREATE_FRACTAL)
			.add_function("InsertCenter", &InsertCenter, grp, "", "", TOOLTIP_INSERT_CENTER);

	//	remeshing
		grp = baseGrp + string("/Remeshing/Polylines");
		reg.add_function("SimplifyPolylines", &SimplifyPolylines, grp, "",
				"mesh # "
				"max curvature angle | default | min=0.0D; max=180.0D; value=5.D # ",
				TOOLTIP_SIMPLIFY_POLYLINES)
			.add_function("SimplifySmoothedPolylines", &SimplifySmoothedPolylines, grp, "",
				"mesh # "
				"max curvature angle | default | min=0.0D; max=180.0D; value=5.D # "
				"smoothing alpha | default | min=0.0D; value=0.9D # "
				"smoothing iterations | default | min=0; value=10 # ",
				TOOLTIP_SIMPLIFY_SMOOTHED_POLYLINES);
			
		grp = baseGrp + string("/Remeshing/Triangulation");
		reg.add_function("ConvertToTriangles", &ConvertToTriangles, grp, "", "", TOOLTIP_CONVERT_TO_TRIANGLES)
			.add_function("TriangleFill", &TriangleFill, grp, "", "", TOOLTIP_TRIANGLE_FILL)
			.add_function("Retriangulate", &Retriangulate, grp, "", "", TOOLTIP_RETRIANGULATE)
			.add_function("AdjustEdgeLength", &AdjustEdgeLength, grp, "", "", TOOLTIP_ADJUST_EDGE_LENGTH)
			.add_function("AdjustEdgeLengthExtended", &AdjustEdgeLengthExtended, grp, "",
				"mesh # "
				"min edge length | default | min=0.0D; value=1.D # "
				"max edge length | default | min=0.0D; value=2.D # "
				"approximation | default | min=0.0D; max=1.0D; value=0.9D # "
				"triangle quality | default | min=0.0D; max=1.0D; value=0.9D # "
				"numIterations | default | min=1; value=10 # "
				"automark boundaries | default |value=true",
				TOOLTIP_ADJUST_EDGE_LENGTH)
			.add_function("AdaptSurfaceToCylinder", &AdaptSurfaceToCylinder, grp, "", "", TOOLTIP_ADAPT_SURFACE_TO_CYLINDER)
			.add_function("ReplaceLowValenceVertices", &ReplaceLowValenceVertices, grp, "",
						  "mesh # max relative height | default | min=0; value = 0.0001D", TOOLTIP_REPLACE_LOW_VALENCE_VERTICES);

		grp = baseGrp + string("/Remeshing/Tetgen");
		reg.add_function("Tetrahedralize", &Tetrahedralize, grp, "", "", TOOLTIP_TETRAHEDRALIZE)
			.add_function("AssignVolumeConstraints", &AssignVolumeConstraints, grp, "", "", TOOLTIP_ASSIGN_VOLUME_CONSTRAINTS)
			.add_function("ClearVolumeConstraints", &ClearVolumeConstraints, grp, "", "", TOOLTIP_CLEAR_VOLUME_CONSTRAINTS)
			.add_function("Retetrahedralize", &Retetrahedralize, grp, "", "", TOOLTIP_RETETRAHEDRALIZE)
			.add_function("Duplicate", &Duplicate, grp, "", "", TOOLTIP_DUPLICATE)
			.add_function("Extrude", &Extrude, grp, "", "", TOOLTIP_EXTRUDE)
			.add_function("ExtrudeCylinders", &ExtrudeCylinders, grp, "", "", TOOLTIP_EXTRUDE_CYLINDERS)
			.add_function("CreateShrinkGeometry", &CreateShrinkGeometry, grp, "", "", TOOLTIP_CREATE_SHRINK_GEOMETRY)
		    .add_function("ExtrudeFacesWithTets", &ExtrudeFacesWithTets, grp, "", "", TOOLTIP_EXTRUDE_FACES_WITH_TETS);

	//	topology
		grp = baseGrp + string("/Remeshing");
		reg.add_function("EraseSelectedElements", &EraseSelectedElements, grp, "", "", TOOLTIP_ERASE_SELECTED_ELEMENTS)
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
		
		grp = baseGrp + string("/Remeshing/Remove Doubles");
		reg.add_function("RemoveDoubles", &RemoveDoubles, grp, "", "", TOOLTIP_REMOVE_DOUBLES)
			.add_function("RemoveDoubleEdges", &RemoveDoubleEdges, grp, "", "", TOOLTIP_REMOVE_DOUBLE_EDGES)
			.add_function("RemoveDoubleFaces", &RemoveDoubleFaces, grp, "", "", TOOLTIP_REMOVE_DOUBLE_FACES);

		grp = baseGrp + string("/Remeshing/Merge Vertices");
		reg.add_function("MergeAtFirst", &MergeAtFirst, grp, "", "", TOOLTIP_MERGE_AT_FIRST)
			.add_function("MergeAtCenter", &MergeAtCenter, grp, "", "", TOOLTIP_MERGE_AT_CENTER)
			.add_function("MergeAtLast", &MergeAtLast, grp, "", "", TOOLTIP_MERGE_AT_LAST);
	}
	UG_REGISTRY_CATCH_THROW(baseGrp);
}

} // end namespace promesh
}// end of namespace

