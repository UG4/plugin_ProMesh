/*
 * Copyright (c) 2015:  G-CSC, Goethe University Frankfurt
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

void RegisterMeshingTools(ProMeshRegistry& reg, string baseGrp)
{
	try{
		string grp;
		grp = baseGrp + "/Grid Generation";
		reg.add_function("CloneMesh", &CloneMesh, grp, "",
				"mesh", TOOLTIP_CLONE_MESH, "", RT_NO_PROMESH)
			.add_function("CopySelection", &CopySelection, grp, "",
				"srcMesh # targetMesh", TOOLTIP_COPY_SELECTION, "", RT_NO_PROMESH);

		grp = baseGrp + "/Grid Generation/Basic Elements";
		reg.add_function("CreateVertex", &CreateVertex, grp, "",
				"mesh # position # subset index", TOOLTIP_CREATE_VERTEX)
			.add_function("CreateEdge", &CreateEdge, grp, "",
				"mesh # subset index", TOOLTIP_CREATE_EDGE)
			.add_function("CreateFace", &CreateFace, grp, "",
				"mesh # subset index", TOOLTIP_CREATE_FACE)
			.add_function("CreateVolume", &CreateVolume, grp, "",
				"mesh # subset index", TOOLTIP_CREATE_VOLUME);

		grp = baseGrp + "/Grid Generation/Geometries";
		reg.add_function("CreatePlane", &CreatePlane, grp, "",
				"mesh # upLeft # upRight # lowLeft # lowRight # subset index # fill",
				TOOLTIP_CREATE_PLANE)
			.add_function("CreateCircle", &CreateCircle, grp, "",
				"mesh # center # radius || value=1 # num rim vertices || value=12 # subset index # fill",
				TOOLTIP_CREATE_CIRCLE)
			.add_function("CreateBox", &CreateBox, grp, "",
				"mesh # boxMin # boxMax # subset index # create volume", TOOLTIP_CREATE_BOX)
			.add_function("CreateSphere", &CreateSphere, grp, "",
				"mesh # center # radius # num refinements # subset index", TOOLTIP_CREATE_SPHERE)
			.add_function("CreateTetrahedron", &CreateTetrahedron, grp, "",
				"mesh # subset index # create volume", TOOLTIP_CREATE_TETRAHEDRON)
			.add_function("CreatePyramid", &CreatePyramid, grp, "",
				"mesh # subset index # create volume", TOOLTIP_CREATE_PYRAMID)
			.add_function("CreatePrism", &CreatePrism, grp, "",
				"mesh # subset index # create volume", TOOLTIP_CREATE_PRISM);

	//	layer meshing
		grp = baseGrp + "/Raster Layers";
		reg.add_function("MeshLayers", &MeshLayers, grp, "",
				"mesh # layers", TOOLTIP_MESH_LAYERS, "", RT_NO_PROMESH)
			.add_function("MeshLayerBoundaries", &MeshLayerBoundaries, grp, "",
				"mesh # layers", TOOLTIP_MESH_LAYER_BOUNDARIES, "", RT_NO_PROMESH)
			.add_function("ExtrudeLayers", &ExtrudeLayers, grp, "",
				"mesh #"
				"layers #"
				"allow for tets and pyras || value = true",
				TOOLTIP_EXTRUDE_LAYERS, "", RT_NO_PROMESH)
			.add_function("ExtrudeLayersAndAddProjector", &ExtrudeLayersAndAddProjector, grp, "",
				"mesh #"
				"layers #"
				"allow for tets and pyras || value = true",
				TOOLTIP_EXTRUDE_LAYERS_AND_ADD_PROJECTOR, "", RT_NO_PROMESH)
			.add_function("ProjectToLayer", &ProjectToLayer, grp, "",
				"mesh #"
				"layers #"
				"layer index",
				TOOLTIP_PROJECT_TO_LAYER, "", RT_NO_PROMESH)
			.add_function("ProjectToTopLayer", &ProjectToTopLayer, grp, "",
				"mesh #"
				"layers",
				TOOLTIP_PROJECT_TO_TOP_LAYER, "", RT_NO_PROMESH);
			
	//	refinement
		grp = baseGrp + "/Remeshing/Refinement";
		reg.add_function("Refine", &Refine, grp, "",
				"mesh # strict subset inheritance # use snap points", TOOLTIP_REFINE)
			.add_function("HangingNodeRefine", &HangingNodeRefine, grp, "",
				"mesh # strict subset inheritance # anisotropic", TOOLTIP_HANGING_NODE_REFINE)
			.add_function("RefineSmooth", &RefineSmooth, grp, "",
				"mesh # strict subset inheritance", TOOLTIP_REFINE_SMOOTH)
			.add_function("InsertCenter", &InsertCenter, grp, "",
				"", TOOLTIP_INSERT_CENTER);

	//	remeshing
		grp = baseGrp + "/Remeshing/Polylines";
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
			
		grp = baseGrp + "/Remeshing/Triangulation";
		reg.add_function("ConvertToTriangles", &ConvertToTriangles, grp, "",
				"mesh", TOOLTIP_CONVERT_TO_TRIANGLES)
			.add_function("TriangleFill", &TriangleFill, grp, "",
				"mesh # quality generation # min angle # subset index", TOOLTIP_TRIANGLE_FILL)
			.add_function("Retriangulate", &Retriangulate, grp, "",
				"mesh # min angle", TOOLTIP_RETRIANGULATE)
			.add_function("AdjustEdgeLength", &AdjustEdgeLength, grp, "",
				"mesh # "
				"min edge length | default | min=0.0D; value=1.D # "
				"max edge length | default | min=0.0D; value=2.D # "
				"num iterations | default | min=1; value=10 # "
				"adaptive | default | value=true # "
				"automark boundaries | default | value=true",
				TOOLTIP_ADJUST_EDGE_LENGTH)
			.add_function("AdjustEdgeLengthExtended", &AdjustEdgeLengthExtended, grp, "",
				"mesh # "
				"min edge length | default | min=0.0D; value=1.D # "
				"max edge length | default | min=0.0D; value=2.D # "
				"approximation | default | min=0.0D; max=1.0D; value=0.9D # "
				"triangle quality | default | min=0.0D; max=1.0D; value=0.9D # "
				"num iterations | default | min=1; value=10 # "
				"automark boundaries | default |value=true",
				TOOLTIP_ADJUST_EDGE_LENGTH)
			.add_function("AdaptSurfaceToCylinder", &AdaptSurfaceToCylinder, grp, "",
				"mesh # radius # threshold", TOOLTIP_ADAPT_SURFACE_TO_CYLINDER)
			.add_function("ReplaceValence3Vertices", &ReplaceValence3Vertices, grp, "",
				"mesh # max relative height | default | min=0; value = 0.0001D", TOOLTIP_REPLACE_VALENCE_3_VERTICES)
			.add_function("ReplaceLowValenceVertices", &ReplaceLowValenceVertices, grp, "",
				"mesh # max relative height | default | min=0; value = 0.0001D", TOOLTIP_REPLACE_LOW_VALENCE_VERTICES);

		grp = baseGrp + "/Remeshing/Quadrilaterals";
		reg.add_function("ConvertToQuadrilaterals", &ConvertToQuadrilaterals, grp, "",
				"mesh", TOOLTIP_CONVERT_TO_QUADRILATERALS);

		grp = baseGrp + "/Remeshing/Tetrahedra";
		reg.add_function("ConvertToTetrahedra",
				static_cast<void (*)(Mesh*)>(&ConvertToTetrahedra), grp, "",
				"mesh", TOOLTIP_CONVERT_TO_TETRAHEDRA)
			.add_function("Tetrahedralize", &Tetrahedralize, grp, "",
				"mesh # quality # preserve outer # preserve all # "
				"separate volumes # append subsets at end # verbosity", TOOLTIP_TETRAHEDRALIZE)
			.add_function("AssignVolumeConstraints", &AssignVolumeConstraints, grp, "",
				"mesh # volume constraint", TOOLTIP_ASSIGN_VOLUME_CONSTRAINTS)
			.add_function("ClearVolumeConstraints", &ClearVolumeConstraints, grp, "",
				"mesh", TOOLTIP_CLEAR_VOLUME_CONSTRAINTS)
			.add_function("Retetrahedralize", &Retetrahedralize, grp, "",
				"mesh # quality # preserve outer # preserve all # "
				"apply volume constraint # verbosity", TOOLTIP_RETETRAHEDRALIZE)
		    .add_function("ExtrudeFacesWithTets", &ExtrudeFacesWithTets, grp, "",
		    	"mesh # from subset # to subset # factor", TOOLTIP_EXTRUDE_FACES_WITH_TETS);


		grp = baseGrp + "/Remeshing/Extrusion";
		reg.add_function("Extrude", &ExtrudeAndMove, grp, "",
				"mesh #"
				"total direction | default | value=[0,0,1] #"
				"num steps | default | min=1;value=1 #"
				"create faces | default | value=true #"
				"create volumes | default | value=true",
				TOOLTIP_EXTRUDE_AND_MOVE)
			.add_function("ExtrudeAndMove", &ExtrudeAndMove, grp, "",
				"mesh #"
				"total direction | default | value=[0,0,1] #"
				"num steps | default | min=1;value=1 #"
				"create faces | default | value=true #"
				"create volumes | default | value=true",
				TOOLTIP_EXTRUDE_AND_MOVE)
			.add_function("ExtrudeAndScale", &ExtrudeAndScale, grp, "",
				"mesh #"
				"total scale | default | value=2D #"
				"scale around pivot #"
				"num steps | default | min=1;value=1 #"
				"create faces | default | value=true #"
				"create volumes | default | value=true",
				TOOLTIP_EXTRUDE_AND_SCALE)
			.add_function("ExtrudeAlongNormal", &ExtrudeAlongNormal, grp, "",
				"mesh #"
				"total length | default | value=1D #"
				"num steps | default | min=1;value=1 #"
				"create faces | default | value=true #"
				"create volumes | default | value=true",
				TOOLTIP_EXTRUDE_ALONG_NORMAL)
			.add_function("ExtrudeCylinders", &ExtrudeCylinders, grp, "",
				"mesh #"
				"height | default | value=1D #"
				"radius | default | value=1D #"
				"snap threshold | default | min=0D;value=0.001D ",
				TOOLTIP_EXTRUDE_CYLINDERS);

	//	topology
		grp = baseGrp + "/Remeshing";
		reg.add_function("EraseSelectedElements", &EraseSelectedElements, grp, "",
				"mesh # erase unused vertices # erase unused edges # erase unused faces",
				TOOLTIP_ERASE_SELECTED_ELEMENTS)
			.add_function("Duplicate", &Duplicate, grp, "",
				"mesh # offset # deselect old # select new", TOOLTIP_DUPLICATE)
			.add_function("PlaneCut", &PlaneCut, grp, "",
				"mesh # plane center # plane normal", TOOLTIP_PLANE_CUT)
			.add_function("CreateShrinkGeometry", &CreateShrinkGeometry, grp, "",
				"mesh # scale", TOOLTIP_CREATE_SHRINK_GEOMETRY);

		grp = baseGrp + "/Remeshing/Edge Operations";
		reg.add_function("CollapseEdge", &CollapseEdge, grp, "",
				"mesh", TOOLTIP_COLLAPSE_EDGE)
			.add_function("SplitEdge", &SplitEdge, grp, "",
				"mesh", TOOLTIP_SPLIT_EDGE)
			.add_function("SwapEdge", &SwapEdge, grp, "",
				"mesh", TOOLTIP_SWAP_EDGE);

		grp = baseGrp + "/Remeshing/Orientation";
		reg.add_function("AdjustEdgeOrientation", &AdjustEdgeOrientation, grp, "",
				"mesh", TOOLTIP_ADJUST_EDGE_ORIENTATION)
			.add_function("FixFaceOrientation", &FixFaceOrientation, grp, "",
				"mesh", TOOLTIP_FIX_FACE_ORIENTATION)
			.add_function("FixFaceSubsetOrientations", &FixFaceSubsetOrientations, grp, "",
				"mesh", TOOLTIP_FIX_FACE_SUBSET_ORIENTATIONS)
			.add_function("FixVolumeOrientation", &FixVolumeOrientation, grp, "",
				"mesh", TOOLTIP_FIX_VOLUME_ORIENTATION)
			.add_function("InvertFaceOrientation", &InvertFaceOrientation, grp, "",
				"mesh", TOOLTIP_INVERT_FACE_ORIENTATION);

		grp = baseGrp + "/Remeshing/Resolve Intersections";
		reg.add_function("ResolveSelfIntersections", &ResolveSelfIntersections, grp, "",
			"mesh # snap threshold", TOOLTIP_RESOLVE_SELF_INTERSECTIONS);

		grp = baseGrp + "/Remeshing/Resolve Intersections/Advanced";
		reg.add_function("ProjectVerticesToCloseEdges", &ProjectVerticesToCloseEdges, grp, "",
				"mesh # snap threshold", TOOLTIP_PROJECT_VERTICES_TO_CLOSE_EDGES)
			.add_function("ProjectVerticesToCloseFaces", &ProjectVerticesToCloseFaces, grp, "",
				"mesh # snap threshold", TOOLTIP_PROJECT_VERTICES_TO_CLOSE_FACES)
			.add_function("IntersectCloseEdges", &IntersectCloseEdges, grp, "",
				"mesh # snap threshold", TOOLTIP_INTERSECT_CLOSE_EDGES)
			.add_function("ResolveEdgeIntersection", &ResolveEdgeIntersection, grp, "",
				"mesh # snap threshold", TOOLTIP_RESOLVE_EDGE_INTERSECTIONS) 
			.add_function("ResolveTriangleIntersections", &ResolveTriangleIntersections, grp, "",
				"mesh # snap threshold", TOOLTIP_RESOLVE_TRIANGLE_INTERSECTIONS);
		
		grp = baseGrp + "/Remeshing/Remove Doubles";
		reg.add_function("RemoveDoubles", &RemoveDoubles, grp, "",
				"mesh # threshold", TOOLTIP_REMOVE_DOUBLES)
			.add_function("RemoveDoubleEdges", &RemoveDoubleEdges, grp, "",
				"mesh", TOOLTIP_REMOVE_DOUBLE_EDGES)
			.add_function("RemoveDoubleFaces", &RemoveDoubleFaces, grp, "",
				"mesh", TOOLTIP_REMOVE_DOUBLE_FACES);

		grp = baseGrp + "/Remeshing/Merge Vertices";
		reg.add_function("MergeAtFirst", &MergeAtFirst, grp, "",
				"mesh", TOOLTIP_MERGE_AT_FIRST)
			.add_function("MergeAtCenter", &MergeAtCenter, grp, "",
				"mesh", TOOLTIP_MERGE_AT_CENTER)
			.add_function("MergeAtLast", &MergeAtLast, grp, "",
				"mesh", TOOLTIP_MERGE_AT_LAST);

		grp = baseGrp + "/Remeshing/Boolean Operations";
		reg.add_function("FaceUnion", &CSGFaceUnion, grp, "",
				"mesh#"
				 "subset 1  | default | value=0#"
				 "subset 2  | default | value=1#"
				 "snap threshold | default | value=0.000001D;min=0;step=0.000001D", TOOLTIP_CSG_FACE_UNION)
			.add_function("FaceIntersection", &CSGFaceIntersection, grp, "",
				"mesh#"
				 "subset 1  | default | value=0#"
				 "subset 2  | default | value=1#"
				 "snap threshold | default | value=0.000001D;min=0;step=0.000001D", TOOLTIP_CSG_FACE_INTERSECTION)
			.add_function("FaceDifference", &CSGFaceDifference, grp, "",
				 "mesh#"
				 "subset 1  | default | value=0#"
				 "subset 2  | default | value=1#"
				 "snap threshold | default | value=0.000001D;min=0;step=0.000001D", TOOLTIP_CSG_FACE_DIFFERENCE);
	}
	UG_REGISTRY_CATCH_THROW(baseGrp);
}

} // end namespace promesh
}// end of namespace

