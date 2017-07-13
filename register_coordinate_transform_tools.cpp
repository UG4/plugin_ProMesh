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

#include "registration_routines.h"
#include "registry/registry.h"
#include "bridge/util.h"
#include "tooltips.h"
#include "tools/coordinate_transform_tools.h"
#include "tools/new_tools.h"

using namespace std;
using namespace ug::bridge;

namespace ug{
namespace promesh{

void RegisterCoordinateTransformTools(ProMeshRegistry& reg, string baseGrp)
{
	try{
		string grp = baseGrp + string("/Coordinate Transform/Pivot");
		reg.add_function("SetPivot", &SetPivot, grp, "",
				"mesh # position", TOOLTIP_SET_PIVOT) 
			.add_function("SetPivotToSelectionCenter", &SetPivotToSelectionCenter, grp, "",
				"mesh", TOOLTIP_SET_PIVOT_TO_SELECTION_CENTER)
			.add_function("SetPivotToMeshCenter", &SetPivotToMeshCenter, grp, "",
				"mesh", TOOLTIP_SET_PIVOT_TO_MESH_CENTER);

		grp = baseGrp + string("/Coordinate Transform/Move");
		reg.add_function("Move", &Move, grp, "",
				"mesh # offset", TOOLTIP_MOVE)
			.add_function("MoveMeshTo", &MoveMeshTo, grp, "",
				"mesh # new position", TOOLTIP_MOVE_MESH_TO) 
			.add_function("MoveSelectionTo", &MoveSelectionTo, grp, "",
				"mesh # center", TOOLTIP_MOVE_SELECTION_TO)
			.add_function("MoveAlongNormal", &MoveAlongNormal, grp, "",
				"mesh # offset | default | value=0.1D", TOOLTIP_MOVE_ALONG_NORMAL)
			.add_function("MoveVerticesAlongEdges",
				&MoveVerticesAlongEdges,
				grp, "", "mesh # relative offset | default | value=0.5D",
				TOOLTIP_MOVE_VERTICES_ALONG_EDGES)
			.add_function("MoveVerticesToEdgeLength",
				&MoveVerticesToEdgeLength,
				grp, "", "mesh # edge length | default | value=1.D",
				TOOLTIP_MOVE_VERTICES_TO_EDGE_LENGTH);
			

		grp = baseGrp + string("/Coordinate Transform/Scale");
		reg.add_function("ScaleAroundCenter", &ScaleAroundCenter, grp, "",
				"mesh # scale | default | value=[1,1,1]", TOOLTIP_SCALE_AROUND_CENTER) 
			.add_function("ScaleAroundPivot", &ScaleAroundPivot, grp, "",
				"mesh # scale | default | value=[1,1,1]", TOOLTIP_SCALE_AROUND_PIVOT)
			.add_function("ScaleAroundPoint", &ScaleAroundPoint, grp, "",
				"mesh # scale | default | value=[1,1,1] # point", TOOLTIP_SCALE_AROUND_POINT);

		grp = baseGrp + string("/Coordinate Transform/Rotate");
		reg.add_function("RotateAroundCenter", &RotateAroundCenter, grp, "",
				"mesh # rotation", TOOLTIP_ROTATE_AROUND_CENTER) 
			.add_function("RotateAroundPivot", &RotateAroundPivot, grp, "",
				"mesh # rotation", TOOLTIP_ROTATE_AROUND_PIVOT);

		grp = baseGrp + string("/Coordinate Transform");
		reg.add_function("Mirror", &Mirror, grp, "",
				"mesh # axis || value=[1, 0, 0] # origin", TOOLTIP_MIRROR);

		grp = baseGrp + string("/Coordinate Transform/Smoothing");
		reg.add_function("LaplacianSmooth", &LaplacianSmooth, grp, "",
				"mesh # alpha | default | value=0.25D # num iterations | default | value=10",
				TOOLTIP_LAPLACIAN_SMOOTH)  
			// .add_function("WeightedEdgeSmooth", &WeightedEdgeSmooth, grp, "",
			// 			  "mesh#"
			// 			  "alpha | default | min=0D; max=1D; value=0.25D#"
			// 			  "num iterations | default | min=0; value=10",
			// 			  TOOLTIP_WEIGHTED_EDGE_SMOOTH)  
			// .add_function("WeightedFaceSmooth", &WeightedFaceSmooth, grp, "",
			// 			  "mesh#"
			// 			  "alpha | default | min=0D; max=1D; value=0.25D#"
			// 			  "num iterations | default | min=0; value=10",
			// 			  TOOLTIP_WEIGHTED_FACE_SMOOTH)  
			// .add_function("WeightedNormalSmooth", &WeightedNormalSmooth, grp, "",
			// 			  "mesh#"
			// 			  "alpha | default | min=0D; max=1D; value=0.25D#"
			// 			  "dot threshold | default | min=-1D; max=1D; value=-1D#"
			// 			  "num iterations | default | min=0; value=10",
			// 			  TOOLTIP_WEIGHTED_NORMAL_SMOOTH)
			.add_function("TangentialSmooth", &TangentialSmooth, grp, "",
				"mesh # alpha | default | value=0.25D # num iterations | default | value=10",
				TOOLTIP_TANGENTIAL_SMOOTH)
			.add_function("SlopeSmooth", &SlopeSmooth, grp, "",
				"mesh#"
				"alpha | default | min=0D; max=1D; value=0.25D#"
				"num iterations | default | min=0; value=10",
				TOOLTIP_SLOPE_SMOOTH);


		grp = baseGrp + string("/Coordinate Transform");
		reg.add_function("GetSelectionCenter", &GetSelectionCenter, grp, "",
				"mesh # centerOut", TOOLTIP_GET_SELECTION_CENTER, "", RT_NO_PROMESH)
			.add_function("ProjectToPlane", &ProjectToPlane, grp, "",
				"mesh # plane point # plane normal || value=[0,0,1]",
				TOOLTIP_RPOJECT_TO_PLANE)
			.add_function("ConeTransform", &ConeTransform, grp, "",
				"mesh # base # axis || value=[0,0,1] # scaleAtTip | default | value=1", TOOLTIP_CONE_TRANSFORM)
			.add_function("FlattenBentQuadrilaterals", &FlattenBentQuadrilaterals, grp, "",
				"mesh # step size | default | value=0.01D # num iterations | default | value=100",
				TOOLTIP_FLATTEN_BENT_QUADRILATERALS)
			.add_function("SnapVerticesToVertices", &SnapVerticesToVertices, grp, "",
				"mesh # targetMesh", TOOLTIP_SNAP_VERTICES_TO_VERTICES, "", RT_NO_PROMESH);

		grp = baseGrp + string("/Coordinate Transform/Subdivision Projection");
		reg.add_function("ProjectToLimitPiecewiseLoop", &ProjectToLimitPLoop, grp, "",
				"mesh", TOOLTIP_PROJECT_TO_LIMIT_PLOOP)  
			.add_function("ProjectToLimitSmoothBoundary", &ProjectToLimitSmoothBoundary, grp, "",
				"mesh", TOOLTIP_PROJECT_TO_LIMIT_SMOOTH_BOUNDARY);
	}
	UG_REGISTRY_CATCH_THROW(baseGrp);
}

} // end namespace promesh
}// end of namespace

