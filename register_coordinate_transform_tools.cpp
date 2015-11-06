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

void RegisterCoordinateTransformTools(Registry& reg, string baseGrp)
{
	try{
		string grp = baseGrp + string("/Coordinate Transform");
		reg.add_function("GetSelectionCenter", &GetSelectionCenter, grp, "",
				"mesh # centerOut", TOOLTIP_GET_SELECTION_CENTER)
			.add_function("SetSelectionCenter", &SetSelectionCenter, grp, "",
				"mesh # center", TOOLTIP_SET_SELECTION_CENTER)
			.add_function("Move", &Move, grp, "",
				"mesh # offset", TOOLTIP_MOVE)
			.add_function("MoveMeshTo", &MoveMeshTo, grp, "",
				"mesh # new position", TOOLTIP_MOVE_MESH_TO) 
			.add_function("MoveAlongNormal", &MoveAlongNormal, grp, "",
				"mesh # offset # use precalculated normals", TOOLTIP_MOVE_ALONG_NORMAL)

			.add_function("ScaleAroundCenter", &ScaleAroundCenter, grp, "",
				"mesh # scale", TOOLTIP_SCALE_AROUND_CENTER) 
			.add_function("ScaleAroundPivot", &ScaleAroundPivot, grp, "",
				"mesh # scale", TOOLTIP_SCALE_AROUND_PIVOT)
			.add_function("ScaleAroundPoint", &ScaleAroundPoint, grp, "",
				"mesh # scale # point", TOOLTIP_SCALE_AROUND_POINT)

			.add_function("RotateAroundCenter", &RotateAroundCenter, grp, "",
				"mesh # rotation", TOOLTIP_ROTATE_AROUND_CENTER) 
			.add_function("RotateAroundPivot", &RotateAroundPivot, grp, "",
				"mesh # rotation", TOOLTIP_ROTATE_AROUND_PIVOT)  

			.add_function("ConeTransform", &ConeTransform, grp, "",
				"mesh # base # axis # scaleAtTip", TOOLTIP_CONE_TRANSFORM)  

			.add_function("LaplacianSmooth", &LaplacianSmooth, grp, "",
				"mesh # alpha # num iterations", TOOLTIP_LAPLACIAN_SMOOTH)  
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
			.add_function("WeightedNormalSmooth", &WeightedNormalSmooth, grp, "",
						  "mesh#"
						  "alpha | default | min=0D; max=1D; value=0.25D#"
						  "dot threshold | default | min=-1D; max=1D; value=-1D#"
						  "num iterations | default | min=0; value=10",
						  TOOLTIP_WEIGHTED_NORMAL_SMOOTH)
			.add_function("SlopeSmooth", &SlopeSmooth, grp, "",
						  "mesh#"
						  "alpha | default | min=0D; max=1D; value=0.25D#"
						  "num iterations | default | min=0; value=10",
						  TOOLTIP_SLOPE_SMOOTH)
			.add_function("TangentialSmooth", &TangentialSmooth, grp, "",
				"mesh # alpha # num iterations", TOOLTIP_TANGENTIAL_SMOOTH) 

			.add_function("ProjectToLimitPLoop", &ProjectToLimitPLoop, grp, "",
				"mesh", TOOLTIP_PROJECT_TO_LIMIT_PLOOP)  
			.add_function("ProjectToLimitSmoothBoundary", &ProjectToLimitSmoothBoundary, grp, "",
				"mesh", TOOLTIP_PROJECT_TO_LIMIT_SMOOTH_BOUNDARY) 

			.add_function("SetPivot", &SetPivot, grp, "",
				"mesh # position", TOOLTIP_SET_PIVOT) 
			.add_function("SetPivotToSelectionCenter", &SetPivotToSelectionCenter, grp, "",
				"mesh", TOOLTIP_SET_PIVOT_TO_SELECTION_CENTER)
			.add_function("SetPivotToMeshCenter", &SetPivotToMeshCenter, grp, "",
				"mesh", TOOLTIP_SET_PIVOT_TO_MESH_CENTER)

			.add_function("FlattenBentQuadrilaterals", &FlattenBentQuadrilaterals, grp, "",
				"mesh # stepSize # num iterations", TOOLTIP_FLATTEN_BENT_QUADRILATERALS)
			
			.add_function("SnapVerticesToVertices", &SnapVerticesToVertices, grp, "",
				"mesh # targetMesh", TOOLTIP_SNAP_VERTICES_TO_VERTICES);
	}
	UG_REGISTRY_CATCH_THROW(baseGrp);
}

} // end namespace promesh
}// end of namespace

