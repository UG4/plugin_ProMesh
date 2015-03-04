// created by Sebastian Reiter
// s.b.reiter@gmail.com

#include "registration_routines.h"
#include "registry/registry.h"
#include "bridge/util.h"
#include "tooltips.h"
#include "tools/coordinate_transform_tools.h"

using namespace std;
using namespace ug::bridge;

namespace ug{
namespace promesh{

void RegisterCoordinateTransformTools(Registry& reg, string baseGrp)
{
	try{
		string grp = baseGrp + string("/Coordinate Transform");
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
			.add_function("TangentialSmooth", &TangentialSmooth, grp, "", "", TOOLTIP_TANGENTIAL_SMOOTH) 
			.add_function("ProjectToLimitPLoop", &ProjectToLimitPLoop, grp, "", "", TOOLTIP_PROJECT_TO_LIMIT_PLOOP)  
			.add_function("ProjectToLimitSmoothBoundary", &ProjectToLimitSmoothBoundary, grp, "", "", TOOLTIP_PROJECT_TO_LIMIT_SMOOTH_BOUNDARY) 
			.add_function("SetPivot", &SetPivot, grp, "", "", TOOLTIP_SET_PIVOT) 
			.add_function("SetPivotToSelectionCenter", &SetPivotToSelectionCenter, grp, "", "", TOOLTIP_SET_PIVOT_TO_SELECTION_CENTER)
			.add_function("SetPivotToMeshCenter", &SetPivotToMeshCenter, grp, "", "", TOOLTIP_SET_PIVOT_TO_MESH_CENTER)
			.add_function("FlattenBentQuadrilaterals", &FlattenBentQuadrilaterals, grp, "", "", TOOLTIP_FLATTEN_BENT_QUADRILATERALS)
			.add_function("SnapVerticesToVertices", &SnapVerticesToVertices, grp, "", "mesh # targetMesh", TOOLTIP_SNAP_VERTICES_TO_VERTICES);
	}
	UG_REGISTRY_CATCH_THROW(baseGrp);
}

} // end namespace promesh
}// end of namespace

