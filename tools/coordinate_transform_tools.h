/*
 * Copyright (c) 2013-2017:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__PROMESH__coordinate_transform_tools__
#define __H__PROMESH__coordinate_transform_tools__

#include <vector>
#include "../mesh.h"
#include "lib_grid/algorithms/subdivision/subdivision_loop.h"
#include "lib_grid/algorithms/selection_util.h"
#include "lib_grid/algorithms/smoothing/manifold_smoothing.h"
#include "lib_grid/algorithms/trees/kd_tree_static.h"
#include "lib_grid/callbacks/callbacks.h"

#define	TOOLTIP_GET_SELECTION_CENTER ""
#define TOOLTIP_MOVE_SELECTION_TO "Moves the selected elements so that the new selection center will lie at the specified position."
#define TOOLTIP_MOVE_ALONG_NORMAL "Moves selected vertices along their normal by the specified offset."
#define TOOLTIP_SCALE_AROUND_CENTER "Scales the selected elements around the center of the current selection."
#define TOOLTIP_SCALE_AROUND_PIVOT "Scales the selected elements around the pivot of the mesh."
#define TOOLTIP_SCALE_AROUND_POINT "Scales the selected geometry around the given point"
#define TOOLTIP_ROTATE_AROUND_CENTER "Rotates the selected elements around the center of the current selection."
#define TOOLTIP_ROTATE_AROUND_PIVOT "Rotates the selected elements around the pivot of the mesh."
#define	TOOLTIP_COORDINATES "Coordinates of the center of the current selection"
#define	TOOLTIP_MOVE "Moves selected vertices."
#define TOOLTIP_MOVE_MESH_TO "Moves the active mesh and its pivot, so that the pivot will be located on the specified position."
#define	TOOLTIP_NORMAL_MOVE "Moves selected vertices along their normal."
#define	TOOLTIP_MOVE_VERTICES_ALONG_EDGES "Moves selected vertices along selected edges by an offset relative to those selected edges. If a selected vertex is connected to multiple selected edges, the new position will be averaged between the individual offsets."
#define	TOOLTIP_MOVE_VERTICES_TO_EDGE_LENGTH "Moves selected vertices along selected edges so that those edges obtain the specified length. If a selected vertex is connected to multiple selected edges, the new position will be averaged between the individual new positions."
#define	TOOLTIP_SCALE "Scales the coordinates of the selected vertices around their center."
#define	TOOLTIP_ROTATE "Rotates the geometry by the given degrees around its center."
#define TOOLTIP_MIRROR "Mirrors the geometry at the given point along the given axis"
#define	TOOLTIP_TRANSFORM "Transforms the vertices with the given matrix"
#define	TOOLTIP_CONE_TRANSFORM "Transforms the vertices with the given cone transformation"
#define	TOOLTIP_LAPLACIAN_SMOOTH "Smoothes vertices in a grid."
#define	TOOLTIP_WEIGHTED_EDGE_SMOOTH "Smoothes vertices along edges in a grid with special weights for non-smoothed vertices."
#define	TOOLTIP_WEIGHTED_FACE_SMOOTH "Smoothes vertices across faces in a grid with special weights for non-smoothed vertices."
#define	TOOLTIP_WEIGHTED_NORMAL_SMOOTH "The higher the dot-product between an outgoing edge and the vertex normal, the higher the influence of that edge during smoothing of that vertex."
#define	TOOLTIP_SLOPE_SMOOTH "Smoothes the grid so that the geometry is linearized along the path of steepest descent."
#define	TOOLTIP_TANGENTIAL_SMOOTH "Smoothes vertices on a manifold."
#define	TOOLTIP_RPOJECT_TO_PLANE "Projects all selected elements to the specified plane"
#define	TOOLTIP_PROJECT_TO_LIMIT_PLOOP "Projects all vertices in the grid to their limit positions as defined by the piecewise loop scheme."
#define	TOOLTIP_PROJECT_TO_LIMIT_SMOOTH_BOUNDARY "Projects all boundary-vertices in the grid to their limit positions as defined by the b-spline subdivision scheme."
#define	TOOLTIP_SET_PIVOT "Sets the pivot point of the selected object."
#define	TOOLTIP_SET_PIVOT_TO_SELECTION_CENTER "Sets the pivot to the center of the current selection."
#define	TOOLTIP_SET_PIVOT_TO_MESH_CENTER "Sets the pivot to the center of the active mesh."
#define	TOOLTIP_FLATTEN_BENT_QUADRILATERALS "Flattens bent quadrilaterals using an iterative flattening scheme"
#define	TOOLTIP_APPLY_HEIGHT_FIELD "Calculates z-values of all nodes in terms of their x and y values." 
#define TOOLTIP_SNAP_VERTICES_TO_VERTICES "Snaps the selected vertices of the specified mesh to the selected vertices of the target mesh"

namespace ug{
namespace promesh{

/// \addtogroup promesh
/// \{
bool GetSelectionCenter(Mesh* obj, vector3& centerOut);

bool MoveSelectionTo(Mesh* obj, const vector3& center);

void Move(Mesh* obj, const vector3& offset);

void MoveMeshTo(Mesh* obj, const vector3& newPos);

void MoveAlongNormal(Mesh* obj, number offset);

void MoveVerticesAlongEdges(Mesh* obj, number relVal);

void MoveVerticesToEdgeLength(Mesh* obj, number edgeLen);

void ScaleAroundCenter(Mesh* obj, const vector3& scale);

void ScaleAroundPivot(Mesh* obj, const vector3& scale);

void ScaleAroundPoint(Mesh* obj, const vector3& scale, const vector3& point);

void RotateAroundCenter(Mesh* obj, const vector3& rot);

void RotateAroundPivot(Mesh* obj, const vector3& rot);

void Mirror(Mesh* obj, const vector3& axis, const vector3& origin);

void ConeTransform(	Mesh* obj,
					const vector3& base,
					const vector3& axis,
					number scaleAtTip);

void LaplacianSmooth(Mesh* obj, number alpha, int numIterations);

void WeightedEdgeSmooth(Mesh* obj, number alpha, int numIterations);

void WeightedFaceSmooth(Mesh* obj, number alpha, int numIterations);

void WeightedNormalSmooth(	Mesh* obj,
							number alpha,
							number dotThreshold,
							int numIterations);

void SlopeSmooth(Mesh* obj, number alpha, int numIterations);

void TangentialSmooth(Mesh* obj, number alpha, int numIterations);

void ProjectToLimitPLoop(Mesh* obj);

void ProjectToLimitSmoothBoundary(Mesh* obj);

void SetPivot(Mesh* obj, const vector3& pos);

void SetPivotToSelectionCenter(Mesh* obj);

void SetPivotToMeshCenter(Mesh* obj);

void FlattenBentQuadrilaterals(Mesh* obj, number stepSize, int numIterations);

void ProjectToPlane(	Mesh* obj,
						const vector3& planeCenter,
						const vector3& planeNormal);

void SnapVerticesToVertices(Mesh* obj, Mesh* targetMesh);

/// \}

}}// end of namespace

#endif
