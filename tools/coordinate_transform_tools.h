/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
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
#define TOOLTIP_ROTATE_AROUND_CENTER "Rotates the selected elements around the center of the current selection."
#define TOOLTIP_ROTATE_AROUND_PIVOT "Rotates the selected elements around the pivot of the mesh."
#define	TOOLTIP_COORDINATES "Coordinates of the center of the current selection"
#define	TOOLTIP_MOVE "Moves selected vertices."
#define TOOLTIP_MOVE_MESH_TO "Moves the active mesh and its pivot, so that the pivot will be located on the specified position."
#define	TOOLTIP_NORMAL_MOVE "Moves selected vertices along their normal."
#define	TOOLTIP_SCALE "Scales the coordinates of the selected vertices around their center."
#define	TOOLTIP_ROTATE "Rotates the geometry by the given degrees around its center."
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
inline bool GetSelectionCenter(Mesh* obj, vector3& centerOut)
{
	return CalculateCenter(centerOut, obj->selector(), obj->position_accessor());
}


inline bool MoveSelectionTo(Mesh* obj, const vector3& center)
{
	vector3 oldCenter;
	if(GetSelectionCenter(obj, oldCenter)){
		vector3 d;
		VecSubtract(d, center, oldCenter);
		TranslateSelection(obj->selector(), d, obj->position_accessor());
		return true;
	}
	return false;
}


inline void Move(Mesh* obj, const vector3& offset)
{
	TranslateSelection(obj->selector(), offset, obj->position_accessor());
}

inline void MoveMeshTo(Mesh* obj, const vector3& newPos)
{
	vector3 offset;
	VecSubtract(offset, newPos, obj->pivot());
	MoveVertices(obj->grid().begin<Vertex>(), obj->grid().end<Vertex>(),
				 obj->position_accessor(), offset);
	obj->set_pivot(newPos);
}

inline void MoveAlongNormal(Mesh* obj, number offset, bool usePrecalculatedNormals)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	Mesh::position_accessor_t& aaPos = obj->position_accessor();
	Mesh::normal_accessor_t& aaNorm = obj->normal_accessor();

	Grid::face_traits::secure_container faces;
	std::vector<Vertex*> vrts;
	CollectVerticesTouchingSelection(vrts, sel);

	for(std::vector<Vertex*>::iterator iter = vrts.begin(); iter != vrts.end(); ++iter)
	{
		vector3& p = aaPos[*iter];
	//	calculate vertex normal by averaging face normals
		vector3 n(0, 0, 0);
		grid.associated_elements(faces, *iter);
		if(usePrecalculatedNormals){
			for(size_t i = 0; i < faces.size(); ++i)
				VecAdd(n, n, aaNorm[faces[i]]);
		}
		else{
			for(size_t i = 0; i < faces.size(); ++i){
				vector3 fn;
				CalculateNormal(fn, faces[i], aaPos);
				VecAdd(n, n, fn);
			}
		}

		VecNormalize(n, n);

		VecScaleAdd(p, 1., p, offset, n);
	}
}

inline void MoveAlongNormal(Mesh* obj, number offset)
{
	MoveAlongNormal(obj, offset, true);
}

inline void ScaleAroundCenter(Mesh* obj, const vector3& scale)
{
	Mesh::position_accessor_t& aaPos = obj->position_accessor();

	std::vector<Vertex*> vrts;
	CollectVerticesTouchingSelection(vrts, obj->selector());
	vector3 center = CalculateCenter(vrts.begin(), vrts.end(), aaPos);

	for(std::vector<Vertex*>::iterator iter = vrts.begin(); iter != vrts.end(); ++iter)
	{
		ug::vector3& v = aaPos[*iter];
		VecSubtract(v, v, center);
		v.x() *= scale.x();
		v.y() *= scale.y();
		v.z() *= scale.z();
		VecAdd(v, v, center);
	}
}

inline void ScaleAroundPivot(Mesh* obj, const vector3& scale)
{
	Mesh::position_accessor_t& aaPos = obj->position_accessor();

	std::vector<Vertex*> vrts;
	CollectVerticesTouchingSelection(vrts, obj->selector());
	vector3 center = obj->pivot();

	for(std::vector<Vertex*>::iterator iter = vrts.begin(); iter != vrts.end(); ++iter)
	{
		vector3& v = aaPos[*iter];
		VecSubtract(v, v, center);
		v.x() *= scale.x();
		v.y() *= scale.y();
		v.z() *= scale.z();
		VecAdd(v, v, center);
	}
}


inline void RotateAroundCenter(Mesh* obj, const vector3& rot)
{
	Mesh::position_accessor_t& aaPos = obj->position_accessor();

	std::vector<Vertex*> vrts;
	CollectVerticesTouchingSelection(vrts, obj->selector());
	vector3 center = CalculateCenter(vrts.begin(), vrts.end(), aaPos);

//	todo: combine rotation and transform matrix and use ugs build in methods.
//	rotation matrix
	matrix33 matRot;
	MatRotationYawPitchRoll(matRot, rot.x(), rot.y(), rot.z());

	for(std::vector<Vertex*>::iterator iter = vrts.begin(); iter != vrts.end(); ++iter)
	{
		vector3 v = aaPos[*iter];
		VecSubtract(v, v, center);
		MatVecMult(aaPos[*iter], matRot, v);
		VecAdd(aaPos[*iter], aaPos[*iter], center);
	}
}


inline void RotateAroundPivot(Mesh* obj, const vector3& rot)
{
	Mesh::position_accessor_t& aaPos = obj->position_accessor();

	std::vector<Vertex*> vrts;
	CollectVerticesTouchingSelection(vrts, obj->selector());
	vector3 center = obj->pivot();

//	todo: combine rotation and transform matrix and use ugs build in methods.
//	rotation matrix
	matrix33 matRot;
	MatRotationYawPitchRoll(matRot, rot.x(), rot.y(), rot.z());

	for(std::vector<Vertex*>::iterator iter = vrts.begin(); iter != vrts.end(); ++iter)
	{
		vector3 v = aaPos[*iter];
		VecSubtract(v, v, center);
		MatVecMult(aaPos[*iter], matRot, v);
		VecAdd(aaPos[*iter], aaPos[*iter], center);
	}
}


inline void ConeTransform(Mesh* obj, const vector3& base, const vector3& axis,
					number scaleAtTip)
{
	Mesh::position_accessor_t& aaPos = obj->position_accessor();

	std::vector<Vertex*> vrts;
	CollectVerticesTouchingSelection(vrts, obj->selector());

	for(std::vector<Vertex*>::iterator iter = vrts.begin(); iter != vrts.end(); ++iter)
	{
		vector3 v = aaPos[*iter];
		vector3 proj(0, 0, 0);

	//	project the vertex onto the ray from src along dir
		number s = ProjectPointToRay(proj, v, base, axis);

	//	s determines the scaling amount:
	//	s = 0 => scale with 1
	//	s = 1 => scale with scaleAtTip
	//	else interpolate linear
		number scaling = (1. - s) + s * scaleAtTip;

		vector3 vtmp;
		VecSubtract(vtmp, v, proj);
		VecScale(vtmp, vtmp, scaling);
		VecAdd(aaPos[*iter], proj, vtmp);
	}
}


inline void LaplacianSmooth(Mesh* obj, number alpha, int numIterations)
{
	std::vector<Vertex*> vrts;
	CollectVerticesTouchingSelection(vrts, obj->selector());

	ug::LaplacianSmooth(obj->grid(), vrts.begin(), vrts.end(),
						obj->position_accessor(), alpha, numIterations);
}

inline void WeightedEdgeSmooth(Mesh* obj, number alpha, int numIterations)
{
	Selector& sel = obj->selector();
	ug::WeightedEdgeSmooth(obj->grid(), sel.begin<Vertex>(), sel.end<Vertex>(),
					   obj->position_accessor(), alpha, numIterations, IsSelected(sel));
}

inline void WeightedFaceSmooth(Mesh* obj, number alpha, int numIterations)
{
	Selector& sel = obj->selector();
	ug::WeightedFaceSmooth(obj->grid(), sel.begin<Vertex>(), sel.end<Vertex>(),
					   obj->position_accessor(), alpha, numIterations, IsSelected(sel));
}

inline void WeightedNormalSmooth(Mesh* obj, number alpha, number dotThreshold,
								 int numIterations)
{
	std::vector<Vertex*> vrts;
	CollectVerticesTouchingSelection(vrts, obj->selector());

	ug::WeightedNormalSmooth(obj->grid(), vrts.begin(), vrts.end(),
							 obj->position_accessor(), alpha, dotThreshold,
							 numIterations);
}

inline void SlopeSmooth(Mesh* obj, number alpha, int numIterations)
{
	std::vector<Vertex*> vrts;
	CollectVerticesTouchingSelection(vrts, obj->selector());

	ug::SlopeSmooth(obj->grid(), vrts.begin(), vrts.end(),
					   obj->position_accessor(), alpha, vector3(0, 0, 1),
					   numIterations);
}

inline void TangentialSmooth(Mesh* obj, number alpha, int numIterations)
{
	std::vector<Vertex*> vrts;
	CollectVerticesTouchingSelection(vrts, obj->selector());

	ug::TangentialSmooth(obj->grid(), vrts.begin(), vrts.end(),
						 obj->position_accessor(), alpha, numIterations);
}

inline void ProjectToLimitPLoop(Mesh* obj)
{
	ProjectToLimitPLoop(obj->grid(), obj->position_attachment(),
						obj->position_attachment());
}


inline void ProjectToLimitSmoothBoundary(Mesh* obj)
{
	ProjectToLimitSubdivBoundary(obj->grid(), obj->position_attachment(),
								 obj->position_attachment());
}


inline void SetPivot(Mesh* obj, const vector3& pos)
{
	obj->set_pivot(pos);
}


inline void SetPivotToSelectionCenter(Mesh* obj)
{
	vector3 center;
	CalculateCenter(center, obj->selector(), obj->position_accessor());
	obj->set_pivot(center);
}

inline void SetPivotToMeshCenter(Mesh* obj)
{
	obj->set_pivot(CalculateCenter(obj->grid().begin<Vertex>(),
								   obj->grid().end<Vertex>(),
								   obj->position_accessor()));
}

inline void FlattenBentQuadrilaterals(Mesh* obj, number stepSize, int numIterations)
{
	Selector& sel = obj->selector();
	Mesh::position_accessor_t& aaPos = obj->position_accessor();

	for(int i = 0; i < numIterations; ++i){
	//	iterate over all quadrilaterals
		for(QuadrilateralIterator iter = sel.begin<Quadrilateral>();
			iter != sel.end<Quadrilateral>(); ++iter)
		{
			Quadrilateral* q = *iter;

		//	find the plane that minimizes the distance to the corner points
		//	of the quadrilateral.
			vector3 p, n;
			vector3 points[4] =	{aaPos[q->vertex(0)], aaPos[q->vertex(1)],
								 aaPos[q->vertex(2)], aaPos[q->vertex(3)]};

			if(FindClosestPlane(p, n, points, 4)){
			//	find the projections of the corners and move the corner
			//	vertices a little into that direction.
				for(size_t i = 0; i < 4; ++i){
					vector3 proj;
					ProjectPointToPlane(proj, points[i], p, n);
					VecScaleAdd(aaPos[q->vertex(i)], 1. - stepSize, points[i],
								stepSize, proj);
				}
			}
		}
	}
}


inline void ProjectToPlane(Mesh* obj, const vector3& planeCenter, const vector3& planeNormal)
{
	Mesh::position_accessor_t& aaPos = obj->position_accessor();

	std::vector<Vertex*> vrts;
	CollectVerticesTouchingSelection(vrts, obj->selector());

	for(std::vector<Vertex*>::iterator iter = vrts.begin(); iter != vrts.end(); ++iter)
	{
		vector3 v = aaPos[*iter];
		ProjectPointToPlane(aaPos[*iter], v, planeCenter, planeNormal);
	}
}


inline void SnapVerticesToVertices(Mesh* obj, Mesh* targetMesh)
{
	KDTreeStatic<APosition> tree;
	Grid& targetGrid = targetMesh->grid();
	Selector& targetSel = targetMesh->selector();

	tree.create_from_grid(targetGrid,
						  targetSel.begin<Vertex>(),
						  targetSel.end<Vertex>(),
						  targetMesh->position_attachment(),
						  20, 10);

	Mesh::position_accessor_t	aaPos = obj->position_accessor();
	Mesh::position_accessor_t	aaTargetPos = targetMesh->position_accessor();

	Selector& srcSel = obj->selector();
	std::vector<Vertex*> closeVrts;

	for(VertexIterator viter = srcSel.begin<Vertex>();
		viter != srcSel.end<Vertex>(); ++viter)
	{
		Vertex* vrt = *viter;
		tree.get_neighbourhood(closeVrts, aaPos[vrt], 1);
		if(!closeVrts.empty()){
			aaPos[vrt] = aaTargetPos[closeVrts[0]];
		}
	}
}

/// \}

}}// end of namespace

#endif
