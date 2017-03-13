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

#include "coordinate_transform_tools.h"
#include "lib_grid/algorithms/orientation_util.h"

using namespace std;

namespace ug{
namespace promesh{

bool GetSelectionCenter(Mesh* obj, vector3& centerOut)
{
	return CalculateCenter(centerOut, obj->selector(), obj->position_accessor());
}


bool MoveSelectionTo(Mesh* obj, const vector3& center)
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


void Move(Mesh* obj, const vector3& offset)
{
	TranslateSelection(obj->selector(), offset, obj->position_accessor());
}

void MoveMeshTo(Mesh* obj, const vector3& newPos)
{
	vector3 offset;
	VecSubtract(offset, newPos, obj->pivot());
	MoveVertices(obj->grid().begin<Vertex>(), obj->grid().end<Vertex>(),
				 obj->position_accessor(), offset);
	obj->set_pivot(newPos);
}

void MoveAlongNormal(Mesh* obj, number offset)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	Mesh::position_accessor_t& aaPos = obj->position_accessor();

	Grid::face_traits::secure_container faces;
	vector<Vertex*> vrts;
	CollectVerticesTouchingSelection(vrts, sel);

	for(vector<Vertex*>::iterator iter = vrts.begin(); iter != vrts.end(); ++iter)
	{
		vector3& p = aaPos[*iter];
	//	calculate vertex normal by averaging face normals
		vector3 n(0, 0, 0);
		grid.associated_elements(faces, *iter);
		for(size_t i = 0; i < faces.size(); ++i){
			vector3 fn;
			CalculateNormal(fn, faces[i], aaPos);
			VecAdd(n, n, fn);
		}

		VecNormalize(n, n);

		VecScaleAdd(p, 1., p, offset, n);
	}
}

void MoveVerticesAlongEdges(Mesh* obj, number relVal)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	Mesh::position_accessor_t& aaPos = obj->position_accessor();
	Grid::edge_traits::secure_container edges;

//	calculate offsets (do not apply, so that we correctly support cases in which
//	both corners of an edge are selected)
	vector<vector3>	offsets;
	offsets.reserve(sel.num<Vertex>());
	lg_for_each(Vertex, vrt, sel){
		vector3 offset(0, 0, 0);
		number numSelEdges = 0;
		grid.associated_elements(edges, vrt);
		for_each_in_vec(Edge* e, edges){
			if(!sel.is_selected(e))
				continue;
			vector3 dir;
			if(e->vertex(0) == vrt)
				VecSubtract(dir, aaPos[e->vertex(1)], aaPos[e->vertex(0)]);
			else
				VecSubtract(dir, aaPos[e->vertex(0)], aaPos[e->vertex(1)]);

			VecScaleAppend(offset, relVal, dir);
			++numSelEdges;
		}end_for;

		if(numSelEdges > 0)
			offset *= 1. / numSelEdges;
		offsets.push_back(offset);
	}lg_end_for;

//	apply offsets
	size_t index = 0;
	lg_for_each(Vertex, vrt, sel){
		VecAppend(aaPos[vrt], offsets[index++]);
	}lg_end_for;
}

void ScaleAroundCenter(Mesh* obj, const vector3& scale)
{
	Mesh::position_accessor_t& aaPos = obj->position_accessor();

	vector<Vertex*> vrts;
	CollectVerticesTouchingSelection(vrts, obj->selector());
	vector3 center = CalculateCenter(vrts.begin(), vrts.end(), aaPos);

	for(vector<Vertex*>::iterator iter = vrts.begin(); iter != vrts.end(); ++iter)
	{
		ug::vector3& v = aaPos[*iter];
		VecSubtract(v, v, center);
		v.x() *= scale.x();
		v.y() *= scale.y();
		v.z() *= scale.z();
		VecAdd(v, v, center);
	}
}

void ScaleAroundPivot(Mesh* obj, const vector3& scale)
{
	Mesh::position_accessor_t& aaPos = obj->position_accessor();

	vector<Vertex*> vrts;
	CollectVerticesTouchingSelection(vrts, obj->selector());
	vector3 center = obj->pivot();

	for(vector<Vertex*>::iterator iter = vrts.begin(); iter != vrts.end(); ++iter)
	{
		vector3& v = aaPos[*iter];
		VecSubtract(v, v, center);
		v.x() *= scale.x();
		v.y() *= scale.y();
		v.z() *= scale.z();
		VecAdd(v, v, center);
	}
}

void ScaleAroundPoint(Mesh* obj, const vector3& scale, const vector3& point)
{
	Mesh::position_accessor_t& aaPos = obj->position_accessor();

	vector<Vertex*> vrts;
	CollectVerticesTouchingSelection(vrts, obj->selector());

	for(vector<Vertex*>::iterator iter = vrts.begin(); iter != vrts.end(); ++iter)
	{
		vector3& v = aaPos[*iter];
		VecSubtract(v, v, point);
		v.x() *= scale.x();
		v.y() *= scale.y();
		v.z() *= scale.z();
		VecAdd(v, v, point);
	}
}

/**	call this method for edges, faces, and volumes after you performed mirroring.
 * The method mirrors all elements of a grid which are connected to at least one
 * of the given vertices.
 * \WARNING	The method may produce invalid results if only a subset of vertices of
 * 			a given element was mirrored.
 * \NOTE	The method uses Grid::mark*/
template <class elem_t, class vrt_iter_t>
static void FixMirrorOrientation(Grid& g, vrt_iter_t vrtsBegin, vrt_iter_t vrtsEnd)
{
	g.begin_marking();
	typename Grid::traits<elem_t>::secure_container elems;
	for(vrt_iter_t ivrt = vrtsBegin; ivrt != vrtsEnd; ++ivrt){
		Vertex* vrt = *ivrt;
		g.associated_elements(elems, vrt);

		for(size_t i = 0; i < elems.size(); ++i){
			elem_t* e = elems[i];
			if(!g.is_marked(e)){
				g.mark(e);
				g.flip_orientation(e);
			}
		}
	}
	g.end_marking();
}

void Mirror(Mesh* obj, const vector3& axisOrig, const vector3& origin)
{
	UG_COND_THROW(VecLengthSq(axisOrig) == 0, "Invalid axis specified in 'Mirror'. "
			"Please choose an axis which differs from (0, 0, 0).");

	vector3 axis;
	VecNormalize(axis, axisOrig);

	Mesh::position_accessor_t& aaPos = obj->position_accessor();
	Grid& g = obj->grid();

	vector<Vertex*> vrts;
	CollectVerticesTouchingSelection(vrts, obj->selector());

	for(vector<Vertex*>::iterator ivrt = vrts.begin();
		ivrt != vrts.end(); ++ivrt)
	{
		Vertex* vrt = *ivrt;
		vector3 p;
		ProjectPointToPlane(p, aaPos[vrt], origin, axis);
		VecSubtract(p, p, aaPos[vrt]);
		p *= 2;
		VecAdd(aaPos[vrt], aaPos[vrt], p);
	}

	FixMirrorOrientation<Edge>(g, vrts.begin(), vrts.end());
	FixMirrorOrientation<Face>(g, vrts.begin(), vrts.end());
	FixMirrorOrientation<Volume>(g, vrts.begin(), vrts.end());
}

void RotateAroundCenter(Mesh* obj, const vector3& rotRad)
{
	const vector3 rot(	deg_to_rad(rotRad.x()),
						deg_to_rad(rotRad.y()),
						deg_to_rad(rotRad.z()));

	Mesh::position_accessor_t& aaPos = obj->position_accessor();

	vector<Vertex*> vrts;
	CollectVerticesTouchingSelection(vrts, obj->selector());
	vector3 center = CalculateCenter(vrts.begin(), vrts.end(), aaPos);

//	todo: combine rotation and transform matrix and use ugs build in methods.
//	rotation matrix
	matrix33 matRot;
	MatRotationYawPitchRoll(matRot, rot.x(), rot.y(), rot.z());

	for(vector<Vertex*>::iterator iter = vrts.begin(); iter != vrts.end(); ++iter)
	{
		vector3 v = aaPos[*iter];
		VecSubtract(v, v, center);
		MatVecMult(aaPos[*iter], matRot, v);
		VecAdd(aaPos[*iter], aaPos[*iter], center);
	}
}


void RotateAroundPivot(Mesh* obj, const vector3& rotRad)
{
	const vector3 rot(	deg_to_rad(rotRad.x()),
						deg_to_rad(rotRad.y()),
						deg_to_rad(rotRad.z()));

	Mesh::position_accessor_t& aaPos = obj->position_accessor();

	vector<Vertex*> vrts;
	CollectVerticesTouchingSelection(vrts, obj->selector());
	vector3 center = obj->pivot();

//	todo: combine rotation and transform matrix and use ugs build in methods.
//	rotation matrix
	matrix33 matRot;
	MatRotationYawPitchRoll(matRot, rot.x(), rot.y(), rot.z());

	for(vector<Vertex*>::iterator iter = vrts.begin(); iter != vrts.end(); ++iter)
	{
		vector3 v = aaPos[*iter];
		VecSubtract(v, v, center);
		MatVecMult(aaPos[*iter], matRot, v);
		VecAdd(aaPos[*iter], aaPos[*iter], center);
	}
}


void ConeTransform(Mesh* obj, const vector3& base, const vector3& axis,
					number scaleAtTip)
{
	Mesh::position_accessor_t& aaPos = obj->position_accessor();

	vector<Vertex*> vrts;
	CollectVerticesTouchingSelection(vrts, obj->selector());

	for(vector<Vertex*>::iterator iter = vrts.begin(); iter != vrts.end(); ++iter)
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


void LaplacianSmooth(Mesh* obj, number alpha, int numIterations)
{
	vector<Vertex*> vrts;
	CollectVerticesTouchingSelection(vrts, obj->selector());

	ug::LaplacianSmooth(obj->grid(), vrts.begin(), vrts.end(),
						obj->position_accessor(), alpha, numIterations);
}

void WeightedEdgeSmooth(Mesh* obj, number alpha, int numIterations)
{
	Selector& sel = obj->selector();
	ug::WeightedEdgeSmooth(obj->grid(), sel.begin<Vertex>(), sel.end<Vertex>(),
					   obj->position_accessor(), alpha, numIterations, IsSelected(sel));
}

void WeightedFaceSmooth(Mesh* obj, number alpha, int numIterations)
{
	Selector& sel = obj->selector();
	ug::WeightedFaceSmooth(obj->grid(), sel.begin<Vertex>(), sel.end<Vertex>(),
					   obj->position_accessor(), alpha, numIterations, IsSelected(sel));
}

void WeightedNormalSmooth(Mesh* obj, number alpha, number dotThreshold,
								 int numIterations)
{
	vector<Vertex*> vrts;
	CollectVerticesTouchingSelection(vrts, obj->selector());

	ug::WeightedNormalSmooth(obj->grid(), vrts.begin(), vrts.end(),
							 obj->position_accessor(), alpha, dotThreshold,
							 numIterations);
}

void SlopeSmooth(Mesh* obj, number alpha, int numIterations)
{
	vector<Vertex*> vrts;
	CollectVerticesTouchingSelection(vrts, obj->selector());

	ug::SlopeSmooth(obj->grid(), vrts.begin(), vrts.end(),
					   obj->position_accessor(), alpha, vector3(0, 0, 1),
					   numIterations);
}

void TangentialSmooth(Mesh* obj, number alpha, int numIterations)
{
	vector<Vertex*> vrts;
	CollectVerticesTouchingSelection(vrts, obj->selector());

	ug::TangentialSmooth(obj->grid(), vrts.begin(), vrts.end(),
						 obj->position_accessor(), alpha, numIterations);
}

void ProjectToLimitPLoop(Mesh* obj)
{
	ProjectToLimitPLoop(obj->grid(), obj->position_attachment(),
						obj->position_attachment());
}


void ProjectToLimitSmoothBoundary(Mesh* obj)
{
	ProjectToLimitSubdivBoundary(obj->grid(), obj->position_attachment(),
								 obj->position_attachment());
}


void SetPivot(Mesh* obj, const vector3& pos)
{
	obj->set_pivot(pos);
}


void SetPivotToSelectionCenter(Mesh* obj)
{
	vector3 center;
	CalculateCenter(center, obj->selector(), obj->position_accessor());
	obj->set_pivot(center);
}

void SetPivotToMeshCenter(Mesh* obj)
{
	obj->set_pivot(CalculateCenter(obj->grid().begin<Vertex>(),
								   obj->grid().end<Vertex>(),
								   obj->position_accessor()));
}

void FlattenBentQuadrilaterals(Mesh* obj, number stepSize, int numIterations)
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


void ProjectToPlane(Mesh* obj, const vector3& planeCenter, const vector3& planeNormal)
{
	Mesh::position_accessor_t& aaPos = obj->position_accessor();

	vector<Vertex*> vrts;
	CollectVerticesTouchingSelection(vrts, obj->selector());

	for(vector<Vertex*>::iterator iter = vrts.begin(); iter != vrts.end(); ++iter)
	{
		vector3 v = aaPos[*iter];
		ProjectPointToPlane(aaPos[*iter], v, planeCenter, planeNormal);
	}
}


void SnapVerticesToVertices(Mesh* obj, Mesh* targetMesh)
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
	vector<Vertex*> closeVrts;

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

}}//	end of namespace
