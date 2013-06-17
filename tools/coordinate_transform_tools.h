// created by Sebastian Reiter
// s.b.reiter@gmail.com
// Jun 17, 2013 (d,m,y)

#ifndef __H__PROMESH__coordinate_transform_tools__
#define __H__PROMESH__coordinate_transform_tools__

#include <vector>
#include "../mesh_object.h"
#include "lib_grid/algorithms/subdivision/subdivision_loop.h"

namespace ug{
namespace promesh{

bool GetSelectionCenter(MeshObject* obj, vector3& centerOut)
{
	return CalculateCenter(centerOut, obj->get_selector(), obj->position_accessor());
}


bool SetSelectionCenter(MeshObject* obj, const vector3& center)
{
	vector3 oldCenter;
	if(GetSelectionCenter(obj, oldCenter)){
		vector3 d;
		VecSubtract(d, center, oldCenter);
		TranslateSelection(obj->get_selector(), d, obj->position_accessor());
		return true;
	}
	return false;
}


void Move(MeshObject* obj, const vector3& offset)
{
	TranslateSelection(obj->get_selector(), offset, obj->position_accessor());
}


void MoveAlongNormal(MeshObject* obj, number offset, bool usePrecalculatedNormals)
{
	Grid& grid = obj->get_grid();
	Selector& sel = obj->get_selector();
	MeshObject::position_accessor_t& aaPos = obj->position_accessor();
	MeshObject::normal_accessor_t& aaNorm = obj->normal_accessor();

	Grid::face_traits::secure_container faces;
	std::vector<VertexBase*> vrts;
	CollectVerticesTouchingSelection(vrts, sel);

	for(std::vector<VertexBase*>::iterator iter = vrts.begin(); iter != vrts.end(); ++iter)
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

void ScaleAroundCenter(MeshObject* obj, const vector3& scale)
{
	MeshObject::position_accessor_t& aaPos = obj->position_accessor();

	std::vector<VertexBase*> vrts;
	CollectVerticesTouchingSelection(vrts, obj->get_selector());
	vector3 center;
	CalculateCenter(vrts.begin(), vrts.end(), aaPos);

	for(std::vector<VertexBase*>::iterator iter = vrts.begin(); iter != vrts.end(); ++iter)
	{
		ug::vector3& v = aaPos[*iter];
		VecSubtract(v, v, center);
		v.x *= scale.x;
		v.y *= scale.y;
		v.z *= scale.z;
		VecAdd(v, v, center);
	}
}

void ScaleAroundPivot(MeshObject* obj, const vector3& scale)
{
	MeshObject::position_accessor_t& aaPos = obj->position_accessor();

	std::vector<VertexBase*> vrts;
	CollectVerticesTouchingSelection(vrts, obj->get_selector());
	vector3 center = obj->pivot();

	for(std::vector<VertexBase*>::iterator iter = vrts.begin(); iter != vrts.end(); ++iter)
	{
		vector3& v = aaPos[*iter];
		VecSubtract(v, v, center);
		v.x *= scale.x;
		v.y *= scale.y;
		v.z *= scale.z;
		VecAdd(v, v, center);
	}
}


void RotateAroundCenter(MeshObject* obj, const vector3& rot)
{
	MeshObject::position_accessor_t& aaPos = obj->position_accessor();

	std::vector<VertexBase*> vrts;
	CollectVerticesTouchingSelection(vrts, obj->get_selector());
	vector3 center;
	CalculateCenter(vrts.begin(), vrts.end(), aaPos);

//	todo: combine rotation and transform matrix and use ugs build in methods.
//	rotation matrix
	matrix33 matRot;
	MatRotationYawPitchRoll(matRot, rot.x, rot.y, rot.z);

	for(std::vector<VertexBase*>::iterator iter = vrts.begin(); iter != vrts.end(); ++iter)
	{
		vector3 v = aaPos[*iter];
		VecSubtract(v, v, center);
		MatVecMult(aaPos[*iter], matRot, v);
		VecAdd(aaPos[*iter], aaPos[*iter], center);
	}
}


void RotateAroundPivot(MeshObject* obj, const vector3& rot)
{
	MeshObject::position_accessor_t& aaPos = obj->position_accessor();

	std::vector<VertexBase*> vrts;
	CollectVerticesTouchingSelection(vrts, obj->get_selector());
	vector3 center = obj->pivot();

//	todo: combine rotation and transform matrix and use ugs build in methods.
//	rotation matrix
	matrix33 matRot;
	MatRotationYawPitchRoll(matRot, rot.x, rot.y, rot.z);

	for(std::vector<VertexBase*>::iterator iter = vrts.begin(); iter != vrts.end(); ++iter)
	{
		vector3 v = aaPos[*iter];
		VecSubtract(v, v, center);
		MatVecMult(aaPos[*iter], matRot, v);
		VecAdd(aaPos[*iter], aaPos[*iter], center);
	}
}


void ConeTransform(MeshObject* obj, const vector3& base, const vector3& axis,
					number scaleAtTip)
{
	MeshObject::position_accessor_t& aaPos = obj->position_accessor();

	std::vector<VertexBase*> vrts;
	CollectVerticesTouchingSelection(vrts, obj->get_selector());

	for(std::vector<VertexBase*>::iterator iter = vrts.begin(); iter != vrts.end(); ++iter)
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


void LaplacianSmooth(MeshObject* obj, number alpha, int numIterations)
{
	std::vector<VertexBase*> vrts;
	CollectVerticesTouchingSelection(vrts, obj->get_selector());

	ug::LaplacianSmooth(obj->get_grid(), vrts.begin(), vrts.end(),
						obj->position_accessor(), alpha, numIterations);
}


void ProjectToLimitPLoop(MeshObject* obj)
{
	ProjectToLimitPLoop(obj->get_grid(), obj->position_attachment(),
						obj->position_attachment());
}


void ProjectToLimitSmoothBoundary(MeshObject* obj)
{
	ProjectToLimitSubdivBoundary(obj->get_grid(), obj->position_attachment(),
								 obj->position_attachment());
}


void SetPivot(MeshObject* obj, const vector3& pos)
{
	obj->set_pivot(pos);
}


void SetPivotToCenter(MeshObject* obj)
{
	vector3 center;
	CalculateCenter(center, obj->get_selector(), obj->position_accessor());
	obj->set_pivot(center);

}


void FlattenBentQuadrilaterals(MeshObject* obj, number stepSize, int numIterations)
{
	Selector& sel = obj->get_selector();
	MeshObject::position_accessor_t& aaPos = obj->position_accessor();

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


void ProjectToPlane(MeshObject* obj, const vector3& planeCenter, const vector3& planeNormal)
{
	MeshObject::position_accessor_t& aaPos = obj->position_accessor();

	std::vector<VertexBase*> vrts;
	CollectVerticesTouchingSelection(vrts, obj->get_selector());

	for(std::vector<VertexBase*>::iterator iter = vrts.begin(); iter != vrts.end(); ++iter)
	{
		vector3 v = aaPos[*iter];
		ProjectPointToPlane(aaPos[*iter], v, planeCenter, planeNormal);
	}
}


}}// end of namespace

#endif
