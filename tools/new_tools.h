// created by Sebastian Reiter
// s.b.reiter@gmail.com
// Jun 27, 2013 (d,m,y)

#ifndef __H__UG_PROMESH__new_tools__
#define __H__UG_PROMESH__new_tools__

#include "../mesh_object.h"
#include "lib_grid/algorithms/geom_obj_util/misc_util.h"

namespace ug{
namespace promesh{

class Box{
	public:
		void set_min(const vector3& val)	{min = val;}
		void set_max(const vector3& val)	{max = val;}
		const vector3& get_min() const		{return min;}
		const vector3& get_max() const		{return max;}

		vector3 min;
		vector3 max;
};


////////////////////////////////////////////////////////////////////////////////
//	INFO
inline SmartPtr<Box> GetBoundingBox(MeshObject* obj)
{
	SmartPtr<Box> box = make_sp(new Box);
	CalculateBoundingBox(box->min, box->max, obj->get_grid().vertices_begin(),
						 obj->get_grid().vertices_end(), obj->position_accessor());
	return box;
}

////////////////////////////////////////////////////////////////////////////////
//	COORDINATE TRANSFORM
inline void ScaleAroundPoint(MeshObject* obj, const vector3& scale, const vector3& point)
{
	MeshObject::position_accessor_t& aaPos = obj->position_accessor();

	std::vector<Vertex*> vrts;
	CollectVerticesTouchingSelection(vrts, obj->get_selector());

	for(std::vector<Vertex*>::iterator iter = vrts.begin(); iter != vrts.end(); ++iter)
	{
		vector3& v = aaPos[*iter];
		VecSubtract(v, v, point);
		v.x() *= scale.x();
		v.y() *= scale.y();
		v.z() *= scale.z();
		VecAdd(v, v, point);
	}
}



////////////////////////////////////////////////////////////////////////////////
//	SELECTION
///	Selects elements whose center lie in a box
template <class TElem>
inline void SelectElementsInBox(MeshObject* obj, const vector3& min, const vector3& max)
{
	MeshObject::position_accessor_t& aaPos = obj->position_accessor();

	Grid& grid = obj->get_grid();
	Selector& sel = obj->get_selector();

	for(typename Grid::traits<TElem>::iterator iter = grid.begin<TElem>();
		iter != grid.end<TElem>(); ++iter)
	{
		if(BoxBoundProbe(CalculateCenter(*iter, aaPos), min, max))
			sel.select(*iter);
	}
}


////////////////////////////////////////////////////////////////////////////////
//	SELECTION
///	Selects elements whose center lie in a cylinder
template <class TElem>
inline void SelectElementsInCylinder(MeshObject* obj, const vector3& cylBase,
						 	  const vector3& cylTop, number radius)
{
	MeshObject::position_accessor_t& aaPos = obj->position_accessor();

	Grid& grid = obj->get_grid();
	Selector& sel = obj->get_selector();

	vector3 from = cylBase;
	vector3 dir;
	VecSubtract(dir, cylTop, cylBase);

	for(typename Grid::traits<TElem>::iterator iter = grid.begin<TElem>();
		iter != grid.end<TElem>(); ++iter)
	{
		vector3 c = CalculateCenter(*iter, aaPos);
		vector3 p;
		number s = ProjectPointToRay(p, c, from, dir);
		if((s > -SMALL) && (s < (1. + SMALL))){
			if(VecDistanceSq(p, c) <= sq(radius))
				sel.select(*iter);
		}
	}
}

}}// end of namespace

#endif
