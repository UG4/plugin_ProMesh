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

#ifndef __H__UG_PROMESH__new_tools__
#define __H__UG_PROMESH__new_tools__

#include "../mesh.h"
#include "lib_grid/algorithms/geom_obj_util/misc_util.h"
#include "lib_grid/algorithms/selection_util.h"

namespace ug{
namespace promesh{

/// \addtogroup promesh
/// \{
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
inline SmartPtr<Box> GetBoundingBox(Mesh* obj)
{
	SmartPtr<Box> box = make_sp(new Box);
	CalculateBoundingBox(box->min, box->max, obj->grid().vertices_begin(),
						 obj->grid().vertices_end(), obj->position_accessor());
	return box;
}

////////////////////////////////////////////////////////////////////////////////
//	COORDINATE TRANSFORM
inline void ScaleAroundPoint(Mesh* obj, const vector3& scale, const vector3& point)
{
	Mesh::position_accessor_t& aaPos = obj->position_accessor();

	std::vector<Vertex*> vrts;
	CollectVerticesTouchingSelection(vrts, obj->selector());

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
inline void SelectElementsInBox(Mesh* obj, const vector3& min, const vector3& max)
{
	Mesh::position_accessor_t& aaPos = obj->position_accessor();

	Grid& grid = obj->grid();
	Selector& sel = obj->selector();

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
inline void SelectElementsInCylinder(Mesh* obj, const vector3& cylBase,
						 	  const vector3& cylTop, number radius)
{
	Mesh::position_accessor_t& aaPos = obj->position_accessor();

	Grid& grid = obj->grid();
	Selector& sel = obj->selector();

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
/// \}
}}// end of namespace

#endif
