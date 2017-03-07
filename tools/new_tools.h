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

		void global_to_local (const vector3& glob, vector3& locOut)
		{
			for(size_t i = 0; i < 3; ++i){
				number d = max[i] - min[i];
				if(d != 0)
					locOut[i] = (glob[i] - min[i]) / d;
				else
					locOut[i] = 0;
			}
		}

		void local_to_global (const vector3& loc, vector3& globOut)
		{
			for(size_t i = 0; i < 3; ++i){
				number d = max[i] - min[i];
				globOut[i] = loc[i] * d + min[i];
			}
		}

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

/// \}

}}// end of namespace

#endif
