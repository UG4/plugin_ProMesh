/*
 * Copyright (c) 2014-2017:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG_measure_tools
#define __H__UG_measure_tools

#include "../mesh.h"
#include "lib_grid/algorithms/volume_calculation.h"

#define TOOLTIP_MEASURE_GRID_LENGTH "Measures the length of all edges of a grid"
#define TOOLTIP_MEASURE_GRID_AREA "Measures the area of all faces of a grid"
#define TOOLTIP_MEASURE_GRID_VOLUME "Measures the volume of all volume elements of a grid"

#define TOOLTIP_MEASURE_SUBSET_LENGTH "Measures the length of all edges of the given subset"
#define TOOLTIP_MEASURE_SUBSET_AREA "Measures the area of all faces of the given subset"
#define TOOLTIP_MEASURE_SUBSET_VOLUME "Measures the volume of all volume elements of the given subset"

#define TOOLTIP_MEASURE_SELECTION_LENGTH "Measures the length of all edges of the current selection"
#define TOOLTIP_MEASURE_SELECTION_AREA "Measures the area of all faces of the current selection"
#define TOOLTIP_MEASURE_SELECTION_VOLUME "Measures the volume of all volume elements of the current selection"

#define TOOLTIP_PROJECTED_DISTANCE "Measures the distance of the projection of two selected vertices into the plane with the specified normal."

namespace ug{
namespace promesh{

/// \addtogroup promesh
/// \{
inline number MeasureGridLength(Mesh* obj)
{
	return CalculateVolume(obj->grid().begin<Edge>(),
						   obj->grid().end<Edge>(),
						   obj->position_accessor());
}

inline number MeasureGridArea(Mesh* obj)
{
	return CalculateVolume(obj->grid().begin<Face>(),
						   obj->grid().end<Face>(),
						   obj->position_accessor());
}

inline number MeasureGridVolume(Mesh* obj)
{
	return CalculateVolume(obj->grid().begin<Volume>(),
						   obj->grid().end<Volume>(),
						   obj->position_accessor());
}



inline number MeasureSubsetLength(Mesh* obj, int subsetInd)
{
	return CalculateVolume(obj->subset_handler().begin<Edge>(subsetInd),
						   obj->subset_handler().end<Edge>(subsetInd),
						   obj->position_accessor());
}

inline number MeasureSubsetArea(Mesh* obj, int subsetInd)
{
	return CalculateVolume(obj->subset_handler().begin<Face>(subsetInd),
						   obj->subset_handler().end<Face>(subsetInd),
						   obj->position_accessor());
}

inline number MeasureSubsetVolume(Mesh* obj, int subsetInd)
{
	return CalculateVolume(obj->subset_handler().begin<Volume>(subsetInd),
						   obj->subset_handler().end<Volume>(subsetInd),
						   obj->position_accessor());
}



inline number MeasureSelectionLength(Mesh* obj)
{
	return CalculateVolume(obj->selector().begin<Edge>(),
						   obj->selector().end<Edge>(),
						   obj->position_accessor());
}

inline number MeasureSelectionArea(Mesh* obj)
{
	return CalculateVolume(obj->selector().begin<Face>(),
						   obj->selector().end<Face>(),
						   obj->position_accessor());
}

inline number MeasureSelectionVolume(Mesh* obj)
{
	return CalculateVolume(obj->selector().begin<Volume>(),
						   obj->selector().end<Volume>(),
						   obj->position_accessor());
}

number ProjectedDistance(Mesh* obj, const vector3& projectionNormal);

/// \}
}//	end of namespace promesh
}//	end of namespace ug

#endif	//__H__measure_tools
