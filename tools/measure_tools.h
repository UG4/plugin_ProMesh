// created by Sebastian Reiter
// s.b.reiter@gmail.com

#ifndef __H__UG_measure_tools
#define __H__UG_measure_tools

#include "../mesh.h"
#include "lib_grid/algorithms/volume_calculation.h"

namespace ug{
namespace promesh{

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

}//	end of namespace promesh
}//	end of namespace ug

#endif	//__H__measure_tools
