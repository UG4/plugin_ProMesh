/*
 * Copyright (c) 2017:  G-CSC, Goethe University Frankfurt
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

#include "quality_tools.h"
#include "lib_grid/algorithms/quality_util.h"
#include "common/util/table.h"

namespace ug {
namespace promesh {

void PrintFaceAspectRatios (Mesh* msh)
{
	Selector& sel = msh->selector();
	AspectRatioInfo ari;
	if(sel.num<Face>() > 0){
		UG_LOG("Face Aspect Ratios (" << sel.num<Face>() << " faces):\n\n");
		ari = GetFaceAspectRatioInfo (sel.faces_begin(), sel.faces_end(),
									  msh->position_accessor());
	}
	else{
		Grid& g = msh->grid();
		UG_LOG("Face Aspect Ratios (" << g.num<Face>() << " faces):\n\n");
		ari = GetFaceAspectRatioInfo (g.faces_begin(), g.faces_end(),
									  msh->position_accessor());
	}

	UG_LOG(ari.to_string() << std::endl);
}

void PrintFaceAspectRatioHistogram (Mesh* msh, int numHistoSecs)
{
	if(numHistoSecs < 1){
		UG_LOG("Can't create histogram with " << numHistoSecs << " sections.\n");
	}

	Selector& sel = msh->selector();
	std::vector<int> histo;
	if(sel.num<Face>() > 0){
		UG_LOG("Face Aspect Ratio Histogram (" << sel.num<Face>() << " faces):\n\n");
		GetFaceAspectRatioHistogram (histo, sel.faces_begin(), sel.faces_end(),
	                             	 numHistoSecs, msh->position_accessor());
	}
	else{
		Grid& g = msh->grid();
		UG_LOG("Face Aspect Ratio Histogram (" << g.num<Face>() << " faces):\n\n");
		GetFaceAspectRatioHistogram (histo, g.faces_begin(), g.faces_end(),
	                             	 numHistoSecs, msh->position_accessor());
	}

	const number stepSize = 1. / (number) numHistoSecs;
	StringStreamTable t;

	for(size_t i = 0; i < histo.size(); ++i){
		t(0, i) << (number)i * stepSize << " - " << (number)(i+1) * stepSize;
		t(1, i) << histo[i];
	}

	UG_LOG(t.to_string() << std::endl);
}

}//	end of namespace
}//	end of namespace
