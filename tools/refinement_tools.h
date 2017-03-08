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

#ifndef __H__UG__refinement_tools__
#define __H__UG__refinement_tools__

#include "../mesh.h"
#include "lib_grid/algorithms/geom_obj_util/edge_util.h"
#include "lib_grid/refinement/regular_refinement.h"
#include "lib_grid/refinement/hanging_node_refiner_grid.h"
#include "lib_grid/refinement/projectors/sphere_projector.h"
#include "lib_grid/refinement/projectors/subdivision_projector.h"
#include "lib_grid/callbacks/callbacks.h"

#define	TOOLTIP_REFINE "Refines selected elements and builds a regular closure."
#define TOOLTIP_REFINE_WITH_SNAP_POINTS "Refines selected elements so that new edges are built between midpoints of selected edges and selected vertices, if possible."
#define	TOOLTIP_HANGING_NODE_REFINE "Refines selected elements using hanging nodes"
#define	TOOLTIP_REFINE_SMOOTH "Refines selected elements using piecewise smooth refinement."
#define	TOOLTIP_REFINE_SMOOTH_BOUNDARY_2D "Refines selected elements using smooth subdivision rules on the boundary edges."
#define	TOOLTIP_FRACTURED_MEDIA_REFINE "Refines selected elements using hanging nodes. Fractures are refined anisotropic."
#define	TOOLTIP_CREATE_FRACTAL "Refines the whole geometry using a fractal-refinement scheme-"
#define	TOOLTIP_INSERT_CENTER "Inserts a central vertex in all selected elements."

namespace ug{
namespace promesh{

/// \addtogroup promesh
/// \{

void Refine(Mesh* obj);
void RefineWithSnapPoints(Mesh* obj);
void Refine(Mesh* obj, bool strictSubsetInheritance, bool useSnapPoints);

void HangingNodeRefine(Mesh* obj, bool anisotropic);
void HangingNodeRefine(Mesh* obj, bool strictSubsetInheritance, bool anisotropic);

void RefineSmooth(Mesh* obj);
void RefineSmooth(Mesh* obj, bool strictSubsetInheritance);

void InsertCenter(Mesh* obj);
void InsertCenter(Mesh* obj, bool strictSubsetInheritance);

/// \}

}}// end of namespace

#endif
