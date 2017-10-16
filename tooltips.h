/*
 * Copyright (c) 2014-2015:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter, Anahita Werner
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


#ifndef __H__PROMESH_TOOLTIPS
#define __H__PROMESH_TOOLTIPS

//	camera tools
#define	TOOLTIP_CENTER_OBJECT "Centers the current object."
#define	TOOLTIP_SLIDER_TEST "Tests the DoubleSlider."
#define	TOOLTIP_CENTER_SELECTION "Centers the current selection."
#define TOOLTIP_TOP_VIEW "View the current scene from the top."
#define	TOOLTIP_FRONT_VIEW "View the current scene from the front."
#define	TOOLTIP_SIDE_VIEW "View the current scene from the side."
#define	TOOLTIP_HIDE_SELECTED_ELEMENTS "Hides all currently selected elements."
#define	TOOLTIP_UNHIDE_ELEMENTS "Unhides all hidden elements."


//info tools
#define	TOOLTIP_PRINT_SELECTION_CENTER "Calculates and prints the position of the center of the current selection."
#define	TOOLTIP_PRINT_GEOMETRY_INFO "Prints info about the current geometry"
#define	TOOLTIP_PRINT_FACE_QUALITY "Prints the quality of selected faces"
#define	TOOLTIP_PRINT_SELECTION_INFO "Prints the quantities of selected elements"
#define	TOOLTIP_PRINT_SELECTION_CONTAINING_SUBSETS "Prints subset indices of all subsets, which contain a selected element."
#define	TOOLTIP_PRINT_VERTEX_DISTANCE "Prints the min and max distance of vertices of selected elements."
#define	TOOLTIP_PRINT_LEAST_SQUARES_PLANE "Prints the position and normal of the least squares fitting plane"


//fracture tools
#define	TOOLTIP_EXPAND_LAYERS_2D "Expands a 1d layer to a 2d layer by introducing quadrilaterals."
#define	TOOLTIP_EXPAND_LAYERS_3D "Expands a 2d layer to a 3d layer by introducing prisms."
#define	TOOLTIP_FRAC_TO_LAYER "Enhances a 2d fracture to a 3d fracture."

//new tools
#define TOOLTIP_GET_BOUNDING_BOX "Returns the bounding box of the specified mesh."

//	all changes should go above the endif statement
#endif
