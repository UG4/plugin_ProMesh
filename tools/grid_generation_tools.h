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

#ifndef __H__UG__grid_generation_tools__
#define __H__UG__grid_generation_tools__

#include <vector>
#include "../mesh.h"
#include "lib_grid/algorithms/selection_util.h"
#include "lib_grid/algorithms/grid_generation/icosahedron.h"
#include "lib_grid/algorithms/remeshing/delaunay_triangulation.h"

#define	TOOLTIP_NEW_OBJECT "Creates a new empty object."
#define TOOLTIP_NEW_CSG_OBJECT "Creates a new csg object."
#define	TOOLTIP_MERGE_OBJECTS "Merges the selected objects into a new one."
#define TOOLTIP_CLONE_MESH "Creates a new mesh and copies all content from the given mesh into the new instance."
#define	TOOLTIP_COPY_SELECTION "Copies the selected elements to a new mesh."
#define	TOOLTIP_CREATE_VERTEX "Creates a new vertex"
#define	TOOLTIP_CREATE_EDGE "Creates an edge between two selected vertices."
#define	TOOLTIP_CREATE_FACE "Creates a face between selected vertices."
#define	TOOLTIP_CREATE_VOLUME "Creates a volume between selected vertices."
#define	TOOLTIP_CREATE_PLANE "Creates a plane."
#define	TOOLTIP_CREATE_CIRCLE "Creates a circle."
#define	TOOLTIP_CREATE_BOX "Creates a box."
#define	TOOLTIP_CREATE_SPHERE "Creates a sphere."
#define	TOOLTIP_CREATE_TETRAHEDRON "Creates a tetrahedron."
#define	TOOLTIP_CREATE_PYRAMID "Creates a pyramid."
#define	TOOLTIP_CREATE_PRISM "Creates a prism."
#define	TOOLTIP_CREATE_DUALGRID "creates the dual grid consisting of control volumes as used in the finite volume method" 
#define TOOLTIP_CREATE_TKD "Creates a tetrakaidecahedral cell"
#define TOOLTIP_CREATE_TKD_WITH_OUTER_LAYER "Creates a tetrakaidecahedral cell with a surrounding layer"

namespace ug{
namespace promesh{

/// \addtogroup promesh
/// \{
SmartPtr<Mesh> CloneMesh(Mesh* mesh);

void CopySelection(Mesh* srcMesh, Mesh* destMesh);

Vertex* CreateVertex(Mesh* obj, const vector3& pos, int subsetInd);

Edge* CreateEdge(Mesh* obj, int subsetInd);

Face* CreateFace(Mesh* obj, int subsetInd);

Volume* CreateVolume(Mesh* obj, int subsetInd);

void CreatePlane(	Mesh* obj,
					const vector3& upLeft,
					const vector3& upRight,
					const vector3& lowLeft,
					const vector3& lowRight,
					int subsetInd,
				 	bool fill);

void CreatePlane(	Mesh* obj,
					number width,
					number height,
					const vector3& center,
					int subsetInd,
					bool fill);

void CreateCircle(	Mesh* obj,
					const vector3& center,
					number radius,
				  	int numRimVertices,
				  	int subsetInd,
					bool fill);

void CreateBox(	Mesh* obj,
				const vector3& boxMin,
				const vector3& boxMax,
				int subsetInd,
				bool fill);

void CreateSphere(	Mesh* obj,
					const vector3& center,
					number radius,
					int numRefinements,
					int subsetInd);

void CreateTetrahedron(Mesh* obj, int subsetInd, bool fill);

void CreatePyramid(Mesh* obj, int subsetInd, bool fill);

void CreatePrism(Mesh* obj, int subsetInd, bool fill);

void CreateTKD(Mesh* obj, int subsetInd, number a, number w, number h);

void CreateTKDWithOuterLayer(	Mesh* obj,
                             	int innerSubsetInd,
                             	int outerSubsetInd,
                             	number a,
                             	number w,
                             	number h,
                             	number d);

/// \}
}}// end of namespace

#endif
