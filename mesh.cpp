/*
 * Copyright (c) 2014-2015:  G-CSC, Goethe University Frankfurt
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

#include "mesh.h"

namespace ug{
namespace promesh{

Mesh::Mesh() : m_creaseHandler(SHE_VERTEX | SHE_EDGE | SHE_FACE)
{
	init();
}

Mesh::Mesh(const Mesh& mesh)
{
	init();
	m_grid = mesh.m_grid;
	// m_selector = mesh.m_selector;
	// m_subsetHandler = mesh.m_subsetHandler;
	// m_creaseHandler = mesh.m_creaseHandler;
	// m_pivot = mesh.m_pivot;
}

void Mesh::init()
{
	m_grid.attach_to_vertices(aPosition, true);
	m_aaPos.access(m_grid, aPosition);
	m_grid.attach_to_faces(aNormal, true);
	m_aaNorm.access(m_grid, aNormal);
	m_subsetHandler.assign_grid(m_grid);
	m_subsetHandler.enable_strict_inheritance(false);
	m_creaseHandler.assign_grid(m_grid);
	m_creaseHandler.subset_info(REM_CREASE).name = "crease";
	m_creaseHandler.subset_info(REM_FIXED).name = "fixed";
	m_selector.assign_grid(m_grid);
	m_pivot = vector3(0, 0, 0);
	m_geometry = make_sp(new Geometry<3, 3>(m_grid, aPosition));
	m_projectionHandler.set_subset_handler(&m_subsetHandler);
	m_projectionHandler.set_geometry(m_geometry);
}

Vertex*	Mesh::
create_vertex(const vector3& p)
{
	Vertex* v = *grid().create<RegularVertex>();
	set_position(v, p);
	return v;
}

Edge* Mesh::
create_edge(Vertex* v0, Vertex* v1)
{
	return *grid().create<RegularEdge>(EdgeDescriptor(v0, v1));
}

Face* Mesh::
create_triangle(Vertex* v0, Vertex* v1, Vertex* v2)
{
	return *grid().create<Triangle>(TriangleDescriptor(v0, v1, v2));
}

Face* Mesh::
create_quadrilateral(Vertex* v0, Vertex* v1, Vertex* v2, Vertex* v3)
{
	return *grid().create<Quadrilateral>(QuadrilateralDescriptor(v0, v1, v2, v3));
}

Volume*	Mesh::
create_tetrahedron(Vertex* v0, Vertex* v1, Vertex* v2, Vertex* v3)
{
	return *grid().create<Tetrahedron>(TetrahedronDescriptor(v0, v1, v2, v3));
}

Volume*	Mesh::
create_pyramid(Vertex* v0, Vertex* v1, Vertex* v2, Vertex* v3, Vertex* v4)
{
	return *grid().create<Pyramid>(PyramidDescriptor(v0, v1, v2, v3, v4));
}

Volume*	Mesh::
create_prism(Vertex* v0, Vertex* v1, Vertex* v2,
			 Vertex* v3, Vertex* v4, Vertex* v5)
{
	return *grid().create<Prism>(PrismDescriptor(v0, v1, v2, v3, v4, v5));
}

Volume*	Mesh::
create_hexahedron(Vertex* v0, Vertex* v1, Vertex* v2, Vertex* v3,
				  Vertex* v4, Vertex* v5, Vertex* v6, Vertex* v7)
{
	return *grid().create<Hexahedron>(HexahedronDescriptor(v0, v1, v2, v3, v4, v5, v6, v7));
}

Volume*	Mesh::
create_octahedron(Vertex* v0, Vertex* v1, Vertex* v2,
				  Vertex* v3, Vertex* v4, Vertex* v5)
{
	return *grid().create<Octahedron>(OctahedronDescriptor(v0, v1, v2, v3, v4, v5));
}


}//	end of namespace	
}//	end of namespace
