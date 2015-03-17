// created by Sebastian Reiter
// s.b.reiter@gmail.com

#include "mesh.h"

namespace ug{
namespace promesh{

Mesh::Mesh() : m_creaseHandler(SHE_VERTEX | SHE_EDGE)
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
	m_subsetHandler.enable_strict_inheritance(true);
	m_creaseHandler.assign_grid(m_grid);
	m_creaseHandler.subset_info(REM_CREASE).name = "crease";
	m_creaseHandler.subset_info(REM_FIXED).name = "fixed";
	m_selector.assign_grid(m_grid);
	m_pivot = vector3(0, 0, 0);
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
