/*
 * Copyright (c) 2015:  G-CSC, Goethe University Frankfurt
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

#include "registration_routines.h"
#include "mesh.h"
#include "registry/registry.h"
#include "bridge/util.h"
#include "tooltips.h"

using namespace std;
using namespace ug::bridge;

#define TOOLTIP_MESH	"The Mesh class stores a Grid, SubsetHandler and Selector. Nearly all algorithms in ProMesh operate on Meshes"
#define TOOLTIP_ITERATOR	"Iterators are used to iterate over the elements (vertices, edges, faces, volumes) of a Mesh"

namespace ug{
namespace promesh{

///	helper function to register mesh-element-iterators
template <class TElem>
static void RegisterElementIterators(ProMeshRegistry& reg, string name, string grp)
{
	typedef ElementIterator<TElem> iter_t;
	reg.add_class_<iter_t>(name, grp, TOOLTIP_ITERATOR)
		.add_method("clone", &iter_t::clone)
		.add_method("assign", &iter_t::assign)
		.add_method("value", &iter_t::value)
		.add_method("advance", &iter_t::advance)
		.add_method("equal", &iter_t::equal)
		.add_method("unequal", &iter_t::unequal)
		.set_construct_as_smart_pointer(true);
}

template <class TElem>
static void RegisterMeshIteratorMethods(ExportedClass<Mesh>& cls, string elemName)
{
	 cls.add_method(elemName + "_begin", &Mesh::begin<TElem>)
		.add_method(elemName + "_selection_begin", &Mesh::selection_begin<TElem>)
		.add_method(elemName + "_subset_begin", &Mesh::subset_begin<TElem>)
		.add_method(elemName + "_end", &Mesh::end<TElem>)
		.add_method(elemName + "_selection_end", &Mesh::selection_end<TElem>)
		.add_method(elemName + "_subset_end", &Mesh::subset_end<TElem>);
}

void RegisterMesh(ProMeshRegistry& reg, string baseGrp)
{
	try{
		string grp = baseGrp;
		RegisterElementIterators<Vertex>(reg, "VertexIterator", grp);
		RegisterElementIterators<Edge>(reg, "EdgeIterator", grp);
		RegisterElementIterators<Face>(reg, "FaceIterator", grp);
		RegisterElementIterators<Volume>(reg, "VolumeIterator", grp);

		ExportedClass<Mesh>& meshCls = reg.add_class_<Mesh>("Mesh", grp, TOOLTIP_MESH);
		meshCls.add_constructor()
				.add_method("crease_handler", &Mesh::crease_handler)
				.add_method("create_vertex", &Mesh::create_vertex)
				.add_method("create_edge", &Mesh::create_edge)
				.add_method("create_triangle", &Mesh::create_triangle)
				.add_method("create_quadrilateral", &Mesh::create_quadrilateral)
				.add_method("create_hexahedron", &Mesh::create_hexahedron)
				.add_method("create_octahedron", &Mesh::create_octahedron)
				.add_method("create_prism", &Mesh::create_prism)
				.add_method("create_pyramid", &Mesh::create_pyramid)
				.add_method("create_tetrahedron", &Mesh::create_tetrahedron)
				.add_method("grid", &Mesh::grid)
				.add_method("pivot", &Mesh::pivot)
				.add_method("position", &Mesh::position)
				.add_method("position_attachment", &Mesh::position_attachment)
				.add_method("selector", &Mesh::selector)
				.add_method("set_pivot", &Mesh::set_pivot)
				.add_method("set_position", &Mesh::set_position)
				.add_method("subset_handler", &Mesh::subset_handler)
				.set_construct_as_smart_pointer(true);

		RegisterMeshIteratorMethods<Vertex>(meshCls, "vertex");
		RegisterMeshIteratorMethods<Edge>(meshCls, "edge");
		RegisterMeshIteratorMethods<Face>(meshCls, "face");
		RegisterMeshIteratorMethods<Triangle>(meshCls, "triangle");
		RegisterMeshIteratorMethods<Quadrilateral>(meshCls, "quadrilateral");
		RegisterMeshIteratorMethods<Volume>(meshCls, "volume");
		RegisterMeshIteratorMethods<Tetrahedron>(meshCls, "tetrahedron");
		RegisterMeshIteratorMethods<Pyramid>(meshCls, "pyramid");
		RegisterMeshIteratorMethods<Prism>(meshCls, "prism");
		RegisterMeshIteratorMethods<Hexahedron>(meshCls, "hexahedron");
		RegisterMeshIteratorMethods<Octahedron>(meshCls, "octahedron");
	}
	UG_REGISTRY_CATCH_THROW(baseGrp);
}

} // end namespace promesh
}// end of namespace

