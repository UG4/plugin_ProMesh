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

#include "grid_generation_tools.h"
#include "lib_grid/algorithms/orientation_util.h"

namespace ug{
namespace promesh{

SmartPtr<Mesh> CloneMesh(Mesh* mesh)
{
	return make_sp(new Mesh(*mesh));
}

/** Make sure that aNewVrt is attached to srcMesh->grid() and contains a
 * pointer to a valid vertex in destMesh for each selected vertex in srcMesh.
 * Also make sure that all vertices belonging to a selected element have been
 * selected, too.*/
template <class TElem>
void CopySelectedElements(Mesh* srcMesh, Mesh* destMesh, AVertex aNewVrt)
{
	Grid& srcGrid						= srcMesh->grid();
	Selector& srcSel					= srcMesh->selector();
	SubsetHandler& srcSH				= srcMesh->subset_handler();
	SubsetHandler& srcCreaseSH			= srcMesh->crease_handler();

	Grid& destGrid						= destMesh->grid();
	SubsetHandler& destSH				= destMesh->subset_handler();
	SubsetHandler& destCreaseSH			= destMesh->crease_handler();

	Grid::VertexAttachmentAccessor<AVertex> aaNewVrt(srcGrid, aNewVrt);

	CustomVertexGroup vrts;
	typedef typename Grid::traits<TElem>::iterator iter_t;

	for(iter_t eiter = srcSel.begin<TElem>();
		eiter != srcSel.end<TElem>(); ++eiter)
	{
		TElem* e = *eiter;
		vrts.resize(e->num_vertices());
		for(size_t iv = 0; iv < e->num_vertices(); ++iv)
			vrts.set_vertex(iv, aaNewVrt[e->vertex(iv)]);

		TElem* ne = *destGrid.create_by_cloning(e, vrts);
		destSH.assign_subset(ne, srcSH.get_subset_index(e));
		
		if(TElem::dim < 2)
			destCreaseSH.assign_subset(ne, srcCreaseSH.get_subset_index(e));
	}
}

void CopySelection(Mesh* srcMesh, Mesh* destMesh){
	Grid& srcGrid						= srcMesh->grid();
	Selector& srcSel					= srcMesh->selector();
	SubsetHandler& srcSH				= srcMesh->subset_handler();
	SubsetHandler& srcCreaseSH			= srcMesh->crease_handler();
	Mesh::position_accessor_t aaPosSrc	= srcMesh->position_accessor();

	Grid& destGrid						= destMesh->grid();
	SubsetHandler& destSH				= destMesh->subset_handler();
	SubsetHandler& destCreaseSH			= destMesh->crease_handler();
	Mesh::position_accessor_t aaPosDest	= destMesh->position_accessor();

	AVertex aNewVrt;
	srcGrid.attach_to_vertices(aNewVrt);
	Grid::VertexAttachmentAccessor<AVertex> aaNewVrt(srcGrid, aNewVrt);

	for(int si = destSH.num_subsets(); si < srcSH.num_subsets(); ++si)
		destSH.subset_info(si) = srcSH.subset_info(si);

	SelectAssociatedGridObjects(srcSel);
//	create new vertices in destGrid
	for(VertexIterator viter = srcSel.begin<Vertex>();
		viter != srcSel.end<Vertex>(); ++viter)
	{
		Vertex* v = *viter;
		Vertex* nv = *destGrid.create_by_cloning(v);
		aaNewVrt[v] = nv;
		aaPosDest[nv] = aaPosSrc[v];
		destSH.assign_subset(nv, srcSH.get_subset_index(v));
		destCreaseSH.assign_subset(nv, srcCreaseSH.get_subset_index(v));
	}

	CopySelectedElements<Edge>(srcMesh, destMesh, aNewVrt);
	CopySelectedElements<Face>(srcMesh, destMesh, aNewVrt);
	CopySelectedElements<Volume>(srcMesh, destMesh, aNewVrt);

	srcGrid.detach_from_vertices(aNewVrt);
}


Vertex* CreateVertex(Mesh* obj, const vector3& pos, int subsetInd)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	SubsetHandler& sh = obj->subset_handler();
	Mesh::position_accessor_t& aaPos = obj->position_accessor();

//	build a new vertex
	RegularVertex* vrt = *grid.create<RegularVertex>();

	if(vrt){
		aaPos[vrt] = pos;
		sel.clear();
		sel.select(vrt);
		sh.assign_subset(vrt, subsetInd);
	}

	return vrt;
}


Edge* CreateEdge(Mesh* obj, int subsetInd)
{
	using namespace std;
	ug::Selector& sel = obj->selector();
	ug::Grid& grid = obj->grid();
	ug::SubsetHandler& sh = obj->subset_handler();

	size_t numVrts = sel.num<ug::Vertex>();
	vector<ug::Vertex*> vrts;
	vrts.reserve(numVrts);
	vrts.assign(sel.begin<ug::Vertex>(), sel.end<ug::Vertex>());

	ug::RegularEdge* e = NULL;
	switch(numVrts){
		case 2:{//	create edge
				if(!grid.get_edge(vrts[0], vrts[1]))
					e = *grid.create<ug::RegularEdge>(ug::EdgeDescriptor(vrts[0], vrts[1]));
				else{
					UG_LOG("Can't create edge: RegularEdge already exists.\n");
				}
			}break;

		default:
			UG_LOG("Can't create edge: Bad number of vertices. 2 are required.\n");
			break;
	}

//todo: use the currently marked subset.
	if(e)
		sh.assign_subset(e, subsetInd);

	return e;
}


Face* CreateFace(Mesh* obj, int subsetInd)
{
	using namespace std;
	Selector& sel = obj->selector();
	Grid& grid = obj->grid();
	SubsetHandler& sh = obj->subset_handler();

	size_t numVrts = sel.num<Vertex>();

	if(numVrts < 3 || numVrts > 4){
		UG_LOG("Bad number of vertices! Can't create a face element from " << numVrts << " vertices.");
		return NULL;
	}

	vector<Vertex*> vrts;
	vrts.reserve(numVrts);
	vrts.assign(sel.begin<Vertex>(), sel.end<Vertex>());

	FaceDescriptor fd(numVrts);
	for(size_t i = 0; i < numVrts; ++i)
		fd.set_vertex(i, vrts[i]);

	if(grid.get_face(fd)){
		UG_LOG("A face connecting the given vertices already exists. Won't create a new one!");
		return NULL;
	}

	Face* f = NULL;
	switch(numVrts){
		case 3:{//	create triangle
			f = *grid.create<Triangle>(TriangleDescriptor(vrts[0], vrts[1], vrts[2]));
		}break;

		case 4:{//	create quadrilateral
			f = *grid.create<Quadrilateral>(QuadrilateralDescriptor(vrts[0], vrts[1],
																			vrts[2], vrts[3]));
		}break;

		default:
			UG_LOG("Can't create face: Bad number of vertices. 3 or 4 are supported.\n");
			break;
	}

	if(f)
		sh.assign_subset(f, subsetInd);

	return f;
}


Volume* CreateVolume(Mesh* obj, int subsetInd)
{
	using namespace std;
	Selector& sel = obj->selector();
	Grid& grid = obj->grid();
	SubsetHandler& sh = obj->subset_handler();

	size_t numVrts = sel.num<Vertex>();

	if(numVrts < 4 || numVrts > 8){
		UG_LOG("Bad number of vertices! Can't create a volume element from " << numVrts << " vertices.");
		return NULL;
	}

	vector<Vertex*> vrts;
	vrts.reserve(numVrts);
	vrts.assign(sel.begin<Vertex>(), sel.end<Vertex>());

	VolumeDescriptor vd(numVrts);
	for(size_t i = 0; i < numVrts; ++i)
		vd.set_vertex(i, vrts[i]);

	if(grid.get_volume(vd)){
		UG_LOG("A volume connecting the given vertices already exists. Won't create a new one!");
		return NULL;
	}

	Volume* v = NULL;
	switch(numVrts){
		case 4:{//	create tetrahedron
			v = *grid.create<Tetrahedron>(TetrahedronDescriptor(vrts[0], vrts[1],
																vrts[2], vrts[3]));
		}break;

		case 5:{//	create pyramid
			v = *grid.create<Pyramid>(PyramidDescriptor(vrts[0], vrts[1],
														vrts[2], vrts[3], vrts[4]));
		}break;

		case 6:{//	create prism
			v = *grid.create<Prism>(PrismDescriptor(vrts[0], vrts[1], vrts[2],
													vrts[3], vrts[4], vrts[5]));
		}break;

		case 8:{//	create hexahedron
			v = *grid.create<Hexahedron>(HexahedronDescriptor(vrts[0], vrts[1], vrts[2], vrts[3],
															  vrts[4], vrts[5], vrts[6], vrts[7]));
		}break;

		default:
			UG_LOG("Can't create volume: Bad number of vertices. 4, 5, 6, and 8 are supported.\n");
			break;
	}

	if(v){
	//	check and fix orientation
		Grid::VertexAttachmentAccessor<APosition> aaPos(grid, aPosition);
		if(!CheckOrientation(v, aaPos)){
			grid.flip_orientation(v);
		}

		sh.assign_subset(v, subsetInd);
	}

	return v;
}


void CreatePlane(	Mesh* obj,
					const vector3& upLeft,
					const vector3& upRight,
					const vector3& lowLeft,
					const vector3& lowRight,
					int subsetInd,
				 	bool fill)
{
	Grid& grid = obj->grid();
	Mesh::position_accessor_t& aaPos = obj->position_accessor();
	Selector& sel = obj->selector();
	SubsetHandler& sh = obj->subset_handler();

	sel.clear();
	bool autoselEnabled = sel.autoselection_enabled();
	sel.enable_autoselection(true);

//	create the vertices
	Vertex* vrts[4];
	for(size_t i = 0; i < 4; ++i)
		vrts[i] = *grid.create<RegularVertex>();

	aaPos[vrts[0]] = upLeft;
	aaPos[vrts[1]] = lowLeft;
	aaPos[vrts[2]] = lowRight;
	aaPos[vrts[3]] = upRight;

//	create the plane
	if(fill)
		grid.create<Quadrilateral>(QuadrilateralDescriptor(vrts[0], vrts[1], vrts[2], vrts[3]));
	else{
		for(size_t i = 0; i < 4; ++i){
			int i0 = i;
			int i1 = (i + 1) % 4;
			grid.create<RegularEdge>(EdgeDescriptor(vrts[i0], vrts[i1]));
		}
	}

//	assign subset
	sh.assign_subset(sel.begin<Vertex>(), sel.end<Vertex>(), subsetInd);
	sh.assign_subset(sel.begin<Edge>(), sel.end<Edge>(), subsetInd);
	sh.assign_subset(sel.begin<Face>(), sel.end<Face>(), subsetInd);

//	restore selector
	sel.enable_autoselection(autoselEnabled);
}

void CreatePlane(	Mesh* obj,
					number width,
					number height,
					const vector3& center,
					int subsetInd,
					bool fill)
{
	const number wh = 0.5 * width;
	const number hh = 0.5 * height;
	const number cx = center.x();
	const number cy = center.y();
	const number cz = center.z();

	const vector3 upLeft (cx - wh, cy + hh, cz);
	const vector3 upRight (cx + wh, cy + hh, cz);
	const vector3 lowLeft (cx - wh, cy - hh, cz);
	const vector3 lowRight (cx + wh, cy - hh, cz);

	CreatePlane(obj, upLeft, upRight, lowLeft, lowRight, subsetInd, fill);
}


void CreateCircle(	Mesh* obj,
					const vector3& center,
					number radius,
				  	int numRimVertices,
				  	int subsetInd,
					bool fill)
{
	Grid& grid = obj->grid();
	Mesh::position_accessor_t& aaPos = obj->position_accessor();
	Selector& sel = obj->selector();
	SubsetHandler& sh = obj->subset_handler();

	sel.clear();
	bool autoselEnabled = sel.autoselection_enabled();
	sel.enable_autoselection(true);

//	create the vertices
//	create one upfront, the others in a loop
	Vertex* centerVrt = NULL;
	if(fill){
		centerVrt = *grid.create<RegularVertex>();
		aaPos[centerVrt] = center;
	}
	Vertex* firstVrt = *grid.create<RegularVertex>();
	aaPos[firstVrt] = vector3(0, radius, 0);
	VecAdd(aaPos[firstVrt], aaPos[firstVrt], center);
	Vertex* lastVrt = firstVrt;
	for(int i = 1; i < numRimVertices; ++i){
	//	create a new vertex
		number ia = (float)i / (float)numRimVertices;
		Vertex* vNew = *grid.create<RegularVertex>();
		aaPos[vNew] = vector3(sin(2. * PI * ia), cos(2. * PI * ia), 0);
		VecScale(aaPos[vNew], aaPos[vNew], radius);
		VecAdd(aaPos[vNew], aaPos[vNew], center);

	//	create a new triangle or a new edge
		if(fill)
			grid.create<Triangle>(TriangleDescriptor(centerVrt, vNew, lastVrt));
		else
			grid.create<RegularEdge>(EdgeDescriptor(lastVrt, vNew));

	//	prepare the next iteration
		lastVrt = vNew;
	}

//	one triangle / edge is still missing
	if(fill)
		grid.create<Triangle>(TriangleDescriptor(centerVrt, firstVrt, lastVrt));
	else
		grid.create<RegularEdge>(EdgeDescriptor(lastVrt, firstVrt));

//	assign subset
	sh.assign_subset(sel.begin<Vertex>(), sel.end<Vertex>(), subsetInd);
	sh.assign_subset(sel.begin<Edge>(), sel.end<Edge>(), subsetInd);
	sh.assign_subset(sel.begin<Face>(), sel.end<Face>(), subsetInd);

//	restore selector
	sel.enable_autoselection(autoselEnabled);
}


void CreateBox(	Mesh* obj,
				const vector3& boxMin,
				const vector3& boxMax,
				int subsetInd,
				bool fill)
{
	Grid& grid = obj->grid();
	Mesh::position_accessor_t& aaPos = obj->position_accessor();
	Selector& sel = obj->selector();
	SubsetHandler& sh = obj->subset_handler();

	sel.clear();
	bool autoselEnabled = sel.autoselection_enabled();
	sel.enable_autoselection(true);

//	create the vertices
	Vertex* vrts[8];
	for(size_t i = 0; i < 8; ++i)
		vrts[i] = *grid.create<RegularVertex>();

	aaPos[vrts[0]] = vector3(boxMin.x(), boxMin.y(), boxMin.z());
	aaPos[vrts[1]] = vector3(boxMax.x(), boxMin.y(), boxMin.z());
	aaPos[vrts[2]] = vector3(boxMax.x(), boxMax.y(), boxMin.z());
	aaPos[vrts[3]] = vector3(boxMin.x(), boxMax.y(), boxMin.z());
	aaPos[vrts[4]] = vector3(boxMin.x(), boxMin.y(), boxMax.z());
	aaPos[vrts[5]] = vector3(boxMax.x(), boxMin.y(), boxMax.z());
	aaPos[vrts[6]] = vector3(boxMax.x(), boxMax.y(), boxMax.z());
	aaPos[vrts[7]] = vector3(boxMin.x(), boxMax.y(), boxMax.z());

//	we'll use a hexahedron and let it create all faces
	Hexahedron hexa(vrts[0], vrts[1], vrts[2], vrts[3],
					vrts[4],vrts[5], vrts[6], vrts[7]);

//	now create the faces and register them at the grid
	for(size_t i = 0; i < hexa.num_faces(); ++i)
		grid.register_element(hexa.create_face(i));

//	if a volume shall be created, do so now
	if(fill)
		grid.create_by_cloning(&hexa, hexa, NULL);

//	assign subset
	sh.assign_subset(sel.begin<Vertex>(), sel.end<Vertex>(), subsetInd);
	sh.assign_subset(sel.begin<Edge>(), sel.end<Edge>(), subsetInd);
	sh.assign_subset(sel.begin<Face>(), sel.end<Face>(), subsetInd);
	sh.assign_subset(sel.begin<Volume>(), sel.end<Volume>(), subsetInd);

//	restore selector
	sel.enable_autoselection(autoselEnabled);
}


void CreateSphere(	Mesh* obj,
					const vector3& center,
					number radius,
					int numRefinements,
					int subsetInd)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	SubsetHandler& sh = obj->subset_handler();

//	create the sphere
	GenerateIcosphere(grid, center, radius, numRefinements, obj->position_attachment(), &sel);

//	assign subset
	sh.assign_subset(sel.begin<Vertex>(), sel.end<Vertex>(), subsetInd);
	sh.assign_subset(sel.begin<Edge>(), sel.end<Edge>(), subsetInd);
	sh.assign_subset(sel.begin<Face>(), sel.end<Face>(), subsetInd);
}


void CreateTetrahedron(Mesh* obj, int subsetInd, bool fill)
{
	Grid& grid = obj->grid();
	Mesh::position_accessor_t& aaPos = obj->position_accessor();
	Selector& sel = obj->selector();
	SubsetHandler& sh = obj->subset_handler();

	sel.clear();
	bool autoselEnabled = sel.autoselection_enabled();
	sel.enable_autoselection(true);

//	create the vertices
	Vertex* vrts[4];
	for(size_t i = 0; i < 4; ++i)
		vrts[i] = *grid.create<RegularVertex>();

	aaPos[vrts[0]] = vector3(0, 0, 0);
	aaPos[vrts[1]] = vector3(1, 0, 0);
	aaPos[vrts[2]] = vector3(0, 1, 0);
	aaPos[vrts[3]] = vector3(0, 0, 1);

//	we'll use a hexahedron and let it create all faces
	Tetrahedron tet(vrts[0], vrts[1], vrts[2], vrts[3]);

//	now create the faces and register them at the grid
	for(size_t i = 0; i < tet.num_faces(); ++i)
		grid.register_element(tet.create_face(i));

//	if a volume shall be created, do so now
	if(fill)
		grid.create_by_cloning(&tet, tet, NULL);

//	assign subset
	sh.assign_subset(sel.begin<Vertex>(), sel.end<Vertex>(), subsetInd);
	sh.assign_subset(sel.begin<Edge>(), sel.end<Edge>(), subsetInd);
	sh.assign_subset(sel.begin<Face>(), sel.end<Face>(), subsetInd);
	sh.assign_subset(sel.begin<Volume>(), sel.end<Volume>(), subsetInd);

//	restore selector
	sel.enable_autoselection(autoselEnabled);
}


void CreatePyramid(Mesh* obj, int subsetInd, bool fill)
{
	Grid& grid = obj->grid();
	Mesh::position_accessor_t& aaPos = obj->position_accessor();
	Selector& sel = obj->selector();
	SubsetHandler& sh = obj->subset_handler();

	sel.clear();
	bool autoselEnabled = sel.autoselection_enabled();
	sel.enable_autoselection(true);

//	create the vertices
	Vertex* vrts[5];
	for(size_t i = 0; i < 5; ++i)
		vrts[i] = *grid.create<RegularVertex>();

	aaPos[vrts[0]] = vector3(-1, -1, -1);
	aaPos[vrts[1]] = vector3(1, -1, -1);
	aaPos[vrts[2]] = vector3(1, 1, -1);
	aaPos[vrts[3]] = vector3(-1, 1, -1);
	aaPos[vrts[4]] = vector3(0, 0, 1);

//	we'll use a hexahedron and let it create all faces
	Pyramid pyra(vrts[0], vrts[1], vrts[2], vrts[3], vrts[4]);

//	now create the faces and register them at the grid
	for(size_t i = 0; i < pyra.num_faces(); ++i)
		grid.register_element(pyra.create_face(i));

//	if a volume shall be created, do so now
	if(fill)
		grid.create_by_cloning(&pyra, pyra, NULL);

//	assign subset
	sh.assign_subset(sel.begin<Vertex>(), sel.end<Vertex>(), subsetInd);
	sh.assign_subset(sel.begin<Edge>(), sel.end<Edge>(), subsetInd);
	sh.assign_subset(sel.begin<Face>(), sel.end<Face>(), subsetInd);
	sh.assign_subset(sel.begin<Volume>(), sel.end<Volume>(), subsetInd);

//	restore selector
	sel.enable_autoselection(autoselEnabled);
}


void CreatePrism(Mesh* obj, int subsetInd, bool fill)
{
	Grid& grid = obj->grid();
	Mesh::position_accessor_t& aaPos = obj->position_accessor();
	Selector& sel = obj->selector();
	SubsetHandler& sh = obj->subset_handler();

	sel.clear();
	bool autoselEnabled = sel.autoselection_enabled();
	sel.enable_autoselection(true);

//	create the vertices
	Vertex* vrts[6];
	for(size_t i = 0; i < 6; ++i)
		vrts[i] = *grid.create<RegularVertex>();

	aaPos[vrts[0]] = vector3(-1, -1, -1);
	aaPos[vrts[1]] = vector3(1, -1, -1);
	aaPos[vrts[2]] = vector3(1, 1, -1);
	aaPos[vrts[3]] = vector3(-1, -1, 1);
	aaPos[vrts[4]] = vector3(1, -1, 1);
	aaPos[vrts[5]] = vector3(1, 1, 1);

//	we'll use a hexahedron and let it create all faces
	Prism prism(vrts[0], vrts[1], vrts[2], vrts[3], vrts[4], vrts[5]);

//	now create the faces and register them at the grid
	for(size_t i = 0; i < prism.num_faces(); ++i)
		grid.register_element(prism.create_face(i));

//	if a volume shall be created, do so now
	if(fill)
		grid.create_by_cloning(&prism, prism, NULL);

//	assign subset
	sh.assign_subset(sel.begin<Vertex>(), sel.end<Vertex>(), subsetInd);
	sh.assign_subset(sel.begin<Edge>(), sel.end<Edge>(), subsetInd);
	sh.assign_subset(sel.begin<Face>(), sel.end<Face>(), subsetInd);
	sh.assign_subset(sel.begin<Volume>(), sel.end<Volume>(), subsetInd);

//	restore selector
	sel.enable_autoselection(autoselEnabled);

}

}}//	end of namespace
