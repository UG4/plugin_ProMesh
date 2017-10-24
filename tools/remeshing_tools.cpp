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

#include "remeshing_tools.h"
#include "lib_grid/algorithms/subset_util.h"
#include "../EigenForUG4/Eigen/Dense"

namespace ug{
namespace promesh{

void SimplifyPolylines(Mesh* m, number curvatureThreshold){
	Grid& grid = m->grid();
	Selector& sel = m->selector();
	Mesh::position_accessor_t& aaPos = m->position_accessor();

	SimplifyPolylines(grid, sel.begin<Edge>(), sel.end<Edge>(), curvatureThreshold, aaPos);
}

void SimplifySmoothedPolylines(
			Mesh* m,
			number curvatureThreshold,
			number smoothingAlpha,
			int smoothingIterations)
{
	Grid& grid = m->grid();
	Selector& sel = m->selector();
	Mesh::position_accessor_t& aaPos = m->position_accessor();

	SimplifySmoothedPolylines(grid, sel.begin<Edge>(), sel.end<Edge>(),
							  curvatureThreshold, aaPos,
							  smoothingAlpha, smoothingIterations);
}

void ConvertToTriangles(Mesh* obj)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	Mesh::position_accessor_t& aaPos = obj->position_accessor();

	Triangulate(grid, sel.begin<Quadrilateral>(),
				sel.end<Quadrilateral>(), &aaPos);
}

void ConvertToQuadrilaterals(Mesh* obj)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();

	if(sel.num<Face>() > 0){
		ReplaceByQuadrilaterals_FaceBased(
				grid,
				sel.begin<Face>(),
				sel.end<Face>(),
				obj->position_accessor());
	}

	if(sel.num<Edge>() > 0){
		ReplaceByQuadrilaterals_EdgeBased(
				grid,
				sel.begin<Edge>(),
				sel.end<Edge>(),
				obj->position_accessor());
	}
}

void ExtrudeFacesWithTets(Mesh* obj, int fromSi, int toSi, const number factor)
{
	using namespace std;
	Grid& grid = obj->grid();
	Mesh::position_accessor_t& aaPos = obj->position_accessor();
	Selector& sel = obj->selector();
	SubsetHandler& sh = obj->subset_handler();

	sel.clear();
	bool autoselEnabled = sel.autoselection_enabled();
	sel.enable_autoselection(true);

	vector<Face*> faces;
	faces.assign(sh.begin<Face>(fromSi), sh.end<Face>(fromSi));
	for (size_t i = 0; i < faces.size(); i++)
	{
		Face* f = faces.at(i);
		vector3 p1, p2, p3;
		vector3 pnormal, center, top;
		number scale;
		if (f->num_vertices() == 3)
		{
			p1 = aaPos[f->vertex(0)];
			p2 = aaPos[f->vertex(1)];
			p3 = aaPos[f->vertex(2)];
			CalculateTriangleNormal(pnormal, p1, p2, p3);
			RegularVertex* n = *grid.create<RegularVertex>();
			center = CalculateCenter(f, aaPos);
			scale = FaceArea(f, aaPos);
			pnormal *= scale * factor;
			VecAdd(top, center, pnormal);
			aaPos[n] = top;
			Tetrahedron tet(f->vertex(0), f->vertex(1), f->vertex(2), n);
			// assumption: the first face in the tet is the base one
			for (size_t i = 1; i < tet.num_faces(); ++i)
			{
				grid.register_element(tet.create_face(i));
			}
		}
		else
		{
			UG_LOG("Face with more than 3 vertices. Doesn't work yet.");
		}
	}
	sh.assign_subset(sel.begin<Vertex>(), sel.end<Vertex>(), toSi);
	sh.assign_subset(sel.begin<Edge>(), sel.end<Edge>(), toSi);
	sh.assign_subset(sel.begin<Face>(), sel.end<Face>(), toSi);
	sh.assign_subset(sel.begin<Volume>(), sel.end<Volume>(), toSi);

	//	restore selector
	sel.enable_autoselection(autoselEnabled);
	// remove doubles
}

void TriangleFill(Mesh* obj, bool qualityGeneration, number minAngle, int si)
{
	if(minAngle < 0)
		minAngle = 0;
	if(minAngle > 30){
		UG_LOG("WARNING in TriangleFill: Restricting minAngle to 30, since the "
				"algorithm may not terminate for minAngle > 30.\n")
	}

	Grid& grid = obj->grid();
	SubsetHandler& sh = obj->subset_handler();
	Selector& sel = obj->selector();

//	if no edges are selected, nothing can be triangulated
	if(sel.num<Edge>() < 3){
		UG_LOG("ERROR in TriangleFill: A closed outer edge-chain has to be selected.\n");
		return;
	}

//	before triangulating, we'll make sure that no double-edges exist
//	in the current selection.
	RemoveDuplicates(grid, sel.begin<Edge>(), sel.end<Edge>());

	Mesh::position_accessor_t& aaPos = obj->position_accessor();
	AInt aInt;
	grid.attach_to_vertices(aInt);

//	we don't want to select new edges. This would be a problem for
//	delaunay constraints.
	bool autoselEnabled = sel.autoselection_enabled();
	sel.enable_autoselection(false);

//	Collect all new faces in this selector.
	FaceSelector faceSel(grid);
	faceSel.enable_autoselection(true);

	if(!TriangleFill_SweepLine(grid, sel.edges_begin(),
							sel.edges_end(), obj->position_attachment(), aInt,
							&sh, si))
	{
		UG_LOG("TriangleFill_SweepLine failed.\n");

	// ONLY FOR DEBUGGING - BEGIN
	/*
		static int fileCounter = 1;
		string filenamePrefix = "/Users/sreiter/Desktop/failed_sweeplines/failed_sweepline_";
		//string filenamePrefix = "C:/sweep_errors/failed_sweepline_";
		stringstream ss2d, ss3d;
		ss2d << filenamePrefix << "2d_" << fileCounter << ".obj";
		ss3d << filenamePrefix << "3d_" << fileCounter << ".obj";
		++fileCounter;
		//UG_LOG("TriangleFill_SweepLine failed!\n");
		UG_LOG("Saving failed geometries to " << ss2d.str() << " and " << ss3d.str() << endl);
		SaveGridToFile(grid, ss3d.str().c_str(), obj->position_attachment());
	//	perform transformation to 2d and save that too.
		std::vector<vector3> vrts;
		for(VertexIterator iter = grid.vertices_begin();
			iter != grid.vertices_end(); ++iter)
		{
			vrts.push_back(aaPos[*iter]);
		}
		std::vector<vector2> vrts2d(vrts.size());
		TransformPointSetTo2D(&vrts2d.front(), &vrts.front(),
							  vrts.size());

		size_t counter = 0;
		for(VertexIterator iter = grid.vertices_begin();
			iter != grid.vertices_end(); ++iter, ++counter)
		{
			aaPos[*iter] = vector3(vrts2d[counter].x(), vrts2d[counter].y(), 0);
		}

		SaveGridToFile(grid, ss2d.str().c_str(), obj->position_attachment());

		counter = 0;
		for(VertexIterator iter = grid.vertices_begin();
			iter != grid.vertices_end(); ++iter, ++counter)
		{
			aaPos[*iter] = vector3(vrts[counter].x(), vrts[counter].y(), 0);
		}
	*/
	// ONLY FOR DEBUGGING - END

	}

	grid.detach_from_vertices(aInt);

	if(qualityGeneration){
		QualityGridGeneration(grid, faceSel.begin(), faceSel.end(),
					 aaPos, minAngle, IsSelected(sel));
	}

	sel.clear ();
	sel.select (faceSel.begin(), faceSel.end());
	CloseSelection (sel);
	CopySubsetIndicesToSides (sh, sel, true);

	sel.enable_autoselection(autoselEnabled);
}

void Retriangulate(Mesh* obj, number minAngle)
{
	Grid& g = obj->grid();
	Selector& sel = obj->selector();
	SubsetHandler& creases = obj->crease_handler();

	Mesh::position_accessor_t& aaPos = obj->position_accessor();

	QualityGridGeneration(g, sel.begin<Triangle>(), sel.end<Triangle>(),
						  aaPos, minAngle, IsNotInSubset(creases, -1));

	CloseSelection (sel);
	CopySubsetIndicesToSides (obj->subset_handler(), sel, true);
}

void AdjustEdgeLength(
			Mesh* obj,
			number minEdgeLen,
			number maxEdgeLen,
			int numIterations,
			bool adaptive,
			bool automarkBoundaries)
{
	Grid& grid = obj->grid();
	SubsetHandler& shCrease = obj->crease_handler();

	if(automarkBoundaries){
		for(EdgeIterator iter = grid.begin<Edge>();
			iter != grid.end<Edge>(); ++iter)
		{
			if(IsBoundaryEdge2D(grid, *iter))
				shCrease.assign_subset(*iter, REM_CREASE);
		}
	}

	AdjustEdgeLength(grid, shCrease, minEdgeLen, maxEdgeLen,
						 numIterations, true, adaptive);

	CopySubsetIndicesToSides (obj->subset_handler(), true);
}

void AdjustEdgeLengthExtended(
			Mesh* obj,
			number minEdgeLen,
			number maxEdgeLen,
			number approximation,
			number triQuality,
			int numIterations,
			bool automarkBoundaries)
{
	Grid& grid = obj->grid();
	SubsetHandler& shCrease = obj->crease_handler();

	if(automarkBoundaries){
		for(EdgeIterator iter = grid.begin<Edge>();
			iter != grid.end<Edge>(); ++iter)
		{
			if(IsBoundaryEdge2D(grid, *iter))
				shCrease.assign_subset(*iter, REM_CREASE);
		}
	}

	AdjustEdgeLengthDesc desc;
	desc.minEdgeLen = minEdgeLen;
	desc.maxEdgeLen = maxEdgeLen;
	desc.approximation = approximation;
	desc.triQuality = triQuality;
	AdjustEdgeLength(grid, shCrease, desc, numIterations);

	CopySubsetIndicesToSides (obj->subset_handler(), true);
}

void AdaptSurfaceToCylinder(Mesh* obj, number radius, number threshold)
{
	using namespace std;
	Grid& g = obj->grid();
	Selector& sel = obj->selector();
	Mesh::position_accessor_t& aaPos = obj->position_accessor();

//	store all source-vertices in a list
	vector<Vertex*> vrts;
	vrts.assign(sel.begin<Vertex>(), sel.end<Vertex>());

//	iterate over selected vertices
	for(vector<Vertex*>::iterator iter = vrts.begin();
		iter != vrts.end(); ++iter)
	{
		Vertex* vrt = *iter;
		vector3 n;
		CalculateVertexNormal(n, g, vrt, aaPos);

		if(!ug::AdaptSurfaceGridToCylinder(sel, g, vrt, n, radius, threshold))
		{
			UG_LOG("AdaptSurfaceGridToCylinder failed for the vertex at " << aaPos[vrt] << "\n");
		}
	}
}

void ConvertToTetrahedra(Mesh* obj)
{
	Selector& sel = obj->selector();
	ConvertToTetrahedra(obj->grid(),
						sel.begin<Volume>(),
						sel.end<Volume>());

}

void Tetrahedralize(
			Mesh* obj,
			number quality,
			bool preserveOuter,
			bool preserveAll,
			bool separateVolumes,
			bool appendSubsetsAtEnd,
			int verbosity)
{
	Grid& grid = obj->grid();
	SubsetHandler& sh = obj->subset_handler();
	UG_LOG("tetrahedralizing using 'tetgen' by Hang Si... ");
UG_LOG("<dbg> \n");
UG_LOG("<dbg> running ug::Tetrahedralize\n");
	ug::Tetrahedralize(grid, sh, quality, preserveOuter, preserveAll,
					   obj->position_attachment(), verbosity);
	UG_LOG("done. Created " << grid.num<Tetrahedron>() << " tetrahedra.\n");

	int oldNumSubsets = sh.num_subsets();
	if(separateVolumes){
		SeparateSubsetsByLowerDimSubsets<Volume>(grid, sh,
														 appendSubsetsAtEnd);
	}
	else if(appendSubsetsAtEnd){
	//todo:	only assign newly generated tetrahedrons.
		sh.assign_subset(grid.begin<Tetrahedron>(),
						 grid.end<Tetrahedron>(), sh.num_subsets());
	}

//	assign a subset name
	for(int i = oldNumSubsets; i < sh.num_subsets(); ++i)
		sh.subset_info(i).name = "tetrahedra";

	CopySubsetIndicesToSides (sh, true);
}

void AssignVolumeConstraints(Mesh* obj, number volConstraint)
{
	Selector& sel = obj->selector();
	Mesh::volume_constraint_accessor_t& aaVolCon = obj->volume_constraint_accessor();

	for(Selector::traits<Volume>::iterator iter = sel.begin<Volume>();
		iter != sel.end<Volume>(); ++iter)
	{
		aaVolCon[*iter] = volConstraint;
	}
}

void ClearVolumeConstraints(Mesh* obj)
{
	obj->clear_volume_constraints();
}

void Retetrahedralize(
			Mesh* obj,
			number quality,
			bool preserveOuter,
			bool preserveAll,
			bool applyVolumeConstraint,
			int verbosity)
{
	UG_LOG("retetrahedralizing using 'tetgen' by Hang Si... ");
	ug::Retetrahedralize(obj->grid(),
					obj->subset_handler(),
					obj->volume_constraint_attachment(),
					quality,
					preserveOuter, preserveAll,
					obj->position_attachment(),
					applyVolumeConstraint,
					verbosity);
	
	CopySubsetIndicesToSides (obj->subset_handler(), true);

	UG_LOG("done.\n");
}

void Duplicate(Mesh* obj, const vector3& offset, bool deselectOld, bool selectNew)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	ug::Duplicate(grid, sel, offset, obj->position_attachment(), deselectOld, selectNew);
}

void ExtrudeAndMove(
			Mesh* obj,
			const vector3& totalDir,
			int numSteps,
			bool createFaces,
			bool createVolumes)
{
	using namespace std;
	if(numSteps < 1)
		return;

	vector3 stepDir;
	VecScale(stepDir, totalDir, 1./(float)numSteps);

	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	SubsetHandler& sh = obj->subset_handler();

	SelectAssociatedGridObjects(sel);

	vector<Vertex*> vrts;
	vrts.assign(sel.vertices_begin(), sel.vertices_end());
	vector<Edge*> edges;
	edges.assign(sel.edges_begin(), sel.edges_end());
	vector<Face*> faces;
	faces.assign(sel.faces_begin(), sel.faces_end());

	uint extrusionOptions = 0;
	if(createFaces)
		extrusionOptions |= EO_CREATE_FACES;
	if(createVolumes)
		extrusionOptions |= EO_CREATE_VOLUMES;

//	we use sel to collect the newly created volumes
	bool strictInheritanceEnabled = sh.strict_inheritance_enabled();
	sh.enable_strict_inheritance(false);
//	mark all elements that were already in the selector.
	for(int i = 0; i < numSteps; ++i)
	{
		Extrude(grid, &vrts, &edges, &faces, stepDir,
					extrusionOptions, obj->position_attachment());
	}

	sh.enable_strict_inheritance(strictInheritanceEnabled);

//	select faces, edges and vertices from the new top-layer.
	sel.clear<Vertex>();
	sel.clear<Edge>();
	sel.clear<Face>();
	sel.clear<Volume>();
	sel.select(vrts.begin(), vrts.end());
	sel.select(edges.begin(), edges.end());
	sel.select(faces.begin(), faces.end());
}

void ExtrudeAndScale(
			Mesh* obj,
			number totalScale,
			bool scaleAroundPivot,
			int numSteps,
			bool createFaces,
			bool createVolumes)
{
	using namespace std;
	if(numSteps < 1)
		return;

	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	SubsetHandler& sh = obj->subset_handler();

	SelectAssociatedGridObjects(sel);

	vector<Vertex*> vrts;
	vrts.assign(sel.vertices_begin(), sel.vertices_end());
	vector<Edge*> edges;
	edges.assign(sel.edges_begin(), sel.edges_end());
	vector<Face*> faces;
	faces.assign(sel.faces_begin(), sel.faces_end());

	uint extrusionOptions = 0;
	if(createFaces)
		extrusionOptions |= EO_CREATE_FACES;
	if(createVolumes)
		extrusionOptions |= EO_CREATE_VOLUMES;

	Mesh::position_accessor_t aaPos = obj->position_accessor();
	vector3 center;
	if(scaleAroundPivot)
		center = obj->pivot();
	else
		center = CalculateCenter(vrts.begin(), vrts.end(), aaPos);

	vector<vector3>	from;
	vector<vector3>	stepOffsets;
	from.reserve(vrts.size());
	stepOffsets.reserve(vrts.size());
	for(size_t ivrt = 0; ivrt < vrts.size(); ++ivrt){
		from.push_back(aaPos[vrts[ivrt]]);
		vector3 dir;
		VecSubtract(dir, from.back(), center);
		dir *= totalScale;
	//	create a directional vector relative to from
		dir += center;
		VecSubtract(dir, dir, from.back());
		dir *= 1. / (number)numSteps;
		stepOffsets.push_back(dir);
	}

//	we use sel to collect the newly created volumes
	bool strictInheritanceEnabled = sh.strict_inheritance_enabled();
	sh.enable_strict_inheritance(false);
//	mark all elements that were already in the selector.
	for(int i = 0; i < numSteps; ++i)
	{
		Extrude(grid, &vrts, &edges, &faces, vector3(0, 0, 0),
					extrusionOptions, obj->position_attachment());

		for(size_t ivrt = 0; ivrt < vrts.size(); ++ivrt){
			VecScaleAdd(aaPos[vrts[ivrt]], 1, from[ivrt], i+1, stepOffsets[ivrt]);
		}
	}

	sh.enable_strict_inheritance(strictInheritanceEnabled);

//	select faces, edges and vertices from the new top-layer.
	sel.clear<Vertex>();
	sel.clear<Edge>();
	sel.clear<Face>();
	sel.clear<Volume>();
	sel.select(vrts.begin(), vrts.end());
	sel.select(edges.begin(), edges.end());
	sel.select(faces.begin(), faces.end());
}


void ExtrudeAlongNormal(
			Mesh* obj,
			number totalLength,
			int numSteps,
			bool createFaces,
			bool createVolumes)
{
	using namespace std;
	if(numSteps < 1)
		return;

	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	SubsetHandler& sh = obj->subset_handler();

	SelectAssociatedGridObjects(sel);

	vector<Vertex*> vrts;
	vrts.assign(sel.vertices_begin(), sel.vertices_end());
	vector<Edge*> edges;
	edges.assign(sel.edges_begin(), sel.edges_end());
	vector<Face*> faces;
	faces.assign(sel.faces_begin(), sel.faces_end());

	uint extrusionOptions = 0;
	if(createFaces)
		extrusionOptions |= EO_CREATE_FACES;
	if(createVolumes)
		extrusionOptions |= EO_CREATE_VOLUMES;

	Mesh::position_accessor_t aaPos = obj->position_accessor();

	vector<vector3>	from;
	vector<vector3>	stepOffsets;
	from.reserve(vrts.size());
	stepOffsets.reserve(vrts.size());

	Grid::face_traits::secure_container	assFaces;
	Grid::edge_traits::secure_container	assEdges;

	for(size_t ivrt = 0; ivrt < vrts.size(); ++ivrt){
		Vertex* vrt = vrts[ivrt];
		from.push_back(aaPos[vrt]);
		vector3 n(0, 0, 0);
		grid.associated_elements(assFaces, vrt);
		if(assFaces.size() > 0){
			vector3 nSel(0, 0, 0);
			for(size_t iface = 0; iface < assFaces.size(); ++iface){
				Face* f = assFaces[iface];
				vector3 fn;
				CalculateNormal(fn, f, aaPos);
				n += fn;
				if(sel.is_selected(f))
					nSel += fn;
			}
			if(VecLengthSq(nSel) > 0)
				n = nSel;
		}
		else{
			grid.associated_elements(assEdges, vrt);
			vector3 nSel(0, 0, 0);
			for(size_t iedge = 0; iedge < assEdges.size(); ++iedge){
				Edge* e = assEdges[iedge];
				vector3 tmpN;
				VecSubtract(tmpN, aaPos[e->vertex(1)], aaPos[e->vertex(0)]);
				number t = tmpN.x();
				tmpN.x() = tmpN.y();
				tmpN.y() = -t;
				tmpN.z() = 0;
				VecNormalize(tmpN, tmpN);
				n += tmpN;
				if(sel.is_selected(e))
					nSel += tmpN;
			}
			if(VecLengthSq(nSel) > 0)
				n = nSel;
		}

		VecNormalize(n, n);

		n *= totalLength / (number)numSteps;
		stepOffsets.push_back(n);
	}

//	we use sel to collect the newly created volumes
	bool strictInheritanceEnabled = sh.strict_inheritance_enabled();
	sh.enable_strict_inheritance(false);
//	mark all elements that were already in the selector.
	for(int i = 0; i < numSteps; ++i)
	{
		Extrude(grid, &vrts, &edges, &faces, vector3(0, 0, 0),
					extrusionOptions, obj->position_attachment());

		for(size_t ivrt = 0; ivrt < vrts.size(); ++ivrt){
			VecScaleAdd(aaPos[vrts[ivrt]], 1, from[ivrt], i+1, stepOffsets[ivrt]);
		}
	}

	sh.enable_strict_inheritance(strictInheritanceEnabled);

//	select faces, edges and vertices from the new top-layer.
	sel.clear<Vertex>();
	sel.clear<Edge>();
	sel.clear<Face>();
	sel.clear<Volume>();
	sel.select(vrts.begin(), vrts.end());
	sel.select(edges.begin(), edges.end());
	sel.select(faces.begin(), faces.end());
}


void ExtrudeToThickness(Mesh* obj,
                        number thickness,
                        int numSteps,
						bool createFaces,
						bool createVolumes)
{
	typedef Eigen::Matrix<number, 2, 1>	EVector2;
	typedef Eigen::Matrix<number, 3, 1>	EVector3;
	typedef Eigen::Matrix<number, Eigen::Dynamic, 1>				EVectorX;
	typedef Eigen::Matrix<number, Eigen::Dynamic, Eigen::Dynamic> 	EMatrixXY;

	using namespace std;
	if(numSteps < 1)
		return;

	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	SubsetHandler& sh = obj->subset_handler();

	SelectAssociatedGridObjects(sel);

	vector<Vertex*> vrts;
	vrts.assign(sel.vertices_begin(), sel.vertices_end());
	vector<Edge*> edges;
	edges.assign(sel.edges_begin(), sel.edges_end());
	vector<Face*> faces;
	faces.assign(sel.faces_begin(), sel.faces_end());

	uint extrusionOptions = 0;
	if(createFaces)
		extrusionOptions |= EO_CREATE_FACES;
	if(createVolumes)
		extrusionOptions |= EO_CREATE_VOLUMES;

	Mesh::position_accessor_t aaPos = obj->position_accessor();

	vector<vector3>	from;
	vector<vector3>	stepOffsets;
	vector<Edge*> selEdges;
	vector<Face*> selFaces;

	from.reserve(vrts.size());
	stepOffsets.reserve(vrts.size());

	Grid::volume_traits::secure_container	assVols;
	Grid::face_traits::secure_container	assFaces;
	Grid::edge_traits::secure_container	assEdges;

	EMatrixXY A;
	EVectorX rhs, np;

	bool extrude2d = false;
	bool extrude3d = false;


	for(size_t ivrt = 0; ivrt < vrts.size(); ++ivrt){
		Vertex* vrt = vrts[ivrt];
		const vector3& p = aaPos[vrt];
		const EVector3 ep (p[0], p[1], p[2]);

		from.push_back(p);

		grid.associated_elements(assFaces, vrt);
		selFaces.clear();
		for(size_t iface = 0; iface < assFaces.size(); ++iface){
			if(sel.is_selected(assFaces[iface]))
				selFaces.push_back(assFaces[iface]);
		}

		if(selFaces.size() > 0){
		//	do a 3d extrusion from selected faces
			extrude3d = true;
			A.resize (selFaces.size(), 3);
			rhs.resize (selFaces.size());

			for(size_t iface = 0; iface < selFaces.size(); ++iface){
				Face* f = selFaces[iface];
				vector3 fn;
				CalculateNormal(fn, f, aaPos);

			//	make sure that the normal does not point inside an existing volume
				grid.associated_elements(assVols, f);
				if(assVols.size() == 1){
					vector3 tmp = CalculateCenter (assVols[0], aaPos);
					tmp -= CalculateCenter (f, aaPos);
					if(VecDot(fn, tmp) > 0)
						fn *= -1.f;
				}

				EVector3 en (fn[0], fn[1], fn[2]);

				A.row(iface) = en;
				rhs[iface] = thickness;
			}

		//	solve with singular value decomposition. Result is solution in least-squares sense.
			Eigen::JacobiSVD<EMatrixXY> svd (A, Eigen::ComputeThinU | Eigen::ComputeThinV);
		//	adjustment of this threshold is important!
			svd.setThreshold (0.00000001);
			np = svd.solve (rhs);

			np /= (number)numSteps;
			stepOffsets.push_back (vector3(np[0], np[1], np[2]));
		}

		else {
		//	if there are no associated faces for this vertex, we try to perform
		//	a 2d extrusion of associated edges
			grid.associated_elements(assEdges, vrt);
			selEdges.clear();
			for(size_t iedge = 0; iedge < assEdges.size(); ++iedge){
				if(sel.is_selected(assEdges[iedge]))
					selEdges.push_back(assEdges[iedge]);
			}

			if(selEdges.size() > 0){
			//	do a 2d extrusion from selected edges
				extrude2d = true;
				A.resize (selEdges.size(), 2);
				rhs.resize (selEdges.size());

				for(size_t iedge = 0; iedge < selEdges.size(); ++iedge){
					Edge* e = selEdges[iedge];
					vector3 n;
					VecSubtract(n, aaPos[e->vertex(1)], aaPos[e->vertex(0)]);
					number t = n.x();
					n.x() = n.y();
					n.y() = -t;
					n.z() = 0;
					VecNormalize(n, n);

				//	make sure that the normal does not point inside an existing face
					grid.associated_elements(assFaces, e);
					if(assFaces.size() == 1){
						vector3 tmp = CalculateCenter (assFaces[0], aaPos);
						tmp -= CalculateCenter (e, aaPos);
						if(VecDot(n, tmp) > 0)
							n *= -1.f;
					}

					EVector2 en (n[0], n[1]);

					A.row(iedge) = en;
					rhs[iedge] = thickness;
				}

			//	solve with singular value decomposition. Result is solution in least-squares sense.
				Eigen::JacobiSVD<EMatrixXY> svd (A, Eigen::ComputeThinU | Eigen::ComputeThinV);
			//	adjustment of this threshold is important!
				svd.setThreshold (0.00000001);
				np = svd.solve (rhs);

				np /= (number)numSteps;
				stepOffsets.push_back (vector3(np[0], np[1], 0));
			}
			else {
				stepOffsets.push_back (vector3(0, 0, 0));
			}
		}
	}

//	we use sel to collect the newly created volumes
	bool strictInheritanceEnabled = sh.strict_inheritance_enabled();
	sh.enable_strict_inheritance(false);
//	mark all elements that were already in the selector.
	for(int i = 0; i < numSteps; ++i)
	{
		Extrude(grid, &vrts, &edges, &faces, vector3(0, 0, 0),
					extrusionOptions, obj->position_attachment());

		for(size_t ivrt = 0; ivrt < vrts.size(); ++ivrt){
			VecScaleAdd(aaPos[vrts[ivrt]], 1, from[ivrt], i+1, stepOffsets[ivrt]);
		}
	}

	sh.enable_strict_inheritance(strictInheritanceEnabled);

//	select faces, edges and vertices from the new top-layer.
	sel.clear<Vertex>();
	sel.clear<Edge>();
	sel.clear<Face>();
	sel.clear<Volume>();
	sel.select(vrts.begin(), vrts.end());
	sel.select(edges.begin(), edges.end());
	if(extrude3d)
		sel.select(faces.begin(), faces.end());
}


void ExtrudeCylinders(Mesh* obj, number height, number radius, number snapThreshold)
{
	using namespace std;
	Grid& g = obj->grid();
	Selector& sel = obj->selector();
	SubsetHandler& sh = obj->subset_handler();
	Mesh::position_accessor_t& aaPos = obj->position_accessor();

//	store all source-vertices in a list
	vector<Vertex*> vrts;
	vrts.assign(sel.begin<Vertex>(), sel.end<Vertex>());

//	iterate over selected vertices
	for(vector<Vertex*>::iterator iter = vrts.begin();
		iter != vrts.end(); ++iter)
	{
		Vertex* vrt = *iter;
		vector3 n;
		CalculateVertexNormal(n, g, vrt, aaPos);

		int numSubs = sh.num_subsets();
		if(!ExtrudeCylinder(g, sh, vrt, n, height, radius, snapThreshold,
							aaPos, numSubs, numSubs + 1))
		{
			UG_LOG("Cylinder-Extrude failed for the vertex at " << aaPos[vrt] << "\n");
		}
	}
}

void CreateShrinkGeometry(Mesh* obj, number scale)
{
	using namespace std;

	Grid& g = obj->grid();
	vector<Volume*>	vols(g.begin<Volume>(), g.end<Volume>());
	vector<Face*>	faces(g.begin<Face>(), g.end<Face>());
	vector<Edge*>	edges(g.begin<Edge>(), g.end<Edge>());

	CreateShrinkElements(obj, scale, vols.begin(), vols.end());
	CreateShrinkElements(obj, scale, faces.begin(), faces.end());
	CreateShrinkElements(obj, scale, edges.begin(), edges.end());
}

void ReplaceLowValenceVertices(Mesh* obj, number maxSquaredHeightToBaseAreaRatio)
{
	Selector& sel = obj->selector();
	ug::ReplaceLowValenceVertices(obj->grid(), sel.begin<Vertex>(), sel.end<Vertex>(),
						   maxSquaredHeightToBaseAreaRatio, obj->position_accessor());
}

void ReplaceValence3Vertices(Mesh* obj, number maxSquaredHeightToBaseAreaRatio)
{
	Selector& sel = obj->selector();
	ug::ReplaceValence3Vertices(obj->grid(), sel.begin<Vertex>(), sel.end<Vertex>(),
						   maxSquaredHeightToBaseAreaRatio, obj->position_accessor());
}

void MeshLayers(Mesh* m, const RasterLayers& layers)
{
	MeshLayers(m->grid(), layers, m->position_accessor(), &m->subset_handler());
}

void MeshLayerBoundaries(Mesh* m, const RasterLayers& layers)
{
	MeshLayerBoundaries(m->grid(), layers, m->position_accessor(), &m->subset_handler());
}


void ProjectToLayer(Mesh* obj, RasterLayers& layers, int layerIndex){
	ProjectToLayer(obj->grid(), layers, layerIndex, obj->position_accessor());
}

void ProjectToTopLayer(Mesh* obj, RasterLayers& layers){
	UG_COND_THROW(layers.size() == 0, "No top-layer exists in ProjectToTopLayer");
	ProjectToLayer(obj->grid(), layers, layers.size() - 1, obj->position_accessor());
}

void ExtrudeLayers(Mesh* obj, RasterLayers& layers, bool allowForTetsAndPyras){
	ExtrudeLayers(obj->grid(), layers, obj->position_accessor(),
				  obj->subset_handler(), allowForTetsAndPyras);
}

void ExtrudeLayersAndAddProjector(
			Mesh* obj,
			SPRasterLayers layers,
			bool allowForTetsAndPyras)
{
	SPRasterLayersProjector proj
		= make_sp(new RasterLayersProjector(obj->geometry(), layers));

	ANumber aRelZ = proj->rel_z_attachment();
	ExtrudeLayers(obj->grid(), *layers, obj->position_accessor(),
				  obj->subset_handler(), allowForTetsAndPyras,
				  &aRelZ);

	obj->projection_handler().set_default_projector(proj);
}

void SnapToHorizontalRaster(Mesh* obj, SPRasterLayers layers)
{
	SnapToHorizontalRaster(obj->grid(), *layers, obj->position_accessor());
}

void CSGFaceOperation(
			Mesh* obj,
			CSGOperation op,
			int subsetIndex0,
			int subsetIndex1,
			number snapThreshold)
{
	using namespace std;
//	both subsets have to be homeomorphic to the sphere!

	Grid& grid = obj->grid();
	SubsetHandler& sh = obj->subset_handler();
	Selector& sel = obj->selector();

	Mesh::position_attachment_t aPos = obj->position_attachment();
	Mesh::position_accessor_t aaPos = obj->position_accessor();

//	in order to create one octree for each subset, we'll copy the
//	triangles of each subset to a fresh grid and create an octree
//	on that grid.

	typedef lg_ntree<3, 3, Face> octree_t;
	octree_t octree[2];
	Grid octreeGrid[2];
	Mesh::position_accessor_t aaPosOT[2];
	int si[2] = {subsetIndex0, subsetIndex1};
	for(int i = 0; i < 2; ++i){
		ug::Triangulate(grid, sh.begin<Quadrilateral>(si[i]),
						sh.end<Quadrilateral>(si[i]),
						&aaPos);
		sel.clear();
		sel.select(sh.begin<Triangle>(si[i]), sh.end<Triangle>(si[i]));
		
		Grid& og = octreeGrid[i];
		og.attach_to_vertices(aPos);
		aaPosOT[i].access(og, aPos);
		CopySelection(sel, og, aaPos, aaPosOT[i]);
		

		octree[i].set_grid(og, aPos);
		octree[i].create_tree(og.begin<Triangle>(), og.end<Triangle>());
	}

	sel.clear();
	SelectSubsetElements<Face>(sel, sh, si[0]);
	SelectSubsetElements<Face>(sel, sh, si[1]);

	ResolveSelfIntersections(obj, snapThreshold);

	
//	remove unnecessary points/edges/faces
	sel.clear();
	std::vector<RayElemIntersectionRecord<Face*> > intersections;

	int remainder[2] = {
		(op == CSG_UNION || op == CSG_DIFFERENCE) ? 1 : 0,
		(op == CSG_UNION) ? 1 : 0
	};

	for(int isub = 0; isub < 2; ++isub){
		const int thisSub = si[isub];
		octree_t& otree = octree[(isub + 1) % 2];

		for(FaceIterator iface = sh.begin<Face>(thisSub);
			iface != sh.end<Face>(thisSub); ++iface)
		{
			Face* f = *iface;
			vector3 center = CalculateCenter(f, aaPos);

		//	check if the center lies inside the other object.
		//	do this by shooting a ray from the center into an arbitrary direction
		//	and count the number of intersections.
			const int maxNumTries = 5;
			for(int itry = 0; itry < maxNumTries; ++itry){
				vector3 dir(urand<number>(-1, 1), urand<number>(-1, 1), urand<number>(-1, 1));
				if(VecLengthSq(dir) == 0)
					continue;
				VecNormalize(dir, dir);

				intersections.clear();

				RayElementIntersections(intersections, otree, center, dir);

				bool repeatTest = false;
				int numPositiveIntersections = 0;
				for(size_t is = 0; is < intersections.size(); ++is){
					RayElemIntersectionRecord<Face*>& rec = intersections[is];
				//	todo:	Check the intersection records to make sure
				//			that the result is reliable. If not, continue
				//			with the next try.
					if(rec.smin == rec.smax){
						if(rec.smin > SMALL)
							++numPositiveIntersections;
						else if(rec.smin > -SMALL){
						//	we most likely have a coplanar face here. In this case
						//	we'll set the number of positive intersections to 1.
							numPositiveIntersections = 1;
							break;
						}
					}
					else{
						numPositiveIntersections = 0;
						repeatTest = true;
						break;
					}
				}

				if(!repeatTest){
					if(numPositiveIntersections % 2 == remainder[isub]){
						sel.select(f);
					}

					break;
				}
			}
		}
	}

	EraseSelectedElements(obj, true, true, false);

	sel.clear();
	SelectSubsetElements<Face>(sel, sh, si[0]);
	SelectSubsetElements<Face>(sel, sh, si[1]);
}


void CSGFaceUnion(
			Mesh* obj,
			int subsetIndex0,
			int subsetIndex1,
			number snapThreshold)
{
	CSGFaceOperation(obj, CSG_UNION, subsetIndex0, subsetIndex1, snapThreshold);
}

void CSGFaceIntersection(
			Mesh* obj,
			int subsetIndex0,
			int subsetIndex1,
			number snapThreshold)
{
	CSGFaceOperation(obj, CSG_INTERSECTION, subsetIndex0, subsetIndex1, snapThreshold);
}

void CSGFaceDifference(
			Mesh* obj,
			int subsetIndex0,
			int subsetIndex1,
			number snapThreshold)
{
	CSGFaceOperation(obj, CSG_DIFFERENCE, subsetIndex0, subsetIndex1, snapThreshold);
}

}}//	end of namespace
