/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__remeshing_tools__
#define __H__UG__remeshing_tools__

#include <vector>
#include "../mesh.h"
#include "lib_grid/algorithms/duplicate.h"
#include "lib_grid/algorithms/remove_duplicates_util.h"
#include "lib_grid/algorithms/extrusion/extrusion.h"
#include "lib_grid/algorithms/grid_generation/horizontal_layers_mesher.h"
#include "lib_grid/algorithms/grid_generation/tetrahedralization.h"
#include "lib_grid/algorithms/grid_generation/triangle_fill_sweep_line.h"
#include "lib_grid/algorithms/remeshing/delaunay_triangulation.h"
#include "lib_grid/algorithms/remeshing/grid_adaption.h"
#include "lib_grid/algorithms/remeshing/edge_length_adjustment.h"
#include "lib_grid/algorithms/remeshing/edge_length_adjustment_extended.h"
#include "lib_grid/algorithms/remeshing/simplification.h"
#include "lib_grid/algorithms/remeshing/simplify_polychain.h"

namespace ug{
namespace promesh{

/// \addtogroup promesh
/// \{

inline void SimplifyPolylines(Mesh* m, number curvatureThreshold){
	Grid& grid = m->grid();
	Selector& sel = m->selector();
	Mesh::position_accessor_t& aaPos = m->position_accessor();

	SimplifyPolylines(grid, sel.begin<Edge>(), sel.end<Edge>(), curvatureThreshold, aaPos);
}

inline void SimplifySmoothedPolylines(Mesh* m, number curvatureThreshold,
									  number smoothingAlpha, int smoothingIterations)
{
	Grid& grid = m->grid();
	Selector& sel = m->selector();
	Mesh::position_accessor_t& aaPos = m->position_accessor();

	SimplifySmoothedPolylines(grid, sel.begin<Edge>(), sel.end<Edge>(),
							  curvatureThreshold, aaPos,
							  smoothingAlpha, smoothingIterations);
}

inline void ConvertToTriangles(Mesh* obj)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	Mesh::position_accessor_t& aaPos = obj->position_accessor();

	Triangulate(grid, sel.begin<Quadrilateral>(),
				sel.end<Quadrilateral>(), &aaPos);
}

inline void ExtrudeFacesWithTets(Mesh* obj, int fromSi, int toSi, const number factor)
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

inline void TriangleFill(Mesh* obj, bool qualityGeneration, number minAngle, int si)
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

	sel.enable_autoselection(autoselEnabled);
}

inline void Retriangulate(Mesh* obj, number minAngle)
{
	Grid& g = obj->grid();
	Selector& sel = obj->selector();
	SubsetHandler& creases = obj->crease_handler();

	Mesh::position_accessor_t& aaPos = obj->position_accessor();

	QualityGridGeneration(g, sel.begin<Triangle>(), sel.end<Triangle>(),
						  aaPos, minAngle, IsNotInSubset(creases, -1));
}

inline void AdjustEdgeLength(Mesh* obj, number minEdgeLen, number maxEdgeLen,
					  int numIterations, bool adaptive, bool automarkBoundaries)
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
}

inline void AdjustEdgeLengthExtended(Mesh* obj, number minEdgeLen, number maxEdgeLen,
								  number approximation, number triQuality,
								  int numIterations, bool automarkBoundaries)
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
}

inline void AdaptSurfaceToCylinder(Mesh* obj, number radius, number threshold)
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

inline void ConvertToTetrahedra(Mesh* obj)
{
	ConvertToTetrahedra(obj->grid(),
						obj->selector().begin<Volume>(),
						obj->selector().end<Volume>());
}

inline void Tetrahedralize(Mesh* obj, number quality, bool preserveOuter, bool preserveAll,
					bool separateVolumes, bool appendSubsetsAtEnd, int verbosity)
{
	Grid& grid = obj->grid();
	SubsetHandler& sh = obj->subset_handler();
	UG_LOG("tetrahedralizing using 'tetgen' by Hang Si... ");
	ug::Tetrahedralize(grid, sh, quality, preserveOuter, preserveAll,
					   obj->position_attachment(), verbosity);
	UG_LOG("done. Created " << grid.num<Tetrahedron>() << " tetrahedrons.\n");

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
		sh.subset_info(i).name = "tetrahedrons";
}

inline void AssignVolumeConstraints(Mesh* obj, number volConstraint)
{
	Selector& sel = obj->selector();
	Mesh::volume_constraint_accessor_t& aaVolCon = obj->volume_constraint_accessor();

	for(Selector::traits<Volume>::iterator iter = sel.begin<Volume>();
		iter != sel.end<Volume>(); ++iter)
	{
		aaVolCon[*iter] = volConstraint;
	}
}

inline void ClearVolumeConstraints(Mesh* obj)
{
	obj->clear_volume_constraints();
}

inline void Retetrahedralize(Mesh* obj, number quality, bool preserveOuter,
					  bool preserveAll, bool applyVolumeConstraint, int verbosity)
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
	UG_LOG("done.\n");
}

inline void Duplicate(Mesh* obj, const vector3& offset, bool deselectOld, bool selectNew)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	ug::Duplicate(grid, sel, offset, obj->position_attachment(), deselectOld, selectNew);
}

inline void ExtrudeAndMove(Mesh* obj, const vector3& totalDir, int numSteps,
			 bool createFaces, bool createVolumes)
{
	using namespace std;
	if(numSteps < 1)
		return;

	vector3 stepDir;
	VecScale(stepDir, totalDir, 1./(float)numSteps);

	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	SubsetHandler& sh = obj->subset_handler();

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

inline void ExtrudeAndScale(Mesh* obj, number totalScale, bool scaleAroundPivot,
							int numSteps, bool createFaces, bool createVolumes)
{
	using namespace std;
	if(numSteps < 1)
		return;

	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	SubsetHandler& sh = obj->subset_handler();

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


inline void ExtrudeAlongNormal(Mesh* obj, number totalLength,
							   int numSteps, bool createFaces, bool createVolumes)
{
	using namespace std;
	if(numSteps < 1)
		return;

	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	SubsetHandler& sh = obj->subset_handler();

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
	Mesh::normal_accessor_t aaNorm = obj->normal_accessor();

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
			for(size_t iface = 0; iface < assFaces.size(); ++iface){
				n += aaNorm[assFaces[iface]];
			}
		}
		else{
			grid.associated_elements(assEdges, vrt);
			for(size_t iedge = 0; iedge < assEdges.size(); ++iedge){
				Edge* e = assEdges[iedge];
				vector3 tmpN;
				VecSubtract(tmpN, aaPos[e->vertex(1)], aaPos[e->vertex(0)]);
				number t = tmpN.x();
				tmpN.x() = tmpN.y();
				tmpN.y() = -t;
				tmpN.z() = 0;
				VecAdd(n, n, tmpN);
			}
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


inline void ExtrudeCylinders(Mesh* obj, number height, number radius, number snapThreshold)
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

/**	For each element of type TElem in obj this method creates a new element with
 * separate corners. The new element will be scaled by 'scale'.
 * All original elements will be deleted before the method terminates.
 * It should thus be called for 'Volumes' first, then for 'Faces' and
 * finally for 'Edges'.*/
template <class TElemIter>
inline void CreateShrinkElements(Mesh* obj, number scale,
						  TElemIter elemsBegin, TElemIter elemsEnd)
{
	using namespace std;
	typedef typename TElemIter::value_type elem_ptr_t;
	Grid& g = obj->grid();
	Mesh::position_accessor_t aaPos = obj->position_accessor();

	vector<elem_ptr_t>	elems(elemsBegin, elemsEnd);

	CustomVertexGroup vrts;
	for(size_t i_elem = 0; i_elem < elems.size(); ++i_elem){
		elem_ptr_t e = elems[i_elem];
		vector3 c = CalculateCenter(e, aaPos);

	//	create new corner-vertices first
		vrts.clear();
		for(size_t i = 0; i < e->num_vertices(); ++i){
			vrts.push_back(*g.create_by_cloning(e->vertex(i), e->vertex(i)));
			VecScaleAdd(aaPos[vrts[i]], scale, aaPos[e->vertex(i)], 1. - scale, c);
		}

	//	now clone the element
		g.create_by_cloning(e, vrts, e);

	//	finally delete the old element
		g.erase(e);
	}
}

/**	For each element in obj this method creates a new element with separate
 * corners. The new element will be scaled by 'scale'.
 * All original elements will be deleted before the method terminates.*/
inline void CreateShrinkGeometry(Mesh* obj, number scale)
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

inline void ReplaceLowValenceVertices(Mesh* obj, number maxSquaredHeightToBaseAreaRatio)
{
	Selector& sel = obj->selector();
	ug::ReplaceLowValenceVertices(obj->grid(), sel.begin<Vertex>(), sel.end<Vertex>(),
						   maxSquaredHeightToBaseAreaRatio, obj->position_accessor());
}

inline void ReplaceValence3Vertices(Mesh* obj, number maxSquaredHeightToBaseAreaRatio)
{
	Selector& sel = obj->selector();
	ug::ReplaceValence3Vertices(obj->grid(), sel.begin<Vertex>(), sel.end<Vertex>(),
						   maxSquaredHeightToBaseAreaRatio, obj->position_accessor());
}

inline void MeshLayers(Mesh* m, const RasterLayers& layers)
{
	MeshLayers(m->grid(), layers, m->position_accessor(), &m->subset_handler());
}

inline void MeshLayerBoundaries(Mesh* m, const RasterLayers& layers)
{
	MeshLayerBoundaries(m->grid(), layers, m->position_accessor(), &m->subset_handler());
}

inline void ExtrudeLayers(Mesh* obj, RasterLayers& layers, bool allowForTetsAndPyras){
	ExtrudeLayers(obj->grid(), layers, obj->position_accessor(),
				  obj->subset_handler(), allowForTetsAndPyras);
}

/// \}

}}// end of namespace

#endif
