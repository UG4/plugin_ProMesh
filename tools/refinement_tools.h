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

#ifndef __H__UG__refinement_tools__
#define __H__UG__refinement_tools__

#include "../mesh.h"
#include "lib_grid/algorithms/refinement/regular_refinement.h"
#include "lib_grid/algorithms/refinement/regular_refinement_new.h"
#include "lib_grid/algorithms/refinement/hanging_node_refiner_grid.h"
#include "lib_grid/algorithms/refinement/projectors/sphere_projector.h"
#include "lib_grid/algorithms/refinement/projectors/subdivision_projector.h"
#include "lib_grid/algorithms/refinement/refinement_projectors_old/loop_subdivision_projectors.h"
#include "lib_grid/algorithms/refinement/refinement_projectors_old/fractal_projector.h"
#include "lib_grid/callbacks/callbacks.h"


namespace ug{
namespace promesh{

/// \addtogroup promesh
/// \{

inline void Refine(Mesh* obj, bool strictSubsetInheritance, bool useSnapPoints)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	SubsetHandler& sh = obj->subset_handler();
	bool siEnabled = sh.strict_inheritance_enabled();
	sh.enable_strict_inheritance(strictSubsetInheritance);

	ProjectionHandler& projector = obj->projection_handler();

//	only for testing...
	// projector.set_projector(1, make_sp(new SphereProjectorNew(obj->geometry(), vector3(0, 0, 0), 1)));
	
	RefineNew(grid, sel, &projector, useSnapPoints);

	sh.enable_strict_inheritance(siEnabled);
}

inline void HangingNodeRefine(Mesh* obj, bool strictSubsetInheritance, bool anisotropic)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	SubsetHandler& sh = obj->subset_handler();
	bool siEnabled = sh.strict_inheritance_enabled();

	sh.enable_strict_inheritance(strictSubsetInheritance);

	HangingNodeRefiner_Grid refiner(grid);
	//refiner.enable_automark_objects_of_higher_dim(true);
	//refiner.enable_node_dependency_order_1(false);

	if(anisotropic){
		refiner.mark(sel.edges_begin(), sel.edges_end(), RM_ANISOTROPIC);
		refiner.mark(sel.faces_begin(), sel.faces_end(), RM_ANISOTROPIC);
		refiner.mark(sel.volumes_begin(), sel.volumes_end(), RM_ANISOTROPIC);
	}
	else{
		refiner.mark(sel.edges_begin(), sel.edges_end(), RM_REFINE);
		refiner.mark(sel.faces_begin(), sel.faces_end(), RM_REFINE);
		refiner.mark(sel.volumes_begin(), sel.volumes_end(), RM_REFINE);
	}

	refiner.refine();

	sh.enable_strict_inheritance(siEnabled);
}

inline void RefineSmooth(Mesh* obj, bool strictSubsetInheritance)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	SubsetHandler& sh = obj->subset_handler();
	bool siEnabled = sh.strict_inheritance_enabled();

	sh.enable_strict_inheritance(strictSubsetInheritance);

//	currently only triangles are supported in smooth refinement.
//	convert all selected quads to triangles first.
	Grid::VertexAttachmentAccessor<APosition> aaPos(grid, aPosition);
	Triangulate(grid, sel.begin<Quadrilateral>(), sel.end<Quadrilateral>(), &aaPos);

	SubdivisionProjector refProj(MakeGeometry3d(grid, aPosition),
								 IsInSubset(obj->crease_handler(), REM_CREASE));
	RefineNew(grid, sel, &refProj);

	sh.enable_strict_inheritance(siEnabled);
}

inline void RefineSmoothBoundary2D(Mesh* obj, bool strictSubsetInheritance)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	SubsetHandler& sh = obj->subset_handler();
	bool siEnabled = sh.strict_inheritance_enabled();

	sh.enable_strict_inheritance(strictSubsetInheritance);

//	since we use a flat hierarchy, a temporary position attachment is required
	APosition aTmpPos;
	grid.attach_to_vertices(aTmpPos);

	SubdivisionLoopBoundaryProjector<APosition>
		refCallbackLoopBnd(grid, aPosition, aTmpPos);

	Refine(grid, sel, &refCallbackLoopBnd);

	sh.enable_strict_inheritance(siEnabled);

//	copy position data of selected vertices
	CopyAttachments(grid, sel.begin<Vertex>(),
						sel.end<Vertex>(),
						aTmpPos, aPosition);

	grid.detach_from_vertices(aTmpPos);
}

inline void CreateFractal(Mesh* obj, size_t numIterations, number scaleFac)
{
	Grid& grid = obj->grid();

//	we'll use a hanging-node refiner
	FractalProjector refCallback(grid, scaleFac);
	HangingNodeRefiner_Grid href(grid);
	href.set_refinement_callback(&refCallback);

//	iterate for the specified number of times
	for(size_t i = 0; i < numIterations; ++i){
		if(grid.num_volumes() > 0){
		//	iterate over all faces and mark them for refinement, if they are boundary faces.
			for(FaceIterator iter = grid.faces_begin();
				iter != grid.faces_end(); ++iter)
			{
				if(IsVolumeBoundaryFace(grid, *iter)){
					href.mark(*iter);
				}

			}
		}
		else if(grid.num_faces() > 0){
		//	markall faces
			href.mark(grid.faces_begin(), grid.faces_end());
		}
		else{
		//	mark all edges
			href.mark(grid.edges_begin(), grid.edges_end());
		}

	//	refine them
		href.refine();

	//	change the scalefac
		refCallback.set_scale_fac(-0.5 * refCallback.get_scale_fac());
		//refCallback.set_scale_fac(refCallback.get_scale_fac() * refCallback.get_scale_fac());
	}
}

inline void InsertCenter(Mesh* obj, bool strictSubsetInheritance)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	SubsetHandler& sh = obj->subset_handler();
	bool siEnabled = sh.strict_inheritance_enabled();

	sh.enable_strict_inheritance(strictSubsetInheritance);

//	access position attachment
	Grid::VertexAttachmentAccessor<APosition> aaPos(grid, aPosition);

	std::vector<Edge*> edges;
	std::vector<Face*> faces;
	std::vector<Volume*> vols;
	vols.assign(sel.begin<Volume>(), sel.end<Volume>());

//todo: support insert center for all selections
	if(grid.num<Volume>()){
		if(sel.num<Face>() > 0){
			UG_LOG("InsertCenter for faces is currently not supported if"
					" volumes are present.\n");
		}
	}
	else
		faces.assign(sel.begin<Face>(), sel.end<Face>());

//todo: support insert center for all selections
	if(grid.num<Face>() > 0){
		if(sel.num<Edge>() > 0){
			UG_LOG("InsertCenter for edges is currently not supported if"
					" faces are present.\n");
		}
	}
	else
		edges.assign(sel.begin<Edge>(), sel.end<Edge>());

//	insert centers
	for(size_t i = 0; i < vols.size(); ++i){
		Volume* vol = vols[i];
		RegularVertex* vrt = *grid.create<RegularVertex>(vol);
		aaPos[vrt] = CalculateCenter(vol, aaPos);
		InsertCenterVertex(grid, vol, vrt, true);
	}

//	insert centers
	for(size_t i = 0; i < faces.size(); ++i){
		Face* f = faces[i];
		RegularVertex* vrt = *grid.create<RegularVertex>(f);
		aaPos[vrt] = CalculateCenter(f, aaPos);
		InsertCenterVertex(grid, f, vrt, true);
	}

//	split edges
	for(size_t i = 0; i < edges.size(); ++i){
		Edge* e = edges[i];
		vector3 center = CalculateCenter(e, aaPos);
		RegularVertex* vrt = SplitEdge<RegularVertex>(grid, e);
		aaPos[vrt] = center;
	}

	sh.enable_strict_inheritance(siEnabled);
}

/// \}

}}// end of namespace

#endif
