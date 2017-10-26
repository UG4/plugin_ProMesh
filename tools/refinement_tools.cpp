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

#include "refinement_tools.h"

using namespace std;

namespace ug{
namespace promesh{

// class TempSelState {
// 	public:
// 		TempSelState (	ISelector& sel,
// 						bool autoselect,
// 						bool selInheritance,
// 						bool strictInheritance) :
// 			m_sel (sel)
// 		{
// 			m_autoselect = m_sel.autoselection_enabled ();
// 			m_selInheritance = m_sel.selection_inheritance_enabled ();
// 			m_strictInheritance = m_sel.strict_inheritance_enabled ();

// 			m_sel.enable_autoselection (autoselect);
// 			m_sel.enable_selection_inheritance (selInheritance);
// 			m_sel.enable_strict_inheritance (strictInheritance);
// 		}

// 		~TempSelState ()
// 		{
// 			m_sel.enable_autoselection (m_autoselect);
// 			m_sel.enable_selection_inheritance (m_selInheritance);
// 			m_sel.enable_strict_inheritance (m_strictInheritance);
// 		}

// 	private:
// 		ISelector&	m_sel;
// 		bool		m_autoselect;
// 		bool		m_selInheritance;
// 		bool		m_strictInheritance;
// };



void Refine(Mesh* obj, bool strictSubsetInheritance, bool useSnapPoints)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	SubsetHandler& sh = obj->subset_handler();
	bool siEnabled = sh.strict_inheritance_enabled();
	sh.enable_strict_inheritance(strictSubsetInheritance);

	ProjectionHandler& projector = obj->projection_handler();

	// TempSelState (sel, false, true, true);

	Refine(grid, sel, &projector, useSnapPoints);

	sh.enable_strict_inheritance(siEnabled);
}

void Refine(Mesh* obj)
{
	Refine(obj, false, false);
}

void RefineWithSnapPoints(Mesh* obj)
{
	Refine(obj, false, true);
}

void RefineWithSnapPointsOrtho(Mesh* obj)
{
	Grid& g = obj->grid();
	Selector& sel = obj->selector();
	Mesh::position_accessor_t aaPos = obj->position_accessor();

	vector<vector3>	projections;
	Grid::face_traits::secure_container	assFaces;

	for(EdgeIterator iedge = sel.begin<Edge>(); iedge != sel.end<Edge>(); ++iedge)
	{
		Edge* e = *iedge;
		const vector3& c0 = aaPos[e->vertex(0)];
		const vector3& c1 = aaPos[e->vertex(1)];
		vector3 dir;
		VecSubtract (dir, c1, c0);

		vector3 projSum (0, 0, 0);
		size_t numProjs = 0;

		g.associated_elements(assFaces, e);
		for(size_t iface = 0; iface < assFaces.size(); ++iface){
			Face* f = assFaces[iface];
			Vertex* selVrt = NULL;
			for(size_t ivrt = 0; ivrt < f->num_vertices(); ++ivrt){
				if(sel.is_selected(f->vertex(ivrt))){
					selVrt = f->vertex(ivrt);
					break;
				}
			}

			if(!selVrt)
				continue;

			vector3 p;
			ProjectPointToRay (p, aaPos[selVrt], c0, dir);
			projSum += p;
			++numProjs;
		}

		if(numProjs == 0){
			projections.push_back (CalculateCenter (e, aaPos));
		}
		else{
			projSum *= (1. / (number)numProjs);
			projections.push_back (projSum);
		}
	}

//	record new vertices in a selector
	VertexSelector selNewVrts(g);
	selNewVrts.enable_autoselection(true);

	Refine(obj, false, true);

//	new edge-vertices have been created in the order of selected edges
	UG_COND_THROW(selNewVrts.num<Vertex>() < projections.size(),
	              "ERROR in RefineWithSnapPointsOrtho: Not enough new vertices created!");

	VertexIterator ivrt = selNewVrts.begin<Vertex>();
	for(size_t i = 0; i < projections.size(); ++i, ++ivrt){
		aaPos[*ivrt] = projections[i];
	}
}


void HangingNodeRefine(Mesh* obj, bool strictSubsetInheritance, bool anisotropic)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	SubsetHandler& sh = obj->subset_handler();
	bool siEnabled = sh.strict_inheritance_enabled();

	sh.enable_strict_inheritance(strictSubsetInheritance);

	HangingNodeRefiner_Grid refiner(grid);

	#ifdef PROMESH_DEBUG_PATH
		string markDbgFile = string(PROMESH_DEBUG_PATH).append("/hnode-marks.ugx");
		UG_LOG("<dbg> SAVING ADJUSTED HNODE MARKS TO " << markDbgFile << endl);
		refiner.set_adjusted_marks_debug_filename(markDbgFile.c_str());
	#endif

	//refiner.enable_automark_objects_of_higher_dim(true);
	//refiner.enable_node_dependency_order_1(false);

	if(anisotropic){
		refiner.mark(sel.edges_begin(), sel.edges_end(), RM_REFINE);
		refiner.mark(sel.faces_begin(), sel.faces_end(), RM_ANISOTROPIC);
		refiner.mark(sel.volumes_begin(), sel.volumes_end(), RM_ANISOTROPIC);
	}
	else{
		refiner.mark(sel.edges_begin(), sel.edges_end(), RM_REFINE);
		refiner.mark(sel.faces_begin(), sel.faces_end(), RM_REFINE);
		refiner.mark(sel.volumes_begin(), sel.volumes_end(), RM_REFINE);
	}


	// TempSelState (sel, false, true, true);
	refiner.refine();

	sh.enable_strict_inheritance(siEnabled);
}

void HangingNodeRefine(Mesh* obj, bool anisotropic)
{
	HangingNodeRefine(obj, false, anisotropic);
}


void RefineSmooth(Mesh* obj, bool strictSubsetInheritance)
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
								 Grid::edge_traits::callback(
								 	IsInSubset(obj->crease_handler(), REM_CREASE)));

	// TempSelState (sel, false, true, true);
	Refine(grid, sel, &refProj);

	sh.enable_strict_inheritance(siEnabled);
}

void RefineSmooth(Mesh* obj)
{
	RefineSmooth(obj, false);
}


void InsertCenter(Mesh* obj, bool strictSubsetInheritance)
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
		RegularVertex* vrt = ug::SplitEdge<RegularVertex>(grid, e);
		aaPos[vrt] = center;
	}

	sh.enable_strict_inheritance(siEnabled);
}

void InsertCenter(Mesh* obj)
{
	InsertCenter(obj, false);
}




class AnisoElemInfo {
public:
	AnisoElemInfo () :
		m_hasShortEdges (false)	{}

	bool is_short_edge (size_t edgeIndex) const	{return m_isShort[edgeIndex];}
	bool has_short_edges () const				{return m_hasShortEdges;}

	template <class TAAPos, class TElem>
	void refresh (TElem* elem, number aspectThreshold, TAAPos aaPos)
	{
		m_tmpLen.clear();
		m_isShort.clear();
		m_hasShortEdges = false;
		number shortest	= numeric_limits<number>::max();
		number longest	= 0;
		EdgeDescriptor ed;
		const int numEdges = (int)elem->num_edges();
		for(int iedge = 0; iedge < numEdges; ++iedge){
			elem->edge_desc(iedge, ed);
			number len = VecDistance(aaPos[ed.vertex(0)], aaPos[ed.vertex(1)]);
			if(len > longest)	longest = len;
			if(len < shortest)	shortest = len;
			m_tmpLen.push_back(len);
		}

		if(longest == 0 || (shortest / longest > aspectThreshold))
			m_isShort.resize(numEdges, false);
		else{
			for(int iedge = 0; iedge < numEdges; ++iedge){
				if(m_tmpLen[iedge] / longest <= aspectThreshold){
					m_isShort.push_back(true);
					m_hasShortEdges = true;
				}
				else
					m_isShort.push_back(false);
			}
		}
	}

private:
	vector<bool>	m_isShort;
	vector<number>	m_tmpLen;
	bool			m_hasShortEdges;
};


template <class TElem>
void RegularizingRefinement_IMPL(Mesh* obj, const number aspectRatio)
{
	typedef TElem	elem_t;
	typedef typename Grid::traits<elem_t>::iterator	iter_t;

	Grid& g = obj->grid();
	Mesh::position_accessor_t aaPos = obj->position_accessor();

	AnisoElemInfo elemInfo;
	HangingNodeRefiner_Grid refiner(g);
	refiner.enable_node_dependency_order_1(false);

	#ifdef PROMESH_DEBUG_PATH
		string markDbgFile = string(PROMESH_DEBUG_PATH).append("/hnode-marks.ugx");
		UG_LOG("<dbg> SAVING ADJUSTED HNODE MARKS TO " << markDbgFile << endl);
		refiner.set_adjusted_marks_debug_filename(markDbgFile.c_str());
	#endif

	for(iter_t ielem = g.begin<elem_t>(); ielem != g.end<elem_t>(); ++ielem)
	{
		elem_t* elem = *ielem;

		elemInfo.refresh(elem, aspectRatio, aaPos);
		int mark = 0;

		if(elemInfo.has_short_edges()){
			const size_t numEdges = elem->num_edges();
			for(size_t iedge = 0; iedge < numEdges; ++iedge){
				if(!elemInfo.is_short_edge(iedge)){
					mark |= 1 << iedge;
				}
			}

			if(elem->is_regular_ref_rule(mark))
				refiner.mark_local(elem, mark);
		}
	}

	refiner.refine();
}


void RegularizingRefinement(Mesh* obj, const number aspectRatio)
{
	Grid& g = obj->grid();
	if(g.num<Volume>() > 0)
		RegularizingRefinement_IMPL<Volume>(obj, aspectRatio);
	else if(g.num<Face>() > 0)
		RegularizingRefinement_IMPL<Face>(obj, aspectRatio);
}

}}//	end of namespace
