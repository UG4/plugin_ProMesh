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

void Refine(Mesh* obj, bool strictSubsetInheritance, bool useSnapPoints)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	SubsetHandler& sh = obj->subset_handler();
	bool siEnabled = sh.strict_inheritance_enabled();
	sh.enable_strict_inheritance(strictSubsetInheritance);

	ProjectionHandler& projector = obj->projection_handler();

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


void HangingNodeRefine(Mesh* obj, bool strictSubsetInheritance, bool anisotropic)
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




class AnisoFaceInfo {
public:
	AnisoFaceInfo () :
		m_hasShortEdges (false)	{}

	bool is_short_edge (size_t edgeIndex) const	{return m_isShort[edgeIndex];}
	bool has_short_edges () const				{return m_hasShortEdges;}

	template <class TAAPos>
	void refresh (Face* f, number aspectThreshold, TAAPos aaPos)
	{
		m_tmpLen.clear();
		m_isShort.clear();
		m_hasShortEdges = false;
		number shortest	= numeric_limits<number>::max();
		number longest	= 0;
		EdgeDescriptor ed;
		const int numEdges = (int)f->num_edges();
		for(int iedge = 0; iedge < numEdges; ++iedge){
			f->edge_desc(iedge, ed);
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



namespace quad_rules {
	bool IsRegularRefMark(int refMark)
	{
		int numMarked = 0;
		bool marks[4];
		for(size_t i = 0; i < 4; ++i){
			int mark = (refMark >> i) & 1;	// 0 or 1
			numMarked += mark;
			marks[i] = mark;
		}

		return		(numMarked == 4)
				||	(numMarked == 2 && (marks[0] == marks[2]));
	}
}

// enum AdvancedHNodeEdgeMarks {
// 	AHN_NONE = 0,
// 	AHN_FULL = 1,
// 	AHN_PARTIAL = 2
// };

// void AdvancedHNodeAdjust(	Grid& g,
// 							vector<Edge*>& refEdges,
// 							vector<Face*>& refFaces,
// 							MultiElementAttachmentAccessor<AInt>& aaMark)
// {

// 	bool adjusting = true;
// 	size_t firstRefFace = 0;

// 	while(adjusting){
// 		adjusting = false;

// 		const size_t numRefFaces = refFaces.size();
// 		for(size_t irefFace = 0; irefFace < numRefFaces; ++irefFace){
// 			Face* refFace = refFaces[irefFace];
// 			const int mark = aaMark[refFace];

// 			}
// 	}
// }


void RegularizingRefinement(Mesh* obj, const number aspectRatio)
{
	Grid& g = obj->grid();
	// Selector& sel = obj->selector();
	Mesh::position_accessor_t aaPos = obj->position_accessor();

	AInt aMark;

	if(g.num<Face>() == 0)
		return;

	// g.attach_to_edges(aMark);
	g.attach_to_faces(aMark);
	g.attach_to_volumes(aMark);

	MultiElementAttachmentAccessor<AInt> aaMark(g, aMark, false, false, true, true);

	AnisoFaceInfo faceInfo;
	HangingNodeRefiner_Grid refiner(g);
	vector<Face*> refFaces;

	for(FaceIterator iface = g.faces_begin(); iface != g.faces_end(); ++iface)
	{
		Face* f = *iface;

		if(f->num_vertices() != 4)
			continue;

		faceInfo.refresh(f, aspectRatio, aaPos);
		int mark = 0;

		if(faceInfo.has_short_edges()){
			const size_t numEdges = f->num_edges();
			for(size_t iedge = 0; iedge < numEdges; ++iedge){
				if(!faceInfo.is_short_edge(iedge)){
					mark |= 1 << iedge;
				}
			}


			if(quad_rules::IsRegularRefMark(mark)){
			//	mark for refinement
			//	don't mark edges, since those are implicitly marked by aaMark
				aaMark[f] = mark;
				refiner.mark(f, RM_ANISOTROPIC);
				for(size_t iedge = 0; iedge < numEdges; ++iedge){
					if(mark & (1 << iedge))
						refiner.mark(g.get_edge(f, iedge), RM_REFINE);
				}
				// refFaces.push_back(f);
			}
		}
	}

	refiner.refine();

	// AdvancedHNodeAdjust(g, refFaces, aMark);
	// AdvancedHNodeRefine(g, refFaces, aMark);


	g.detach_from_faces(aMark);
	g.detach_from_volumes(aMark);	
}

}}//	end of namespace
