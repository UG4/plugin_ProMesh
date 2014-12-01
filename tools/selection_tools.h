// created by Sebastian Reiter
// s.b.reiter@gmail.com
// Jun 18, 2013 (d,m,y)

#ifndef __H__UG__selection_tools__
#define __H__UG__selection_tools__

#include <vector>
#include <stack>
#include "../mesh.h"
#include "common/node_tree/node_tree.h"
#include "lib_grid/algorithms/remeshing/edge_length_adjustment.h"
#include "lib_grid/algorithms/trees/octree.h"
#include "lib_grid/algorithms/mark_util.h"
#include "lib_grid/algorithms/problem_detection_util.h"

namespace ug{
namespace promesh{

inline void ClearSelection(Mesh* obj)
{
	obj->selector().clear();
}

inline void SelectAll(Mesh* obj)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	sel.select(grid.vertices_begin(), grid.vertices_end());
	sel.select(grid.edges_begin(), grid.edges_end());
	sel.select(grid.faces_begin(), grid.faces_end());
	sel.select(grid.volumes_begin(), grid.volumes_end());
}

inline void ExtendSelection(Mesh* obj, int neighborhoodSize)
{
	Selector& sel = obj->selector();
	ExtendSelection(sel, (size_t)neighborhoodSize);
}

template <class TElem>
TElem* SelectElemByCoordinate(Mesh* obj, const vector3& coord)
{
	Grid& grid = obj->grid();
	Mesh::position_accessor_t& aaPos = obj->position_accessor();

	TElem* e = FindClosestByCoordinate<TElem>(coord,
											 grid.begin<TElem>(),
											 grid.end<TElem>(),
											 aaPos);

	if(e)
		obj->selector().select(e);

	return e;
}

template <class TElem>
TElem* SelectElemByCylindricalCoordinate(Mesh* obj, number rho, number phi, number z)
{
	number x = rho * cos(phi);
	number y = rho * sin(phi);

	vector3 coord = vector3(x,y,z);
	return SelectElemByCoordinate<TElem>(obj, coord);
}

inline void SelectSubset(Mesh* obj, int si, bool selVrts, bool selEdges,
				  		 bool selFaces, bool selVols)
{
	SubsetHandler& sh = obj->subset_handler();
	Selector& sel = obj->selector();

	if(si >= 0){
		if(selVrts)
			sel.select(sh.begin<Vertex>(si), sh.end<Vertex>(si));
		if(selEdges)
			sel.select(sh.begin<Edge>(si), sh.end<Edge>(si));
		if(selFaces)
			sel.select(sh.begin<Face>(si), sh.end<Face>(si));
		if(selVols)
			sel.select(sh.begin<Volume>(si), sh.end<Volume>(si));
	}
	else{
		Grid& grid = obj->grid();
	//	subset -1 has to be selected. Those are not directly accessible.
		if(selVrts){
			for(VertexIterator iter = grid.vertices_begin();
				iter != grid.vertices_end(); ++iter)
			{
				if(sh.get_subset_index(*iter) == -1)
					sel.select(*iter);
			}
		}
		if(selEdges){
			for(EdgeIterator iter = grid.edges_begin();
				iter != grid.edges_end(); ++iter)
			{
				if(sh.get_subset_index(*iter) == -1)
					sel.select(*iter);
			}
		}
		if(selFaces){
			for(FaceIterator iter = grid.faces_begin();
				iter != grid.faces_end(); ++iter)
			{
				if(sh.get_subset_index(*iter) == -1)
					sel.select(*iter);
			}
		}
		if(selVols){
			for(VolumeIterator iter = grid.volumes_begin();
				iter != grid.volumes_end(); ++iter)
			{
				if(sh.get_subset_index(*iter) == -1)
					sel.select(*iter);
			}
		}
	}
}

inline void SelectSubsetBoundary(Mesh* obj, int si, bool edgeBnds,
						  bool faceBnds, bool volBnds)
{
	SubsetHandler& sh = obj->subset_handler();
	Selector& sel = obj->selector();

	if(edgeBnds)
		SelectAreaBoundary(sel, sh.begin<Edge>(si), sh.end<Edge>(si));
	if(faceBnds)
		SelectAreaBoundary(sel, sh.begin<Face>(si), sh.end<Face>(si));
	if(volBnds)
		SelectAreaBoundary(sel, sh.begin<Volume>(si), sh.end<Volume>(si));
}

template <class TGeomObj>
static inline void SelectUnassignedElementsHelper(Grid& grid, SubsetHandler& sh, Selector& sel)
{
	typedef typename geometry_traits<TGeomObj>::iterator	iterator;
	for(iterator iter = grid.begin<TGeomObj>(); iter != grid.end<TGeomObj>(); ++iter)
	{
		if(sh.get_subset_index(*iter) == -1){
			sel.select(*iter);
		}
	}
}

inline void SelectUnassignedElements(Mesh* obj, bool selVrts, bool selEdges,
							  bool selFaces, bool selVols)
{
	Grid& grid = obj->grid();
	SubsetHandler& sh = obj->subset_handler();
	Selector& sel = obj->selector();

	if(selVrts)
		SelectUnassignedElementsHelper<Vertex>(grid, sh, sel);
	if(selEdges)
		SelectUnassignedElementsHelper<Edge>(grid, sh, sel);
	if(selFaces)
		SelectUnassignedElementsHelper<Face>(grid, sh, sel);
	if(selVols)
		SelectUnassignedElementsHelper<Volume>(grid, sh, sel);
}

inline void InvertSelection(Mesh* obj, bool invVrts, bool invEdges,
					 bool invFaces, bool invVols)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();

	if(invVrts)
		InvertSelection(sel, grid.begin<Vertex>(),
							grid.end<Vertex>());

	if(invEdges)
		InvertSelection(sel, grid.begin<Edge>(),
							grid.end<Edge>());

	if(invFaces)
		InvertSelection(sel, grid.begin<Face>(),
							grid.end<Face>());

	if(invVols)
		InvertSelection(sel, grid.begin<Volume>(),
							grid.end<Volume>());
}

inline void SelectSelectionBoundary(Mesh* obj)
{
	Selector& sel = obj->selector();

	if(sel.num<Volume>() > 0)
		SelectAreaBoundary(sel, sel.begin<Volume>(), sel.end<Volume>());
	if(sel.num<Face>() > 0)
		SelectAreaBoundary(sel, sel.begin<Face>(), sel.end<Face>());
}

inline void CloseSelection(Mesh* obj)
{
	Selector& sel = obj->selector();
	SelectAssociatedFaces(sel, sel.begin<Volume>(), sel.end<Volume>());
	SelectAssociatedEdges(sel, sel.begin<Face>(), sel.end<Face>());
	SelectAssociatedVertices(sel, sel.begin<Edge>(), sel.end<Edge>());
}


////////////////////////////////////////////////////////////////////////////////
//	VERTICES
inline void SelectBoundaryVertices(Mesh* obj)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	ug::SelectBoundaryElements(sel, grid.begin<Vertex>(), grid.end<Vertex>());
}

inline void SelectInnerVertices(Mesh* obj)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	SelectInnerElements(sel, grid.begin<Vertex>(), grid.end<Vertex>());
}

inline void SelectAssociatedVertices(Mesh* obj)
{
	Selector& sel = obj->selector();
	SelectAssociatedVertices(sel, sel.begin<Edge>(), sel.end<Edge>());
	SelectAssociatedVertices(sel, sel.begin<Face>(), sel.end<Face>());
	SelectAssociatedVertices(sel, sel.begin<Volume>(), sel.end<Volume>());
}

inline void SelectAllVertices(Mesh* obj)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	sel.select(grid.vertices_begin(), grid.vertices_end());
}

inline void DeselectAllVertices(Mesh* obj)
{
	Selector& sel = obj->selector();
	sel.deselect(sel.vertices_begin(), sel.vertices_end());
}

inline void SelectMarkedVertices(Mesh* obj)
{
	obj->selector().select(
		obj->crease_handler().begin<Vertex>(REM_FIXED),
		obj->crease_handler().end<Vertex>(REM_FIXED));
}

inline bool SelectVertexByIndex(Mesh* obj, int index)
{
	Grid& grid = obj->grid();
	int counter = 0;

	for(VertexIterator iter = grid.begin<Vertex>();
		iter != grid.end<Vertex>(); ++iter, ++counter)
	{
		if(counter == index){
			obj->selector().select(*iter);
			return true;
		}
	}
	return false;
}

template <class TElem>
inline size_t SelectUnconnectedVerticesHelper(Grid& grid, Selector& sel)
{
	using namespace ug;
	typename Grid::traits<TElem>::secure_container	elems;

	size_t numUnconnected = 0;
	for(VertexIterator iter = grid.vertices_begin();
		iter != grid.vertices_end(); ++iter)
	{
		if(!sel.is_selected(*iter)){
			grid.associated_elements(elems, *iter);
			if(elems.size() == 0){
				sel.select(*iter);
				numUnconnected++;
			}
		}
	}
	return numUnconnected;
}

inline size_t SelectUnconnectedVertices(Mesh* obj, bool edgeCons, bool faceCons, bool volCons)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();

	size_t numUnconnected = 0;
	if(edgeCons)
		numUnconnected += SelectUnconnectedVerticesHelper<Edge>(grid, sel);
	if(faceCons)
		numUnconnected += SelectUnconnectedVerticesHelper<Face>(grid, sel);
	if(volCons)
		numUnconnected += SelectUnconnectedVerticesHelper<Volume>(grid, sel);
	return numUnconnected;
}

////////////////////////////////////////////////////////////////////////////////
//	EDGES
inline void SelectBoundaryEdges(Mesh* obj)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	SelectBoundaryElements(sel, grid.begin<Edge>(), grid.end<Edge>());
}

inline void SelectInnerEdges(Mesh* obj)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	SelectInnerElements(sel, grid.begin<Edge>(), grid.end<Edge>());
}

inline void SelectNonManifoldEdges(Mesh* obj)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();

//	iterate over all edges and check how many adjacent faces each has.
//	if there are more than 2, the edge will be selected
	for(EdgeIterator iter = grid.edges_begin();
		iter != grid.edges_end(); ++iter)
	{
		if(NumAssociatedFaces(grid, *iter) != 2)
			sel.select(*iter);
	}
}

inline void SelectSmoothEdgePath(Mesh* obj, number maxDeviation, bool stopAtSelectedVrts)
{
	SelectSmoothEdgePath(obj->selector(), maxDeviation, stopAtSelectedVrts);
}

inline void SelectShortEdges(Mesh* obj, number maxLength)
{
	number maxLengthSq = maxLength * maxLength;
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	Mesh::position_accessor_t& aaPos = obj->position_accessor();

	for(EdgeIterator iter = grid.begin<Edge>();
		iter != grid.end<Edge>(); ++iter)
	{
		Edge* e = *iter;
		if(VecDistanceSq(aaPos[e->vertex(0)], aaPos[e->vertex(1)]) < maxLengthSq)
			sel.select(e);
	}
}

inline void SelectLongEdges(Mesh* obj, number minLength)
{
	number minLengthSq = minLength * minLength;
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	Mesh::position_accessor_t& aaPos = obj->position_accessor();

	for(EdgeIterator iter = grid.begin<Edge>();
		iter != grid.end<Edge>(); ++iter)
	{
		Edge* e = *iter;
		if(VecDistanceSq(aaPos[e->vertex(0)], aaPos[e->vertex(1)]) > minLengthSq)
			sel.select(e);
	}
}

inline void SelectCreaseEdges(Mesh* obj, number minAngle)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	SelectCreaseEdges(sel, grid.begin<Edge>(), grid.end<Edge>(),
					  minAngle, obj->position_attachment());
}

inline void SelectLinkedBoundaryEdges(Mesh* obj, bool stopAtSelectedVrts)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();

	if(stopAtSelectedVrts)
		SelectLinkedElements<Edge>(sel, IsOnBoundary(grid), IsNotSelected(sel));
	else
		SelectLinkedElements<Edge>(sel, IsOnBoundary(grid));
}

inline void SelectAssociatedEdges(Mesh* obj)
{
	Selector& sel = obj->selector();
	SelectAssociatedEdges(sel,sel.begin<Face>(), sel.end<Face>());
	SelectAssociatedEdges(sel, sel.begin<Volume>(), sel.end<Volume>());
}

inline void SelectAllEdges(Mesh* obj)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	sel.select(grid.edges_begin(), grid.edges_end());
}

inline void DeselectAllEdges(Mesh* obj)
{
	Selector& sel = obj->selector();
	sel.deselect(sel.edges_begin(), sel.edges_end());
}

inline void SelectMarkedEdges(Mesh* obj)
{
	obj->selector().select(
		obj->crease_handler().begin<Edge>(REM_CREASE),
		obj->crease_handler().end<Edge>(REM_CREASE));
}

inline bool SelectEdgeByIndex(Mesh* obj, int index)
{
	Grid& grid = obj->grid();
	int counter = 0;

	for(EdgeIterator iter = grid.begin<Edge>();
		iter != grid.end<Edge>(); ++iter, ++counter)
	{
		if(counter == index){
			obj->selector().select(*iter);
			return true;
		}
	}
	return false;
}

inline void EdgeSelectionFill(Mesh* obj)
{
	Selector& sel = obj->selector();
	SelectionFill<Edge>(sel, IsSelected(sel));
}

////////////////////////////////////////////////////////////////////////////////
//	FACES
inline void SelectBoundaryFaces(Mesh* obj)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	SelectBoundaryElements(sel, grid.begin<Face>(), grid.end<Face>());
}

inline void SelectInnerFaces(Mesh* obj)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	SelectInnerElements(sel, grid.begin<Face>(), grid.end<Face>());
}

inline void SelectLinkedManifoldFaces(Mesh* obj)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
//	select associated faces of selected elements
	SelectAssociatedFaces(sel, sel.vertices_begin(), sel.vertices_end());
	SelectAssociatedFaces(sel, sel.edges_begin(), sel.edges_end());
	SelectAssociatedFaces(sel, sel.volumes_begin(), sel.volumes_end());

//	push all selected faces to a stack
	std::stack<Face*> stk;
	for(FaceIterator iter = sel.faces_begin(); iter != sel.faces_end(); ++iter){
		stk.push(*iter);
	}

//	while there are faces in the stack, get their associated edges.
//	if those edges are adjacent to exactly two faces, then select the
//	neighboured face and push it to the stack (if it was unselected previously)
	std::vector<Edge*> edges;
	while(!stk.empty()){
		Face* f = stk.top();
		stk.pop();
		CollectEdges(edges, grid, f);

		for(size_t i = 0; i < edges.size(); ++i){
			Face* faces[2];
			if(GetAssociatedFaces(faces, grid, edges[i], 2) == 2){
				for(size_t j = 0; j < 2; ++j){
					if(!sel.is_selected(faces[j])){
						sel.select(faces[j]);
						stk.push(faces[j]);
					}
				}
			}
		}
	}
}

inline void SelectLinkedBoundaryFaces(Mesh* obj, bool stopAtSelectedEdges)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();

	if(stopAtSelectedEdges)
		SelectLinkedElements<Face>(sel, IsOnBoundary(grid), IsNotSelected(sel));
	else
		SelectLinkedElements<Face>(sel, IsOnBoundary(grid));
}

inline void SelectDegenerateFaces(Mesh* obj, number maxHeight)
{
	number maxHeightSq = maxHeight * maxHeight;
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	Mesh::position_accessor_t& aaPos = obj->position_accessor();

	for(FaceIterator iter = grid.begin<Face>();
		iter != grid.end<Face>(); ++iter)
	{
		Face* f = *iter;
	//	iterate over the edges and check which is the longest.
		number maxLenSq = -1;
		size_t bestInd = -1;
		EdgeDescriptor ed;
		for(size_t i = 0; i < f->num_edges(); ++i){
			f->edge_desc(i, ed);
			number lenSq = VecDistanceSq(aaPos[ed.vertex(0)], aaPos[ed.vertex(1)]);
			if(lenSq > maxLenSq){
				maxLenSq = lenSq;
				bestInd = i;
			}
		}

		if(maxLenSq < maxHeightSq)
			sel.select(f);
		else{
		//	project the other vertices to the line and check the height
		//todo: this is not enough for quadrilaterals
			vector3 p;
			vector3& v = aaPos[f->vertex((bestInd + 2) % f->num_vertices())];
			f->edge_desc(bestInd, ed);
			number t;
			if(DistancePointToLine(t, v, aaPos[ed.vertex(0)], aaPos[ed.vertex(1)]) < maxHeight)
				sel.select(f);
		}
	}
}

inline void SelectLinkedFlatFaces(Mesh* obj, number maxDeviationAngle, bool traverseFlipped,
						   bool traverseDegeneratedFaces, bool stopAtSelectedEdges)
{
	Selector& sel = obj->selector();
	if(traverseDegeneratedFaces)
		ug::SelectLinkedFlatAndDegeneratedFaces(sel, maxDeviationAngle,
												traverseFlipped,
												stopAtSelectedEdges);
	else
		ug::SelectLinkedFlatFaces(sel, maxDeviationAngle, traverseFlipped,
								  stopAtSelectedEdges);
}

inline void SelectIntersectingTriangles(Mesh* obj)
{
	using namespace ug::node_tree;

	Grid& grid = obj->grid();
	Selector& sel = obj->selector();

//	make sure that triangles are present.
	if(grid.num<Triangle>() == 0){
		UG_LOG("Given grid does not contain any triangles. Aborting 'Select Intersecting Triangles'.\n");
		return;
	}

	if(grid.num<Quadrilateral>() > 0){
		UG_LOG("WARNING: Quadrilateral intersections are ignored during 'Select Intersecting Triangles'.\n");
	}

//	create an octree
//	sort the triangles of grid into an octree to speed-up projection performance
	SPOctree octree;
	octree = CreateOctree(grid, grid.begin<Triangle>(),
								grid.end<Triangle>(),
								10, 30, false, obj->position_attachment());

	if(!octree.valid()){
		UG_LOG("  Octree creation failed in 'Select Intersecting Triangles'. Aborting.\n");
		return;
	}

//	access the position attachment of the grid
	Mesh::position_accessor_t& aaPos = obj->position_accessor();

//	we'll use a traverser to find the intersections.
	Traverser_IntersectFaces intersectionTraverser;

//	for each face the face itself and direct neighbors shall be ignored
	std::vector<Face*> ignoreList;

	size_t triCount = 0;
//	now iterate over all triangles of the grid and find intersections
	for(TriangleIterator iter = grid.begin<Triangle>();
		iter != grid.end<Triangle>(); ++iter, ++triCount)
	{
		Triangle* t = *iter;

	//	add neighbors and the face itself to the ignore list
		CollectNeighbors(ignoreList, t, grid, NHT_VERTEX_NEIGHBORS);
		intersectionTraverser.clear_ignore_list();
		for(size_t i = 0; i < ignoreList.size(); ++i){
			intersectionTraverser.ignore_element(ignoreList[i]);
		}
		intersectionTraverser.ignore_element(t);

		if(intersectionTraverser.intersect_tri(aaPos[t->vertex(0)], aaPos[t->vertex(1)],
												aaPos[t->vertex(2)], octree))
		{
		//	check the intersecting face:
			/*
			const std::vector<CollisionElementID>& faces =
				intersectionTraverser.get_intersected_element_ids();
			*/
		//	an intersection occurred. Log the index and select the triangle.
			sel.select(t);
		}
	}

}

inline void SelectAssociatedFaces(Mesh* obj)
{
	Selector& sel = obj->selector();
	SelectAssociatedFaces(sel, sel.begin<Volume>(), sel.end<Volume>());
}

inline void SelectAllFaces(Mesh* obj)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	sel.select(grid.faces_begin(), grid.faces_end());
}

inline void DeselectAllFaces(Mesh* obj)
{
	Selector& sel = obj->selector();
	sel.deselect(sel.faces_begin(), sel.faces_end());
}

inline bool SelectFaceByIndex(Mesh* obj, int index)
{
	Grid& grid = obj->grid();
	int counter = 0;

	for(FaceIterator iter = grid.begin<Face>();
		iter != grid.end<Face>(); ++iter, ++counter)
	{
		if(counter == index){
			obj->selector().select(*iter);
			return true;
		}
	}
	return false;
}

inline void SelectFacesByNormal(Mesh* obj, const vector3& refNormal, number maxDeviationAngle)
{
	number dotThreshold = cos(deg_to_rad(maxDeviationAngle));
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	Grid::AttachmentAccessor<Vertex, APosition> aaPos(grid, aPosition);

	vector3 n;
	VecNormalize(n, refNormal);

	for(FaceIterator iter = grid.faces_begin(); iter != grid.faces_end(); ++iter){
		vector3 fn;
		CalculateNormal(fn, *iter, aaPos);
		if(VecDot(fn, n) >= dotThreshold)
			sel.select(*iter);
	}
}

inline void FaceSelectionFill(Mesh* obj)
{
	Selector& sel = obj->selector();
	SelectionFill<Face>(sel, IsSelected(sel));
}

inline size_t SelectBentQuadrilaterals(Mesh* obj, number dotThreshold)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	Grid::AttachmentAccessor<Vertex, APosition> aaPos(grid, aPosition);

//	iterate over all quadrilaterals and search for bent ones
	size_t selCount = 0;
	for(QuadrilateralIterator iter = grid.begin<Quadrilateral>();
		iter != grid.end<Quadrilateral>(); ++iter)
	{
		Quadrilateral* q = *iter;

	//	we'll compare the dot-product of the normals
		vector3 n1, n2;
		CalculateTriangleNormal(n1, aaPos[q->vertex(0)],
								aaPos[q->vertex(1)], aaPos[q->vertex(2)]);
		CalculateTriangleNormal(n2, aaPos[q->vertex(2)],
								aaPos[q->vertex(3)], aaPos[q->vertex(0)]);

		number d1 = VecDot(n1, n2);

		if(d1 < dotThreshold){
			++selCount;
			sel.select(q);
		}
	}

	return selCount;
}


////////////////////////////////////////////////////////////////////////////////
//	VOLUMES
inline void SelectAllVolumes(Mesh* obj)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	sel.select(grid.volumes_begin(), grid.volumes_end());
}

inline void DeselectAllVolumes(Mesh* obj)
{
	Selector& sel = obj->selector();
	sel.deselect(sel.volumes_begin(), sel.volumes_end());
}

inline int SelectUnorientableVolumes(Mesh* obj)
{
	Grid& grid = obj->grid();
	Selector& sel = obj->selector();

	Grid::AttachmentAccessor<Vertex, APosition> aaPos(grid, aPosition);

	int numUnorientable = 0;

//	iterate over all volumes. Check whether the orientation can be determined.
	for(VolumeIterator iter = grid.volumes_begin();
		iter != grid.volumes_end(); ++iter)
	{
		Volume* v = *iter;
	//	get orientation of original volume
		bool bOriented = CheckOrientation(v, aaPos);
	//	flip the volume
		grid.flip_orientation(v);

	//	if orientations match, then the volume can not be oriented.
		if(bOriented == CheckOrientation(v, aaPos)){
			sel.select(v);
			++numUnorientable;
		}

	//	reflip orientation
		grid.flip_orientation(v);
	}

	return numUnorientable;
}

inline int SelectSlivers(Mesh* obj, number thresholdRatio)
{
	std::vector<Tetrahedron*> slivers;
	Grid& g = obj->grid();
	FindSlivers(slivers, g.begin<Tetrahedron>(), g.end<Tetrahedron>(),
			    thresholdRatio, obj->position_accessor());
	obj->selector().select(slivers.begin(), slivers.end());
	return static_cast<int>(slivers.size());
}

inline bool SelectVolumeByIndex(Mesh* obj, int index)
{
	Grid& grid = obj->grid();
	int counter = 0;

	for(VolumeIterator iter = grid.begin<Volume>();
		iter != grid.end<Volume>(); ++iter, ++counter)
	{
		if(counter == index){
			obj->selector().select(*iter);
			return true;
		}
	}
	return false;
}

inline void VolumeSelectionFill(Mesh* obj)
{
	Selector& sel = obj->selector();
	SelectionFill<Volume>(sel, IsSelected(sel));
}

// template <class TVertexIterator, class TAPosition>
// UG_API 
// void MarkCorners(Grid& grid, ISubsetHandler& sh,
// 					TVertexIterator vrtsBegin, TVertexIterator vrtsEnd,
// 					Grid::edge_traits::callback cbPathEdge,
// 					int subsetIndex, number angle,
// 					TAPosition& aPos);

inline void MarkCornersOfMarkedEdges(Mesh* obj, number angle)
{
	Selector& sel = obj->selector();
	MarkCorners(obj->grid(), obj->crease_handler(),
				sel.begin<Vertex>(), sel.end<Vertex>(),
				IsInSubset(obj->crease_handler(), ug::REM_CREASE),
				ug::REM_FIXED, angle, obj->position_attachment());
}


inline void MarkSelection(Mesh* obj)
{
	obj->crease_handler().assign_subset(
		obj->selector().begin<Vertex>(),
		obj->selector().end<Vertex>(),
		ug::REM_FIXED);
	obj->crease_handler().assign_subset(
		obj->selector().begin<Edge>(),
		obj->selector().end<Edge>(),
		ug::REM_CREASE);
}

inline void UnmarkSelection(Mesh* obj)
{
	obj->crease_handler().assign_subset(
		obj->selector().begin<Vertex>(),
		obj->selector().end<Vertex>(),
		ug::REM_NONE);

	obj->crease_handler().assign_subset(
		obj->selector().begin<Edge>(),
		obj->selector().end<Edge>(),
		ug::REM_NONE);
}

}}// end of namespace

#endif
