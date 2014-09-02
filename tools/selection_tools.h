// created by Sebastian Reiter
// s.b.reiter@gmail.com
// Jun 18, 2013 (d,m,y)

#ifndef __H__UG__selection_tools__
#define __H__UG__selection_tools__

#include <vector>
#include <stack>
#include "../mesh_object.h"
#include "common/node_tree/node_tree.h"
#include "lib_grid/algorithms/remeshing/edge_length_adjustment.h"
#include "lib_grid/algorithms/trees/octree.h"

namespace ug{
namespace promesh{

inline void ClearSelection(MeshObject* obj)
{
	obj->get_selector().clear();
}

inline void SelectAll(MeshObject* obj)
{
	Grid& grid = obj->get_grid();
	Selector& sel = obj->get_selector();
	sel.select(grid.vertices_begin(), grid.vertices_end());
	sel.select(grid.edges_begin(), grid.edges_end());
	sel.select(grid.faces_begin(), grid.faces_end());
	sel.select(grid.volumes_begin(), grid.volumes_end());
}

inline void ExtendSelection(MeshObject* obj, int neighborhoodSize)
{
	Selector& sel = obj->get_selector();
	ExtendSelection(sel, (size_t)neighborhoodSize);
}

template <class TElem>
TElem* SelectElemByCoordinate(MeshObject* obj, const vector3& coord)
{
	Grid& grid = obj->get_grid();
	MeshObject::position_accessor_t& aaPos = obj->position_accessor();

	TElem* e = FindClosestByCoordinate<TElem>(coord,
											 grid.begin<TElem>(),
											 grid.end<TElem>(),
											 aaPos);

	if(e)
		obj->get_selector().select(e);

	return e;
}

template <class TElem>
TElem* SelectElemByCylindricalCoordinate(MeshObject* obj, number rho, number phi, number z)
{
	Grid& grid = obj->get_grid();

	number x = rho * cos(phi);
	number y = rho * sin(phi);

	vector3 coord = vector3(x,y,z);
	return SelectElemByCoordinate<TElem>(obj, coord);
}

inline void SelectSubset(MeshObject* obj, int si, bool selVrts, bool selEdges,
				  		 bool selFaces, bool selVols)
{
	SubsetHandler& sh = obj->get_subset_handler();
	Selector& sel = obj->get_selector();

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
		Grid& grid = obj->get_grid();
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

inline void SelectSubsetBoundary(MeshObject* obj, int si, bool edgeBnds,
						  bool faceBnds, bool volBnds)
{
	SubsetHandler& sh = obj->get_subset_handler();
	Selector& sel = obj->get_selector();

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

inline void SelectUnassignedElements(MeshObject* obj, bool selVrts, bool selEdges,
							  bool selFaces, bool selVols)
{
	Grid& grid = obj->get_grid();
	SubsetHandler& sh = obj->get_subset_handler();
	Selector& sel = obj->get_selector();

	if(selVrts)
		SelectUnassignedElementsHelper<Vertex>(grid, sh, sel);
	if(selEdges)
		SelectUnassignedElementsHelper<Edge>(grid, sh, sel);
	if(selFaces)
		SelectUnassignedElementsHelper<Face>(grid, sh, sel);
	if(selVols)
		SelectUnassignedElementsHelper<Volume>(grid, sh, sel);
}

inline void InvertSelection(MeshObject* obj, bool invVrts, bool invEdges,
					 bool invFaces, bool invVols)
{
	Grid& grid = obj->get_grid();
	Selector& sel = obj->get_selector();

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

inline void SelectSelectionBoundary(MeshObject* obj)
{
	Selector& sel = obj->get_selector();

	if(sel.num<Volume>() > 0)
		SelectAreaBoundary(sel, sel.begin<Volume>(), sel.end<Volume>());
	if(sel.num<Face>() > 0)
		SelectAreaBoundary(sel, sel.begin<Face>(), sel.end<Face>());
}

inline void CloseSelection(MeshObject* obj)
{
	Selector& sel = obj->get_selector();
	SelectAssociatedFaces(sel, sel.begin<Volume>(), sel.end<Volume>());
	SelectAssociatedEdges(sel, sel.begin<Face>(), sel.end<Face>());
	SelectAssociatedVertices(sel, sel.begin<Edge>(), sel.end<Edge>());
}


////////////////////////////////////////////////////////////////////////////////
//	VERTICES
inline void SelectBoundaryVertices(MeshObject* obj)
{
	Grid& grid = obj->get_grid();
	Selector& sel = obj->get_selector();
	ug::SelectBoundaryElements(sel, grid.begin<Vertex>(), grid.end<Vertex>());
}

inline void SelectInnerVertices(MeshObject* obj)
{
	Grid& grid = obj->get_grid();
	Selector& sel = obj->get_selector();
	SelectInnerElements(sel, grid.begin<Vertex>(), grid.end<Vertex>());
}

inline void SelectAssociatedVertices(MeshObject* obj)
{
	Selector& sel = obj->get_selector();
	SelectAssociatedVertices(sel, sel.begin<Edge>(), sel.end<Edge>());
	SelectAssociatedVertices(sel, sel.begin<Face>(), sel.end<Face>());
	SelectAssociatedVertices(sel, sel.begin<Volume>(), sel.end<Volume>());
}

inline void SelectAllVertices(MeshObject* obj)
{
	Grid& grid = obj->get_grid();
	Selector& sel = obj->get_selector();
	sel.select(grid.vertices_begin(), grid.vertices_end());
}

inline void DeselectAllVertices(MeshObject* obj)
{
	Selector& sel = obj->get_selector();
	sel.deselect(sel.vertices_begin(), sel.vertices_end());
}

inline void SelectMarkedVertices(MeshObject* obj)
{
	obj->get_selector().select(
		obj->get_crease_handler().begin<Vertex>(REM_FIXED),
		obj->get_crease_handler().end<Vertex>(REM_FIXED));
}

inline bool SelectVertexByIndex(MeshObject* obj, int index)
{
	Grid& grid = obj->get_grid();
	int counter = 0;

	for(VertexIterator iter = grid.begin<Vertex>();
		iter != grid.end<Vertex>(); ++iter, ++counter)
	{
		if(counter == index){
			obj->get_selector().select(*iter);
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

inline size_t SelectUnconnectedVertices(MeshObject* obj, bool edgeCons, bool faceCons, bool volCons)
{
	Grid& grid = obj->get_grid();
	Selector& sel = obj->get_selector();

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
inline void SelectBoundaryEdges(MeshObject* obj)
{
	Grid& grid = obj->get_grid();
	Selector& sel = obj->get_selector();
	SelectBoundaryElements(sel, grid.begin<Edge>(), grid.end<Edge>());
}

inline void SelectInnerEdges(MeshObject* obj)
{
	Grid& grid = obj->get_grid();
	Selector& sel = obj->get_selector();
	SelectInnerElements(sel, grid.begin<Edge>(), grid.end<Edge>());
}

inline void SelectNonManifoldEdges(MeshObject* obj)
{
	Grid& grid = obj->get_grid();
	Selector& sel = obj->get_selector();

//	iterate over all edges and check how many adjacent faces each has.
//	if there are more than 2, the edge will be selected
	for(EdgeIterator iter = grid.edges_begin();
		iter != grid.edges_end(); ++iter)
	{
		if(NumAssociatedFaces(grid, *iter) != 2)
			sel.select(*iter);
	}
}

inline void SelectSmoothEdgePath(MeshObject* obj, number maxDeviation, bool stopAtSelectedVrts)
{
	SelectSmoothEdgePath(obj->get_selector(), maxDeviation, stopAtSelectedVrts);
}

inline void SelectShortEdges(MeshObject* obj, number maxLength)
{
	number maxLengthSq = maxLength * maxLength;
	Grid& grid = obj->get_grid();
	Selector& sel = obj->get_selector();
	MeshObject::position_accessor_t& aaPos = obj->position_accessor();

	for(EdgeIterator iter = grid.begin<Edge>();
		iter != grid.end<Edge>(); ++iter)
	{
		Edge* e = *iter;
		if(VecDistanceSq(aaPos[e->vertex(0)], aaPos[e->vertex(1)]) < maxLengthSq)
			sel.select(e);
	}
}

inline void SelectLongEdges(MeshObject* obj, number minLength)
{
	number minLengthSq = minLength * minLength;
	Grid& grid = obj->get_grid();
	Selector& sel = obj->get_selector();
	MeshObject::position_accessor_t& aaPos = obj->position_accessor();

	for(EdgeIterator iter = grid.begin<Edge>();
		iter != grid.end<Edge>(); ++iter)
	{
		Edge* e = *iter;
		if(VecDistanceSq(aaPos[e->vertex(0)], aaPos[e->vertex(1)]) > minLengthSq)
			sel.select(e);
	}
}

inline void SelectCreaseEdges(MeshObject* obj, number minAngle)
{
	Grid& grid = obj->get_grid();
	Selector& sel = obj->get_selector();
	SelectCreaseEdges(sel, grid.begin<Edge>(), grid.end<Edge>(),
					  minAngle, obj->position_attachment());
}

inline void SelectLinkedBoundaryEdges(MeshObject* obj, bool stopAtSelectedVrts)
{
	Grid& grid = obj->get_grid();
	Selector& sel = obj->get_selector();

	if(stopAtSelectedVrts)
		SelectLinkedElements<Edge>(sel, IsOnBoundary(grid), IsNotSelected(sel));
	else
		SelectLinkedElements<Edge>(sel, IsOnBoundary(grid));
}

inline void SelectAssociatedEdges(MeshObject* obj)
{
	Selector& sel = obj->get_selector();
	SelectAssociatedEdges(sel,sel.begin<Face>(), sel.end<Face>());
	SelectAssociatedEdges(sel, sel.begin<Volume>(), sel.end<Volume>());
}

inline void SelectAllEdges(MeshObject* obj)
{
	Grid& grid = obj->get_grid();
	Selector& sel = obj->get_selector();
	sel.select(grid.edges_begin(), grid.edges_end());
}

inline void DeselectAllEdges(MeshObject* obj)
{
	Selector& sel = obj->get_selector();
	sel.deselect(sel.edges_begin(), sel.edges_end());
}

inline void SelectMarkedEdges(MeshObject* obj)
{
	obj->get_selector().select(
		obj->get_crease_handler().begin<Edge>(REM_CREASE),
		obj->get_crease_handler().end<Edge>(REM_CREASE));
}

inline bool SelectEdgeByIndex(MeshObject* obj, int index)
{
	Grid& grid = obj->get_grid();
	int counter = 0;

	for(EdgeIterator iter = grid.begin<Edge>();
		iter != grid.end<Edge>(); ++iter, ++counter)
	{
		if(counter == index){
			obj->get_selector().select(*iter);
			return true;
		}
	}
	return false;
}

inline void EdgeSelectionFill(MeshObject* obj)
{
	Selector& sel = obj->get_selector();
	SelectionFill<Edge>(sel, IsSelected(sel));
}

////////////////////////////////////////////////////////////////////////////////
//	FACES
inline void SelectBoundaryFaces(MeshObject* obj)
{
	Grid& grid = obj->get_grid();
	Selector& sel = obj->get_selector();
	SelectBoundaryElements(sel, grid.begin<Face>(), grid.end<Face>());
}

inline void SelectInnerFaces(MeshObject* obj)
{
	Grid& grid = obj->get_grid();
	Selector& sel = obj->get_selector();
	SelectInnerElements(sel, grid.begin<Face>(), grid.end<Face>());
}

inline void SelectLinkedManifoldFaces(MeshObject* obj)
{
	Grid& grid = obj->get_grid();
	Selector& sel = obj->get_selector();
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

inline void SelectLinkedBoundaryFaces(MeshObject* obj, bool stopAtSelectedEdges)
{
	Grid& grid = obj->get_grid();
	Selector& sel = obj->get_selector();

	if(stopAtSelectedEdges)
		SelectLinkedElements<Face>(sel, IsOnBoundary(grid), IsNotSelected(sel));
	else
		SelectLinkedElements<Face>(sel, IsOnBoundary(grid));
}

inline void SelectDegenerateFaces(MeshObject* obj, number maxHeight)
{
	number maxHeightSq = maxHeight * maxHeight;
	Grid& grid = obj->get_grid();
	Selector& sel = obj->get_selector();
	MeshObject::position_accessor_t& aaPos = obj->position_accessor();

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

inline void SelectLinkedFlatFaces(MeshObject* obj, number maxDeviationAngle, bool traverseFlipped,
						   bool traverseDegeneratedFaces, bool stopAtSelectedEdges)
{
	Selector& sel = obj->get_selector();
	if(traverseDegeneratedFaces)
		ug::SelectLinkedFlatAndDegeneratedFaces(sel, maxDeviationAngle,
												traverseFlipped,
												stopAtSelectedEdges);
	else
		ug::SelectLinkedFlatFaces(sel, maxDeviationAngle, traverseFlipped,
								  stopAtSelectedEdges);
}

inline void SelectIntersectingTriangles(MeshObject* obj)
{
	using namespace ug::node_tree;

	Grid& grid = obj->get_grid();
	Selector& sel = obj->get_selector();

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
	MeshObject::position_accessor_t& aaPos = obj->position_accessor();

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

inline void SelectAssociatedFaces(MeshObject* obj)
{
	Selector& sel = obj->get_selector();
	SelectAssociatedFaces(sel, sel.begin<Volume>(), sel.end<Volume>());
}

inline void SelectAllFaces(MeshObject* obj)
{
	Grid& grid = obj->get_grid();
	Selector& sel = obj->get_selector();
	sel.select(grid.faces_begin(), grid.faces_end());
}

inline void DeselectAllFaces(MeshObject* obj)
{
	Selector& sel = obj->get_selector();
	sel.deselect(sel.faces_begin(), sel.faces_end());
}

inline bool SelectFaceByIndex(MeshObject* obj, int index)
{
	Grid& grid = obj->get_grid();
	int counter = 0;

	for(FaceIterator iter = grid.begin<Face>();
		iter != grid.end<Face>(); ++iter, ++counter)
	{
		if(counter == index){
			obj->get_selector().select(*iter);
			return true;
		}
	}
	return false;
}

inline void SelectFacesByNormal(MeshObject* obj, const vector3& refNormal, number maxDeviationAngle)
{
	number dotThreshold = cos(deg_to_rad(maxDeviationAngle));
	Grid& grid = obj->get_grid();
	Selector& sel = obj->get_selector();
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

inline void FaceSelectionFill(MeshObject* obj)
{
	Selector& sel = obj->get_selector();
	SelectionFill<Face>(sel, IsSelected(sel));
}

inline size_t SelectBentQuadrilaterals(MeshObject* obj, number dotThreshold)
{
	Grid& grid = obj->get_grid();
	Selector& sel = obj->get_selector();
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
inline void SelectAllVolumes(MeshObject* obj)
{
	Grid& grid = obj->get_grid();
	Selector& sel = obj->get_selector();
	sel.select(grid.volumes_begin(), grid.volumes_end());
}

inline void DeselectAllVolumes(MeshObject* obj)
{
	Selector& sel = obj->get_selector();
	sel.deselect(sel.volumes_begin(), sel.volumes_end());
}

inline int SelectUnorientableVolumes(MeshObject* obj)
{
	Grid& grid = obj->get_grid();
	Selector& sel = obj->get_selector();

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

inline bool SelectVolumeByIndex(MeshObject* obj, int index)
{
	Grid& grid = obj->get_grid();
	int counter = 0;

	for(VolumeIterator iter = grid.begin<Volume>();
		iter != grid.end<Volume>(); ++iter, ++counter)
	{
		if(counter == index){
			obj->get_selector().select(*iter);
			return true;
		}
	}
	return false;
}

inline void VolumeSelectionFill(MeshObject* obj)
{
	Selector& sel = obj->get_selector();
	SelectionFill<Volume>(sel, IsSelected(sel));
}

inline void MarkSelection(MeshObject* obj)
{
	obj->get_crease_handler().assign_subset(
		obj->get_selector().begin<Vertex>(),
		obj->get_selector().end<Vertex>(),
		ug::REM_FIXED);
	obj->get_crease_handler().assign_subset(
		obj->get_selector().begin<Edge>(),
		obj->get_selector().end<Edge>(),
		ug::REM_CREASE);
}

inline void UnmarkSelection(MeshObject* obj)
{
	obj->get_crease_handler().assign_subset(
		obj->get_selector().begin<Vertex>(),
		obj->get_selector().end<Vertex>(),
		ug::REM_NONE);

	obj->get_crease_handler().assign_subset(
		obj->get_selector().begin<Edge>(),
		obj->get_selector().end<Edge>(),
		ug::REM_NONE);
}

}}// end of namespace

#endif
