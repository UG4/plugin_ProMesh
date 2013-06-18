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

void ClearSelection(MeshObject* obj)
{
	obj->get_selector().clear();
}

void SelectAll(MeshObject* obj)
{
	Grid& grid = obj->get_grid();
	Selector& sel = obj->get_selector();
	sel.select(grid.vertices_begin(), grid.vertices_end());
	sel.select(grid.edges_begin(), grid.edges_end());
	sel.select(grid.faces_begin(), grid.faces_end());
	sel.select(grid.volumes_begin(), grid.volumes_end());
}

void ExtendSelection(MeshObject* obj, int neighborhoodSize)
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

void SelectSubset(MeshObject* obj, int si, bool selVrts, bool selEdges,
				  bool selFaces, bool selVols)
{
	SubsetHandler& sh = obj->get_subset_handler();
	Selector& sel = obj->get_selector();

	if(si >= 0){
		if(selVrts)
			sel.select(sh.begin<VertexBase>(si), sh.end<VertexBase>(si));
		if(selEdges)
			sel.select(sh.begin<EdgeBase>(si), sh.end<EdgeBase>(si));
		if(selFaces)
			sel.select(sh.begin<Face>(si), sh.end<Face>(si));
		if(selVols)
			sel.select(sh.begin<Volume>(si), sh.end<Volume>(si));
	}
	else{
		Grid& grid = obj->get_grid();
	//	subset -1 has to be selected. Those are not directly accessible.
		if(selVrts){
			for(VertexBaseIterator iter = grid.vertices_begin();
				iter != grid.vertices_end(); ++iter)
			{
				if(sh.get_subset_index(*iter) == -1)
					sel.select(*iter);
			}
		}
		if(selEdges){
			for(EdgeBaseIterator iter = grid.edges_begin();
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

void SelectSubsetBoundary(MeshObject* obj, int si, bool edgeBnds,
						  bool faceBnds, bool volBnds)
{
	SubsetHandler& sh = obj->get_subset_handler();
	Selector& sel = obj->get_selector();

	if(edgeBnds)
		SelectAreaBoundary(sel, sh.begin<EdgeBase>(si), sh.end<EdgeBase>(si));
	if(faceBnds)
		SelectAreaBoundary(sel, sh.begin<Face>(si), sh.end<Face>(si));
	if(volBnds)
		SelectAreaBoundary(sel, sh.begin<Volume>(si), sh.end<Volume>(si));
}

template <class TGeomObj>
static void SelectUnassignedElementsHelper(Grid& grid, SubsetHandler& sh, Selector& sel)
{
	typedef typename geometry_traits<TGeomObj>::iterator	iterator;
	for(iterator iter = grid.begin<TGeomObj>(); iter != grid.end<TGeomObj>(); ++iter)
	{
		if(sh.get_subset_index(*iter) == -1){
			sel.select(*iter);
		}
	}
}

void SelectUnassignedElements(MeshObject* obj, bool selVrts, bool selEdges,
							  bool selFaces, bool selVols)
{
	Grid& grid = obj->get_grid();
	SubsetHandler& sh = obj->get_subset_handler();
	Selector& sel = obj->get_selector();

	if(selVrts)
		SelectUnassignedElementsHelper<VertexBase>(grid, sh, sel);
	if(selEdges)
		SelectUnassignedElementsHelper<EdgeBase>(grid, sh, sel);
	if(selFaces)
		SelectUnassignedElementsHelper<Face>(grid, sh, sel);
	if(selVols)
		SelectUnassignedElementsHelper<Volume>(grid, sh, sel);
}

void InvertSelection(MeshObject* obj, bool invVrts, bool invEdges,
					 bool invFaces, bool invVols)
{
	Grid& grid = obj->get_grid();
	Selector& sel = obj->get_selector();

	if(invVrts)
		InvertSelection(sel, grid.begin<VertexBase>(),
							grid.end<VertexBase>());

	if(invEdges)
		InvertSelection(sel, grid.begin<EdgeBase>(),
							grid.end<EdgeBase>());

	if(invFaces)
		InvertSelection(sel, grid.begin<Face>(),
							grid.end<Face>());

	if(invVols)
		InvertSelection(sel, grid.begin<Volume>(),
							grid.end<Volume>());
}

void SelectSelectionBoundary(MeshObject* obj)
{
	Selector& sel = obj->get_selector();

	if(sel.num<Volume>() > 0)
		SelectAreaBoundary(sel, sel.begin<Volume>(), sel.end<Volume>());
	if(sel.num<Face>() > 0)
		SelectAreaBoundary(sel, sel.begin<Face>(), sel.end<Face>());
}

void CloseSelection(MeshObject* obj)
{
	Selector& sel = obj->get_selector();
	SelectAssociatedFaces(sel, sel.begin<Volume>(), sel.end<Volume>());
	SelectAssociatedEdges(sel, sel.begin<Face>(), sel.end<Face>());
	SelectAssociatedVertices(sel, sel.begin<EdgeBase>(), sel.end<EdgeBase>());
}


////////////////////////////////////////////////////////////////////////////////
//	VERTICES
void SelectBoundaryVertices(MeshObject* obj)
{
	Grid& grid = obj->get_grid();
	Selector& sel = obj->get_selector();
	ug::SelectBoundaryElements(sel, grid.begin<VertexBase>(), grid.end<VertexBase>());
}

void SelectInnerVertices(MeshObject* obj)
{
	Grid& grid = obj->get_grid();
	Selector& sel = obj->get_selector();
	SelectInnerElements(sel, grid.begin<VertexBase>(), grid.end<VertexBase>());
}

void SelectAssociatedVertices(MeshObject* obj)
{
	Selector& sel = obj->get_selector();
	SelectAssociatedVertices(sel, sel.begin<EdgeBase>(), sel.end<EdgeBase>());
	SelectAssociatedVertices(sel, sel.begin<Face>(), sel.end<Face>());
	SelectAssociatedVertices(sel, sel.begin<Volume>(), sel.end<Volume>());
}

void SelectAllVertices(MeshObject* obj)
{
	Grid& grid = obj->get_grid();
	Selector& sel = obj->get_selector();
	sel.select(grid.vertices_begin(), grid.vertices_end());
}

void DeselectAllVertices(MeshObject* obj)
{
	Selector& sel = obj->get_selector();
	sel.deselect(sel.vertices_begin(), sel.vertices_end());
}

void SelectMarkedVertices(MeshObject* obj)
{
	obj->get_selector().select(
		obj->get_crease_handler().begin<VertexBase>(REM_FIXED),
		obj->get_crease_handler().end<VertexBase>(REM_FIXED));
}

bool SelectVertexByIndex(MeshObject* obj, int index)
{
	Grid& grid = obj->get_grid();
	int counter = 0;

	VertexBaseIterator iter = grid.begin<VertexBase>();
	while(counter < index && iter != grid.end<VertexBase>()){
		++counter;
		++iter;
	}

	if(counter == index){
		obj->get_selector().select(*iter);
		return true;
	}
	return false;
}

template <class TElem>
size_t SelectUnconnectedVerticesHelper(Grid& grid, Selector& sel)
{
	using namespace ug;
	std::vector<TElem*> elems;

	size_t numUnconnected = 0;
	for(VertexBaseIterator iter = grid.vertices_begin();
		iter != grid.vertices_end(); ++iter)
	{
		if(!sel.is_selected(*iter)){
			CollectAssociated(elems, grid, *iter);
			if(elems.size() == 0){
				sel.select(*iter);
				numUnconnected++;
			}
		}
	}
	return numUnconnected;
}

size_t SelectUnconnectedVertices(MeshObject* obj, bool edgeCons, bool faceCons, bool volCons)
{
	Grid& grid = obj->get_grid();
	Selector& sel = obj->get_selector();

	size_t numUnconnected = 0;
	if(edgeCons)
		numUnconnected += SelectUnconnectedVerticesHelper<EdgeBase>(grid, sel);
	if(faceCons)
		numUnconnected += SelectUnconnectedVerticesHelper<Face>(grid, sel);
	if(volCons)
		numUnconnected += SelectUnconnectedVerticesHelper<Volume>(grid, sel);
	return numUnconnected;
}

////////////////////////////////////////////////////////////////////////////////
//	EDGES
void SelectBoundaryEdges(MeshObject* obj)
{
	Grid& grid = obj->get_grid();
	Selector& sel = obj->get_selector();
	SelectBoundaryElements(sel, grid.begin<EdgeBase>(), grid.end<EdgeBase>());
}

void SelectInnerEdges(MeshObject* obj)
{
	Grid& grid = obj->get_grid();
	Selector& sel = obj->get_selector();
	SelectInnerElements(sel, grid.begin<EdgeBase>(), grid.end<EdgeBase>());
}

void SelectNonManifoldEdges(MeshObject* obj)
{
	Grid& grid = obj->get_grid();
	Selector& sel = obj->get_selector();

//	iterate over all edges and check how many adjacent faces each has.
//	if there are more than 2, the edge will be selected
	for(EdgeBaseIterator iter = grid.edges_begin();
		iter != grid.edges_end(); ++iter)
	{
		if(NumAssociatedFaces(grid, *iter) != 2)
			sel.select(*iter);
	}
}

void SelectSmoothEdgePath(MeshObject* obj, number maxDeviation, bool stopAtSelectedVrts)
{
	SelectSmoothEdgePath(obj->get_selector(), maxDeviation, stopAtSelectedVrts);
}

void SelectShortEdges(MeshObject* obj, number maxLength)
{
	number maxLengthSq = maxLength * maxLength;
	Grid& grid = obj->get_grid();
	Selector& sel = obj->get_selector();
	MeshObject::position_accessor_t& aaPos = obj->position_accessor();

	for(EdgeBaseIterator iter = grid.begin<EdgeBase>();
		iter != grid.end<EdgeBase>(); ++iter)
	{
		EdgeBase* e = *iter;
		if(VecDistanceSq(aaPos[e->vertex(0)], aaPos[e->vertex(1)]) < maxLengthSq)
			sel.select(e);
	}
}

void SelectLongEdges(MeshObject* obj, number minLength)
{
	number minLengthSq = minLength * minLength;
	Grid& grid = obj->get_grid();
	Selector& sel = obj->get_selector();
	MeshObject::position_accessor_t& aaPos = obj->position_accessor();

	for(EdgeBaseIterator iter = grid.begin<EdgeBase>();
		iter != grid.end<EdgeBase>(); ++iter)
	{
		EdgeBase* e = *iter;
		if(VecDistanceSq(aaPos[e->vertex(0)], aaPos[e->vertex(1)]) > minLengthSq)
			sel.select(e);
	}
}

void SelectLinkedBoundaryEdges(MeshObject* obj, bool stopAtSelectedVrts)
{
	Grid& grid = obj->get_grid();
	Selector& sel = obj->get_selector();

	if(stopAtSelectedVrts)
		SelectLinkedElements<EdgeBase>(sel, IsOnBoundary(grid), IsNotSelected(sel));
	else
		SelectLinkedElements<EdgeBase>(sel, IsOnBoundary(grid));
}

void SelectAssociatedEdges(MeshObject* obj)
{
	Selector& sel = obj->get_selector();
	SelectAssociatedEdges(sel,sel.begin<Face>(), sel.end<Face>());
	SelectAssociatedEdges(sel, sel.begin<Volume>(), sel.end<Volume>());
}

void SelectAllEdges(MeshObject* obj)
{
	Grid& grid = obj->get_grid();
	Selector& sel = obj->get_selector();
	sel.select(grid.edges_begin(), grid.edges_end());
}

void DeselectAllEdges(MeshObject* obj)
{
	Selector& sel = obj->get_selector();
	sel.deselect(sel.edges_begin(), sel.edges_end());
}

void SelectMarkedEdges(MeshObject* obj)
{
	obj->get_selector().select(
		obj->get_crease_handler().begin<EdgeBase>(REM_CREASE),
		obj->get_crease_handler().end<EdgeBase>(REM_CREASE));
}

bool SelectEdgeByIndex(MeshObject* obj, int index)
{
	Grid& grid = obj->get_grid();
	int counter = 0;

	EdgeBaseIterator iter = grid.begin<EdgeBase>();
	while(counter < index && iter != grid.end<EdgeBase>()){
		++counter;
		++iter;
	}

	if(counter == index){
		obj->get_selector().select(*iter);
		return true;
	}
	return false;
}

void EdgeSelectionFill(MeshObject* obj)
{
	Selector& sel = obj->get_selector();
	SelectionFill<EdgeBase>(sel, IsSelected(sel));
}

////////////////////////////////////////////////////////////////////////////////
//	FACES
void SelectBoundaryFaces(MeshObject* obj)
{
	Grid& grid = obj->get_grid();
	Selector& sel = obj->get_selector();
	SelectBoundaryElements(sel, grid.begin<Face>(), grid.end<Face>());
}

void SelectInnerFaces(MeshObject* obj)
{
	Grid& grid = obj->get_grid();
	Selector& sel = obj->get_selector();
	SelectInnerElements(sel, grid.begin<Face>(), grid.end<Face>());
}

void SelectLinkedManifoldFaces(MeshObject* obj)
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
	std::vector<EdgeBase*> edges;
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

void SelectLinkedBoundaryFaces(MeshObject* obj, bool stopAtSelectedEdges)
{
	Grid& grid = obj->get_grid();
	Selector& sel = obj->get_selector();

	if(stopAtSelectedEdges)
		SelectLinkedElements<Face>(sel, IsOnBoundary(grid), IsNotSelected(sel));
	else
		SelectLinkedElements<Face>(sel, IsOnBoundary(grid));
}

void SelectDegenerateFaces(MeshObject* obj, number maxHeight)
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

void SelectLinkedFlatFaces(MeshObject* obj, number maxDeviationAngle, bool traverseFlipped,
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

void SelectIntersectingTriangles(MeshObject* obj)
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

void SelectAssociatedFaces(MeshObject* obj)
{
	Selector& sel = obj->get_selector();
	SelectAssociatedFaces(sel, sel.begin<Volume>(), sel.end<Volume>());
}

void SelectAllFaces(MeshObject* obj)
{
	Grid& grid = obj->get_grid();
	Selector& sel = obj->get_selector();
	sel.select(grid.faces_begin(), grid.faces_end());
}

void DeselectAllFaces(MeshObject* obj)
{
	Selector& sel = obj->get_selector();
	sel.deselect(sel.faces_begin(), sel.faces_end());
}

bool SelectFaceByIndex(MeshObject* obj, int index)
{
	Grid& grid = obj->get_grid();
	int counter = 0;

	FaceIterator iter = grid.begin<Face>();
	while(counter < index && iter != grid.end<Face>()){
		++counter;
		++iter;
	}

	if(counter == index){
		obj->get_selector().select(*iter);
		return true;
	}
	return false;
}

void FaceSelectionFill(MeshObject* obj)
{
	Selector& sel = obj->get_selector();
	SelectionFill<Face>(sel, IsSelected(sel));
}

size_t SelectBentQuadrilaterals(MeshObject* obj, number dotThreshold)
{
	Grid& grid = obj->get_grid();
	Selector& sel = obj->get_selector();
	Grid::AttachmentAccessor<VertexBase, APosition> aaPos(grid, aPosition);

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
void SelectAllVolumes(MeshObject* obj)
{
	Grid& grid = obj->get_grid();
	Selector& sel = obj->get_selector();
	sel.select(grid.volumes_begin(), grid.volumes_end());
}

void DeselectAllVolumes(MeshObject* obj)
{
	Selector& sel = obj->get_selector();
	sel.deselect(sel.volumes_begin(), sel.volumes_end());
}

int SelectUnorientableVolumes(MeshObject* obj)
{
	Grid& grid = obj->get_grid();
	Selector& sel = obj->get_selector();

	Grid::AttachmentAccessor<VertexBase, APosition> aaPos(grid, aPosition);

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

bool SelectVolumeByIndex(MeshObject* obj, int index)
{
	Grid& grid = obj->get_grid();
	int counter = 0;

	VolumeIterator iter = grid.begin<Volume>();
	while(counter < index && iter != grid.end<Volume>()){
		++counter;
		++iter;
	}

	if(counter == index){
		obj->get_selector().select(*iter);
		return true;
	}
	return false;
}

void VolumeSelectionFill(MeshObject* obj)
{
	Selector& sel = obj->get_selector();
	SelectionFill<Volume>(sel, IsSelected(sel));
}

}}// end of namespace

#endif
