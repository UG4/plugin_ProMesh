// created by Sebastian Reiter
// s.b.reiter@gmail.com
// Jun 20, 2013 (d,m,y)

#ifndef __H__UG__refinement_tools__
#define __H__UG__refinement_tools__

#include "../mesh_object.h"

namespace ug{
namespace promesh{

void Refine(MeshObject* obj, bool strictSubsetInheritance)
{
	Grid& grid = obj->get_grid();
	Selector& sel = obj->get_selector();
	SubsetHandler& sh = obj->get_subset_handler();
	bool siEnabled = sh.strict_inheritance_enabled();
	sh.enable_strict_inheritance(strictSubsetInheritance);

	Refine(grid, sel);

	sh.enable_strict_inheritance(siEnabled);
}

void HangingNodeRefine(MeshObject* obj, bool strictSubsetInheritance, bool anisotropic)
{
	Grid& grid = obj->get_grid();
	Selector& sel = obj->get_selector();
	SubsetHandler& sh = obj->get_subset_handler();
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

void RefineSmooth(MeshObject* obj, bool strictSubsetInheritance)
{
	Grid& grid = obj->get_grid();
	Selector& sel = obj->get_selector();
	SubsetHandler& sh = obj->get_subset_handler();
	bool siEnabled = sh.strict_inheritance_enabled();

	sh.enable_strict_inheritance(strictSubsetInheritance);

//	currently only triangles are supported in smooth refinement.
//	convert all selected quads to triangles first.
	Grid::VertexAttachmentAccessor<APosition> aaPos(grid, aPosition);
	Triangulate(grid, sel.begin<Quadrilateral>(), sel.end<Quadrilateral>(), &aaPos);

//	since we use a flat hierarchy, a temporary position attachment is required
	APosition aTmpPos;
	grid.attach_to_vertices(aTmpPos);

	RefinementCallbackSubdivisionLoop<APosition>
		refCallbackLoop(grid, aPosition, aTmpPos);

	Refine(grid, sel, &refCallbackLoop);

//	copy position data of selected vertices
	CopyAttachments(grid, sel.begin<VertexBase>(),
						sel.end<VertexBase>(),
						aTmpPos, aPosition);


	grid.detach_from_vertices(aTmpPos);

	sh.enable_strict_inheritance(siEnabled);
}

void RefineSmoothBoundary2D(MeshObject* obj, bool strictSubsetInheritance)
{
	Grid& grid = obj->get_grid();
	Selector& sel = obj->get_selector();
	SubsetHandler& sh = obj->get_subset_handler();
	bool siEnabled = sh.strict_inheritance_enabled();

	sh.enable_strict_inheritance(strictSubsetInheritance);

//	since we use a flat hierarchy, a temporary position attachment is required
	APosition aTmpPos;
	grid.attach_to_vertices(aTmpPos);

	RefinementCallbackSubdivBoundary<APosition>
		refCallbackLoopBnd(grid, aPosition, aTmpPos);

	Refine(grid, sel, &refCallbackLoopBnd);

	sh.enable_strict_inheritance(siEnabled);

//	copy position data of selected vertices
	CopyAttachments(grid, sel.begin<VertexBase>(),
						sel.end<VertexBase>(),
						aTmpPos, aPosition);

	grid.detach_from_vertices(aTmpPos);
}

void CreateFractal(MeshObject* obj, size_t numIterations, number scaleFac)
{
	Grid& grid = obj->get_grid();

//	we'll use a hanging-node refiner
	RefinementCallbackFractal refCallback(grid, scaleFac);
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

void InsertCenter(MeshObject* obj, bool strictSubsetInheritance)
{
	Grid& grid = obj->get_grid();
	Selector& sel = obj->get_selector();
	SubsetHandler& sh = obj->get_subset_handler();
	bool siEnabled = sh.strict_inheritance_enabled();

	sh.enable_strict_inheritance(strictSubsetInheritance);

//	access position attachment
	Grid::VertexAttachmentAccessor<APosition> aaPos(grid, aPosition);

	std::vector<EdgeBase*> edges;
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
		if(sel.num<EdgeBase>() > 0){
			UG_LOG("InsertCenter for edges is currently not supported if"
					" faces are present.\n");
		}
	}
	else
		edges.assign(sel.begin<EdgeBase>(), sel.end<EdgeBase>());

//	insert centers
	for(size_t i = 0; i < vols.size(); ++i){
		Volume* vol = vols[i];
		Vertex* vrt = *grid.create<Vertex>(vol);
		aaPos[vrt] = CalculateCenter(vol, aaPos);
		InsertCenterVertex(grid, vol, vrt, true);
	}

//	insert centers
	for(size_t i = 0; i < faces.size(); ++i){
		Face* f = faces[i];
		Vertex* vrt = *grid.create<Vertex>(f);
		aaPos[vrt] = CalculateCenter(f, aaPos);
		InsertCenterVertex(grid, f, vrt, true);
	}

//	split edges
	for(size_t i = 0; i < edges.size(); ++i){
		EdgeBase* e = edges[i];
		vector3 center = CalculateCenter(e, aaPos);
		Vertex* vrt = SplitEdge<Vertex>(grid, e);
		aaPos[vrt] = center;
	}

	sh.enable_strict_inheritance(siEnabled);
}

}}// end of namespace

#endif
