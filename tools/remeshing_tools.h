// created by Sebastian Reiter
// s.b.reiter@gmail.com
// Jun 20, 2013 (d,m,y)

#ifndef __H__UG__remeshing_tools__
#define __H__UG__remeshing_tools__

#include <vector>
#include "../mesh_object.h"
#include "lib_grid/algorithms/remeshing/delaunay_triangulation.h"
#include "lib_grid/algorithms/extrusion/extrusion.h"
#include "lib_grid/algorithms/grid_generation/tetrahedralization.h"
#include "lib_grid/algorithms/remeshing/grid_adaption.h"
#include "lib_grid/algorithms/remeshing/edge_length_adjustment.h"
#include "lib_grid/algorithms/duplicate.h"

namespace ug{
namespace promesh{

void ConvertToTriangles(MeshObject* obj)
{
	Grid& grid = obj->get_grid();
	Selector& sel = obj->get_selector();
	MeshObject::position_accessor_t& aaPos = obj->position_accessor();

	Triangulate(grid, sel.begin<Quadrilateral>(),
				sel.end<Quadrilateral>(), &aaPos);
}

void TriangleFill(MeshObject* obj, bool qualityGeneration, number minAngle, int si)
{
	if(minAngle < 0)
		minAngle = 0;
	if(minAngle > 30){
		UG_LOG("WARNING in TriangleFill: Restricting minAngle to 30, since the "
				"algorithm may not terminate for minAngle > 30.\n")
	}

	Grid& grid = obj->get_grid();
	SubsetHandler& sh = obj->get_subset_handler();
	Selector& sel = obj->get_selector();

//	if no edges are selected, nothing can be triangulated
	if(sel.num<EdgeBase>() < 3){
		UG_LOG("ERROR in TriangleFill: A closed outer edge-chain has to be selected.\n");
		return;
	}

//	before triangulating, we'll make sure that no double-edges exist
//	in the current selection.
	RemoveDoubleEdges(grid, sel.begin<EdgeBase>(), sel.end<EdgeBase>());

	MeshObject::position_accessor_t& aaPos = obj->position_accessor();
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
		for(VertexBaseIterator iter = grid.vertices_begin();
			iter != grid.vertices_end(); ++iter)
		{
			vrts.push_back(aaPos[*iter]);
		}
		std::vector<vector2> vrts2d(vrts.size());
		TransformPointSetTo2D(&vrts2d.front(), &vrts.front(),
							  vrts.size());

		size_t counter = 0;
		for(VertexBaseIterator iter = grid.vertices_begin();
			iter != grid.vertices_end(); ++iter, ++counter)
		{
			aaPos[*iter] = vector3(vrts2d[counter].x, vrts2d[counter].y, 0);
		}

		SaveGridToFile(grid, ss2d.str().c_str(), obj->position_attachment());

		counter = 0;
		for(VertexBaseIterator iter = grid.vertices_begin();
			iter != grid.vertices_end(); ++iter, ++counter)
		{
			aaPos[*iter] = vector3(vrts[counter].x, vrts[counter].y, 0);
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

void Retriangulate(MeshObject* obj, number minAngle)
{
	Grid& g = obj->get_grid();
	Selector& sel = obj->get_selector();
	SubsetHandler& creases = obj->get_crease_handler();

	MeshObject::position_accessor_t& aaPos = obj->position_accessor();

	QualityGridGeneration(g, sel.begin<Triangle>(), sel.end<Triangle>(),
						  aaPos, minAngle, IsNotInSubset(creases, -1));
}

void AdjustEdgeLength(MeshObject* obj, number minEdgeLen, number maxEdgeLen,
					  int numIterations, bool adaptive, bool automarkBoundaries)
{
	Grid& grid = obj->get_grid();
	SubsetHandler& shCrease = obj->get_crease_handler();

	if(automarkBoundaries){
		for(EdgeBaseIterator iter = grid.begin<EdgeBase>();
			iter != grid.end<EdgeBase>(); ++iter)
		{
			if(IsBoundaryEdge2D(grid, *iter))
				shCrease.assign_subset(*iter, REM_CREASE);
		}
	}

	AdjustEdgeLength(grid, shCrease, minEdgeLen, maxEdgeLen,
						 numIterations, true, adaptive);
}

void AdaptSurfaceToCylinder(MeshObject* obj, number radius, number threshold)
{
	using namespace std;
	Grid& g = obj->get_grid();
	Selector& sel = obj->get_selector();
	MeshObject::position_accessor_t& aaPos = obj->position_accessor();

//	store all source-vertices in a list
	vector<VertexBase*> vrts;
	vrts.assign(sel.begin<VertexBase>(), sel.end<VertexBase>());

//	iterate over selected vertices
	for(vector<VertexBase*>::iterator iter = vrts.begin();
		iter != vrts.end(); ++iter)
	{
		VertexBase* vrt = *iter;
		vector3 n;
		CalculateVertexNormal(n, g, vrt, aaPos);

		if(!ug::AdaptSurfaceGridToCylinder(sel, g, vrt, n, radius, threshold))
		{
			UG_LOG("AdaptSurfaceGridToCylinder failed for the vertex at " << aaPos[vrt] << "\n");
		}
	}
}

void Tetrahedralize(MeshObject* obj, number quality, bool preserveOuter, bool preserveAll,
					bool separateVolumes, bool appendSubsetsAtEnd)
{
	Grid& grid = obj->get_grid();
	SubsetHandler& sh = obj->get_subset_handler();
	UG_LOG("tetrahedralizing using 'tetgen' by Hang Si... ");
	ug::Tetrahedralize(grid, sh, quality, preserveOuter, preserveAll, obj->position_attachment());
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

void AssignVolumeConstraints(MeshObject* obj, number volConstraint)
{
	Selector& sel = obj->get_selector();
	MeshObject::volume_constraint_accessor_t& aaVolCon = obj->volume_constraint_accessor();

	for(Selector::traits<Volume>::iterator iter = sel.begin<Volume>();
		iter != sel.end<Volume>(); ++iter)
	{
		aaVolCon[*iter] = volConstraint;
	}
}

void ClearVolumeConstraints(MeshObject* obj)
{
	obj->clear_volume_constraints();
}

void Retetrahedralize(MeshObject* obj, number quality, bool preserveOuter,
					  bool preserveAll, bool applyVolumeConstraint)
{
	UG_LOG("retetrahedralizing using 'tetgen' by Hang Si... ");
	ug::Retetrahedralize(obj->get_grid(),
					obj->get_subset_handler(),
					obj->volume_constraint_attachment(),
					quality,
					preserveOuter, preserveAll,
					obj->position_attachment(),
					applyVolumeConstraint);
	UG_LOG("done.\n");
}

void Duplicate(MeshObject* obj, const vector3& offset, bool deselectOld, bool selectNew)
{
	Grid& grid = obj->get_grid();
	Selector& sel = obj->get_selector();
	ug::Duplicate(grid, sel, offset, obj->position_attachment(), deselectOld, selectNew);
}

void Extrude(MeshObject* obj, const vector3& totalDir, int numSteps, int si,
			 bool createFaces, bool createVolumes)
{
	using namespace std;
	vector3 stepDir;
	VecScale(stepDir, totalDir, 1./(float)numSteps);

	Grid& grid = obj->get_grid();
	Selector& sel = obj->get_selector();
	SubsetHandler& sh = obj->get_subset_handler();

	vector<VertexBase*> vrts;
	vrts.assign(sel.vertices_begin(), sel.vertices_end());
	vector<EdgeBase*> edges;
	edges.assign(sel.edges_begin(), sel.edges_end());
	vector<Face*> faces;
	faces.assign(sel.faces_begin(), sel.faces_end());

	uint extrusionOptions = 0;
	if(createFaces)
		extrusionOptions |= EO_CREATE_FACES;
	if(createVolumes)
		extrusionOptions |= EO_CREATE_VOLUMES;

//	we use sel to collect the newly created volumes
	sel.clear();
	sel.enable_autoselection(true);

//	mark all elements that were already in the selector.
	for(int i = 0; i < numSteps; ++i)
	{
		Extrude(grid, &vrts, &edges, &faces, stepDir,
					extrusionOptions, obj->position_attachment());
	}

	sel.enable_autoselection(false);
	sh.assign_subset(sel.volumes_begin(), sel.volumes_end(), si);
	sh.assign_subset(sel.faces_begin(), sel.faces_end(), si);
	sh.assign_subset(sel.edges_begin(), sel.edges_end(), si);
	sh.assign_subset(sel.vertices_begin(), sel.vertices_end(), si);


//	select faces, edges and vertices from the new top-layer.
	sel.clear<VertexBase>();
	sel.clear<EdgeBase>();
	sel.clear<Face>();
	sel.select(vrts.begin(), vrts.end());
	sel.select(edges.begin(), edges.end());
	sel.select(faces.begin(), faces.end());
}

void ExtrudeCylinders(MeshObject* obj, number height, number radius, number snapThreshold)
{
	using namespace std;
	Grid& g = obj->get_grid();
	Selector& sel = obj->get_selector();
	SubsetHandler& sh = obj->get_subset_handler();
	MeshObject::position_accessor_t& aaPos = obj->position_accessor();

//	store all source-vertices in a list
	vector<VertexBase*> vrts;
	vrts.assign(sel.begin<VertexBase>(), sel.end<VertexBase>());

//	iterate over selected vertices
	for(vector<VertexBase*>::iterator iter = vrts.begin();
		iter != vrts.end(); ++iter)
	{
		VertexBase* vrt = *iter;
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

}}// end of namespace

#endif
