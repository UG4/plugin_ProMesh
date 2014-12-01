// created by Sebastian Reiter
// s.b.reiter@gmail.com
// Jun 24, 2013 (d,m,y)

#ifndef __H__UG_PROMESH__subset_tools__
#define __H__UG_PROMESH__subset_tools__

#include <vector>
#include <queue>
#include "../mesh.h"
#include "lib_grid/algorithms/subset_util.h"
#include "lib_grid/algorithms/grid_statistics.h"

namespace ug{
namespace promesh{

inline void AssignSubset(Mesh* obj, int newIndex)
{
	Selector& sel = obj->selector();
	SubsetHandler& sh = obj->subset_handler();

	sh.assign_subset(sel.begin<Vertex>(),
					 sel.end<Vertex>(), newIndex);
	sh.assign_subset(sel.begin<Edge>(),
					 sel.end<Edge>(), newIndex);
	sh.assign_subset(sel.begin<Face>(),
					 sel.end<Face>(), newIndex);
	sh.assign_subset(sel.begin<Volume>(),
					 sel.end<Volume>(), newIndex);
}

inline void SetSubsetName(Mesh* obj, int si, const char* name)
{
	SubsetHandler& sh = obj->subset_handler();
	sh.subset_info(si).name = name;
}

inline void AssignSubsetColors(Mesh* obj)
{
	AssignSubsetColors(obj->subset_handler());
}

inline void SeparateFacesByEdgeSubsets(Mesh* obj)
{
	Grid& grid = obj->grid();
	SubsetHandler& sh = obj->subset_handler();
	SeparateSubsetsByLowerDimSubsets<Face>(grid, sh);
}

inline void SeparateFacesBySelectedEdges(Mesh* obj)
{
	Grid& grid = obj->grid();
	SubsetHandler& sh = obj->subset_handler();
	Selector& sel = obj->selector();
	SeparateSubsetsByLowerDimSelection<Face>(grid, sh, sel);
}

inline void SeparateVolumesByFaceSubsets(Mesh* obj)
{
	Grid& grid = obj->grid();
	SubsetHandler& sh = obj->subset_handler();
	SeparateSubsetsByLowerDimSubsets<Volume>(grid, sh);
}

inline void SeparateVolumesBySelectedFaces(Mesh* obj)
{
	Grid& grid = obj->grid();
	SubsetHandler& sh = obj->subset_handler();
	Selector& sel = obj->selector();
	SeparateSubsetsByLowerDimSelection<Volume>(grid, sh, sel);
}

inline void SeparateIrregularManifoldSubsets(Mesh* obj)
{
	SubsetHandler& sh = obj->subset_handler();
	for(int i = 0; i < sh.num_subsets(); ++i){
		int firstFree = GetMaxSubsetIndex<Face>(sh) + 1;
		SplitIrregularManifoldSubset(sh, i, firstFree);
	}
}

inline void MoveSubset(Mesh* obj, int oldIndex, int newIndex)
{
	if(newIndex != oldIndex){
		SubsetHandler& sh = obj->subset_handler();
		sh.move_subset(oldIndex, newIndex);
	}
}

inline void SwapSubsets(Mesh* obj, int oldIndex, int newIndex)
{
	SubsetHandler& sh = obj->subset_handler();
	if(newIndex != oldIndex && newIndex < sh.num_subsets()
		&& oldIndex < sh.num_subsets())
	{
		sh.swap_subsets(oldIndex, newIndex);
	}
}

inline void JoinSubsets(Mesh* obj, int target, int si1, int si2, bool eraseUnused)
{
	SubsetHandler& sh = obj->subset_handler();
	sh.join_subsets(target, si1, si2, eraseUnused);
}

inline void EraseSubset(Mesh* obj, int si, bool eraseGeometry)
{
	Grid& grid = obj->grid();
	SubsetHandler& sh = obj->subset_handler();

	if(si < sh.num_subsets())
	{
		if(eraseGeometry){
			grid.erase(sh.begin<Volume>(si), sh.end<Volume>(si));
			grid.erase(sh.begin<Face>(si), sh.end<Face>(si));
			grid.erase(sh.begin<Edge>(si), sh.end<Edge>(si));
			grid.erase(sh.begin<Vertex>(si), sh.end<Vertex>(si));
		}
		sh.erase_subset(si);
	}
}

inline void EraseEmptySubsets(Mesh* obj)
{
	SubsetHandler& sh = obj->subset_handler();
	int i = 0;
	while(i < sh.num_subsets()){
		if(sh.empty(i))
			sh.erase_subset(i);
		else
			++i;
	}
}

inline void AdjustSubsetsForUG3(Mesh* obj, bool keepIntfSubs)
{
	AdjustSubsetsForLgmNg(obj->grid(), obj->subset_handler(), keepIntfSubs);
}

inline void AdjustSubsetsForUG4(Mesh* obj, bool preserveExistingSubsets)
{
	AdjustSubsetsForSimulation(obj->subset_handler(), preserveExistingSubsets);
}

inline void SeparateFaceSubsetsByNormal(Mesh* obj)
{
	ug::SeparateFaceSubsetsByNormal(obj->grid(), obj->subset_handler());
}

inline void SeparateFaceSubsetByNormal(Mesh* obj, int si)
{
	if(si < obj->subset_handler().num_subsets())
		SeparateFaceSubsetsByNormal(obj->grid(), obj->subset_handler(),
									obj->position_attachment(), NULL, si);
}

inline void AssignSubsetsByQuality(Mesh* obj, int numSections)
{
	using namespace std;
	std::vector<number> intervals;
	intervals.push_back(0);
	for(int i = 1; i < numSections; ++i)
		intervals.push_back((number)i / (number)numSections);
	intervals.push_back(1.);

	Grid& grid = obj->grid();
	Selector& sel = obj->selector();
	SubsetHandler& sh = obj->subset_handler();

	ug::AssignSubsetsByQuality(grid, sh, sel.begin<Face>(),
								sel.end<Face>(), intervals);

//	log how many faces were assigned to the different subsets.
//	since potentially only a subset of faces has been considered,
//	we may not simply output the subset sizes.

	UG_LOG("Assigned faces to subsets:\n");
	for(size_t i = 0; i < intervals.size() - 1; ++i){
	//	count the number of selected faces in this section
		size_t counter = 0;
		for(FaceIterator iter = sel.begin<Face>(); iter != sel.end<Face>(); ++iter)
		{
			if(sh.get_subset_index(*iter) == (int)i)
				++counter;
		}

		UG_LOG("  quality " << intervals[i] << " - " << intervals[i+1] << ": \t" << counter << "\n");
	}

	UG_LOG(endl);
}

inline void SeparateDegeneratedBoundaryFaceSubsets(Mesh* obj, number angle)
{
	using namespace std;
	number thresholdDot = cos(deg_to_rad(angle));

	Grid& g = obj->grid();
	SubsetHandler& sh = obj->subset_handler();
	Mesh::position_accessor_t& aaPos = obj->position_accessor();

	vector<Edge*> edges;
	vector<Edge*> edges2;
	vector<Face*> faces;
	vector<Face*> assembledSubset;
	queue<Face*> queFaces;

//	we'll use this vector to check whether we have to assign faces which
//	we assembled to a subset to a new subset or whether it can stay where it
//	is. The first assembled face-pack can always stay in its subset.
	vector<bool> vAssignNewSubset(sh.num_subsets(), false);

//	the index at which we'll add new subsets (increases during the algorithm)
	int newSI = GetMaxSubsetIndex<Face>(sh) + 1;

//	we use marks to mark all processed elements
	g.begin_marking();

//	iterate over all faces and search for a degenerated boundary face.
	for(FaceIterator iter = g.begin<Face>(); iter != g.end<Face>(); ++iter){
		Face* f = *iter;
	//	the face may have been marked during subset assembly below
		if(g.is_marked(f))
			continue;
		g.mark(f);

		if(IsDegenerated(f, aaPos)){
			if(IsVolumeBoundaryFace(g, f)){
			//	the face is a candidate. Get subset index and push it to the queue
				int origSI = sh.get_subset_index(f);
				queFaces.push(f);
				assembledSubset.clear();

				while(!queFaces.empty()){
					Face* curFace = queFaces.front();
					queFaces.pop();
					assembledSubset.push_back(curFace);

				//	check all degenerated neighbor faces
					CollectAssociated(edges, g, curFace);

					vector3 dir(0, 0, 0);
					bool gotOne = false;
				//	the first non-degenerated edge defines the direction of the face
					for(size_t i_edge = 0; i_edge < edges.size(); ++i_edge){
						Edge* e = edges[i_edge];
						if(VecDistanceSq(aaPos[e->vertex(0)], aaPos[e->vertex(1)]) >= SMALL*SMALL){
							VecSubtract(dir, aaPos[e->vertex(1)], aaPos[e->vertex(0)]);
							VecNormalize(dir, dir);
							gotOne = true;
							break;
						}
					}

				//	if we haven't found a non-degenerated edge, we won't continue.
					if(!gotOne)
						continue;

				//	now find associated degenerated faces
					for(size_t i_edge = 0; i_edge < edges.size(); ++i_edge){
						Edge* e = edges[i_edge];

					//	we have to know whether the edge is degenerated or not.
						bool bDegEdge = (VecDistanceSq(aaPos[e->vertex(0)], aaPos[e->vertex(1)]) < SMALL*SMALL);

						CollectAssociated(faces, g, e);

						for(size_t i_face = 0; i_face < faces.size(); ++i_face){
							Face* nbrFace = faces[i_face];
							if(!g.is_marked(nbrFace) && sh.get_subset_index(nbrFace) == origSI){
								if(IsVolumeBoundaryFace(g, nbrFace)){
									if(IsDegenerated(nbrFace, aaPos)){
									//	if the edge was non-degenerated, it is automatically part of the
									//	assembled subset.
										if(!bDegEdge){
											g.mark(nbrFace);
											queFaces.push(nbrFace);
										}
										else{
										//	we have to compare the directions of the faces
											CollectAssociated(edges2, g, nbrFace);

											vector3 dir2(0, 0, 0);
											bool gotOne2 = false;
										//	the first non-degenerated edge defines the direction of the face
											for(size_t i_edge = 0; i_edge < edges2.size(); ++i_edge){
												Edge* e = edges2[i_edge];
												if(VecDistanceSq(aaPos[e->vertex(0)], aaPos[e->vertex(1)]) >= SMALL*SMALL){
													VecSubtract(dir2, aaPos[e->vertex(1)], aaPos[e->vertex(0)]);
													VecNormalize(dir2, dir2);
													gotOne2 = true;
													break;
												}
											}

											if(gotOne2){
											//	we now got both directions. compare the angle.
												if(fabs(VecDot(dir, dir2)) >= thresholdDot){
												//	both are in the same subset
													g.mark(nbrFace);
													queFaces.push(nbrFace);
												}
											}
										}
									}
								}
							}
						}
					}

				}

			//	now add the assembledSubset to its new destination
				if(vAssignNewSubset[origSI]){
					sh.assign_subset(assembledSubset.begin(), assembledSubset.end(), newSI);
					++newSI;
				}
				else{
				//	if more degenerated faces are contained in this subset, they shall be
				//	assigned to another subset.
					vAssignNewSubset[origSI] = true;
				}
			}
		}
	}

	g.end_marking();
}

inline void AssignSubsetsByElementType(Mesh* obj)
{
	SubsetHandler& sh = obj->subset_handler();
	AssignSubsetsByElementType(sh);
	int i = 0;
	while(i < sh.num_subsets()){
		if(sh.empty(i))
			sh.erase_subset(i);
		else
			++i;
	}
}

inline void CopySubsetIndicesToSides(Mesh* obj, bool selectionOnly,
							  bool toUnassignedOnly)
{
	SubsetHandler& sh = obj->subset_handler();

	if(selectionOnly){
		Selector& sel = obj->selector();
		CopySubsetIndicesToSides(sh, sel.begin<Volume>(), sel.end<Volume>(), toUnassignedOnly);
		CopySubsetIndicesToSides(sh, sel.begin<Face>(), sel.end<Face>(), toUnassignedOnly);
		CopySubsetIndicesToSides(sh, sel.begin<Edge>(), sel.end<Edge>(), toUnassignedOnly);
		CopySubsetIndicesToSides(sh, sel.begin<Vertex>(), sel.end<Vertex>(), toUnassignedOnly);
	}
	else{
		for(int i = 0; i < sh.num_subsets(); ++i){
			CopySubsetIndicesToSides(sh, sh.begin<Volume>(i), sh.end<Volume>(i), toUnassignedOnly);
			CopySubsetIndicesToSides(sh, sh.begin<Face>(i), sh.end<Face>(i), toUnassignedOnly);
			CopySubsetIndicesToSides(sh, sh.begin<Edge>(i), sh.end<Edge>(i), toUnassignedOnly);
			CopySubsetIndicesToSides(sh, sh.begin<Vertex>(i), sh.end<Vertex>(i), toUnassignedOnly);
		}
	}
}

}}// end of namespace

#endif
