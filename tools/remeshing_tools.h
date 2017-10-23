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

#ifndef __H__UG__remeshing_tools__
#define __H__UG__remeshing_tools__

#include <vector>
#include "../mesh.h"
#include "topology_tools.h"
#include "lib_grid/algorithms/duplicate.h"
#include "lib_grid/algorithms/grid_util.h"
#include "lib_grid/algorithms/quadrilateral_util.h"
#include "lib_grid/algorithms/remove_duplicates_util.h"
#include "lib_grid/algorithms/extrusion/extrusion.h"
#include "lib_grid/algorithms/geom_obj_util/face_util.h"
#include "lib_grid/algorithms/grid_generation/horizontal_layers_mesher.h"
#include "lib_grid/algorithms/grid_generation/tetrahedralization.h"
#include "lib_grid/algorithms/grid_generation/triangle_fill_sweep_line.h"
#include "lib_grid/algorithms/remeshing/delaunay_triangulation.h"
#include "lib_grid/algorithms/remeshing/grid_adaption.h"
#include "lib_grid/algorithms/remeshing/edge_length_adjustment.h"
#include "lib_grid/algorithms/remeshing/edge_length_adjustment_extended.h"
#include "lib_grid/algorithms/remeshing/simplification.h"
#include "lib_grid/algorithms/remeshing/simplify_polychain.h"
#include "lib_grid/algorithms/space_partitioning/lg_ntree.h"
#include "lib_grid/refinement/projectors/raster_layers_projector.h"

#define TOOLTIP_SIMPLIFY_POLYLINES "Removes vertices from the selected polyline which have a smaller curvature than the specified angle."
#define TOOLTIP_SIMPLIFY_SMOOTHED_POLYLINES "Removes vertices from the selected polyline which have a smaller smoothed curvature than the specified angle."
#define	TOOLTIP_CONVERT_TO_TRIANGLES "Converts selected quadrilaterals to triangles."
#define TOOLTIP_CONVERT_TO_QUADRILATERALS "Converts selected triangles or triangles connected by selected edges to quadrilaterals."
#define	TOOLTIP_TRIANGLE_FILL "Performs triangle fill using the sweep-line algorithm followed by an optional Constrained Delaunay retriangulation."
#define	TOOLTIP_RETRIANGULATE "Inserts vertices as required and performs Constrained Delaunay triangulation."
#define	TOOLTIP_ADJUST_EDGE_LENGTH "Remeshes the active grid so that all edges approximatly have a certain length."
#define	TOOLTIP_ADAPT_SURFACE_TO_CYLINDER "Introduces edges in a grid around a selected vertex which roughly correspond to the intersection of a cylinder with the surface."
#define TOOLTIP_REPLACE_VALENCE_3_VERTICES "Replaces selected valence-3 vertices by triangles, if the curvature of associated triangles is low"
#define TOOLTIP_REPLACE_LOW_VALENCE_VERTICES "Replaces selected valence-3 and valence-4 vertices by triangles, if the curvature of associated triangles is low"
#define	TOOLTIP_CONVERT_TO_TETRAHEDRA "Converts selected volume elements to tetrahedra."
#define	TOOLTIP_TETRAHEDRALIZE "Fills a closed surface with tetrahedra using TetGen."
#define	TOOLTIP_ASSIGN_VOLUME_CONSTRAINTS "Assigns volume constraints to selected tetrahedra."
#define	TOOLTIP_CLEAR_VOLUME_CONSTRAINTS "Clears all assigned volume constraints."
#define	TOOLTIP_RETETRAHEDRALIZE "Given a tetrahedralization and volume constraints, this method adapts the tetrahedra using TetGen."
#define	TOOLTIP_DUPLICATE "Duplicates the selected geometry."
#define	TOOLTIP_EXTRUDE_AND_MOVE "Extrudes selected geometry (vertices, edges, faces) and moves new vertices by the specified offset."
#define	TOOLTIP_EXTRUDE_AND_SCALE "Extrudes selected geometry (vertices, edges, faces) and scales new vertices by the specified scale."
#define	TOOLTIP_EXTRUDE_ALONG_NORMAL "Extrudes selected geometry (vertices, edges, faces) and moves new vertices along their normal."
#define TOOLTIP_EXTRUDE_TO_TICKNESS "Extrudes selected geometry (vertices, edges, faces) and moves new vertices so that created elements have the specified thickness."
#define	TOOLTIP_EXTRUDE_CYLINDERS "Extrudes cylinders around selected points of a 2d manifold."
#define TOOLTIP_CREATE_SHRINK_GEOMETRY "Creates new elements from existing ones, providing each with a unique set of corner vertices. Those corners are scaled towards the center using the given scale-parameter."
#define TOOLTIP_EXTRUDE_FACES_WITH_TETS "Experimental function to create 'plaque'-like geometry based on extruding faces with tetrahedra."

#define TOOLTIP_MESH_LAYERS "Creates triangle/quadrilateral grids for the given raster-layers"
#define TOOLTIP_MESH_LAYER_BOUNDARIES "Creates boundary grids for the given raster-layers"
#define TOOLTIP_PROJECT_TO_LAYER "Projects a (surface-)mesh to the specified raster-layer. Only height values in valid regions are adjusted."
#define TOOLTIP_PROJECT_TO_TOP_LAYER "Projects a (surface-)mesh to the top-layer of the specified raster-stack. Only height values in valid regions are adjusted."
#define TOOLTIP_EXTRUDE_LAYERS "Creates volumes for a given stack of raster-layers and an initial triangulation of the surface."
#define TOOLTIP_EXTRUDE_LAYERS_AND_ADD_PROJECTOR "Creates volumes for a given stack of raster-layers and an initial triangulation of the surface. It also generates a raster-based refinement-projector for the whole geometry."
#define TOOLTIP_SNAP_TO_HORIZONTAL_RASTER "Snaps all vertices of the given (surface) grid horizontally to the closest raster-node."

#define TOOLTIP_CSG_FACE_UNION "Performs a union operation on the geometry of the two specifed subsets. IMPORTANT: Both subsets have to be closed manifolds, i.e., homeomorphic to the sphere."
#define TOOLTIP_CSG_FACE_INTERSECTION "Performs an intersection operation on the geometry of the two specifed subsets. IMPORTANT: Both subsets have to be closed manifolds, i.e., homeomorphic to the sphere."
#define TOOLTIP_CSG_FACE_DIFFERENCE "Performs a difference operation on the geometry of the two specifed subsets. IMPORTANT: Both subsets have to be closed manifolds, i.e., homeomorphic to the sphere."

namespace ug{
namespace promesh{

/// \addtogroup promesh
/// \{

void SimplifyPolylines(Mesh* m, number curvatureThreshold);

void SimplifySmoothedPolylines(
			Mesh* m,
			number curvatureThreshold,
			number smoothingAlpha,
			int smoothingIterations);

void ConvertToTriangles(Mesh* obj);

void ConvertToQuadrilaterals(Mesh* obj);

void ExtrudeFacesWithTets(Mesh* obj, int fromSi, int toSi, const number factor);

void TriangleFill(Mesh* obj, bool qualityGeneration, number minAngle, int si);

void Retriangulate(Mesh* obj, number minAngle);

void AdjustEdgeLength(
			Mesh* obj,
			number minEdgeLen,
			number maxEdgeLen,
			int numIterations,
			bool adaptive,
			bool automarkBoundaries);

void AdjustEdgeLengthExtended(
			Mesh* obj,
			number minEdgeLen,
			number maxEdgeLen,
			number approximation,
			number triQuality,
			int numIterations,
			bool automarkBoundaries);

void AdaptSurfaceToCylinder(Mesh* obj, number radius, number threshold);

void ConvertToTetrahedra(Mesh* obj);

void Tetrahedralize(
			Mesh* obj,
			number quality,
			bool preserveOuter,
			bool preserveAll,
			bool separateVolumes,
			bool appendSubsetsAtEnd,
			int verbosity);

void AssignVolumeConstraints(Mesh* obj, number volConstraint);

void ClearVolumeConstraints(Mesh* obj);

void Retetrahedralize(
			Mesh* obj,
			number quality,
			bool preserveOuter,
			bool preserveAll,
			bool applyVolumeConstraint,
			int verbosity);

void Duplicate(Mesh* obj, const vector3& offset, bool deselectOld, bool selectNew);

void ExtrudeAndMove(
			Mesh* obj,
			const vector3& totalDir,
			int numSteps,
			bool createFaces,
			bool createVolumes);

void ExtrudeAndScale(
			Mesh* obj,
			number totalScale,
			bool scaleAroundPivot,
			int numSteps,
			bool createFaces,
			bool createVolumes);

void ExtrudeAlongNormal(
			Mesh* obj,
			number totalLength,
			int numSteps,
			bool createFaces,
			bool createVolumes);

void ExtrudeToThickness(Mesh* obj,
                        number thickness,
                        int numSteps,
						bool createFaces,
						bool createVolumes);

void ExtrudeCylinders(Mesh* obj, number height, number radius, number snapThreshold);

/**	For each element of type TElem in obj this method creates a new element with
 * separate corners. The new element will be scaled by 'scale'.
 * All original elements will be deleted before the method terminates.
 * It should thus be called for 'Volumes' first, then for 'Faces' and
 * finally for 'Edges'.*/
template <class TElemIter>
void CreateShrinkElements(
			Mesh* obj,
			number scale,
			TElemIter elemsBegin,
			TElemIter elemsEnd)
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
void CreateShrinkGeometry(Mesh* obj, number scale);

void ReplaceLowValenceVertices(Mesh* obj, number maxSquaredHeightToBaseAreaRatio);

void ReplaceValence3Vertices(Mesh* obj, number maxSquaredHeightToBaseAreaRatio);

void MeshLayers(Mesh* m, const RasterLayers& layers);

void MeshLayerBoundaries(Mesh* m, const RasterLayers& layers);

void ProjectToLayer(Mesh* obj, RasterLayers& layers, int layerIndex);

void ProjectToTopLayer(Mesh* obj, RasterLayers& layers);

void ExtrudeLayers(Mesh* obj, RasterLayers& layers, bool allowForTetsAndPyras);

void ExtrudeLayersAndAddProjector(
			Mesh* obj,
			SPRasterLayers layers,
			bool allowForTetsAndPyras);

void SnapToHorizontalRaster(Mesh* obj, SPRasterLayers layers);

enum CSGOperation {
	CSG_UNION,
	CSG_INTERSECTION,
	CSG_DIFFERENCE
};

void CSGFaceOperation(
			Mesh* obj,
			CSGOperation op,
			int subsetIndex0,
			int subsetIndex1,
			number snapThreshold);


void CSGFaceUnion(
			Mesh* obj,
			int subsetIndex0,
			int subsetIndex1,
			number snapThreshold);

void CSGFaceIntersection(
			Mesh* obj,
			int subsetIndex0,
			int subsetIndex1,
			number snapThreshold);

void CSGFaceDifference(
			Mesh* obj,
			int subsetIndex0,
			int subsetIndex1,
			number snapThreshold);


// }
/// \}

}}// end of namespace

#endif
