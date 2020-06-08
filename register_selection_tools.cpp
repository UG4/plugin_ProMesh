/*
 * Copyright (c) 2015:  G-CSC, Goethe University Frankfurt
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

#include "registration_routines.h"
#include "registry/registry.h"
#include "bridge/util.h"
#include "tools/new_tools.h"
#include "tools/selection_tools.h"

using namespace std;
using namespace ug::bridge;

namespace ug{
namespace promesh{

void RegisterSelectionTools(ProMeshRegistry& reg, string baseGrp)
{
	baseGrp.append("/Selection");
	try{
		string grp = baseGrp;
		reg.add_function("SelectAll", &SelectAll, grp, "",
				"mesh", TOOLTIP_SELECT_ALL) 
			.add_function("ClearSelection", &ClearSelection, grp, "",
				"mesh", TOOLTIP_CLEAR_SELECTION)
			.add_function("CloseSelection", &CloseSelection, grp, "",
				"mesh", TOOLTIP_CLOSE_SELECTION, "", RT_DEFAULT, Key_C)
			.add_function("SelectSelectionBoundary", &SelectSelectionBoundary, grp, "",
				"mesh", TOOLTIP_SELECT_SELECTION_BOUNDARY)
			.add_function("SelectElementsByIndexRange", &SelectElementsByIndexRange, grp, "",
			    "mesh#vertex inds#edge inds#face inds#volume inds#clear selection",
			    TOOLTIP_SELECT_ELEMENTS_BY_INDEX_RANGE);

		grp = baseGrp + "/Vertices";
		reg.add_function("SelectAllVertices", &SelectAllVertices, grp, "",
				"mesh", TOOLTIP_SELECT_ALL_VERTICES)
			.add_function("DeselectAllVertices", &DeselectAllVertices, grp, "",
				"mesh", TOOLTIP_DESELECT_ALL_VERTICES)
			.add_function("SelectBoundaryVertices", &SelectBoundaryVertices, grp, "",
				"mesh", TOOLTIP_SELECT_BOUNDARY_VERTICES)
			.add_function("SelectInnerVertices", &SelectInnerVertices, grp, "",
				"mesh", TOOLTIP_SELECT_INNER_VERTICES)
			.add_function("SelectAssociatedVertices", &SelectAssociatedVertices, grp, "",
				"mesh", TOOLTIP_SELECT_ASSOCIATED_VERTICES)
			.add_function("SelectMarkedVertices", &SelectMarkedVertices, grp, "",
				"mesh", TOOLTIP_SELECT_MARKED_VERTICES)
			.add_function("SelectInterfaceVertices", &SelectInterfaceElements<Vertex>, grp, "",
				"mesh # regard selected neighbors only", TOOLTIP_SELECT_INTERFACE_ELEMENTS)
			.add_function("SelectVertexByIndex", &SelectVertexByIndex, grp, "bool",
				"mesh # index", TOOLTIP_SELECT_VERTEX_BY_INDEX)
			.add_function("SelectUnconnectedVertices", &SelectUnconnectedVertices, grp, "size_t",
				"mesh #"
				"unconnected to edges || value=true #"
				"unconnected to faces || value=true #"
				"unconnected to volumes || value=true",
				TOOLTIP_SELECT_UNCONNECTED_VERTICES);
			
		grp = baseGrp + "/Vertices/Kinks";
		reg.add_function("SelectSelectionKinkVertices", &SelectSelectionKinkVertices, grp,
						  "num selected",
						  "mesh #"
						  "threshold angle || min=0D; value=20D; max=180D #"
						  "select darts || value=true",
						  TOOLTIP_SELECT_SELECTION_KINK_VERTICES)
			.add_function("SelectSubsetKinkVertices", &SelectSubsetKinkVertices, grp,
						  "num selected",
						  "mesh #"
						  "subset index || min=0; value=0 #"
						  "threshold angle || min=0D; value=20D; max=180D #"
						  "select darts || value=true",
						  TOOLTIP_SELECT_SUBSET_KINK_VERTICES);

		grp = baseGrp + "/Vertices/Coordinate Based";
		reg.add_function("SelectVertexByCoordinate", &SelectElemByCoordinate<Vertex>, grp, "",
				"mesh # coordinate", TOOLTIP_SELECT_VERTEX_BY_COORDINATE)
			.add_function("SelectVertexByCylindricalCoordinate", &SelectElemByCylindricalCoordinate<Vertex>, grp, "",
				"mesh # rho # phi # z", TOOLTIP_SELECT_VERTEX_BY_CYL_COORDINATE)
			.add_function("SelectVerticesInBox", &SelectElementsInBox<Vertex>, grp, "",
				"mesh # min # max || value=[1,1,1]", TOOLTIP_SELECT_VERTEX_IN_BOX)
			.add_function("SelectVerticesInCylinder", &SelectElementsInCylinder<Vertex>, grp, "",
				"mesh # base # top || value=[0,0,1] # radius || value=1", TOOLTIP_SELECT_VERTEX_IN_CYLINDER);


		grp = baseGrp + "/Edges";
		reg
			.add_function("SelectAllEdges", &SelectAllEdges, grp, "",
				"mesh", TOOLTIP_SELECT_ALL_EDGES)
			.add_function("DeselectAllEdges", &DeselectAllEdges, grp, "",
				"mesh", TOOLTIP_DESELECT_ALL_EDGES)
			.add_function("SelectAssociatedEdges", &SelectAssociatedEdges, grp, "",
				"mesh", TOOLTIP_SELECT_ASSOCIATED_EDGES)
			.add_function("SelectBoundaryEdges", &SelectBoundaryEdges, grp, "",
				"mesh", TOOLTIP_SELECT_BOUNDARY_EDGES)
			.add_function("SelectInnerEdges", &SelectInnerEdges, grp, "",
				"mesh", TOOLTIP_SELECT_INNER_EDGES)
			.add_function("SelectLinkedEdges", &SelectLinkedElements<Edge>, grp, "",
				"mesh", TOOLTIP_SELECT_LINKED_EDGES)
			.add_function("SelectNonManifoldEdges", &SelectNonManifoldEdges, grp, "",
				"mesh", TOOLTIP_SELECT_NON_MANIFOLD_EDGES)
			.add_function("SelectMarkedEdges", &SelectMarkedEdges, grp, "",
				"mesh", TOOLTIP_SELECT_MARKED_EDGES)
			.add_function("EdgeSelectionFill", &EdgeSelectionFill, grp, "",
				"mesh", TOOLTIP_EDGE_SELECTION_FILL)
			.add_function("SelectShortEdges", &SelectShortEdges, grp, "",
				"mesh # max length || value = 0.1D; min=0D;", TOOLTIP_SELECT_SHORT_EDGES)
			.add_function("SelectLongEdges", &SelectLongEdges, grp, "",
				"mesh # min length || value = 1", TOOLTIP_SELECT_LONG_EDGES)
			.add_function("SelectCreaseEdges", &SelectCreaseEdges, grp, "",
				"mesh # min angle || value = 20", TOOLTIP_SELECT_CREASE_EDGES)
			.add_function("SelectLinkedBoundaryEdges", &SelectLinkedBoundaryEdges, grp, "",
				"mesh # stop at selected vertices", TOOLTIP_SELECT_LINKED_BOUNDARY_EDGES)
			.add_function("SelectSmoothEdgePath", &SelectSmoothEdgePath, grp, "",
				"mesh #"
				"max deviation angle || value=20D #"
				"normal weight || value=0.9D; step=0.05D#"
				"stop at selected vertices",
				TOOLTIP_SELECT_SMOOTH_EDGE_PATH)
			.add_function("SelectInterfaceEdges", &SelectInterfaceElements<Edge>, grp, "",
				"mesh # regard selected neighbors only", TOOLTIP_SELECT_INTERFACE_ELEMENTS)
			.add_function("SelectEdgesByDirection", &SelectEdgesByDirection, grp, "",
				"mesh #"
				"direction#"
				"min deviation angle || value = 0D; min = 0D; max = 180D #"
				"max deviation angle || value = 1D; min = 0D; max = 180D #"
				"select flipped edges || value = true",
				TOOLTIP_SELECT_EDGES_BY_DIRECTION)
			.add_function("SelectSubsetEdgesByDirection", &SelectSubsetEdgesByDirection, grp, "",
				"mesh #"
				"subset index#"
				"direction#"
				"min deviation angle || value = 0D; min = 0D; max = 180D #"
				"max deviation angle || value = 1D; min = 0D; max = 180D #"
				"select flipped edges || value = true",
				TOOLTIP_SELECT_SUBSET_EDGES_BY_DIRECTION)
			.add_function("SelectShortPolychains", &SelectShortPolychains, grp, "",
				  "mesh # maxChainLength || min=0D; value=1D # closedChainsOnly",
				  TOOLTIP_SELECT_SHORT_POLYCHAINS)
			.add_function("SelectEdgeByIndex", &SelectEdgeByIndex, grp, "",
				"mesh # index", TOOLTIP_SELECT_EDGE_BY_INDEX);

		grp = baseGrp + "/Edges/Coordinate Based";
			reg.add_function("SelectEdgeByCoordinate", &SelectElemByCoordinate<Edge>, grp, "",
				"mesh # coordinate", TOOLTIP_SELECT_EDGE_BY_COORDINATE)
			.add_function("SelectEdgeByCylindricalCoordinate", &SelectElemByCylindricalCoordinate<Edge>, grp, "",
				"mesh # rho # phi # z", TOOLTIP_SELECT_EDGE_BY_CYL_COORDINATE)
			.add_function("SelectEdgesInBox", &SelectElementsInBox<Edge>, grp, "",
				"mesh # min # max || value=[1,1,1]", TOOLTIP_SELECT_EDGE_IN_BOX)  
			.add_function("SelectEdgesInCylinder", &SelectElementsInCylinder<Edge>, grp, "",
				"mesh # base # top || value=[0,0,1] # radius || value=1", TOOLTIP_SELECT_EDGE_IN_CYLINDER);
			
		grp = baseGrp + "/Faces";
		reg
			.add_function("SelectAllFaces", &SelectAllFaces, grp, "",
				"mesh", TOOLTIP_SELECT_ALL_FACES)
			.add_function("DeselectAllFaces", &DeselectAllFaces, grp, "",
				"mesh", TOOLTIP_DESELECT_ALL_FACES)
			.add_function("SelectAssociatedFaces", &SelectAssociatedFaces, grp, "",
				"mesh", TOOLTIP_SELECT_ASSOCIATED_FACES)
			.add_function("SelectBoundaryFaces", &SelectBoundaryFaces, grp, "",
				"mesh", TOOLTIP_SELECT_BOUNDARY_FACES)
			.add_function("SelectInnerFaces", &SelectInnerFaces, grp, "",
				"mesh", TOOLTIP_SELECT_INNER_FACES)
			.add_function("SelectMarkedFaces", &SelectMarkedFaces, grp, "",
				"mesh", TOOLTIP_SELECT_MARKED_FACES)
			.add_function("FaceSelectionFill", &FaceSelectionFill, grp, "",
				"mesh", TOOLTIP_FACE_SELECTION_FILL)
			.add_function("SelectLinkedFaces", &SelectLinkedElements<Face>, grp, "",
				"mesh", TOOLTIP_SELECT_LINKED_FACES)
			.add_function("SelectLinkedManifoldFaces", &SelectLinkedManifoldFaces, grp, "",
				"mesh", TOOLTIP_SELECT_LINKED_MANIFOLD_FACES)
			.add_function("SelectLinkedBoundaryFaces", &SelectLinkedBoundaryFaces, grp, "",
				"mesh # stop at selected edges", TOOLTIP_SELECT_LINKED_BOUNDARY_FACES)
			.add_function("SelectLinkedFlatFaces", &SelectLinkedFlatFaces, grp, "",
				"mesh #"
				"max deviation angle || value=1 #"
				"ignore orientation || value=true#"
				"traverse degenerated faces#"
				"stop at selected edges", TOOLTIP_SELECT_LINKED_FLAT_FACES)

			.add_function("SelectInterfaceFaces", &SelectInterfaceElements<Face>, grp, "",
				"mesh # regard selected neighbors only", TOOLTIP_SELECT_INTERFACE_ELEMENTS)
			
			.add_function("SelectFacesByNormal",
			    static_cast<void (*)(Mesh*, const vector3&, number)>(&SelectFacesByNormal),
			    grp, "",
				"mesh #"
				"normal #"
				"max deviation angle || value=10",
				TOOLTIP_SELECT_FACES_BY_NORMAL,
				"", RT_NO_PROMESH)

			.add_function("SelectFacesByNormal",
			    static_cast<void (*)(Mesh*, const vector3&, number, number, bool)>(&SelectFacesByNormal),
			    grp, "",
				"mesh #"
				"normal #"
				"min deviation angle || value=0 #"
				"max deviation angle || value=10 #"
				"ignore inner faces	|| value=true",
				TOOLTIP_SELECT_FACES_BY_NORMAL)

			.add_function("SelectAnisotropicFaces", &SelectAnisotropicElements<Face>, grp, "",
				"mesh # min edge ratio || value=0.5D", TOOLTIP_SELECT_ANISOTROPIC_ELEMENTS)
			.add_function("SelectDegenerateFaces", &SelectDegenerateFaces, grp, "",
				"mesh # max height || value=0.01D", TOOLTIP_SELECT_DEGENERATE_FACES)
			.add_function("SelectBentQuadrilaterals", &SelectBentQuadrilaterals, grp, "",
				"mesh # dot threshold || value = 0.95D", TOOLTIP_SELECT_BENT_QUADRILATERALS)
			.add_function("SelectFaceByIndex", &SelectFaceByIndex, grp, "",
				"mesh # index", TOOLTIP_SELECT_FACE_BY_INDEX);

		grp = baseGrp + "/Faces/Coordinate Based";
			reg.add_function("SelectFaceByCoordinate", &SelectElemByCoordinate<Face>, grp, "",
				"mesh # coordinate", TOOLTIP_SELECT_FACE_BY_COORDINATE)
			.add_function("SelectFaceByCylindricalCoordinate", &SelectElemByCylindricalCoordinate<Face>, grp, "",
				"mesh # rho # phi # z", TOOLTIP_SELECT_FACE_BY_CYL_COORDINATE)
			.add_function("SelectFacesInBox", &SelectElementsInBox<Face>, grp, "",
				"mesh # min # max || value=[1,1,1]", TOOLTIP_SELECT_FACE_IN_BOX)
			.add_function("SelectFacesInCylinder", &SelectElementsInCylinder<Face>, grp, "",
				"mesh # base # top || value=[0,0,1] # radius || value=1", TOOLTIP_SELECT_FACE_IN_CYLINDER);

		grp = baseGrp + "/Faces";
			reg.add_function("SelectIntersectingTriangles", &SelectIntersectingTriangles, grp, "",
				"mesh", TOOLTIP_SELECT_INTERSECTING_TRIANGLES);

		grp = baseGrp + "/Volumes";
		reg.add_function("SelectAllVolumes", &SelectAllVolumes, grp, "",
				"mesh", TOOLTIP_SELECT_ALL_VOLUMES)
			.add_function("DeselectAllVolumes", &DeselectAllVolumes, grp, "",
				"mesh", TOOLTIP_DESELECT_ALL_VOLUMES)
			.add_function("VolumeSelectionFill", &VolumeSelectionFill, grp, "",
				"mesh", TOOLTIP_VOLUME_SELECTION_FILL)
			.add_function("SelectLinkedVolumes", &SelectLinkedElements<Volume>, grp, "",
				"mesh", TOOLTIP_SELECT_LINKED_VOLUMES)
			.add_function("SelectSlivers", &SelectSlivers, grp, "",
				"mesh # threshold ratio || min=0; value=0.01D", TOOLTIP_SELECT_SLIVERS)
			.add_function("SelectAnisotropicVolumes", &SelectAnisotropicElements<Volume>, grp, "",
				"mesh # min edge ratio || value=0.5D", TOOLTIP_SELECT_ANISOTROPIC_ELEMENTS)
			.add_function("SelectUnorientableVolumes", &SelectUnorientableVolumes, grp, "",
				"mesh", TOOLTIP_SELECT_UNORIENTABLE_VOLUMES)
			.add_function("SelectVolumeByIndex", &SelectVolumeByIndex, grp, "",
				"mesh # index", TOOLTIP_SELECT_VOLUME_BY_INDEX)
			.add_function("SelectVolumesByType", &SelectVolumesByType, grp, "",
				"mesh #"
				"Hexahedra || value=true #"
				"Octahedra || value=true #"
				"Prisms || value=true #"
				"Pyramids || value=true #"
				"Tetrahedra || value=true #",
				TOOLTIP_SELECT_VOLUMES_BY_TYPE);

		grp = baseGrp + "/Volumes/Coordinate Based";
		reg.add_function("SelectVolumeByCoordinate", &SelectElemByCoordinate<Volume>, grp, "",
				"mesh # coordinate", TOOLTIP_SELECT_VOLUME_BY_COORDINATE)
			.add_function("SelectVolumeByCylindricalCoordinate", &SelectElemByCylindricalCoordinate<Volume>, grp, "",
				"mesh # rho # phi # z", TOOLTIP_SELECT_VOLUME_BY_CYL_COORDINATE)
			.add_function("SelectVolumesInBox", &SelectElementsInBox<Volume>, grp, "",
				"mesh # min # max || value=[1,1,1]", TOOLTIP_SELECT_VOLUME_IN_BOX)
			.add_function("SelectVolumesInCylinder", &SelectElementsInCylinder<Volume>, grp, "",
				"mesh # base # top || value=[0,0,1] # radius || value=1", TOOLTIP_SELECT_VOLUME_IN_CYLINDER);


		grp = baseGrp;
		reg.add_function("InvertSelection", &InvertSelection, grp, "",
				"mesh #"
				"vertices || value=true #"
				"edges || value=true #"
				"faces || value=true #"
				"volumes || value=true", TOOLTIP_INVERT_SELECTION)
			.add_function("ExtendSelection", &ExtendSelection, grp, "",
				"mesh # neighborhood size || value=1; min=0", TOOLTIP_EXTEND_SELECTION)
			.add_function("ExtendSelectionInDirection",
			    &ExtendSelectionInDirection, grp, "",
				"mesh #"
				"neighborhood size || value=1; min=0 #"
				"direction #"
				"minAngle || value=0; min=0 #"
				"maxAngle || value=1; min=0",
				TOOLTIP_EXTEND_SELECTION_IN_DIRECTION);


		grp = baseGrp + "/Subset Based";
		reg.add_function("SelectSubset", &SelectSubset, grp, "",
				"mesh #"
				"subset index #"
				"vertices || value=true #"
				"edges || value=true #"
				"faces || value=true #"
				"volumes || value=true",
				TOOLTIP_SELECT_SUBSET)
			.add_function("SelectSubsetBoundary", &SelectSubsetBoundary, grp, "",
				"mesh #"
				"subset index #"
				"boundaries of edges || value=true #"
				"boundaries of faces || value=true #"
				"boundaries of volumes || value=true",
				TOOLTIP_SELECT_SUBSET_BOUNDARY)
			.add_function("SelectUnassignedElements", &SelectUnassignedElements, grp, "",
				"mesh #"
				"vertices || value=true #"
				"edges || value=true #"
				"faces || value=true #"
				"volumes || value=true", TOOLTIP_SELECT_UNASSIGNED_ELEMENTS)
			.add_function("RestrictSelectionToSubset", &RestrictSelectionToSubset, grp, "",
				"mesh #"
				"subset index", TOOLTIP_RESTRICT_SELECTION_TO_SUBSET);

		
		grp = baseGrp + "/Coordinate Range";
		reg.add_function("SelectElementsBySplitPlane", &SelectElementsBySplitPlane, grp, "",
				"mesh #"
				"select vrts || value=true #"
				"select edges || value=true #"
				"select faces || value=true #"
				"select volumes || value=true #"
				"pivot #"
				"normal",
				TOOLTIP_SELECT_ELEMENTS_BY_SPLIT_PLANE)
			.add_function("SelectElementsByRangeX", &SelectElementsInCoordinateRange<0, 1>,
				grp, "",
				"mesh #"
				"min #"
				"max #"
				"select vrts || value=true #"
				"select edges || value=true #"
				"select faces || value=true #"
				"select volumes || value=true",
				TOOLTIP_SELECT_ELEMENTS_IN_COORDINATE_RANGE)
			.add_function("SelectElementsByRangeY", &SelectElementsInCoordinateRange<1, 1>,
				grp, "",
				"mesh #"
				"min #"
				"max #"
				"select vrts || value=true #"
				"select edges || value=true #"
				"select faces || value=true #"
				"select volumes || value=true",
				TOOLTIP_SELECT_ELEMENTS_IN_COORDINATE_RANGE)
			.add_function("SelectElementsByRangeZ", &SelectElementsInCoordinateRange<2, 1>,
				grp, "",
				"mesh #"
				"min #"
				"max #"
				"select vrts || value=true #"
				"select edges || value=true #"
				"select faces || value=true #"
				"select volumes || value=true",
				TOOLTIP_SELECT_ELEMENTS_IN_COORDINATE_RANGE)
			.add_function("DeselectElementsByRangeX", &SelectElementsInCoordinateRange<0, 0>,
				grp, "",
				"mesh #"
				"min #"
				"max #"
				"deselect vrts || value=true #"
				"deselect edges || value=true #"
				"deselect faces || value=true #"
				"deselect volumes || value=true",
				TOOLTIP_DESELECT_ELEMENTS_IN_COORDINATE_RANGE)
			.add_function("DeselectElementsByRangeY", &SelectElementsInCoordinateRange<1, 0>,
				grp, "",
				"mesh #"
				"min #"
				"max #"
				"deselect vrts || value=true #"
				"deselect edges || value=true #"
				"deselect faces || value=true #"
				"deselect volumes || value=true",
				TOOLTIP_DESELECT_ELEMENTS_IN_COORDINATE_RANGE)
			.add_function("DeselectElementsByRangeZ", &SelectElementsInCoordinateRange<2, 0>,
				grp, "",
				"mesh #"
				"min #"
				"max #"
				"deselect vrts || value=true #"
				"deselect edges || value=true #"
				"deselect faces || value=true #"
				"deselect volumes || value=true",
				TOOLTIP_DESELECT_ELEMENTS_IN_COORDINATE_RANGE);

		grp = baseGrp + "/Marks";
		reg.add_function("ClearMarks", &ClearMarks, grp, "",
				"mesh", TOOLTIP_CLEAR_MARKS)
			.add_function("MarkSelection", &MarkSelection, grp, "",
				"mesh", TOOLTIP_MARK_SELECTION, "", RT_DEFAULT, Key_M)
			.add_function("UnmarkSelection", &UnmarkSelection, grp, "",
				"mesh", TOOLTIP_UNMARK_SELECTION)
			.add_function("MarkCreaseEdges", &MarkCreaseEdges, grp, "",
				"mesh # min angle || value=25D; min=0; max=180; step=1# clear marks",
				TOOLTIP_MARK_CREASE_EDGES)
			.add_function("MarkCornersOfMarkedEdges", &MarkCornersOfMarkedEdges, grp, "",
				"mesh # min angle | default | value=75D",
				 TOOLTIP_MARK_CORNERS_OF_MARKED_EDGES);
	}
	UG_REGISTRY_CATCH_THROW(baseGrp);
}

} // end namespace promesh
}// end of namespace

