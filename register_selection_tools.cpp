// created by Sebastian Reiter
// s.b.reiter@gmail.com

#include "registration_routines.h"
#include "registry/registry.h"
#include "bridge/util.h"
#include "tooltips.h"
#include "tools/selection_tools.h"

using namespace std;
using namespace ug::bridge;

namespace ug{
namespace promesh{

void RegisterSelectionTools(Registry& reg, string baseGrp)
{
	baseGrp.append("/Selection");
	try{
		string grp = baseGrp;
		reg.add_function("ClearSelection", &ClearSelection, grp, "", "", TOOLTIP_CLEAR_SELECTION)
			.add_function("SelectAll", &SelectAll, grp, "", "", TOOLTIP_SELECT_ALL) 
			.add_function("ExtendSelection", &ExtendSelection, grp, "", "", TOOLTIP_EXTEND_SELECTION)
			.add_function("SelectSubset", &SelectSubset, grp, "", "", TOOLTIP_SELECT_SUBSET)
			.add_function("SelectSubsetBoundary", &SelectSubsetBoundary, grp, "", "", TOOLTIP_SELECT_SUBSET_BOUNDARY)
			.add_function("SelectUnassignedElements", &SelectUnassignedElements, grp, "", "", TOOLTIP_SELECT_UNASSIGNED_ELEMENTS)
			.add_function("InvertSelection", &InvertSelection, grp, "", "", TOOLTIP_INVERT_SELECTION)
			.add_function("SelectSelectionBoundary", &SelectSelectionBoundary, grp, "", "", TOOLTIP_SELECT_SELECTION_BOUNDARY)
			.add_function("CloseSelection", &CloseSelection, grp, "", "", TOOLTIP_CLOSE_SELECTION);

		grp = baseGrp + "/Vertices";
		reg.add_function("SelectVertexByCoordinate", &SelectElemByCoordinate<Vertex>, grp,TOOLTIP_SELECT_VERTEX_BY_COORDINATE)
			.add_function("SelectVertexByCylindricalCoordinate", &SelectElemByCylindricalCoordinate<Vertex>, grp,TOOLTIP_SELECT_VERTEX_BY_CYL_COORDINATE)
			.add_function("SelectBoundaryVertices", &SelectBoundaryVertices, grp, "", "", TOOLTIP_SELECT_BOUNDARY_VERTICES)
			.add_function("SelectInnerVertices", &SelectInnerVertices, grp,TOOLTIP_SELECT_INNER_VERTICES)
			.add_function("SelectAssociatedVertices", &SelectAssociatedVertices, grp, "", "", TOOLTIP_SELECT_ASSOCIATED_VERTICES)
			.add_function("SelectAllVertices", &SelectAllVertices, grp, "", "", TOOLTIP_SELECT_ALL_VERTICES)
			.add_function("DeselectAllVertices", &DeselectAllVertices, grp, "", "", TOOLTIP_DESELECT_ALL_VERTICES)
			.add_function("SelectMarkedVertices", &SelectMarkedVertices, grp, "", "", TOOLTIP_SELECT_MARKED_VERTICES)
			.add_function("SelectVertexByIndex", &SelectVertexByIndex, grp, "", "", TOOLTIP_SELECT_VERTEX_BY_INDEX)
			.add_function("SelectUnconnectedVertices", &SelectUnconnectedVertices, grp, "", "", TOOLTIP_SELECT_UNCONNECTED_VERTICES)
			.add_function("SelectSelectionKinkVertices", &SelectSelectionKinkVertices, grp,
						  "num selected",
						  "mesh #"
						  "threshold angle | default | min=0D; value=20D; max=180D #"
						  "select darts | default | value=true",
						  TOOLTIP_SELECT_SELECTION_KINK_VERTICES)
			.add_function("SelectSubsetKinkVertices", &SelectSubsetKinkVertices, grp,
						  "num selected",
						  "mesh #"
						  "subset index | default | min=0; value=0 #"
						  "threshold angle | default | min=0D; value=20D; max=180D #"
						  "select darts | default | value=true",
						  TOOLTIP_SELECT_SUBSET_KINK_VERTICES);

		grp = baseGrp + "/Edges";
		reg.add_function("SelectEdgeByCoordinate", &SelectElemByCoordinate<Edge>, grp, "", "", TOOLTIP_SELECT_EDGE_BY_COORDINATE)
			.add_function("SelectEdgeByCylindricalCoordinate", &SelectElemByCylindricalCoordinate<Edge>, grp, "", "", TOOLTIP_SELECT_EDGE_BY_CYL_COORDINATE)
			.add_function("SelectBoundaryEdges", &SelectBoundaryEdges, grp, "", "", TOOLTIP_SELECT_BOUNDARY_EDGES)
			.add_function("SelectInnerEdges", &SelectInnerEdges, grp, "", "", TOOLTIP_SELECT_INNER_EDGES)
			.add_function("SelectNonManifoldEdges", &SelectNonManifoldEdges, grp, "", "", TOOLTIP_SELECT_NON_MANIFOLD_EDGES)
			.add_function("SelectSmoothEdgePath", &SelectSmoothEdgePath, grp, "", "", TOOLTIP_SELECT_SMOOTH_EDGE_PATH)
			.add_function("SelectShortEdges", &SelectShortEdges, grp, "", "", TOOLTIP_SELECT_SHORT_EDGES)
			.add_function("SelectLongEdges", &SelectLongEdges, grp, "", "", TOOLTIP_SELECT_LONG_EDGES)
			.add_function("SelectCreaseEdges", &SelectCreaseEdges, grp, "", "", TOOLTIP_SELECT_CREASE_EDGES)
			.add_function("SelectLinkedEdges", &SelectLinkedElements<Edge>, grp, "", "", TOOLTIP_SELECT_LINKED_EDGES)
			.add_function("SelectLinkedBoundaryEdges", &SelectLinkedBoundaryEdges, grp, "", "", TOOLTIP_SELECT_LINKED_BOUNDARY_EDGES)
			.add_function("SelectAssociatedEdges", &SelectAssociatedEdges, grp, "", "", TOOLTIP_SELECT_ASSOCIATED_EDGES)
			.add_function("SelectAllEdges", &SelectAllEdges, grp, "", "", TOOLTIP_SELECT_ALL_EDGES)
			.add_function("DeselectAllEdges", &DeselectAllEdges, grp, "", "", TOOLTIP_DESELECT_ALL_EDGES)
			.add_function("SelectMarkedEdges", &SelectMarkedEdges, grp, "", "", TOOLTIP_SELECT_MARKED_EDGES)
			.add_function("SelectEdgeByIndex", &SelectEdgeByIndex, grp, "", "", TOOLTIP_SELECT_EDGE_BY_INDEX)
			.add_function("EdgeSelectionFill", &EdgeSelectionFill, grp, "", "", TOOLTIP_EDGE_SELECTION_FILL)
			.add_function("SelectShortPolychains", &SelectShortPolychains, grp, "",
						  "mesh # maxChainLength | default | min=0D; value=1D # closedChainsOnly",
						  TOOLTIP_SELECT_SHORT_POLYCHAINS);
			
		grp = baseGrp + "/Faces";
		reg.add_function("SelectFaceByCoordinate", &SelectElemByCoordinate<Face>, grp, "", "", TOOLTIP_SELECT_FACE_BY_COORDINATE)
			.add_function("SelectFaceByCylindricalCoordinate", &SelectElemByCylindricalCoordinate<Face>, grp, "", "", TOOLTIP_SELECT_FACE_BY_CYL_COORDINATE)
			.add_function("SelectBoundaryFaces", &SelectBoundaryFaces, grp, "", "", TOOLTIP_SELECT_BOUNDARY_FACES)
			.add_function("SelectInnerFaces", &SelectInnerFaces, grp, "", "", TOOLTIP_SELECT_INNER_FACES)
			.add_function("SelectLinkedFaces", &SelectLinkedElements<Face>, grp, "", "", TOOLTIP_SELECT_LINKED_FACES)
			.add_function("SelectLinkedManifoldFaces", &SelectLinkedManifoldFaces, grp, "", "", TOOLTIP_SELECT_LINKED_MANIFOLD_FACES)
			.add_function("SelectLinkedBoundaryFaces", &SelectLinkedBoundaryFaces, grp, "", "", TOOLTIP_SELECT_LINKED_BOUNDARY_FACES)
			.add_function("SelectDegenerateFaces", &SelectDegenerateFaces, grp, "", "", TOOLTIP_SELECT_DEGENERATE_FACES)
			.add_function("SelectLinkedFlatFaces", &SelectLinkedFlatFaces, grp, "", "", TOOLTIP_SELECT_LINKED_FLAT_FACES)
			.add_function("SelectIntersectingTriangles", &SelectIntersectingTriangles, grp, "", "", TOOLTIP_SELECT_INTERSECTING_TRIANGLES)
			.add_function("SelectAssociatedFaces", &SelectAssociatedFaces, grp, "", "", TOOLTIP_SELECT_ASSOCIATED_FACES)
			.add_function("SelectAllFaces", &SelectAllFaces, grp, "", "", TOOLTIP_SELECT_ALL_FACES)
			.add_function("DeselectAllFaces", &DeselectAllFaces, grp, "", "", TOOLTIP_DESELECT_ALL_FACES)
			.add_function("SelectFaceByIndex", &SelectFaceByIndex, grp, "", "", TOOLTIP_SELECT_FACE_BY_INDEX)
			.add_function("SelectFacesByNormal", &SelectFacesByNormal, grp, "", "", TOOLTIP_SELECT_FACES_BY_NORMAL)
			.add_function("FaceSelectionFill", &FaceSelectionFill, grp, "", "", TOOLTIP_FACE_SELECTION_FILL)
			.add_function("SelectBentQuadrilaterals", &SelectBentQuadrilaterals, grp, "", "", TOOLTIP_SELECT_BENT_QUADRILATERALS);

		grp = baseGrp + string("/Volumes");
		reg.add_function("SelectVolumeByCoordinate", &SelectElemByCoordinate<Volume>, grp, "", "", TOOLTIP_SELECT_VOLUME_BY_COORDINATE)
			.add_function("SelectVolumeByCylindricalCoordinate", &SelectElemByCylindricalCoordinate<Volume>, grp, "", "", TOOLTIP_SELECT_VOLUME_BY_CYL_COORDINATE)
			.add_function("SelectAllVolumes", &SelectAllVolumes, grp, "", "", TOOLTIP_SELECT_ALL_VOLUMES)
			.add_function("DeselectAllVolumes", &DeselectAllVolumes, grp, "", "", TOOLTIP_DESELECT_ALL_VOLUMES)
			.add_function("SelectLinkedVolumes", &SelectLinkedElements<Volume>, grp, "", "", TOOLTIP_SELECT_LINKED_VOLUMES)
			.add_function("SelectUnorientableVolumes", &SelectUnorientableVolumes, grp, "", "", TOOLTIP_SELECT_UNORIENTABLE_VOLUMES)
			.add_function("SelectSlivers", &SelectSlivers, grp, "", "mesh # threshold ratio | default | min=0; val=0.01D", TOOLTIP_SELECT_SLIVERS)
			.add_function("SelectVolumeByIndex", &SelectVolumeByIndex, grp, "", "", TOOLTIP_SELECT_VOLUME_BY_INDEX)
			.add_function("VolumeSelectionFill", &VolumeSelectionFill, grp, "", "", TOOLTIP_VOLUME_SELECTION_FILL);


		grp = baseGrp + string("/Marks");
		reg.add_function("ClearMarks", &ClearMarks, grp, "", "", TOOLTIP_CLEAR_MARKS)
			.add_function("MarkSelection", &MarkSelection, grp, "", "", TOOLTIP_MARK_SELECTION)
			.add_function("UnmarkSelection", &UnmarkSelection, grp, "", "", TOOLTIP_UNMARK_SELECTION)
			.add_function("MarkCornersOfMarkedEdges", &MarkCornersOfMarkedEdges,
						  grp, "", "mesh # min angle | default | value=75D",
						  TOOLTIP_MARK_CORNERS_OF_MARKED_EDGES);
	}
	UG_REGISTRY_CATCH_THROW(baseGrp);
}

} // end namespace promesh
}// end of namespace

