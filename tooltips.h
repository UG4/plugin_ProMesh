
#ifndef __H__PROMESH_TOOLTIPS
#define __H__PROMESH_TOOLTIPS

//	camera tools
#define	TOOLTIP_CENTER_OBJECT "Centers the current object."
#define	TOOLTIP_SLIDER_TEST "Tests the DoubleSlider."
#define	TOOLTIP_CENTER_SELECTION "Centers the current selection."
#define TOOLTIP_TOP_VIEW "View the current scene from the top."
#define	TOOLTIP_FRONT_VIEW "View the current scene from the front."
#define	TOOLTIP_SIDE_VIEW "View the current scene from the side."
#define	TOOLTIP_HIDE_SELECTED_ELEMENTS "Hides all currently selected elements."
#define	TOOLTIP_UNHIDE_ELEMENTS "Unhides all hidden elements."

// grid generation tools
#define	TOOLTIP_NEW_OBJECT "Creates a new empty object."
#define	TOOLTIP_MERGE_OBJECTS "Merges the selected objects into a new one."
#define TOOLTIP_CLONE_MESH "Creates a new mesh and copies all content from the given mesh into the new instance."
#define	TOOLTIP_COPY_SELECTION "Copies the selected elements to a new mesh."
#define	TOOLTIP_CREATE_VERTEX "Creates a new vertex"
#define	TOOLTIP_CREATE_EDGE "Creates an edge between two selected vertices."
#define	TOOLTIP_CREATE_FACE "Creates a face between selected vertices."
#define	TOOLTIP_CREATE_VOLUME "Creates a volume between selected vertices."
#define	TOOLTIP_CREATE_PLANE "Creates a plane."
#define	TOOLTIP_CREATE_CIRCLE "Creates a circle."
#define	TOOLTIP_CREATE_BOX "Creates a box."
#define	TOOLTIP_CREATE_SPHERE "Creates a sphere."
#define	TOOLTIP_CREATE_TETRAHEDRON "Creates a tetrahedron."
#define	TOOLTIP_CREATE_PYRAMID "Creates a pyramid."
#define	TOOLTIP_CREATE_PRISM "Creates a prism."
#define	TOOLTIP_CREATE_DUALGRID "creates the dual grid consisting of control volumes as used in the finite volume method" 

#define TOOLTIP_MESH_LAYERS "Creates triangle/quadrilateral grids for the given raster-layers"
#define TOOLTIP_MESH_LAYER_BOUNDARIES "Creates boundary grids for the given raster-layers"
#define TOOLTIP_EXTRUDE_LAYERS "Creates volumes for a given stack of raster-layers and an initial triangulation of the surface."

//remeshing tools
#define TOOLTIP_SIMPLIFY_POLYLINES "Removes vertices from the selected polyline which have a smaller curvature than the specified angle."
#define TOOLTIP_SIMPLIFY_SMOOTHED_POLYLINES "Removes vertices from the selected polyline which have a smaller smoothed curvature than the specified angle."
#define	TOOLTIP_CONVERT_TO_TRIANGLES "Converts selected quadrilaterals to triangles."
#define	TOOLTIP_TRIANGLE_FILL "Performs triangle fill using the sweep-line algorithm."
#define	TOOLTIP_RETRIANGULATE "Inserts vertices as required and performs Constrained Delaunay triangulation."
#define	TOOLTIP_ADJUST_EDGE_LENGTH "Remeshes the active grid so that all edges approximatly have a certain length."
#define	TOOLTIP_ADAPT_SURFACE_TO_CYLINDER "Introduces edges in a grid around a selected vertex which roughly correspond to the intersection of a cylinder with the surface."
#define TOOLTIP_REPLACE_LOW_VALENCE_VERTICES "Replaces selected valence-3 and valence-4 vertices by triangles, if the curvature of associated triangles is low"
#define	TOOLTIP_TETRAHEDRALIZE "Fills a closed surface with tetrahedrons."
#define	TOOLTIP_ASSIGN_VOLUME_CONSTRAINTS "Assigns volume constraints to selected tetrahedrons."
#define	TOOLTIP_CLEAR_VOLUME_CONSTRAINTS "Clears all assigned volume constraints."
#define	TOOLTIP_RETETRAHEDRALIZE "Given a tetrahedralization and volume constraints, this method adapts the tetrahedrons."
#define	TOOLTIP_DUPLICATE "Duplicates the selected geometry."
#define	TOOLTIP_EXTRUDE "Extrudes selected geometry (vertices, edges, faces)."
#define	TOOLTIP_EXTRUDE_CYLINDERS "Extrudes cylinders around selected points of a 2d manifold."
#define TOOLTIP_CREATE_SHRINK_GEOMETRY "Creates new elements from existing ones, providing each with a unique set of corner vertices. Those corners are scaled towards the center using the given scale-parameter."
#define TOOLTIP_EXTRUDE_FACES_WITH_TETS "Experimental function to create 'plaque'-like geometry based on extruding faces with tetrahedrons."

//refinement tools
#define	TOOLTIP_REFINE "Refines selected elements and builds a regular closure."
#define	TOOLTIP_HANGING_NODE_REFINE "Refines selected elements using hanging nodes"
#define	TOOLTIP_REFINE_SMOOTH "Refines selected elements using piecewise smooth refinement."
#define	TOOLTIP_REFINE_SMOOTH_BOUNDARY_2D "Refines selected elements using smooth subdivision rules on the boundary edges."
#define	TOOLTIP_FRACTURED_MEDIA_REFINE "Refines selected elements using hanging nodes. Fractures are refined anisotropic."
#define	TOOLTIP_CREATE_FRACTAL "Refines the whole geometry using a fractal-refinement scheme-"
#define	TOOLTIP_INSERT_CENTER "Inserts a central vertex in all selected elements."


//topology tools
#define TOOLTIP_RESOLVE_SELF_INTERSECTIONS "Resolves self intersections of faces and edges."
#define	TOOLTIP_RESOLVE_EDGE_INTERSECTIONS "Makes sure that all edge intersections are represented by a vertex."
#define	TOOLTIP_RESOLVE_TRIANGLE_INTERSECTIONS "Makes sure that all triangle intersections are represented by an edge and vertices."
#define	TOOLTIP_PROJECT_VERTICES_TO_CLOSE_EDGES "Projects selected vertices to selected close edges."
#define	TOOLTIP_PROJECT_VERTICES_TO_CLOSE_FACES "Projects selected vertices to selected close faces."
#define	TOOLTIP_INTERSECT_CLOSE_EDGES "Performs intersections between selected close edges."
#define	TOOLTIP_ERASE_SELECTED_ELEMENTS "Erases selected elements and associated unreferenced geometry."
#define	TOOLTIP_REMOVE_DOUBLES "Removes selected vertices that are close to each other"
#define	TOOLTIP_REMOVE_DOUBLE_EDGES "Removes selected duplicates of selected edges."
#define	TOOLTIP_REMOVE_DOUBLE_FACES "Removes selected duplicates of selected faces."
#define	TOOLTIP_MERGE_AT_FIRST "Merges all selected objects into a single vertex at the first selected vertex."
#define	TOOLTIP_MERGE_AT_CENTER "Merges all selected objects into a single vertex at the center of the selection."
#define	TOOLTIP_MERGE_AT_LAST "Merges all selected objects into a single vertex at the last selected vertex."
#define	TOOLTIP_COLLAPSE_EDGE "Collapses the edge and removes adjacent triangles."
#define	TOOLTIP_SPLIT_EDGE "Splits the edge and inserts new triangles."
#define	TOOLTIP_SWAP_EDGE "Swaps selected edges that are adjacent to exactly two triangles."
#define	TOOLTIP_PLANE_CUT "Cuts selected edges along the given plane."
#define	TOOLTIP_ADJUST_EDGE_ORIENTATION "Adjusts the orientation of boundary edges to associated faces."
#define	TOOLTIP_FIX_FACE_ORIENTATION "Tries to change orientation of selected faces so that all neighbouring faces point into the same direction. Only works correctly for manifold selections."
#define	TOOLTIP_FIX_FACE_SUBSET_ORIENTATIONS "Iterates over all subset and tries to fix face orientation for each. Only works correctly for manifold subsets."
#define	TOOLTIP_FIX_VOLUME_ORIENTATION "Changes orientation of selected volumes, so that they are oriented correctly."
#define	TOOLTIP_INVERT_FACE_ORIENTATION "Inverts the orientation of all selected faces."

//info tools
#define	TOOLTIP_PRINT_SELECTION_CENTER "Calculates and prints the position of the center of the current selection."
#define	TOOLTIP_PRINT_GEOMETRY_INFO "Prints info about the current geometry"
#define	TOOLTIP_PRINT_FACE_QUALITY "Prints the quality of selected faces"
#define	TOOLTIP_PRINT_SELECTION_INFO "Prints the quantities of selected elements"
#define	TOOLTIP_PRINT_SELECTION_CONTAINING_SUBSETS "Prints subset indices of all subsets, which contain a selected element."
#define	TOOLTIP_PRINT_VERTEX_DISTANCE "Prints the min and max distance of vertices of selected elements."
#define	TOOLTIP_PRINT_LEAST_SQUARES_PLANE "Prints the position and normal of the least squares fitting plane"

//subset tools
#define TOOLTIP_SET_SUBSET_NAME ""
#define	TOOLTIP_ASSIGN_SUBSET "Assigns the selected elements to a subset."
#define	TOOLTIP_ASSIGN_SUBSET_COLORS "assigns subset colors by a procedural scheme."
#define	TOOLTIP_SEPARATE_FACES_BY_EDGE_SUBSETS "Assigns faces that are surrounded by a set of edge-subsets to a common subset."
#define	TOOLTIP_SEPARATE_FACES_BY_SELECTED_EDGES "Assigns faces that are surrounded by a set of selected edges to a common subset."
#define	TOOLTIP_SEPARATE_VOLUMES_BY_FACE_SUBSETS "Assigns volumes that are surrounded by a set of face-subsets to a common subset."
#define	TOOLTIP_SEPARATE_VOLUMES_BY_SELECTED_FACES "Assigns volumes that are surrounded by a set of selected faces to a common subset."
#define	TOOLTIP_SEPARATE_IRREGULAR_MANIFOLD_SUBSETS "After this algorithm all face-subsets are regular manifolds."
#define	TOOLTIP_MOVE_SUBSET "Moves a subset to another index."
#define	TOOLTIP_SWAP_SUBSETS "Swaps two subsets"
#define	TOOLTIP_JOIN_SUBSETS "Joins two subsets"
#define	TOOLTIP_ERASE_SUBSET "Erases a subset, but not its associated geometry."
#define	TOOLTIP_ERASE_EMPTY_SUBSETS "Erases Subsets, which do not contain any elements at all."
#define	TOOLTIP_ADJUST_SUBSETS_FOR_UG3 "Assigns face and edge indices so that the geometry can be used with ug3."
#define	TOOLTIP_ADJUST_SUBSETS_FOR_UG4 "Adjusts subsets for simulation with ug4."
#define	TOOLTIP_SEPARATE_FACE_SUBSETS_BY_NORMAL "Collects faces of each subset that have a similar normal and assigns them to new subsets."
#define	TOOLTIP_SEPARATE_FACE_SUBSET_BY_NORMAL "Collects faces of a given subset that have a similar normal and assigns them to new subsets."
#define	TOOLTIP_ASSIGN_SUBSETS_BY_QUALITY "Assigns the selected to a subset depending on their quality."
#define	TOOLTIP_SEPARATE_DEGENERATED_BOUNDARY_FACE_SUBSETS "Separates degenerated boundary face subsets at sharp creases."
#define	TOOLTIP_COPY_SUBSET_INDICES_TO_SIDES "Copies subset indices of selected elements to sides of those elements."
#define	TOOLTIP_ASSIGN_SUBSETS_BY_ELEMENT_TYPE "Assigns elemets to subsets based on their concrete type."

//coordinate transform tools
#define	TOOLTIP_GET_SELECTION_CENTER ""
#define TOOLTIP_SET_SELECTION_CENTER ""
#define TOOLTIP_MOVE_ALONG_NORMAL ""
#define TOOLTIP_SCALE_AROUND_CENTER ""
#define TOOLTIP_SCALE_AROUND_PIVOT ""
#define TOOLTIP_ROTATE_AROUND_CENTER ""
#define TOOLTIP_ROTATE_AROUND_PIVOT ""
#define	TOOLTIP_COORDINATES "Coordinates of the center of the current selection"
#define	TOOLTIP_MOVE "Moves selected vertices."
#define TOOLTIP_MOVE_MESH_TO "Moves the active mesh and its pivot, so that the pivot will be located on the specified position."
#define	TOOLTIP_NORMAL_MOVE "Moves selected vertices along their normal."
#define	TOOLTIP_SCALE "Scales the coordinates of the selected vertices around their center."
#define	TOOLTIP_ROTATE "Rotates the geometry by the given degrees around its center."
#define	TOOLTIP_TRANSFORM "Transforms the vertices with the given matrix"
#define	TOOLTIP_CONE_TRANSFORM "Transforms the vertices with the given cone transformation"
#define	TOOLTIP_LAPLACIAN_SMOOTH "Smoothes vertices in a grid."
#define	TOOLTIP_WEIGHTED_EDGE_SMOOTH "Smoothes vertices along edges in a grid with special weights for non-smoothed vertices."
#define	TOOLTIP_WEIGHTED_FACE_SMOOTH "Smoothes vertices across faces in a grid with special weights for non-smoothed vertices."
#define	TOOLTIP_WEIGHTED_NORMAL_SMOOTH "The higher the dot-product between an outgoing edge and the vertex normal, the higher the influence of that edge during smoothing of that vertex."
#define	TOOLTIP_SLOPE_SMOOTH "Smoothes the grid so that the geometry is linearized along the path of steepest descent."
#define	TOOLTIP_TANGENTIAL_SMOOTH "Smoothes vertices on a manifold."
#define	TOOLTIP_RPOJECT_TO_PLANE "Projects all selected elements to the specified plane"
#define	TOOLTIP_PROJECT_TO_LIMIT_PLOOP "Projects all vertices in the grid to their limit positions as defined by the piecewise loop scheme."
#define	TOOLTIP_PROJECT_TO_LIMIT_SMOOTH_BOUNDARY "Projects all boundary-vertices in the grid to their limit positions as defined by the b-spline subdivision scheme."
#define	TOOLTIP_SET_PIVOT "Sets the pivot point of the selected object."
#define	TOOLTIP_SET_PIVOT_TO_SELECTION_CENTER "Sets the pivot to the center of the current selection."
#define	TOOLTIP_SET_PIVOT_TO_MESH_CENTER "Sets the pivot to the center of the active mesh."
#define	TOOLTIP_FLATTEN_BENT_QUADRILATERALS "Flattens bent quadrilaterals using an iterative flattening scheme"
#define	TOOLTIP_APPLY_HEIGHT_FIELD "Calculates z-values of all nodes in terms of their x and y values." 

//	file io
#define TOOLTIP_LOAD_MESH ""
#define TOOLTIP_SAVE_MESH ""
#define TOOLTIP_EXPORT_TO_UG3 ""


//mark tools
#define	TOOLTIP_CLEAR_MARKS "Clears all marks"
#define	TOOLTIP_MARK_CREASE_EDGES "Marks edges whose associated faces have a certain angle as crease-edge."
#define	TOOLTIP_MARK_SELECTION "Marks selected vertices and edges."
#define	TOOLTIP_UNMARK_SELECTION "Unmarks selected elements."
#define TOOLTIP_MARK_CORNERS_OF_MARKED_EDGES "Marks selected vertices as fixed, if they lie at a sharp corner of a marked path or if they are at endpoints or at junctions of marked edges."

//fracture tools
#define	TOOLTIP_EXPAND_LAYERS_2D "Expands a 1d layer to a 2d layer by introducing quadrilaterals."
#define	TOOLTIP_EXPAND_LAYERS_3D "Expands a 2d layer to a 3d layer by introducing prisms."
#define	TOOLTIP_FRAC_TO_LAYER "Enhances a 2d fracture to a 3d fracture."

//info tools
#define TOOLTIP_MEASURE_GRID_LENGTH "Measures the length of all edges of a grid"
#define TOOLTIP_MEASURE_GRID_AREA "Measures the area of all faces of a grid"
#define TOOLTIP_MEASURE_GRID_VOLUME "Measures the volume of all volume elements of a grid"

#define TOOLTIP_MEASURE_SUBSET_LENGTH "Measures the length of all edges of the given subset"
#define TOOLTIP_MEASURE_SUBSET_AREA "Measures the area of all faces of the given subset"
#define TOOLTIP_MEASURE_SUBSET_VOLUME "Measures the volume of all volume elements of the given subset"

#define TOOLTIP_MEASURE_SELECTION_LENGTH "Measures the length of all edges of the current selection"
#define TOOLTIP_MEASURE_SELECTION_AREA "Measures the area of all faces of the current selection"
#define TOOLTIP_MEASURE_SELECTION_VOLUME "Measures the volume of all volume elements of the current selection"

//new tools
#define TOOLTIP_GET_BOUNDING_BOX "Returns the bounding box of the specified mesh."
#define TOOLTIP_SELECT_VERTEX_IN_BOX "Selects all vertices in the given box"
#define TOOLTIP_SELECT_EDGE_IN_BOX "Selects all edges in the given box"
#define TOOLTIP_SELECT_FACE_IN_BOX "Selects all faces in the given box"
#define TOOLTIP_SELECT_VOLUME_IN_BOX "Selects all volumes in the given box"
#define TOOLTIP_SELECT_VERTEX_IN_CYLINDER "Selects all vertices in the given cylinder"
#define TOOLTIP_SELECT_EDGE_IN_CYLINDER "Selects all edges in the given cylinder"
#define TOOLTIP_SELECT_FACE_IN_CYLINDER "Selects all faces in the given cylinder"
#define TOOLTIP_SELECT_VOLUME_IN_CYLINDER "Selects all volumes in the given cylinder"




//	all changes should go above the endif statement
#endif
