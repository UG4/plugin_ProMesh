# Copyright (c) 2013-2017:  G-CSC, Goethe University Frankfurt
# Author: Sebastian Reiter
# 
# This file is part of UG4.
# 
# UG4 is free software: you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License version 3 (as published by the
# Free Software Foundation) with the following additional attribution
# requirements (according to LGPL/GPL v3 §7):
# 
# (1) The following notice must be displayed in the Appropriate Legal Notices
# of covered and combined works: "Based on UG4 (www.ug4.org/license)".
# 
# (2) The following notice must be displayed at a prominent place in the
# terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
# 
# (3) The following bibliography is recommended for citation and must be
# preserved in all covered files:
# "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
#   parallel geometric multigrid solver on hierarchically distributed grids.
#   Computing and visualization in science 16, 4 (2013), 151-164"
# "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
#   flexible software system for simulating pde based models on high performance
#   computers. Computing and visualization in science 16, 4 (2013), 165-179"
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.

################################################################################
# PROMESH
#
# The promesh plugin.
################################################################################

# Name of your plugin, path to your plugin (relative to ug's root path) and sources.
set(pluginName	ProMesh)
set(SOURCES		mesh.cpp
				promesh_plugin.cpp
				register_coordinate_transform_tools.cpp
				register_mesh.cpp
				register_meshing_tools.cpp
				register_selection_tools.cpp
				register_subset_tools.cpp
				tools/coordinate_transform_tools.cpp
				tools/file_io_tools.cpp
				tools/grid_generation_tools.cpp
				tools/measure_tools.cpp
				tools/refinement_tools.cpp
				tools/remeshing_tools.cpp
				tools/selection_tools.cpp
				tools/subset_tools.cpp
				tools/topology_tools.cpp)


################################################################################
# The code below doesn't have to be changed (usually)
################################################################################
cmake_minimum_required(VERSION 2.6)

project(UG_PLUGIN_${pluginName})


option(DEBUG_PATH "Specify a path to which debug output may be written. Only for development builds.")

message(STATUS "Info: Additional options for the ProMesh plugin:")
message(STATUS "Info:   DEBUG_PATH   Specify a path to which debug output may be written. Only for development builds.")

if(DEBUG_PATH)
	message(STATUS "Info:   DEBUG_PATH = \"${DEBUG_PATH}\"")
	ExportDefinitions(-DPROMESH_DEBUG_PATH=\"${DEBUG_PATH}\")
endif(DEBUG_PATH)


# include the definitions and dependencies for ug-plugins.
include(${UG_ROOT_CMAKE_PATH}/ug_plugin_includes.cmake)

if(buildEmbeddedPlugins)
	# add the sources to ug4's sources
	EXPORTSOURCES(${CMAKE_CURRENT_SOURCE_DIR} ${SOURCES})
else(buildEmbeddedPlugins)
	# create a shared library from the sources and link it against ug4.
	add_library(${pluginName} SHARED ${SOURCES})
	target_link_libraries (${pluginName} ug4)
endif(buildEmbeddedPlugins)
