################################################################################
# PROMESH
#
# The promesh plugin.
################################################################################

# Name of your plugin, path to your plugin (relative to ug's root path) and sources.
set(pluginName	ProMesh)
set(pluginPath	plugins/experimental/promesh)
set(SOURCES		mesh.cpp
				promesh_plugin.cpp)


################################################################################
# The code below doesn't have to be changed (usually)
################################################################################
cmake_minimum_required(VERSION 2.6)

project(UG_PLUGIN_${pluginName})

# include the definitions and dependencies for ug-plugins.
include("../../../cmake/ug_plugin_includes.cmake")

if(buildEmbeddedPlugins)
	# add the sources to ug4's sources
	EXPORTSOURCES(${pluginPath} ${SOURCES})
else(buildEmbeddedPlugins)
	# create a shared library from the sources and link it against ug4.
	add_library(${pluginName} SHARED ${SOURCES})
	target_link_libraries (${pluginName} ug4)
endif(buildEmbeddedPlugins)

