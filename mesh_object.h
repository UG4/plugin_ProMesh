/*
 * Copyright (c) 2013-2014:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__mesh_object__
#define __H__UG__mesh_object__

#include "lib_grid/grid/grid.h"
#include "lib_grid/common_attachments.h"
#include "lib_grid/selector.h"
#include "lib_grid/subset_handler.h"
#include "lib_grid/algorithms/remeshing/edge_length_adjustment.h"

namespace ug{
namespace promesh{


class MeshObject
{
	public:
		typedef APosition position_attachment_t;
		typedef Grid::VertexAttachmentAccessor<position_attachment_t>	position_accessor_t;

		typedef ANormal normal_attachment_t;
		typedef Grid::FaceAttachmentAccessor<normal_attachment_t>	normal_accessor_t;

		typedef ANumber volume_constraint_attachment_t;
		typedef Grid::VolumeAttachmentAccessor<volume_constraint_attachment_t>	volume_constraint_accessor_t;

		MeshObject() : m_creaseHandler(SHE_VERTEX | SHE_EDGE)
		{
			m_grid.attach_to_vertices(aPosition);
			m_aaPos.access(m_grid, aPosition);
			m_grid.attach_to_faces(aNormal);
			m_aaNorm.access(m_grid, aNormal);
			m_subsetHandler.assign_grid(m_grid);
			m_subsetHandler.enable_strict_inheritance(true);
			m_creaseHandler.assign_grid(m_grid);
			m_creaseHandler.subset_info(REM_CREASE).name = "crease";
			m_creaseHandler.subset_info(REM_FIXED).name = "fixed";
			m_selector.assign_grid(m_grid);
			m_pivot = vector3(0, 0, 0);
		}

		virtual ~MeshObject()	{}

		Grid&	get_grid()					{return m_grid;}
		SubsetHandler& get_subset_handler()	{return m_subsetHandler;}
		SubsetHandler& get_crease_handler()	{return m_creaseHandler;}
		Selector& get_selector()			{return m_selector;}

		Grid&			grid()				{return m_grid;}
		SubsetHandler&	subset_handler()	{return m_subsetHandler;}
		SubsetHandler&	crease_handler()	{return m_creaseHandler;}
		Selector&		selector()			{return m_selector;}

	//	pivot
		void set_pivot(const vector3& pivot)	{m_pivot = pivot;}
		const vector3& get_pivot() const		{return m_pivot;}

		vector3& pivot()						{return m_pivot;}
		const vector3& pivot() const			{return m_pivot;}

	///	returns accessor to vertex positions
		position_accessor_t& position_accessor()		{return m_aaPos;}
		position_attachment_t& position_attachment()	{return aPosition;}

	///	returns accessor to face normals
		normal_accessor_t& normal_accessor()		{return m_aaNorm;}
		normal_attachment_t& normal_attachment()	{return aNormal;}

	///	returns accessor to volume constraints.
		volume_constraint_accessor_t& volume_constraint_accessor()
		{
			volume_constraints_required();
			return m_aaVolumeConstraint;
		}

	///	returns the volume constraint attachment
		volume_constraint_attachment_t& volume_constraint_attachment()
		{
			volume_constraints_required();
			return m_aVolumeConstraint;
		}

	///	clears the volume constraints (removes the attachment)
		void clear_volume_constraints()
		{
			if(m_aaVolumeConstraint.valid()){
				m_grid.detach_from_volumes(m_aVolumeConstraint);
				m_aaVolumeConstraint.invalidate();
			}
		}

	protected:
		void volume_constraints_required()
		{
			if(!m_aaVolumeConstraint.valid()){
				m_grid.attach_to_volumes(m_aVolumeConstraint);
				m_aaVolumeConstraint.access(m_grid, m_aVolumeConstraint);
			}
		}

	protected:
		Grid				m_grid;
		SubsetHandler		m_subsetHandler;
		SubsetHandler		m_creaseHandler;
		Selector			m_selector;
		position_accessor_t	m_aaPos;
		normal_accessor_t	m_aaNorm;
		vector3				m_pivot;
		volume_constraint_attachment_t		m_aVolumeConstraint;
		volume_constraint_accessor_t		m_aaVolumeConstraint;

};


}}// end of namespace

#endif
