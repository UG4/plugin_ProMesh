/*
 * Copyright (c) 2014-2015:  G-CSC, Goethe University Frankfurt
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
#include "lib_grid/grid/geometry.h"
#include "lib_grid/common_attachments.h"
#include "lib_grid/selector.h"
#include "lib_grid/subset_handler.h"
#include "lib_grid/refinement/projectors/projection_handler.h"
#include "lib_grid/algorithms/remeshing/edge_length_adjustment.h"

namespace ug{
namespace promesh{


template <class TElem>
class ElementIterator{
	public:
		typedef typename Grid::traits<TElem>::iterator	iterator_t;

		ElementIterator(iterator_t i) : m_iter(i)	{}

		SmartPtr<ElementIterator<TElem> >
		clone()
		{return make_sp(new ElementIterator(m_iter));}
		
		void	assign(ElementIterator& i)	{m_iter = i.m_iter;}
		TElem*	value()						{return *m_iter;}
		void	advance()					{++m_iter;}
		bool	equal(ElementIterator& i)	{return (m_iter == i.m_iter);}
		bool	unequal(ElementIterator& i)	{return (m_iter != i.m_iter);}

	private:
		iterator_t	m_iter;
};


/// \addtogroup promesh
/// \{
class Mesh
{
	public:
		typedef APosition 	position_attachment_t;
		typedef Grid::VertexAttachmentAccessor<position_attachment_t>	position_accessor_t;

		typedef ANormal normal_attachment_t;
		typedef Grid::FaceAttachmentAccessor<normal_attachment_t>	normal_accessor_t;

		typedef ANumber volume_constraint_attachment_t;
		typedef Grid::VolumeAttachmentAccessor<volume_constraint_attachment_t>	volume_constraint_accessor_t;

		typedef ElementIterator<Vertex>	vertex_iter_t;
		typedef ElementIterator<Edge>	edge_iter_t;
		typedef ElementIterator<Face>	face_iter_t;
		typedef ElementIterator<Volume>	volume_iter_t;

		typedef SmartPtr<vertex_iter_t>	sp_vertex_iter_t;
		typedef SmartPtr<edge_iter_t>	sp_edge_iter_t;
		typedef SmartPtr<face_iter_t>	sp_face_iter_t;
		typedef SmartPtr<volume_iter_t>	sp_volume_iter_t;

		Mesh();
		Mesh(const Mesh& m);

		virtual ~Mesh()	{}

		Grid&			grid()						{return m_grid;}
		SubsetHandler&	subset_handler()			{return m_subsetHandler;}
		SubsetHandler&	crease_handler()			{return m_creaseHandler;}
		Selector&		selector()					{return m_selector;}
		ProjectionHandler&	projection_handler()	{return m_projectionHandler;}

		SmartPtr<Geometry<3, 3> >		geometry()			{return m_geometry;}
		ConstSmartPtr<Geometry<3, 3> >	geometry() const	{return m_geometry;}

	//	pivot
		void set_pivot(const vector3& pivot)	{m_pivot = pivot;}
		vector3& pivot()						{return m_pivot;}

	///	returns accessor to vertex positions
		position_accessor_t& position_accessor()		{return m_aaPos;}
		position_attachment_t& position_attachment()	{return aPosition;}

		void set_position(Vertex* v, const vector3& p)	{m_aaPos[v] = p;}
		vector3& position(Vertex* v)					{return m_aaPos[v];}

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

	///	element creation and deletion
		Vertex*	create_vertex(const vector3& p);
		Edge*	create_edge(Vertex* v0, Vertex* v1);
		Face*	create_triangle(Vertex* v0, Vertex* v1, Vertex* v2);
		Face*	create_quadrilateral(Vertex* v0, Vertex* v1, Vertex* v2, Vertex* v3);
		Volume*	create_tetrahedron(Vertex* v0, Vertex* v1, Vertex* v2, Vertex* v3);
		Volume*	create_pyramid(Vertex* v0, Vertex* v1, Vertex* v2,
							   Vertex* v3, Vertex* v4);
		Volume*	create_prism(Vertex* v0, Vertex* v1, Vertex* v2,
							 Vertex* v3, Vertex* v4, Vertex* v5);
		Volume*	create_hexahedron(Vertex* v0, Vertex* v1, Vertex* v2, Vertex* v3,
							 	  Vertex* v4, Vertex* v5, Vertex* v6, Vertex* v7);
		Volume*	create_octahedron(Vertex* v0, Vertex* v1, Vertex* v2,
							 	  Vertex* v3, Vertex* v4, Vertex* v5);

	//	iterators
	///	returns an iterator to the first element of the specified type
	/**	\note	This iterator is intended for use through bindings only! If you're
	 *			actually writing C++ code, please use the iterators of the
	 *			Grid/Selector/SubsetHandler classes directly.
	 * \{ */
		template <class TElem>
		SmartPtr<ElementIterator<typename TElem::grid_base_object> >
		begin()
		{return make_sp(new ElementIterator<typename TElem::grid_base_object>(grid().begin<TElem>()));}

		template <class TElem>
		SmartPtr<ElementIterator<typename TElem::grid_base_object> >
		selection_begin()
		{return make_sp(new ElementIterator<typename TElem::grid_base_object>(selector().begin<TElem>()));}

		template <class TElem>
		SmartPtr<ElementIterator<typename TElem::grid_base_object> >
		subset_begin(int si)
		{return make_sp(new ElementIterator<typename TElem::grid_base_object>(subset_handler().begin<TElem>(si)));}
	/** \} */

	///	returns an iterator to the position behind the last element of the specified type
	/**	\note	This iterator is intended for use through bindings only! If you're
	 *			actually writing C++ code, please use the iterators of the
	 *			Grid/Selector/SubsetHandler classes directly.
	 * \{ */
		template <class TElem>
		SmartPtr<ElementIterator<typename TElem::grid_base_object> >
		end()
		{return make_sp(new ElementIterator<typename TElem::grid_base_object>(grid().end<TElem>()));}

		template <class TElem>
		SmartPtr<ElementIterator<typename TElem::grid_base_object> >
		selection_end()
		{return make_sp(new ElementIterator<typename TElem::grid_base_object>(selector().end<TElem>()));}

		template <class TElem>
		SmartPtr<ElementIterator<typename TElem::grid_base_object> >
		subset_end(int si)
		{return make_sp(new ElementIterator<typename TElem::grid_base_object>(subset_handler().end<TElem>(si)));}
	/** \} */

	protected:
		void volume_constraints_required()
		{
			if(!m_aaVolumeConstraint.valid()){
				m_grid.attach_to_volumes_dv(m_aVolumeConstraint, -1, true);
				m_aaVolumeConstraint.access(m_grid, m_aVolumeConstraint);
			}
		}

	protected:
		void init();

		Grid				m_grid;
		SubsetHandler		m_subsetHandler;
		SubsetHandler		m_creaseHandler;
		Selector			m_selector;
		ProjectionHandler	m_projectionHandler;
		position_accessor_t	m_aaPos;
		normal_accessor_t	m_aaNorm;
		vector3				m_pivot;
		volume_constraint_attachment_t		m_aVolumeConstraint;
		volume_constraint_accessor_t		m_aaVolumeConstraint;
		SmartPtr<Geometry<3, 3> >			m_geometry;

};

/// \}
}}// end of namespace

#endif
