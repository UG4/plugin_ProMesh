// created by Sebastian Reiter
// s.b.reiter@gmail.com
// Jun 14, 2013 (d,m,y)

#ifndef __H__UG__mesh_object__
#define __H__UG__mesh_object__

#include "lib_grid/grid/grid.h"
#include "lib_grid/common_attachments.h"
#include "lib_grid/selector.h"
#include "lib_grid/subset_handler.h"


namespace ug{
namespace promesh{


class MeshObject
{
	public:
		typedef APosition position_attachment_t;
		typedef Grid::VertexAttachmentAccessor<APosition>	position_accessor_t;

		MeshObject()
		{
			m_grid.attach_to_vertices(aPosition);
			m_aaPos.access(m_grid, aPosition);
			m_grid.attach_to_faces(aNormal);
			m_subsetHandler.assign_grid(m_grid);
			m_subsetHandler.enable_strict_inheritance(true);
			m_creaseHandler.assign_grid(m_grid);
			m_selector.assign_grid(m_grid);
			m_pivot = vector3(0, 0, 0);
		}

		virtual ~MeshObject()	{}

		Grid&	get_grid()					{return m_grid;}
		SubsetHandler& get_subset_handler()	{return m_subsetHandler;}
		SubsetHandler& get_crease_handler()	{return m_creaseHandler;}
		Selector& get_selector()			{return m_selector;}

	//	pivot
		void set_pivot(const vector3& pivot)	{m_pivot = pivot;}
		const vector3& get_pivot() const		{return m_pivot;}

		vector3& pivot()						{return m_pivot;}
		const vector3& pivot() const			{return m_pivot;}

		position_accessor_t& position_accessor()	{return m_aaPos;}
		position_attachment_t& position_attachment()	{return aPosition;}

	protected:
		Grid				m_grid;
		SubsetHandler		m_subsetHandler;
		SubsetHandler		m_creaseHandler;
		Selector			m_selector;
		position_accessor_t	m_aaPos;
		vector3				m_pivot;

};


}}// end of namespace

#endif
