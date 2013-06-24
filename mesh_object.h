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
		typedef Grid::VertexAttachmentAccessor<position_attachment_t>	position_accessor_t;

		typedef ANormal normal_attachment_t;
		typedef Grid::FaceAttachmentAccessor<normal_attachment_t>	normal_accessor_t;

		typedef ANumber volume_constraint_attachment_t;
		typedef Grid::VolumeAttachmentAccessor<volume_constraint_attachment_t>	volume_constraint_accessor_t;

		MeshObject()
		{
			m_grid.attach_to_vertices(aPosition);
			m_aaPos.access(m_grid, aPosition);
			m_grid.attach_to_faces(aNormal);
			m_aaNorm.access(m_grid, aNormal);
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
