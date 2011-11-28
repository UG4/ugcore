/*
 * dof_storage_manager.h
 *
 *  Created on: 28.11.2011
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__DOF_MANAGER__STORAGE_MANAGER__
#define __H__UG__LIB_DISC__DOF_MANAGER__STORAGE_MANAGER__

#include "lib_grid/lg_base.h"

namespace ug{

/// accesses storage for DoFs and handles grid attachment process
class DoFStorageManager
{
	public:
	///	type of DoF attachment
		typedef ug::Attachment<size_t> ADoF;

	///	type of accessor
		typedef Grid::AttachmentAccessor<VertexBase, ADoF>
				vertex_attachment_accessor_type;

	///	type of accessor
		typedef Grid::AttachmentAccessor<EdgeBase, ADoF>
				edge_attachment_accessor_type;

	///	type of accessor
		typedef Grid::AttachmentAccessor<Face, ADoF>
				face_attachment_accessor_type;

	///	type of accessor
		typedef Grid::AttachmentAccessor<Volume, ADoF>
				volume_attachment_accessor_type;

	///	type of needed attachment
		enum AttachmentType
		{
			DSM_VERTEX = 1 << 0,
			DSM_EDGE = 1 << 1,
			DSM_FACE = 1 << 2,
			DSM_VOLUME = 1 << 3,
			DSM_ALL = DSM_VERTEX | DSM_EDGE | DSM_FACE | DSM_VOLUME
		};

	public:
	///	Constructor
		DoFStorageManager() : m_pSH(NULL) {}

	/// attach indices
		void update_attachments(AttachmentType type);

	///	returns the associated grid
		Grid* grid() {return m_pSH->get_assigned_grid();}

	///	sets the underlying subset handler
		void set_subset_handler(ISubsetHandler* pSH);

	///	returns the underlying subset handler
		ISubsetHandler* subset_handler() {return m_pSH;}

	///	returns the attachment accessor
		vertex_attachment_accessor_type& vertex_att_acc() 	{return m_aaIndexVRT;}
		edge_attachment_accessor_type& edge_att_acc()		{return m_aaIndexEDGE;}
		face_attachment_accessor_type& face_att_acc()		{return m_aaIndexFACE;}
		volume_attachment_accessor_type& volume_att_acc() 	{return m_aaIndexVOL;}

	/// clear all dofs
		void clear();

	/// destructor
		~DoFStorageManager() {clear();};

	protected:
	/// subset handler
		ISubsetHandler* m_pSH;

	///	Attachment Accessor
		vertex_attachment_accessor_type m_aaIndexVRT;
		edge_attachment_accessor_type m_aaIndexEDGE;
		face_attachment_accessor_type m_aaIndexFACE;
		volume_attachment_accessor_type m_aaIndexVOL;

	///	Attachment (for vertices)
		ADoF m_aIndex;
};


} // end namespace ug

#endif /* __H__UG__LIB_DISC__DOF_MANAGER__STORAGE_MANAGER__ */
