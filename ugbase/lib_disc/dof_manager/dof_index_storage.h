/*
 * dof_index_storage.h
 *
 *  Created on: 29.11.2011
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__DOF_MANAGER__DOF_INDEX_STORAGE__
#define __H__UG__LIB_DISC__DOF_MANAGER__DOF_INDEX_STORAGE__

#include "dof_distribution_info.h"

namespace ug{

class DoFIndexStorage : public DoFDistributionInfoProvider
{
	public:
		DoFIndexStorage(SmartPtr<MultiGrid> spMG,
		                ConstSmartPtr<DoFDistributionInfo> spDDInfo);

		~DoFIndexStorage();

	public:
		///	returns the multigrid
		SmartPtr<MultiGrid> multi_grid() {return m_spMG;}
		ConstSmartPtr<MultiGrid> multi_grid() const {return m_spMG;}

		///	returns first algebra index of a geometric object
		/// \{
			   size_t& obj_index(GridObject* obj);
		inline size_t& obj_index(Vertex* vrt) 	{return m_aaIndexVRT[vrt];}
		inline size_t& obj_index(Edge* ed) 		{return m_aaIndexEDGE[ed];}
		inline size_t& obj_index(Face* face)     	{return m_aaIndexFACE[face];}
		inline size_t& obj_index(Volume* vol)     	{return m_aaIndexVOL[vol];}
		/// \}

		///	const access to first algebra index of a geometric object
		/// \{
			   const size_t& obj_index(GridObject* obj) const;
		inline const size_t& obj_index(Vertex* vrt) const {return m_aaIndexVRT[vrt];}
		inline const size_t& obj_index(Edge* ed)    const {return m_aaIndexEDGE[ed];}
		inline const size_t& obj_index(Face* face)      const {return m_aaIndexFACE[face];}
		inline const size_t& obj_index(Volume* vol)     const {return m_aaIndexVOL[vol];}
		/// \}

	protected:
		/// initializes the attachments
		void init_attachments();

		/// removes the attachments
		void clear_attachments();

	protected:
		///	Multi Grid
		SmartPtr<MultiGrid> m_spMG;

		///	Attachment type
		typedef ug::Attachment<size_t> ADoF;
		ADoF m_aIndex;

		///	Attachment Accessors
		///	\{
		typedef Grid::AttachmentAccessor<Vertex, ADoF> vertex_attachment_accessor_type;
		typedef Grid::AttachmentAccessor<Edge, ADoF> edge_attachment_accessor_type;
		typedef Grid::AttachmentAccessor<Face, ADoF> face_attachment_accessor_type;
		typedef Grid::AttachmentAccessor<Volume, ADoF> volume_attachment_accessor_type;
		/// \}

		///	Attachments
		///	\{
		vertex_attachment_accessor_type m_aaIndexVRT;
		edge_attachment_accessor_type m_aaIndexEDGE;
		face_attachment_accessor_type m_aaIndexFACE;
		volume_attachment_accessor_type m_aaIndexVOL;
		///	\}
};

} // end namespace ug


#endif /* __H__UG__LIB_DISC__DOF_MANAGER__DOF_INDEX_STORAGE__ */
