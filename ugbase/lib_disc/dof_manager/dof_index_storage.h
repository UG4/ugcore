/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel
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
		[[nodiscard]] ConstSmartPtr<MultiGrid> multi_grid() const {return m_spMG;}

		///	returns first algebra index of a geometric object
		/// \{
			   size_t& obj_index(GridObject* obj);
		inline size_t& obj_index(Vertex* vrt) {return m_aaIndexVRT[vrt];}
		inline size_t& obj_index(Edge* ed) {return m_aaIndexEDGE[ed];}
		inline size_t& obj_index(Face* face) {return m_aaIndexFACE[face];}
		inline size_t& obj_index(Volume* vol) {return m_aaIndexVOL[vol];}
		/// \}

		///	const access to first algebra index of a geometric object
		/// \{
			   const size_t& obj_index(GridObject* obj) const;
		inline const size_t& obj_index(Vertex* vrt) const {return m_aaIndexVRT[vrt];}
		inline const size_t& obj_index(Edge* ed) const {return m_aaIndexEDGE[ed];}
		inline const size_t& obj_index(Face* face) const {return m_aaIndexFACE[face];}
		inline const size_t& obj_index(Volume* vol) const {return m_aaIndexVOL[vol];}
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
		using ADoF = Attachment<size_t>;
		ADoF m_aIndex;

		///	Attachment Accessors
		///	\{
		using vertex_attachment_accessor_type = Grid::AttachmentAccessor<Vertex, ADoF>;
		using edge_attachment_accessor_type = Grid::AttachmentAccessor<Edge, ADoF>;
		using face_attachment_accessor_type = Grid::AttachmentAccessor<Face, ADoF>;
		using volume_attachment_accessor_type = Grid::AttachmentAccessor<Volume, ADoF>;
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


#endif