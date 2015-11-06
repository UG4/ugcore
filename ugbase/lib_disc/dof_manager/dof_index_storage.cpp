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

#include "dof_index_storage.h"

namespace ug{

DoFIndexStorage::
DoFIndexStorage(SmartPtr<MultiGrid> spMG,
                ConstSmartPtr<DoFDistributionInfo> spDDInfo)
:	DoFDistributionInfoProvider(spDDInfo),
 	m_spMG(spMG)
{
	init_attachments();
}

DoFIndexStorage::
~DoFIndexStorage()
{
	clear_attachments();
}

void DoFIndexStorage::init_attachments()
{
//	attach DoFs to vertices
	if(max_dofs(VERTEX)) {
		multi_grid()->attach_to_dv<Vertex>(m_aIndex, (size_t)-1);
		m_aaIndexVRT.access(*multi_grid(), m_aIndex);
	}
	if(max_dofs(EDGE)) {
		multi_grid()->attach_to_dv<Edge>(m_aIndex, (size_t)-1);
		m_aaIndexEDGE.access(*multi_grid(), m_aIndex);
	}
	if(max_dofs(FACE)) {
		multi_grid()->attach_to_dv<Face>(m_aIndex, (size_t)-1);
		m_aaIndexFACE.access(*multi_grid(), m_aIndex);
	}
	if(max_dofs(VOLUME)) {
		multi_grid()->attach_to_dv<Volume>(m_aIndex, (size_t)-1);
		m_aaIndexVOL.access(*multi_grid(), m_aIndex);
	}
}

void DoFIndexStorage::clear_attachments()
{
//	detach DoFs
	if(m_aaIndexVRT.valid()) multi_grid()->detach_from<Vertex>(m_aIndex);
	if(m_aaIndexEDGE.valid()) multi_grid()->detach_from<Edge>(m_aIndex);
	if(m_aaIndexFACE.valid()) multi_grid()->detach_from<Face>(m_aIndex);
	if(m_aaIndexVOL.valid()) multi_grid()->detach_from<Volume>(m_aIndex);

	m_aaIndexVRT.invalidate();
	m_aaIndexEDGE.invalidate();
	m_aaIndexFACE.invalidate();
	m_aaIndexVOL.invalidate();
}

size_t& DoFIndexStorage::obj_index(GridObject* obj)
{
	switch(obj->base_object_id())
	{
		case VERTEX: return obj_index(static_cast<Vertex*>(obj));
		case EDGE:   return obj_index(static_cast<Edge*>(obj));
		case FACE:   return obj_index(static_cast<Face*>(obj));
		case VOLUME: return obj_index(static_cast<Volume*>(obj));
		default: UG_THROW("Base Object type not found.");
	}
}

const size_t& DoFIndexStorage::obj_index(GridObject* obj) const
{
	switch(obj->base_object_id())
	{
		case VERTEX: return obj_index(static_cast<Vertex*>(obj));
		case EDGE:   return obj_index(static_cast<Edge*>(obj));
		case FACE:   return obj_index(static_cast<Face*>(obj));
		case VOLUME: return obj_index(static_cast<Volume*>(obj));
		default: UG_THROW("Base Object type not found.");
	}
};

} // end namespace ug
