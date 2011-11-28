/*
 * dof_storage_manager.cpp
 *
 *  Created on: 28.11.2011
 *      Author: andreasvogel
 */

#include "dof_storage_manager.h"

namespace ug{

void DoFStorageManager::set_subset_handler(ISubsetHandler* pSH)
{
	if(pSH == NULL)
		UG_THROW_FATAL("Null pointer passed as ISubsetHandler.");

//	nothing to do if equal
	if(pSH == m_pSH) return;

//	clear first  if not equal
	if(m_pSH != NULL) clear();

//	remember new one
	m_pSH = pSH;
}

void DoFStorageManager::clear()
{
	if(grid() == NULL)
		UG_THROW_FATAL("No Grid in SubsetHandler used in DoFStorageManager.");

//	detach DoFs
	if(m_aaIndexVRT.valid()) grid()->detach_from<VertexBase>(m_aIndex);
	if(m_aaIndexEDGE.valid()) grid()->detach_from<EdgeBase>(m_aIndex);
	if(m_aaIndexFACE.valid()) grid()->detach_from<Face>(m_aIndex);
	if(m_aaIndexVOL.valid()) grid()->detach_from<Volume>(m_aIndex);

	m_aaIndexVRT.invalidate();
	m_aaIndexEDGE.invalidate();
	m_aaIndexFACE.invalidate();
	m_aaIndexVOL.invalidate();
}

void DoFStorageManager::update_attachments(AttachmentType type)
{
	if(grid() == NULL)
		UG_THROW_FATAL("No Grid in SubsetHandler used in DoFStorageManager.");

//	attach DoFs to vertices
	if(type & DSM_VERTEX) {
		grid()->attach_to<VertexBase>(m_aIndex);
		m_aaIndexVRT.access(*grid(), m_aIndex);
	}
	if(type & DSM_EDGE) {
		grid()->attach_to<EdgeBase>(m_aIndex);
		m_aaIndexEDGE.access(*grid(), m_aIndex);
	}
	if(type & DSM_FACE) {
		grid()->attach_to<Face>(m_aIndex);
		m_aaIndexFACE.access(*grid(), m_aIndex);
	}
	if(type & DSM_VOLUME) {
		grid()->attach_to<Volume>(m_aIndex);
		m_aaIndexVOL.access(*grid(), m_aIndex);
	}

}

} // end namespace ug
