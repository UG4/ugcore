
#include "data_items.h"
#include <algorithm>

namespace ug{

//////////////////////////////
// Data Possibility Item
//////////////////////////////

bool
DataPossibilityItem::
delete_data_export(DataExportItem* exportItem)
{
	std::vector<DataExportItem*>::iterator iter;
	iter = find(m_vCreatedDataExports.begin(), m_vCreatedDataExports.end(), exportItem);
	if(iter == m_vCreatedDataExports.end())
	{
		UG_ASSERT(0, "DataExportItem must exist in list, but does not. Cannot remove from list.");
		return false;
	}
	m_vCreatedDataExports.erase(iter);
	return true;
}

bool
DataPossibilityItem::
clear_slot(size_t slot) {
	// if not linkedm return false
	if(slot >= m_numSlots) return false;

	// remove from all possibilities slots this possibility is linked to as an export possibility
	for(size_t i = 0; i < m_vLinkedPosItems.size(); ++i)
	{
		if(!(m_vLinkedPosItems[i].item)->clear_slot(m_vLinkedPosItems[i].slot)) return false;
	}

	// remove from linkPosItem list of plugged in possibility
	LinkedPossibility linkedPos; linkedPos.item = this; linkedPos.slot = slot;
	if(!m_vSlotPosItems[slot]->remove_linked_possibility(linkedPos))
	{
		UG_LOG("DataPossibilityItem::clear_slot: Error while clearing slot.\n");
		return false;
	}

	// delete slot
	m_vSlotPosItems[slot] = NULL;

	// TODO: Now not instantiable anymore. Should we remove all created exports ???

	return true;
};


void
DataPossibilityItem::
add_linked_possibility(LinkedPossibility& linkedPos)
{
	m_vLinkedPosItems.push_back(linkedPos);
}

bool
DataPossibilityItem::
remove_linked_possibility(LinkedPossibility& linkedPos)
{
	for(size_t i = 0; i < m_vLinkedPosItems.size(); ++i)
	{
		if(m_vLinkedPosItems[i].item == linkedPos.item)
		{
			if(m_vLinkedPosItems[i].slot == linkedPos.slot)
			{
					m_vLinkedPosItems.erase(m_vLinkedPosItems.begin() + i);
					return true;
			}
		}
	}
	return false;
}


}
