
#include "element_data_items.h"

namespace ug{

bool
DataPossibilityItem::
delete_data_export(DataExportItem* exportItem)
{
	std::vector<DataExportItem*>::iterator iter;
	iter = find(m_createdDataExports.begin(), m_createdDataExports.end(), exportItem);
	if(iter == m_createdDataExports.end())
	{
		UG_ASSERT(0, "DataExportItem must exist in list, but does not. Cannot remove from list.");
		return false;
	}
	m_createdDataExports.erase(iter);

	return true;
}

bool
DataPossibilityItem::
clear_slot(std::size_t slot) {
	UG_ASSERT(slot < m_num_slots, "Slot does not exist.");

	// remove from all possibilities slots this possibility is linked to.
	for(std::size_t i = 0; i < m_linkedPosItems.size(); ++i)
	{
		if((m_linkedPosItems[i].item)->clear_slot(m_linkedPosItems[i].slot) != true) return false;
	}

	// remove from linkPosItem list of plugged in possibility
	LinkedPossibility linkedPos; linkedPos.item = this; linkedPos.slot = slot;
	if(m_slotPosItems[slot]->remove_linked_possibility(linkedPos) != true)
	{
		UG_LOG("DataPossibilityItem::clear_slot: Error while clearing slot.\n");
		return false;
	}

	// delete slot
	m_slotPosItems[slot] = NULL;

	// TODO: Now not instantiable anymore. Should we remove all created exports ???

	return true;
};


void
DataPossibilityItem::
add_linked_possibility(LinkedPossibility& linkedPos)
{
	m_linkedPosItems.push_back(linkedPos);
}

bool
DataPossibilityItem::
remove_linked_possibility(LinkedPossibility& linkedPos)
{
	for(std::size_t i = 0; i < m_linkedPosItems.size(); ++i)
	{
		if(m_linkedPosItems[i].item == linkedPos.item)
		{
			if(m_linkedPosItems[i].slot == linkedPos.slot)
			{
					m_linkedPosItems.erase(m_linkedPosItems.begin() + i);
					return true;
			}
		}
	}
	return false;
}


}
