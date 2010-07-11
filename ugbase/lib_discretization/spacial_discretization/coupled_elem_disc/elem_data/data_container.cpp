/*
 * elementdata.cpp
 *
 *  Created on: 30.06.2009
 *      Author: andreasvogel
 */

#include "data_container.h"

namespace ug{

//////////////////////////////
// Data Container
//////////////////////////////

DataContainer::
~DataContainer()
{
	UG_DLOG(LIB_DISC_LINKER, 2, "DataContainer::~DataContainer: Deleting Data Container.\n");

	// delete export items
	for(size_t i = 0; i < m_vExportItemList.size(); ++i)
	{
		UG_DLOG(LIB_DISC_LINKER, 2, "DataContainer::~DataContainer: Deleting ExportItem '" << m_vExportItemList[i]->name() << "' [" << i << "].\n");
		delete m_vExportItemList[i];
	}
	m_vExportItemList.clear();
}

////////////////////////////////
// register
////////////////////////////////

bool
DataContainer::
register_item(DataImportItem& ImportItem)
{
	std::vector<DataImportItem*>::iterator iter;
	iter = find(m_vImportItemList.begin(), m_vImportItemList.end(), &ImportItem);
	if(iter != m_vImportItemList.end())
	{
		UG_LOG("Container already contains DataExportItem. Can not add.\n");
		return false;
	}

	m_vImportItemList.push_back(&ImportItem);
	return true;
}

bool
DataContainer::
register_item(DataExportItem& ExportItem)
{
	std::vector<DataExportItem*>::iterator iter;
	iter = find(m_vExportItemList.begin(), m_vExportItemList.end(), &ExportItem);
	if(iter != m_vExportItemList.end())
	{
		UG_LOG("Container already contains DataExportItem. Can not add.\n");
		return false;
	}
	m_vExportItemList.push_back(&ExportItem);
	return true;
}

bool
DataContainer::
register_item(DataPossibilityItem& PossibilityItem)
{
	std::vector<DataPossibilityItem*>::iterator iter;
	iter = find(m_vPossibilityItemList.begin(), m_vPossibilityItemList.end(), &PossibilityItem);
	if(iter != m_vPossibilityItemList.end())
	{
		UG_LOG("Container already contains DataPossibilityItem. Can not add.\n");
		return false;
	}
	m_vPossibilityItemList.push_back(&PossibilityItem);
	return true;
}

////////////////////////////////
// unregister
////////////////////////////////

bool
DataContainer::
unregister_item(DataImportItem& ImportItem)
{
	std::vector<DataImportItem*>::iterator iter;
	iter = find(m_vImportItemList.begin(), m_vImportItemList.end(), &ImportItem);
	if(iter == m_vImportItemList.end())
	{
		UG_LOG("DataContainer::unregister_item: Container does not ImportItem. Can not remove it.\n");
		return false;
	}

	m_vImportItemList.erase(iter);
	return true;
}

bool
DataContainer::
unregister_item(DataExportItem& ExportItem)
{
	std::vector<DataExportItem*>::iterator iter;
	iter = find(m_vExportItemList.begin(), m_vExportItemList.end(), &ExportItem);
	if(iter == m_vExportItemList.end())
	{
		UG_LOG("DataContainer::unregister_export: Container does not contain DataExportItem. Can not remove it.\n");
		return false;
	}

	m_vExportItemList.erase(iter);
	return true;
}

bool
DataContainer::
unregister_item(DataPossibilityItem& PossibilityItem)
{
	std::vector<DataPossibilityItem*>::iterator iter;
	iter = find(m_vPossibilityItemList.begin(), m_vPossibilityItemList.end(), &PossibilityItem);
	if(iter == m_vPossibilityItemList.end())
	{
		UG_LOG("DataContainer::unregister_item: Container does not PossibilityItem. Can not remove it.\n");
		return false;
	}

	UG_DLOG(LIB_DISC_LINKER, 2, "DataContainer::unregister_item: Deleting " << PossibilityItem.num_created_data_export_items()<<" Data Export created by Data Possibility " << PossibilityItem.name() << ".\n");
	for(size_t i = 0; i < PossibilityItem.num_created_data_export_items(); ++i)
	{
		unregister_item(*PossibilityItem.get_created_data_export_item(i));
	}

	m_vPossibilityItemList.erase(iter);
	return true;
}

////////////////////////////////
// linking
////////////////////////////////
bool
DataContainer::
link_possibility(size_t nr_pos1, size_t slot, size_t nr_pos2)
{
	// check if linking is possible
	if(nr_pos1 >= m_vPossibilityItemList.size()) return false;
	if(nr_pos2 >= m_vPossibilityItemList.size()) return false;

	return link_possibility(*m_vPossibilityItemList[nr_pos1], slot, *m_vPossibilityItemList[nr_pos2]);
}

bool
DataContainer::
link_possibility(DataPossibilityItem& posItem1, size_t slot, DataPossibilityItem& posItem2)
{
	std::vector<DataPossibilityItem*>::iterator posIter;

	posIter = find(m_vPossibilityItemList.begin(), m_vPossibilityItemList.end(), &posItem1);
	if(posIter == m_vPossibilityItemList.end())
	{
		UG_LOG("DataContainer::link: Container does not contain DataExportItem. Can not link.\n");
		return false;
	}
	posIter = find(m_vPossibilityItemList.begin(), m_vPossibilityItemList.end(), &posItem2);
	if(posIter == m_vPossibilityItemList.end())
	{
		UG_LOG("DataContainer::link: Container does not contain DataExportItem. Can not link.\n");
		return false;
	}

	if(!(slot < posItem1.num_slots()))
	{
		UG_LOG("DataContainer::link: Data slot not found.\n");
		return false;
	}

	// check, that import is not already linked to other export
	if(posItem1.is_linked(slot) == true)
	{
		UG_LOG("DataContainer::link: Slot already used. Invalid operation.\n");
		return false;
	}

	if(posItem1.link(&posItem2, slot) != true)
	{
		UG_LOG("DataContainer::link: Error while registering possibility '" << posItem1.name() << "' at linker '" << posItem2.name() <<"'. Can not link.\n");
		return false;
	}

	return true;
}

bool
DataContainer::
link(size_t nr_imp, size_t nr_exp)
{
	// check if linking is possible
	if(nr_imp >= m_vImportItemList.size())
		{UG_LOG("Import with number " << nr_imp << " does not exist.\n"); return false;}
	if(nr_exp >= m_vExportItemList.size())
		{UG_LOG("Export with number " << nr_exp << " does not exist.\n"); return false;}

	return link(*m_vImportItemList[nr_imp], *m_vExportItemList[nr_exp]);
}

bool
DataContainer::
link(size_t nr_exp1, size_t slot, size_t nr_exp2)
{
	// check if linking is possible
	if(nr_exp1 >= m_vExportItemList.size()) return false;
	if(nr_exp2 >= m_vExportItemList.size()) return false;

	return link(*m_vExportItemList[nr_exp1], slot, *m_vExportItemList[nr_exp2]);
}

bool
DataContainer::
link(DataExportItem& Export1Item, size_t slot, DataExportItem& Export2Item)
{
	std::vector<DataExportItem*>::iterator expIter;

	expIter = find(m_vExportItemList.begin(), m_vExportItemList.end(), &Export1Item);
	if(expIter == m_vExportItemList.end())
	{
		UG_LOG("DataContainer::link: Container does not contain DataExportItem. Can not link.\n");
		return false;
	}
	expIter = find(m_vExportItemList.begin(), m_vExportItemList.end(), &Export2Item);
	if(expIter == m_vExportItemList.end())
	{
		UG_LOG("DataContainer::link: Container does not contain DataExportItem. Can not link.\n");
		return false;
	}

	if(!(slot < Export1Item.num_slots()))
	{
		UG_LOG("DataContainer::link: Data Import slot not found.\n");
		return false;
	}

	// check, that import is not already linked to other export
	if(Export1Item.is_linked(slot) == true)
	{
		UG_LOG("DataContainer::link: DataImport already registered to an DataExport. Invalid operation.\n");
		return false;
	}

	if(Export1Item.link(&Export2Item, slot) != true)
	{
		UG_LOG("DataContainer::link: Error while registering export '" << Export1Item.name() << "' at linker '" << Export2Item.name() <<"'. Can not link.\n");
		return false;
	}

	return true;
}

bool
DataContainer::
link(DataImportItem& ImportItem, DataExportItem& ExportItem)
{
	// check, that import is not already linked to other export
	if(ImportItem.is_linked() == true)
	{
		UG_LOG("DataContainer::link: DataImport already registered to an DataExport. Invalid operation.\n");
		return false;
	}

	std::vector<DataImportItem*>::iterator importIter;
	std::vector<DataExportItem*>::iterator expIter;

	// Find requested Import in Import list
	importIter = find(m_vImportItemList.begin(), m_vImportItemList.end(), &ImportItem);
	if(importIter == m_vImportItemList.end())
	{
		UG_LOG("DataContainer::link: Container does not contain DataImportItem. Can not link.\n");
		return false;
	}

	// Find requested ExportPossibility in ExportPossibility list
	expIter = find(m_vExportItemList.begin(), m_vExportItemList.end(), &ExportItem);
	if(expIter == m_vExportItemList.end())
	{
		UG_LOG("DataContainer::link: Container does not contain DataExportItem. Can not link.\n");
		return false;
	}

	if(ImportItem.link_data_export(*expIter) != true)
	{
		UG_LOG("DataContainer::link: Error while registering export '" << (*expIter)->name() << "' at import '" << ImportItem.name() <<"'. Can not link.\n");
		return false;
	}

	return true;
}

DataExportItem*
DataContainer::
create_export(size_t nr_pos)
{
	// check if creation is possible
	if(nr_pos >= m_vPossibilityItemList.size()) return false;

	return create_export(* m_vPossibilityItemList[nr_pos]);
}

DataExportItem*
DataContainer::
create_export(DataPossibilityItem& PossibilityItem)
{
	// Find requested ExportPossibility in ExportPossibility list
	std::vector<DataPossibilityItem*>::iterator possibilityIter;
	possibilityIter = find(m_vPossibilityItemList.begin(), m_vPossibilityItemList.end(), &PossibilityItem);
	if(possibilityIter == m_vPossibilityItemList.end())
	{
		UG_LOG("DataContainer::create_export: Container does not contain DataExportPossiblityItem. Can not create Data Export.\n");
		return NULL;
	}

	// create new DataExport
	DataExportItem* expItem = PossibilityItem.create_data_export();
	if(expItem == NULL)
	{
		UG_LOG("DataContainer::create_export: Cannot create ExportItem (Not enough memory?).\n");
		return NULL;
	}

	// add export to list
	m_vExportItemList.push_back(expItem);

	// handle recursive creation of dependent exports if linker
	if(create_export_recursive(PossibilityItem, *expItem) != true)
	{
		UG_LOG("DataContainer::create_export: Error while creating recursive exports. Can not create Data Export.\n");
		return NULL;
	}

	UG_DLOG(LIB_DISC_LINKER, 2, "DataContainer::create_export: ExportItem '" << expItem->name() << "' created.\n");
	return expItem;
}

bool
DataContainer::
create_export_recursive(DataPossibilityItem& PossibilityItem, DataExportItem& exp)
{
	for(size_t i = 0; i < PossibilityItem.num_slots(); ++i)
	{
		UG_ASSERT(PossibilityItem.is_linked(i) == true, "Can only create recursively from Possibility if all slots are linked.");
		DataPossibilityItem* slotPos = PossibilityItem.get_linked_possibility(i);

		UG_DLOG(LIB_DISC_LINKER, 2, "DataContainer::create_export_recursive: Trying to create Export from Possibility '" << slotPos->name() << "'.\n");
		// Find requested ExportPossibility in ExportPossibility list
		std::vector<DataPossibilityItem*>::iterator possibilityIter;
		possibilityIter = find(m_vPossibilityItemList.begin(), m_vPossibilityItemList.end(), slotPos);
		if(possibilityIter == m_vPossibilityItemList.end())
		{
			UG_LOG("DataContainer::create_export_recursive: Container does not contain DataExportPossiblityItem. Can not create Data Export.\n");
			return false;
		}

		// create new DataExport
		DataExportItem* expItem = slotPos->create_data_export();
		if(expItem == NULL)
		{
			UG_LOG("DataContainer::create_export_recursive: Cannot create ExportItem (Not enough memory?).\n");
			return false;
		}

		std::vector<DataExportItem*>::iterator expIter;
		expIter = find(m_vExportItemList.begin(), m_vExportItemList.end(), &exp);
		if(expIter == m_vExportItemList.end())
		{
			UG_LOG("DataContainer::create_export_recursive: Container does not contain ExportItem.\n");
			return false;
		}

		// add export to list
		m_vExportItemList.insert(expIter, expItem);
		UG_DLOG(LIB_DISC_LINKER, 2, "DataContainer::create_export_recursive: ExportItem '" << expItem->name() << "' created.\n");

		// link directly to slot
		if(link(exp, i, *expItem) != true)
		{
			UG_LOG("DataContainer::create_export_recursive: Cannot link to recursive Data Export. Can not create Data Export.\n");
			return false;
		}
		UG_DLOG(LIB_DISC_LINKER, 2, "DataContainer::create_export_recursive: ExportItem '" << expItem->name() << "' linked to ExportItem '" << exp.name() << "', slot " << i << ".\n");

		// handle recursive creation of dependent exports if linker
		if(create_export_recursive(*slotPos, *expItem) != true)
		{
			UG_LOG("DataContainer::create_export_recursive: Error while creating recursive exports. Can not create Data Export.\n");
			return false;
		}
	}

	return true;
}


////////////////////////////////
// identify
////////////////////////////////

bool
DataContainer::
identify_exports()
{
	std::vector<DataExportItem*>::iterator iter, iterCopy;

	for(iter = m_vExportItemList.begin(); iter != m_vExportItemList.end(); ++iter)
	{
		UG_DLOG(LIB_DISC_LINKER, 3, "DataContainer::identify_exports: Checking for copies of: '" << (*iter)->name() << "' [" << find_nr(*iter) << "].\n");

		for(iterCopy = iter, ++iterCopy; iterCopy != m_vExportItemList.end(); )
		{
			UG_DLOG(LIB_DISC_LINKER, 3, "DataContainer::identify_exports:  -- comparing to : '" << (*iterCopy)->name() << "' [" << find_nr(*iterCopy) << "].\n");

			// if export items are equal
			if((*iter)->equal(*(*iterCopy)) == true)
			{
				// get export item and decrese iterator by 1 to start correctly in next search
				DataExportItem* expCopy = *iterCopy;
				UG_ASSERT(expCopy != NULL, "Export must exist, but does not. Internal error.");
				DataExportItem* exp = *iter;
				UG_ASSERT(exp != NULL, "Export must exist, but does not. Internal error.");

				UG_DLOG(LIB_DISC_LINKER, 3, "DataContainer::identify_exports:  ---- Equal exports detected: '" << exp->name() << "' [" << find_nr(exp) << "] == '"<< expCopy->name() << "' [" << find_nr(expCopy) << "]  .\n");

				// link all imports to the original
				for(size_t i = 0; i < expCopy->num_data_imports(); ++i)
				{
					// get data import i
					DataImportItem* imp = expCopy->get_data_import(i);
					UG_ASSERT(imp != NULL, "Import must exist, but does not. Internal error.");
					UG_DLOG(LIB_DISC_LINKER, 3, "DataContainer::identify_exports:  -------- Relink import nr " << i <<" '" << imp->name() << "'.\n");

					// remove data export
					if(imp->clear_data_export() != true)
					{
						UG_LOG("DataContainer::identify_exports:  Cannot clear Import. Aborting.\n");
						return false;
					}

					// add to original
					if(exp->add_data_import( imp ) != true)
					{
						UG_LOG("DataContainer::identify_exports:  Cannot add Export. Aborting.\n");
						return false;
					}
				}

				// remove export item from container and set iterator to next element
				iterCopy = m_vExportItemList.erase(iterCopy);
			}
			else
			{
				// if not equal continue with next element
				++iterCopy;
			}
		}
	}

	UG_DLOG(LIB_DISC_LINKER, 2, "DataContainer::identify_exports: Exports identified.\n");
	return true;
}

bool
DataContainer::
clear_identification()
{
	std::vector<DataExportItem*>::reverse_iterator iter;

	for(iter = m_vExportItemList.rbegin(); iter != m_vExportItemList.rend(); ++iter)
	{
		// get export
		DataExportItem* exp = *iter;
		UG_ASSERT(exp != NULL, "Export must exist, but does not. Internal error.");

		size_t num_imports;
		// if this export has more than on import registered
		if((num_imports = exp->num_data_imports()) > 1)
		{
			// loop over registered imports, except first one
			for(size_t i = num_imports - 1; i != 0; --i)
			{
				// get data import
				DataImportItem* imp = exp->get_data_import(i);
				UG_ASSERT(imp != NULL, "Import must exist, but does not. Internal error.");

				// get data export possibility item
				DataPossibilityItem* expPos = exp->get_possibility_item();
				UG_ASSERT(expPos != NULL, "Possibility must exist, but does not. Internal error.");

				// remove data export from import
				if(imp->clear_data_export() != true)
				{
					UG_LOG("DataContainer::clear_identification:  Cannot clear Export. Aborting.\n");
					return false;
				}

				DataExportItem* new_exp;
				// create new export
				if((new_exp = create_export(*expPos)) == NULL)
				{
					UG_LOG("DataContainer::clear_identification:  Cannot create new Export. Aborting.\n");
					return false;
				}

				// create new linkage
				if(link(*imp, *new_exp) != true)
				{
					UG_LOG("DataContainer::clear_identification:  Cannot link to new Export. Aborting.\n");
					return false;
				}
			}
		}
	}

	UG_DLOG(LIB_DISC_LINKER, 2, "DataContainer::identify_exports: Export identification cleared.\n");
	return true;
}

////////////////////////////////
// compute
////////////////////////////////

void
DataContainer::
compute(bool compute_derivatives)
{
	for(size_t i=0; i < m_vExportItemList.size(); ++i)
	{
		UG_DLOG(LIB_DISC_LINKER, 2, "DataContainer::compute: Computing Export '" << m_vExportItemList[i]->name() << "' [" << i << "] ...");
		m_vExportItemList[i]->compute(compute_derivatives);
		UG_DLOG(LIB_DISC_LINKER, 2, "done.\n");
	}
}

////////////////////////////////
// output
////////////////////////////////

bool
DataContainer::
print_export_possibilities(DataContainerInfoType type)
{
	// Export Possibilities
	UG_LOG("\n   Export Possibilities: " << m_vPossibilityItemList.size() << "\n");
	for(size_t i=0; i < m_vPossibilityItemList.size(); ++i)
	{
		DataPossibilityItem* pos = m_vPossibilityItemList[i];
		UG_LOG("   # [" << i << "]  <- " << pos->name() << " ");

		if(type == DCI_TYPES)
		{
			UG_LOG(",  [DataType: " << pos->data_type()->name());
			UG_LOG(",  PositionType: " << pos->position_type()->name() << "]");
		}

		if(type == DCI_LINKS)
		{
			std::string output = std::string("   # [i]  <- ") + pos->name() + std::string(" | -<");
			std::string offset;

			if(pos->num_slots() == 0)
			{
				UG_LOG("\n");
			}

			for(size_t slot = 0; slot < pos->num_slots(); ++slot)
			{
				if(slot==0) offset = " ";
				else offset = std::string(output.size() - 3, ' ');

				UG_LOG(offset << "| -<");

				if(pos->is_linked(slot) == true)
				{
					if(print_export_possibilities( pos->get_linked_possibility(slot), output) != true) return false;
				}
				else
				{
					UG_LOG("\n");
				}
			}
		}
	}
	return true;
}

bool
DataContainer::
print_export_possibilities(const DataPossibilityItem* pos, std::string output)
{
	std::string offset;

	UG_LOG(" <- " << pos->name() << " [ ] ");
	output = output + std::string(" <- ") + pos->name() + std::string(" [ ] ");

	if(pos->num_slots() == 0)
	{
		UG_LOG("\n");
		return true;
	}

	for(size_t i = 0; i < pos->num_slots(); ++i)
	{
		if(i==0) offset = " ";
		else offset = std::string(output.size() + 1, ' ');

		UG_LOG(offset << "| -<");

		if(pos->is_linked(i) == true)
		{
			if(print_export_possibilities(pos->get_linked_possibility(i), output) != true) return false;
		}
	}

	return true;
}

bool
DataContainer::
print_linkage()
{
	std::string output;

	UG_LOG("\n   Linkage: \n");
	for(size_t i=0; i < m_vImportItemList.size(); ++i)
	{
		DataImportItem* imp = m_vImportItemList[i];

		// print import
		UG_LOG("   # [" << i << "] " <<imp->name() << " -<");
		output = std::string("   # [i] ") + imp->name() + std::string(" -<");

		if(imp->is_linked() == true)
		{
			// print export
			const DataExportItem* exp = imp->get_data_export();

			if(print_linkage(exp, output) != true) return false;
		}
		else
		{
			UG_LOG("\n");
		}
	}
	UG_LOG("\n");
	return true;
}

bool
DataContainer::
print_linkage(const DataExportItem* exp, std::string output)
{
	std::string offset;

	UG_LOG(" <- " << exp->name() << " [" << find_nr(const_cast<DataExportItem*>(exp)) << "] ");
	output = output + std::string(" <- ") + exp->name() + std::string(" [") + std::string(" ") + std::string("] ");

	if(exp->num_slots() == 0)
	{
		UG_LOG("\n");
		return true;
	}

	for(size_t i = 0; i < exp->num_slots(); ++i)
	{
		if(i==0) offset = " ";
		else offset = std::string(output.size() + 1, ' ');

		UG_LOG(offset << "| -<");

		if(exp->is_linked(i) == true)
		{
			if(print_linkage(exp->get_data_export(i), output) != true) return false;
		}
	}

	return true;
}


bool
DataContainer::
print_imports(DataContainerInfoType type)
{
	UG_LOG("\n   Imports: " << m_vImportItemList.size() << "\n");
	for(size_t i=0; i < m_vImportItemList.size(); ++i)
	{
		UG_LOG("   # [" << i << "] " << m_vImportItemList[i]->name() << " -<   ");

		if(type == DCI_IPS)
		{
			UG_LOG(" [ at  ");
			if(m_vImportItemList[i]->print_positions() != true) return false;
			UG_LOG(" ]");
		}

		if(type == DCI_TYPES)
		{
			UG_LOG(",  [DataType: " << m_vImportItemList[i]->data_type()->name());
			UG_LOG(",  PositionType: " << m_vImportItemList[i]->position_type()->name() << "]");
		}

		if(type == DCI_VALUES)
		{
			UG_LOG("\n");
			if(m_vImportItemList[i]->print_info("             - ") != true) return false;
		}
		UG_LOG("\n");
	}
	return true;
}

bool
DataContainer::
print_exports(DataContainerInfoType type)
{
	UG_LOG("\n   Exports: " << m_vExportItemList.size() << "\n");
	for(size_t i=0; i < m_vExportItemList.size(); ++i)
	{
		UG_LOG("   # [" << i << "]  <- " << m_vExportItemList[i]->name() << "  ");

		if(type == DCI_IPS)
		{
			UG_LOG(" [ at  ");
			if(m_vExportItemList[i]->print_positions() != true) return false;
			UG_LOG(" ]");
		}

		if(type == DCI_TYPES)
		{
			UG_LOG(",  [DataType: " << m_vExportItemList[i]->data_type()->name());
			UG_LOG(",  PositionType: " << m_vExportItemList[i]->position_type()->name() << "]");
		}

		if(type == DCI_VALUES)
		{
			UG_LOG("\n");
			if(m_vExportItemList[i]->print_info("             - ") != true) return false;
		}

		UG_LOG("\n");
	}
	return true;
}

bool
DataContainer::
print_linker_imports(DataContainerInfoType type)
{
	if(type < 0 || type > 3) return false;

	// Imports
	UG_LOG("\n   Linker (Imports): " << m_vExportItemList.size() << "\n");
	for(size_t l=0; l < m_vExportItemList.size(); ++l)
	{
		UG_LOG("   # [" << l << "] " << m_vExportItemList[l]->name() << " [" << m_vExportItemList[l]->num_slots() <<" Slot(s)]");

		for(size_t i = 0; i < m_vExportItemList[l]->num_slots(); ++i)
		{
			UG_LOG("\n           | " <<  m_vExportItemList[l]->slot_name(i) << " -<   ");
		}
		UG_LOG("\n");
	}
	return true;
}


size_t
DataContainer::
find_nr(DataExportItem* exp)
{
	for(size_t i = 0; i < m_vExportItemList.size(); ++i)
	{
		if(m_vExportItemList[i] == exp) return i;
	}

	UG_ASSERT(0, "Export Item '" << exp->name() << "' not found in list.");
	return (size_t)-1;
}

}
