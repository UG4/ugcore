

#ifndef __H__LIB_DISCRETIZATION__ELEMENT_DATA_CONTAINER__
#define __H__LIB_DISCRETIZATION__ELEMENT_DATA_CONTAINER__

#include <vector>
#include "common/common.h"

#include "lib_discretization/domain_discretization/disc_coupling/element_data_items.h"


namespace ug{

enum DataContainerInfoType {
	DCI_NONE = 0,
	DCI_IPS,
	DCI_TYPES,
	DCI_VALUES,
	DCI_LINKS,
	NUM_DCI
};

class DataContainer {
	public:
		// import items (registered by user)
		bool register_item(DataImportItem& ImportItem);

		// export possibilites or linker possiblity items (registered by user)
		bool register_item(DataPossibilityItem& PossibilityItem);

	//protected:
		// linker items (registered by container)
		bool register_item(DataExportItem& ExportItem);

	public:
		// unregister import
		bool unregister_item(DataImportItem& ImportItem);

		// unregister export possibility or linker possibility
		bool unregister_item(DataPossibilityItem& PossibilityItem);

	//protected:
		// unregister (unregistered by container)
		bool unregister_item(DataExportItem& ExportItem);

	public:
		bool link_possibility(DataPossibilityItem& posItem1, std::size_t slot, DataPossibilityItem& posItem2);
		bool link_possibility(std::size_t nr_pos1, std::size_t slot, std::size_t nr_pos2);

		// create an instance of the Data Possibility (called DataExport)
		DataExportItem* create_export(DataPossibilityItem& PossibilityItem);
		DataExportItem* create_export(std::size_t nr_pos);
		bool create_export_recursive(DataPossibilityItem& PossibilityItem, DataExportItem& expItem);

		// linking
		bool link(DataImportItem& ImportItem, DataExportItem& ExportItem);
		bool link(DataExportItem& Export1Item, std::size_t slot, DataExportItem& Export2Item);

		// linking by number
		bool link(std::size_t nr_imp, std::size_t nr_exp);
		bool link(std::size_t nr_exp1, std::size_t slot, std::size_t nr_exp2);

		// identify exports: Those from same ExportPossibility that have the same positions
		// Afterwards only one export of each type exist and all imports are linked to this one,
		// that have previously been linked to one of the copies
		bool identify_exports();

		// resets the identification, i.e. afterwards every import is linked to on export
		// and number of exports == number of imports
		// ATTENTION: Values are not computed after wards. A compute() must be called in order to rely on valid values
		bool clear_identification();

		// compute values of all exports
		// TODO: should be return a bool ?
		void compute(bool compute_derivatives);

		// print informations
		bool print_export_possibilities(DataContainerInfoType type);
		bool print_linkage();

	protected:
		bool print_export_possibilities(const DataPossibilityItem* exp, std::string output);
		bool print_linkage(const DataExportItem* exp, std::string output);

	public:
		bool print_imports(DataContainerInfoType type = DCI_NONE);
		bool print_exports(DataContainerInfoType type = DCI_NONE);
		bool print_linker_imports(DataContainerInfoType type = DCI_NONE);

		// destructor
		~DataContainer();

	protected:
		uint find_nr(DataExportItem* exp);

	protected:
		std::vector<DataPossibilityItem*> m_possibilityItemList;

		std::vector<DataImportItem*> m_importItemList;

		std::vector<DataExportItem*> m_exportItemList;
};

}

#endif /*__H__LIB_DISCRETIZATION__ELEMENT_DATA_CONTAINER__*/
