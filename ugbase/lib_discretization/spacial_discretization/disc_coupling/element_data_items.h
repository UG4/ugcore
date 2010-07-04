
#ifndef __H__LIB_DISCRETIZATION__ELEMENT_DATA_ITEMS__
#define __H__LIB_DISCRETIZATION__ELEMENT_DATA_ITEMS__

#include <typeinfo>
#include <vector>
#include <string>

#include "common/common.h"
#include "lib_algebra/local_matrix_vector/flex_local_matrix_vector.h"

namespace ug{

// predeclaration
class DataPossibilityItem;
class DataImportItem;
class DataExportItem;

class DataItem {
	public:
		DataItem(std::string name, const std::type_info* dataType, const std::type_info* posType) :
			m_name(name), m_dataType(dataType), m_positionType(posType) {};

		const std::string name() const {return m_name;};
		const std::type_info* data_type() const {return m_dataType;};
		const std::type_info* position_type() const {return m_positionType;};

		virtual ~DataItem(){};

	private:
		std::string m_name;
		const std::type_info* m_dataType;
		const std::type_info* m_positionType;
};

// Factory to create Exports
class DataPossibilityItem : public DataItem {
	public:
		DataPossibilityItem(std::string name, std::size_t num_slots, const std::type_info* dataType, const std::type_info* posType) :
			DataItem(name, dataType, posType), m_num_slots(num_slots)
			{m_slotPosItems.resize(m_num_slots); m_linkedPosItems.clear(); m_createdDataExports.clear();};

	public:
		// links a Possibility to slot 'slot'
		virtual bool link(DataPossibilityItem* posItem, std::size_t slot) = 0;

	public:
		// create data export
		virtual DataExportItem* create_data_export() = 0;

		// delete reference to created export
		bool delete_data_export(DataExportItem* exportItem);

		// number of exports created by this factory
		std::size_t num_created_data_export_items() const {	return m_createdDataExports.size();}

		// the i'th created export
		DataExportItem* get_created_data_export_item(std::size_t i) const {return m_createdDataExports[i];}

	public:
		// number of slots this possibility needs to be linked to
		std::size_t num_slots() const {return m_num_slots;};

		// returns if slot 'slot' is already linked
		bool is_linked(std::size_t slot) const {
			UG_ASSERT(slot < m_num_slots, "Slot does not exist.");
			return (m_slotPosItems[slot] != NULL);};

		// returns if all slots are linked
		bool is_linked() const {
			for(std::size_t i = 0; i < num_slots(); ++i){if(is_linked(i) == false) return false;};
			return true;}

		// returns if this possibility is ready to create an export
		bool is_instantiable() const {return is_linked();};

		// clear slot 'slot'
		bool clear_slot(std::size_t slot);

		// returns the possibility linked to slot 'slot', NULL if not linked
		DataPossibilityItem* get_linked_possibility(std::size_t slot) const {
			UG_ASSERT(slot < m_num_slots, "Slot does not exist.");
			return m_slotPosItems[slot];};

	public:
		struct LinkedPossibility
		{
			DataPossibilityItem* item;
			std::size_t slot;
		};

		void add_linked_possibility(LinkedPossibility& linkedPos);
		bool remove_linked_possibility(LinkedPossibility& linkedPos);

	protected:
		// number of slots of this possibility
		std::size_t m_num_slots;

		// linked possibilities to the slots
		std::vector<DataPossibilityItem*> m_slotPosItems;

		// possibilities, this export is linked to as an export
		std::vector<LinkedPossibility> m_linkedPosItems;

		// created exports of this possibility
		std::vector<DataExportItem*> m_createdDataExports;
};

// A Data Export is linked to several imports/linkers providing the same values to all of them.
class DataExportItem : public DataItem {
	friend class DataImportItem;
	friend class DataContainer;

	public:
		DataExportItem(std::string name, const std::type_info* dataType, const std::type_info* posType, DataPossibilityItem* possibility) :
			DataItem(name, dataType, posType),
			m_num_sys(0), m_possibilityItem(possibility)
			{m_importList.clear(); m_sys.clear(); m_num_sh.clear();};

		// number of system unknowns this export depends on
		inline std::size_t num_sys() const {return m_num_sys;};

		// number of unknowns the data depends on
		inline std::size_t num_sh(std::size_t loc_sys) const {
			UG_ASSERT(loc_sys < m_num_sys, "Accessing system, that is not registered, sys = " << loc_sys << ", num_sys = " << m_num_sys);
			UG_ASSERT(loc_sys < m_sys.size(), "Accessing system, that is not registered, sys = " << loc_sys << ", m_sys.size() = " <<  m_sys.size());
			return m_num_sh[loc_sys];};

		// id of system this exports depends on
		inline std::size_t sys(std::size_t loc_sys) const{
			UG_ASSERT(loc_sys < m_num_sys, "Accessing system, that is not registered, sys = " << loc_sys << ", num_sys = " << m_num_sys);
			UG_ASSERT(loc_sys < m_sys.size(), "Accessing system, that is not registered, sys = " << loc_sys << ", m_sys.size() = " <<  m_sys.size());
			return m_sys[loc_sys];}

		inline bool depends_on_sys(std::size_t glob_sys, std::size_t& loc_sys) const
		{
			for(std::size_t k = 0; k < num_sys(); ++k)
			{
				if(sys(k) == glob_sys)
				{
					UG_DLOG(LIB_DISC_LINKER, 2, "DataExport::derivatives: Depending on system " << glob_sys << ".\n");
					loc_sys = k;
					return true;
				}
			}
			return false;
		}

	public:
		// compute the values
		virtual void compute(bool compute_derivatives) = 0;


	// Export side
	public:
		// add data import
		virtual bool add_data_import(DataImportItem* importItem) = 0;

		// remove data import
		virtual bool remove_data_import(DataImportItem* importItem) = 0;

		// number of registered imports
		std::size_t num_data_imports() const {return m_importList.size();};

		// get registered import
		DataImportItem* get_data_import(std::size_t i) {return m_importList[i];};
		const DataImportItem* get_data_import(std::size_t i) const {return m_importList[i];};

	protected:
		// get underlying data export possibility
		DataPossibilityItem* get_possibility_item() {return m_possibilityItem;}
		const DataPossibilityItem* get_possibility_item() const {return m_possibilityItem;}

	// Import (Linker) side
	public:
		// number of data exports linked by this linker
		virtual std::size_t num_slots() const = 0;

		// name of import i
		virtual std::string slot_name(std::size_t slot) const = 0;

		// add a Data Export number i
		virtual bool link(DataExportItem* exportItem, std::size_t slot) = 0;

		// remove Data Export number i
		virtual bool clear_slot(std::size_t slot) = 0;

		// get registered export of slot i
		virtual const DataExportItem* get_data_export(std::size_t slot) const = 0;

		// return if an export is set at slot i
		virtual bool is_linked(std::size_t slot) const = 0;

		// return if all exports are set
		virtual bool is_linked() const = 0;

	// informations
	public:
		// returns, if two export items are equal
		virtual bool equal(const DataExportItem& v) const = 0;

		// print positions
		virtual bool print_positions() const = 0;

		// print values
		virtual bool print_values() const = 0;

		// print derivatives
		virtual bool print_derivatives(std::string offset) const = 0;

		// print general informations
		virtual bool print_info(std::string offset) const = 0;

	protected:
		// number of system this export depends on
		std::size_t m_num_sys;

		// id of system this export depends on
		std::vector<std::size_t> m_sys;

		// number of unknowns the export depends on
		std::vector<std::size_t> m_num_sh;

		// Imports linked to this export
		std::vector<DataImportItem*> m_importList;

		// underlying
		DataPossibilityItem* m_possibilityItem;
};

// A Data Import is linked to exactly on data export.
class DataImportItem : public DataItem{
	friend class DataExportItem;

	public:
		DataImportItem(std::string name, const std::type_info* dataType, const std::type_info* posType) :
			DataItem(name, dataType, posType),
			m_export(NULL)
			{};

		// number of systems the data depends on (num_sys)
		inline std::size_t num_sys() const {
			UG_ASSERT(m_export != NULL, "No Export set.");
			return m_export->num_sys();};

		// number of systems the data depends on (num_sys)
		inline std::size_t sys(std::size_t s) const {
			UG_ASSERT(m_export != NULL, "No Export set.");
			return m_export->sys(s);};

		// number of unknowns the data depends on per system (num_sh)
		inline std::size_t num_sh(std::size_t s) const {
			UG_ASSERT(m_export != NULL, "No Export set.");
			return m_export->num_sh(s);};

		// number of unknowns the data depends on for all systems (num_sh)
		inline std::size_t num_sh() const {
			UG_ASSERT(m_export != NULL, "No Export set.");
			std::size_t t = 0;
			for(std::size_t s = 0; s < num_sys(); ++s) t += num_sh(s);
			return t;};

		inline bool depends_on_sys(std::size_t glob_sys, std::size_t& loc_sys) const {
			UG_ASSERT(m_export != NULL, "No Export set.");
			return m_export->depends_on_sys(glob_sys, loc_sys);
		}

		virtual bool add_offdiagonal(FlexLocalMatrix& J, size_t s, number s_a) = 0;

	public:
		// link this import to an export
		virtual bool link_data_export(DataExportItem* exportItem) = 0;

		// remove export
		virtual bool clear_data_export() = 0;

		// returns if export is set
		bool is_linked() const {return m_export != NULL;};

		// return export
		const DataExportItem* get_data_export() const {return m_export;};

	public:
		// print positions
		virtual bool print_positions() const = 0;

		// print values
		virtual bool print_values() const = 0;

		// print derivatives
		virtual bool print_derivatives(std::string offset) const = 0;

		// print general informations
		virtual bool print_info(std::string offset) const = 0;

	protected:
		// Data Export, to which this import is linked (may be a Data Linker)
		DataExportItem* m_export;
};


} // namespace ug
#endif /*__H__LIB_DISCRETIZATION__ELEMENT_DATA_ITEMS__*/
