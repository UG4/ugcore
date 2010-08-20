
#ifndef __H__LIB_DISCRETIZATION__ELEMENT_DATA_ITEMS__
#define __H__LIB_DISCRETIZATION__ELEMENT_DATA_ITEMS__

#include <typeinfo>
#include <vector>
#include <string>

#include "common/common.h"
#include "lib_discretization/common/local_algebra.h"

namespace ug{

// predeclaration
class DataPossibilityItem;
class DataImportItem;
class DataExportItem;

class DataItem {
	public:
		DataItem(std::string name, const std::type_info* dataType) :
			m_name(name), m_pDataType(dataType) {};

		/// name
		const std::string name() const {return m_name;};

		/// data type
		const std::type_info* data_type() const {return m_pDataType;};

		/// destructor
		virtual ~DataItem(){};

	private:
		std::string m_name;
		const std::type_info* m_pDataType;
};

// Factory to create Exports, a possibility can have slots that must be filled by other possbilities
class DataPossibilityItem : public DataItem {
	public:
		DataPossibilityItem(std::string name, size_t num_slots, const std::type_info* dataType) :
			DataItem(name, dataType), m_numSlots(num_slots)
			{m_vSlotPosItems.resize(m_numSlots, NULL); m_vLinkedPosItems.clear(); m_vCreatedDataExports.clear();};

	public:
		// create data export
		virtual DataExportItem* create_data_export() = 0;

		// delete reference to created export
		bool delete_data_export(DataExportItem* exportItem);

		// number of exports created by this factory
		size_t num_created_data_export_items() const {	return m_vCreatedDataExports.size();}

		// the i'th created export
		DataExportItem* get_created_data_export_item(size_t i) {return m_vCreatedDataExports[i];}
		const DataExportItem* get_created_data_export_item(size_t i) const {return m_vCreatedDataExports[i];}

	public:
		// number of slots this possibility needs to be linked to
		size_t num_slots() const {return m_numSlots;};

		// name of slot
		virtual std::string slot_name(size_t slot) const {return "";};

		// links a Possibility to slot 'slot'
		virtual bool link(DataPossibilityItem* posItem, size_t slot) {return false;};

		// clear slot 'slot'
		bool clear_slot(size_t slot);

		// returns if slot 'slot' is already linked
		bool is_linked(size_t slot) const
			{if(slot >= m_numSlots) return false; return (m_vSlotPosItems[slot] != NULL);};

		// returns if all slots are linked
		bool is_linked() const
			{for(size_t i = 0; i < num_slots(); ++i){if(!is_linked(i)) return false;}; return true;}

		// returns if this possibility is ready to create an export
		bool is_instantiable() const {return is_linked();};

		// returns the possibility linked to slot 'slot', NULL if not linked
		DataPossibilityItem* get_linked_possibility(size_t slot) const
			{if(slot >= m_numSlots) return NULL;return m_vSlotPosItems[slot];};

	protected:
		struct LinkedPossibility
		{
			DataPossibilityItem* item;
			size_t slot;
		};

		// remember all (possibilities+slots) this possibility is linked to as an export possibility
		void add_linked_possibility(LinkedPossibility& linkedPos);

		// remove this possibility from the list
		bool remove_linked_possibility(LinkedPossibility& linkedPos);

	protected:
		// created exports of this possibility
		std::vector<DataExportItem*> m_vCreatedDataExports;

		// number of slots of this possibility
		size_t m_numSlots;

		// linked export possibilities that are linked to the slots of this possibility
		std::vector<DataPossibilityItem*> m_vSlotPosItems;

		// possibilities, this export is linked to as an export
		std::vector<LinkedPossibility> m_vLinkedPosItems;
};

// A Data Export is linked to several imports/linkers providing the same values to all of them.
class DataExportItem : public DataItem {
	friend class DataImportItem;
	friend class DataContainer;

	public:
		DataExportItem(std::string name, const std::type_info* dataType, DataPossibilityItem* possibility) :
			DataItem(name, dataType),
			m_numSys(0), m_pPossibilityItem(possibility)
			{m_vImportList.clear(); m_vSysId.clear(); m_vNumFct.clear();};

		// number of system unknowns this export depends on
		inline size_t num_sys() const {return m_numSys;};

		// number of fct the data depends on
		inline size_t num_fct(size_t loc_sys) const {
			UG_ASSERT(loc_sys < m_numSys, "Accessing system, that is not registered, sys = " << loc_sys << ", num_sys = " << m_numSys);
			UG_ASSERT(loc_sys < m_vSysId.size(), "Accessing system, that is not registered, sys = " << loc_sys << ", m_vSysId.size() = " <<  m_vSysId.size());
			return m_vNumFct[loc_sys];};

		// number of unknowns per fct the data depends on
		inline size_t num_dofs(size_t fct, size_t loc_sys) const {
			UG_ASSERT(loc_sys < m_numSys, "Accessing system, that is not registered, sys = " << loc_sys << ", num_sys = " << m_numSys);
			UG_ASSERT(loc_sys < m_vSysId.size(), "Accessing system, that is not registered, sys = " << loc_sys << ", m_vSysId.size() = " <<  m_vSysId.size());
			return m_vvDoFsPerFct[loc_sys][fct];};

		// number of unknowns per system the data depends on
		inline size_t num_dofs(size_t loc_sys) const {
			size_t sum = 0;
			for(size_t fct = 0; fct < num_fct(loc_sys); ++fct) sum += num_dofs(fct, loc_sys);
			return sum;
		}

		// number of unknowns this export depends on
		inline size_t num_dofs() const {
			size_t sum = 0;
			for(size_t s = 0; s < num_sys(); ++s) sum += num_dofs(s);
			return sum;
		}

		// id of system this exports depends on
		inline size_t sys_id(size_t loc_sys) const{
			UG_ASSERT(loc_sys < m_numSys, "Accessing system, that is not registered, sys = " << loc_sys << ", num_sys = " << m_numSys);
			UG_ASSERT(loc_sys < m_vSysId.size(), "Accessing system, that is not registered, sys = " << loc_sys << ", m_vSysId.size() = " <<  m_vSysId.size());
			return m_vSysId[loc_sys];}

		// returns false if export does not depend on sys_id
		// returns true if export depends on sys_id and the loc_sys corresponding to sys_id
		inline bool depends_on_sys(size_t sys_id, size_t& loc_sys) const
		{
			for(size_t k = 0; k < num_sys(); ++k)
			{
				if(this->sys_id(k) == sys_id)
				{
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
		size_t num_data_imports() const {return m_vImportList.size();};

		// get registered import
		DataImportItem* get_data_import(size_t i) {return m_vImportList[i];};
		const DataImportItem* get_data_import(size_t i) const {return m_vImportList[i];};

	protected:
		// get underlying data export possibility
		DataPossibilityItem* get_possibility_item() {return m_pPossibilityItem;}
		const DataPossibilityItem* get_possibility_item() const {return m_pPossibilityItem;}

	// Slot side (an export can itself depend on exports. Thus it has internally some Imports, called slots)
	// As default an export item does not have any slots. In case a derived class has, this functions have to be overwritten.
	public:
		// number of data exports linked by this linker
		virtual size_t num_slots() const {return 0;}

		// name of slot
		virtual std::string slot_name(size_t slot) const {return "";};

		// add a Data Export number i
		virtual bool link(DataExportItem* exportItem, size_t slot) {return false;}

		// remove Data Export number i
		virtual bool clear_slot(size_t slot) {return false;}

		// get registered export of slot i
		virtual const DataExportItem* get_data_export(size_t slot) const {return NULL;}

		// return if an export is set at slot i
		virtual bool is_linked(size_t slot) const {return false;}

		// return if all exports are set
		virtual bool is_linked() const {return false;}

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
		// sets the number of systems the export depends on
		virtual bool set_num_sys(size_t num_sys)
		{
			m_numSys = num_sys;
			m_vSysId.resize(m_numSys);
			m_vNumFct.resize(m_numSys, 0);
			m_vvDoFsPerFct.resize(m_numSys);
			return true;
		}

		// sets the system id for the local system
		virtual bool set_sys_id(size_t sys_id, size_t loc_sys = 0)
		{
			if(loc_sys >= m_numSys) return false;
			m_vSysId[loc_sys] = sys_id;
			return true;
		}

		// sets the number of shape functions for the local system
		virtual bool set_num_fct(size_t num_fct, size_t loc_sys = 0)
		{
			if(loc_sys >= m_numSys)
				{UG_LOG("Local system does not exist for this export.\n"); return false;}
			m_vNumFct[loc_sys] = num_fct;
			m_vvDoFsPerFct[loc_sys].resize(num_fct, 0);
			return true;
		}

		// sets the number of dofs for a functions for the local system
		virtual bool set_num_dofs(size_t fct, size_t num_dofs, size_t loc_sys = 0)
		{
			if(loc_sys >= m_numSys)
				{UG_LOG("Local system does not exist for this export.\n"); return false;}
			if(fct >= m_vNumFct[loc_sys])
				{UG_LOG("Function does not exist for this export and system "<<loc_sys<< ".\n"); return false;}

			m_vvDoFsPerFct[loc_sys][fct] = num_dofs;
			return true;
		}

	protected:
		// number of system this export depends on
		size_t m_numSys;

		// ids of systems this export depends on
		std::vector<size_t> m_vSysId;

		// number of unknowns the export depends on
		std::vector<size_t> m_vNumFct;

		// number of dofs the export depends on for each function
		std::vector<std::vector<size_t> > m_vvDoFsPerFct;

		// Imports linked to this export
		std::vector<DataImportItem*> m_vImportList;

		// underlying possibility this export is created from
		DataPossibilityItem* m_pPossibilityItem;
};

// A Data Import is linked to exactly on data export.
class DataImportItem : public DataItem{
	friend class DataExportItem;

	public:
		DataImportItem(std::string name, const std::type_info* dataType) :
			DataItem(name, dataType),
			m_numFct(0), m_pExport(NULL)
			{m_vDoFsPerFct.clear();};

		// number of systems the data depends on (num_sys)
		inline size_t num_sys() const {
			UG_ASSERT(m_pExport != NULL, "No Export set.");
			return m_pExport->num_sys();};

		// id of system this import depends on
		inline size_t sys_id(size_t loc_sys) const {
			UG_ASSERT(m_pExport != NULL, "No Export set.");
			return m_pExport->sys_id(loc_sys);};

		// number of functions the data depends on per system
		inline size_t num_fct(size_t loc_sys) const {
			UG_ASSERT(m_pExport != NULL, "No Export set.");
			return m_pExport->num_fct(loc_sys);};

		// number of functions the data depends on for all systems
		inline size_t num_fct() const {
			UG_ASSERT(m_pExport != NULL, "No Export set.");
			size_t t = 0;
			for(size_t s = 0; s < num_sys(); ++s) t += num_fct(s);
			return t;};

		// number of unknowns the data depends on per function and system
		inline size_t num_dofs(size_t fct, size_t loc_sys) const {
			UG_ASSERT(m_pExport != NULL, "No Export set.");
			return m_pExport->num_dofs(fct, loc_sys);};

		// number of unknowns the data depends on for all systems
		inline size_t num_dofs(size_t fct) const {
			UG_ASSERT(m_pExport != NULL, "No Export set.");
			size_t t = 0;
			for(size_t s = 0; s < num_sys(); ++s) t += num_dofs(fct, s);
			return t;};

		// number of fct the data depends on
		inline size_t num_eq_fct() const {return m_numFct;};

		// number of unknowns per fct the data depends on
		inline size_t num_eq_dofs(size_t fct) const {
			UG_ASSERT(fct < num_eq_fct(), "Invalid index.");
			return m_vDoFsPerFct[fct];};

		// returns false if import does not depend on sys_id
		// returns true if import depends on sys_id and the loc_sys corresponding to sys_id
		inline bool depends_on_sys(size_t sys_id, size_t& loc_sys) const {
			UG_ASSERT(m_pExport != NULL, "No Export set.");
			return m_pExport->depends_on_sys(sys_id, loc_sys);
		}

		// sets the number of shape functions for the local system
		virtual bool set_num_eq_fct(size_t num_fct)
		{
			m_numFct = num_fct;
			m_vDoFsPerFct.resize(m_numFct);
			return true;
		}

		// sets the number of dofs for a functions for the local system
		virtual bool set_num_eq_dofs(size_t fct, size_t num_dofs)
		{
			if(fct >= m_numFct)
				{UG_LOG("Function does not exist for this import.\n"); return false;}

			m_vDoFsPerFct[fct] = num_dofs;
			return true;
		}

		// add of diagonal couplings to local Matrix J for coupling with system loc_sys
		virtual bool add_offdiagonal(LocalMatrixBase& J, size_t loc_sys, number s_a) = 0;

	public:
		// link this import to an export
		virtual bool link_data_export(DataExportItem* exportItem) = 0;

		// remove export
		virtual bool clear_data_export() = 0;

		// returns if export is set
		bool is_linked() const {return m_pExport != NULL;};

		// return export
		const DataExportItem* get_data_export() const {return m_pExport;};

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
		// number of functions for equations
		size_t m_numFct;

		// number of dofs per function for equations
		std::vector<size_t> m_vDoFsPerFct;

		// Data Export, to which this import is linked (may be a Data Linker)
		DataExportItem* m_pExport;
};


} // namespace ug
#endif /*__H__LIB_DISCRETIZATION__ELEMENT_DATA_ITEMS__*/
