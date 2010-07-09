
#ifndef __H__LIB_DISCRETIZATION__ELEMENT_DATA_CLASS_EXPORT__
#define __H__LIB_DISCRETIZATION__ELEMENT_DATA_CLASS_EXPORT__

#include <typeinfo>
#include <string>
#include <vector>

#include "common/common.h"

#include "element_data_import_export.h"

#include "../coupled_elem_disc_interface.h"

namespace ug{
/*
template <typename TDataType, typename TPositionType, typename TAlgebra>
class DataExportingClass
{
public:
	typedef TDataType data_type;
	typedef TPositionType position_type;
	typedef TAlgebra algebra_type;
	typedef typename TAlgebra::vector_type::local_vector_type local_vector_type;

	virtual void export1(std::vector<data_type>&, std::vector<std::vector<data_type> >&, const std::vector<position_type>&, const local_vector_type&, bool)
	{ UG_ASSERT(0, "When this function is called an export tries to use a function, that does not exist, falling back to this Standard implementation.");};

	virtual ~DataExportingClass(){};
};
*/


template <typename TDataType, typename TPositionType, typename TAlgebra>
class DataClassExportPossibility : public DataPossibilityItem
{
	public:
		typedef TDataType data_type;
		typedef TPositionType position_type;
		typedef TAlgebra algebra_type;
		typedef typename TAlgebra::vector_type::local_vector_type local_vector_type;

	public:
		DataClassExportPossibility(std::string name, ICoupledElemDisc<TAlgebra>* Class, size_t nrExport) :
			DataPossibilityItem(name, 0, &typeid(TDataType), &typeid(TPositionType)),
			m_sys(0), m_num_sh(0), m_u(NULL), m_nrExport(nrExport), m_pExportingClass(Class)
		{m_createdDataExports.clear();};

	public:
		// create a data export from this possibility
		virtual DataExportItem* create_data_export();

		// links a Possibility to slot 'slot'
		virtual bool link(DataPossibilityItem* posItem, size_t slot) {return false;};

	public:
		// set the number of unknowns this possibility depends on of this possibility and all exports created from this possibility
		bool set_num_sh(size_t sys, size_t num_sh);

		// set local solution and global indices
		bool set_local_solution(const local_vector_type& u);

		virtual ~DataClassExportPossibility();

	protected:
		size_t m_sys;
		size_t m_num_sh;

		const local_vector_type* m_u;
		size_t m_nrExport;
		ICoupledElemDisc<TAlgebra>* m_pExportingClass;
};

template <typename TDataType, typename TPositionType, typename TAlgebra>
class DataClassExport : public DataExport<TDataType, TPositionType>{
	friend class DataImport<TDataType,TPositionType>;
	friend class DataClassExportPossibility<TDataType,TPositionType,TAlgebra>;

	public:
		typedef TDataType data_type;
		typedef TPositionType position_type;
		typedef TAlgebra algebra_type;
		typedef typename TAlgebra::vector_type::local_vector_type local_vector_type;

	protected:
		// Only Data Possibility can create an instance
		DataClassExport(std::string name, DataPossibilityItem* possibility, ICoupledElemDisc<TAlgebra>* expClass, size_t nrExport) 	:
			DataExport<TDataType, TPositionType>(name, possibility),
			m_u(NULL), m_nrExport(nrExport), m_pExportingClass(expClass)
			{};

		// set number of unknowns, this export depends on
		bool set_num_sh(size_t sys, size_t num_sh);

		// set local solution and global indices (element local)
		bool set_local_solution(const local_vector_type& u);

	public:
		// compute
		virtual void compute(bool compute_derivatives);

	public:
		// number of data exports linked by this linker
		virtual size_t num_slots() const {return 0;}

		// add a Data Export number i
		virtual bool link(DataExportItem* exportItem, size_t slot) {return false;}

		// remove Data Export number i
		virtual bool clear_slot(size_t slot) {return false;}

		// name of slot
		virtual std::string slot_name(size_t slot) const {return "";};

		// get registered export of slot i
		virtual const DataExportItem* get_data_export(size_t slot) const {return NULL;}

		// return if an export is set at slot i
		virtual bool is_linked(size_t slot) const {return false;}

		// return if all exports are set
		virtual bool is_linked() const {return false;}

	protected:
		// current local solution (size: (0, ... , num_sh-1))
		const local_vector_type* m_u;

		// evaluation function of this export
		size_t m_nrExport;
		ICoupledElemDisc<TAlgebra>* m_pExportingClass;
};

}

#include "element_data_class_export_impl.h"

#endif

