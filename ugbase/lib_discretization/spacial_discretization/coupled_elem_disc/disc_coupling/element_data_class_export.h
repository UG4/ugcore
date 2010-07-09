
#ifndef __H__LIB_DISCRETIZATION__ELEMENT_DATA_CLASS_EXPORT__
#define __H__LIB_DISCRETIZATION__ELEMENT_DATA_CLASS_EXPORT__

#include <typeinfo>
#include <string>
#include <vector>

#include "common/common.h"

#include "element_data_import_export.h"

#include "../../elem_disc/elem_disc_interface.h"

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

	protected:
		//typedef void (IElemDisc<TAlgebra>::*EvalFunction)(std::vector<data_type>&, std::vector<std::vector<data_type> >&, const std::vector<position_type>&, const local_vector_type&, bool);
		typedef void (IElemDisc<TAlgebra>::*EvalFunction)(std::vector<data_type>&, std::vector<std::vector<data_type> >&, const std::vector<position_type>&, const local_vector_type&, bool);

	public:
		DataClassExportPossibility(std::string name, EvalFunction func = NULL, IElemDisc<TAlgebra>* Class = NULL) :
			DataPossibilityItem(name, 0, &typeid(TDataType), &typeid(TPositionType)),
			m_sys(0), m_num_sh(0), m_evalFunction(func), m_u(NULL), m_ExportingClass(Class)
		{m_createdDataExports.clear();};

	public:
		// create a data export from this possibility
		virtual DataExportItem* create_data_export();

		// links a Possibility to slot 'slot'
		virtual bool link(DataPossibilityItem* posItem, std::size_t slot) {return false;};

	public:
		// set the eval function of this possibility and all exports created from this possibility
		bool set_eval_function(EvalFunction func, IElemDisc<TAlgebra>* Class);

		// set the number of unknowns this possibility depends on of this possibility and all exports created from this possibility
		bool set_num_sh(std::size_t sys, std::size_t num_sh);

		// set local solution and global indices
		bool set_local_solution(const local_vector_type& u);

		virtual ~DataClassExportPossibility();

	protected:
		std::size_t m_sys;
		std::size_t m_num_sh;

		EvalFunction m_evalFunction;
		const local_vector_type* m_u;
		IElemDisc<TAlgebra>* m_ExportingClass;
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
		typedef void (IElemDisc<TAlgebra>::*EvalFunction)(std::vector<data_type>&, std::vector<std::vector<data_type> >&, const std::vector<position_type>&, const local_vector_type&, bool);

	protected:
		// Only Data Possibility can create an instance
		DataClassExport(std::string name, DataPossibilityItem* possibility, EvalFunction func, IElemDisc<TAlgebra>* expClass) 	:
			DataExport<TDataType, TPositionType>(name, possibility),
			m_u(NULL), m_evalFunction(func), m_ExportingClass(expClass)
			{};

		// set evaluation function
		bool set_eval_function(EvalFunction func, IElemDisc<TAlgebra>* Class);

		// set number of unknowns, this export depends on
		bool set_num_sh(std::size_t sys, std::size_t num_sh);

		// set local solution and global indices (element local)
		bool set_local_solution(const local_vector_type& u);

	public:
		// compute
		virtual void compute(bool compute_derivatives);

	public:
		// number of data exports linked by this linker
		virtual std::size_t num_slots() const {return 0;}

		// add a Data Export number i
		virtual bool link(DataExportItem* exportItem, std::size_t slot) {return false;}

		// remove Data Export number i
		virtual bool clear_slot(std::size_t slot) {return false;}

		// name of slot
		virtual std::string slot_name(std::size_t slot) const {return "";};

		// get registered export of slot i
		virtual const DataExportItem* get_data_export(std::size_t slot) const {return NULL;}

		// return if an export is set at slot i
		virtual bool is_linked(std::size_t slot) const {return false;}

		// return if all exports are set
		virtual bool is_linked() const {return false;}

	protected:
		// current local solution (size: (0, ... , num_sh-1))
		const local_vector_type* m_u;

		// evaluation function of this export
		EvalFunction m_evalFunction;
		IElemDisc<TAlgebra>* m_ExportingClass;
};

}

#include "element_data_class_export_impl.h"

#endif

