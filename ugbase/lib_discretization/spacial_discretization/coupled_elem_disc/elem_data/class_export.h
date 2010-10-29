
#ifndef __H__LIB_DISCRETIZATION__ELEMENT_DATA_CLASS_EXPORT__
#define __H__LIB_DISCRETIZATION__ELEMENT_DATA_CLASS_EXPORT__

#include <typeinfo>
#include <string>
#include <vector>

#include "common/common.h"

#include "import_export.h"

#include "../coupled_elem_disc_interface.h"
#include "lib_discretization/common/local_algebra.h"

namespace ug{

// Export Possibility of a CoupledElemDisc.
// It has no slots.
template <typename TDataType, typename TAlgebra>
class DataClassExportPossibility : public DataPossibilityItem
{
	public:
		typedef TDataType data_type;
		typedef TAlgebra algebra_type;
		typedef LocalVector<typename TAlgebra::vector_type::value_type> local_vector_type;

	public:
		DataClassExportPossibility(std::string name, ICoupledElemDisc<TAlgebra>* Class, size_t nrExport) :
			DataPossibilityItem(name, 0, &typeid(TDataType)),
			m_sysId(0), m_numFct(0), m_pSolution(NULL), m_nrExport(nrExport), m_pExportingClass(Class)
			{m_vCreatedDataExports.clear(); m_vNumDoFsPerFct.clear();};

	public:
		// create a data export from this possibility
		virtual DataExportItem* create_data_export();

	public:
		// set the number of unknowns this possibility depends on of this possibility and all exports created from this possibility
		bool set_num_fct(size_t num_fct);

		// set number of dofs per fct
		bool set_num_dofs(size_t fct, size_t num_dofs);

		// set the system id
		bool set_sys_id(size_t sys_id);

		// set local solution and global indices
		bool set_local_solution(const local_vector_type& u);

		virtual ~DataClassExportPossibility();

	protected:
		size_t m_sysId;
		size_t m_numFct;
		std::vector<size_t> m_vNumDoFsPerFct;
		const local_vector_type* m_pSolution;

		size_t m_nrExport;
		ICoupledElemDisc<TAlgebra>* m_pExportingClass;
};

template <typename TDataType, typename TAlgebra>
class DataClassExport : public DataExport<TDataType>{
	friend class DataImport<TDataType>;
	friend class DataClassExportPossibility<TDataType,TAlgebra>;

	public:
		typedef TDataType data_type;
		typedef TAlgebra algebra_type;
		typedef LocalVector<typename TAlgebra::vector_type::value_type> local_vector_type;

	protected:
		// Only Data Possibility can create an instance
		DataClassExport(std::string name, DataPossibilityItem* possibility, ICoupledElemDisc<TAlgebra>* expClass, size_t nrExport) 	:
			DataExport<TDataType>(name, possibility),
			m_pSolution(NULL), m_nrExport(nrExport), m_pExportingClass(expClass)
			{DataExport<TDataType>::set_num_sys(1);};

		// set the number of unknowns this possibility depends on of this possibility and all exports created from this possibility
		bool set_num_fct(size_t num_fct);

		// set number of dofs per fct
		bool set_num_dofs(size_t fct, size_t num_dofs);

		// set sys_id, this export depends on
		bool set_sys_id(size_t sys_id);

		// set local solution and global indices (element local)
		bool set_local_solution(const local_vector_type& u);

	public:
		// compute
		virtual void compute(bool compute_derivatives);

	protected:
		// current local solution (size: (0, ... , num_sh-1))
		const SubFunctionMap* m_pSubFunctionMap;
		local_vector_type* m_pSolution;

		// evaluation function of this export
		size_t m_nrExport;
		ICoupledElemDisc<TAlgebra>* m_pExportingClass;
};

}

#include "class_export_impl.h"

#endif

