
#ifndef __H__LIB_DISCRETIZATION__ELEMENT_DATA_CLASS_EXPORT__
#define __H__LIB_DISCRETIZATION__ELEMENT_DATA_CLASS_EXPORT__

#include <typeinfo>
#include <string>
#include <vector>

#include "common/common.h"

#include "import_export.h"

#include "../coupled_elem_disc_interface.h"

namespace ug{

// Export Possibility of a CoupledElemDisc.
// It has no slots.
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
			m_sysId(0), m_numSh(0), m_pSolution(NULL), m_nrExport(nrExport), m_pExportingClass(Class)
			{m_vCreatedDataExports.clear();};

	public:
		// create a data export from this possibility
		virtual DataExportItem* create_data_export();

	public:
		// set the number of unknowns this possibility depends on of this possibility and all exports created from this possibility
		bool set_num_sh(size_t num_sh);

		// set the system id
		bool set_sys_id(size_t sys_id);

		// set local solution and global indices
		bool set_local_solution(const local_vector_type& u);

		virtual ~DataClassExportPossibility();

	protected:
		size_t m_sysId;
		size_t m_numSh;
		const local_vector_type* m_pSolution;

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
			m_pSolution(NULL), m_nrExport(nrExport), m_pExportingClass(expClass)
			{};

		// set number of unknowns, this export depends on
		bool set_num_sh(size_t num_sh);

		// set sys_id, this export depends on
		bool set_sys_id(size_t sys_id);

		// set local solution and global indices (element local)
		bool set_local_solution(const local_vector_type& u);

	public:
		// compute
		virtual void compute(bool compute_derivatives);

	protected:
		// current local solution (size: (0, ... , num_sh-1))
		const local_vector_type* m_pSolution;

		// evaluation function of this export
		size_t m_nrExport;
		ICoupledElemDisc<TAlgebra>* m_pExportingClass;
};

}

#include "class_export_impl.h"

#endif

