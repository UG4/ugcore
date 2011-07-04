/*
 * data_evaluator.h
 *
 *  Created on: 17.12.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__DATA_EVALUATOR__
#define __H__UG__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__DATA_EVALUATOR__

#include "user_data.h"
#include "../elem_disc/elem_disc_interface.h"
#include "data_import_export.h"

namespace ug{

class DataEvaluator
{
	public:
	///	Type of local vector
		typedef IElemDisc::local_vector_type local_vector_type;

	///	Type of local matrix
		typedef IElemDisc::local_matrix_type local_matrix_type;

	///	Type of local indices
		typedef IElemDisc::local_index_type local_index_type;


	public:
	///	sets the elem discs to evaluate
		bool set_elem_discs(const std::vector<IElemDisc*>& vElemDisc,
		                    const FunctionPattern& fctPat,
		                    bool bNonRegularGrid);

	///	Mapping between local functions of ElemDisc i and common FunctionGroup
		const FunctionIndexMapping& map(size_t i) const
		{
			UG_ASSERT(i < m_vMap.size(), "Wrong index");
			return m_vMap[i];
		}

	///	Function group of all needed functions
		const FunctionGroup& fct_group() const {return m_commonFctGroup;}

	///	clears imports and ip data
		void clear();

	//	this function tries to add the last entry of vTryingToAdd to the eval
	//	data
		bool add_data_to_eval_data(std::vector<IIPData*>& vEvalData,
		                           std::vector<IIPData*>& vTryingToAdd);

		bool extract_imports_and_ipdata();

		bool use_hanging() const {return m_bUseHanging;}

		////////////////////////////////////////////
		// Regular assembling
		////////////////////////////////////////////

		bool set_time_dependent(bool bTimeDep, number time = 0.0,
		                        LocalVectorTimeSeries* locTimeSeries = NULL);

		bool set_non_regular_grid(bool bNonRegularGrid);

		template <typename TElem>
		bool prepare_elem_loop(local_index_type& ind, number time = 0.0);

		template <typename TElem>
		bool prepare_element(TElem* elem,
		                     local_vector_type& u,
		                     const local_index_type& ind);

		bool compute_elem_data(local_vector_type & u,
		                       local_index_type& ind,
		                       bool computeDeriv = false);

		bool assemble_JA(local_matrix_type& A, local_vector_type& u);
		bool assemble_JM(local_matrix_type& M, local_vector_type& u);

		bool assemble_A(local_vector_type& d, local_vector_type& u);
		bool assemble_M(local_vector_type& d, local_vector_type& u);

		bool assemble_rhs(local_vector_type& rhs);

		bool finish_element_loop();

		////////////////////////////////////////////
		// Coupling
		////////////////////////////////////////////

		bool compute_lin_defect_JA(local_vector_type& u, local_index_type& ind);

		bool compute_lin_defect_JM(local_vector_type& u, local_index_type& ind);

		bool add_coupl_JA(local_matrix_type& J, local_index_type& indRow);

		bool add_coupl_JM(local_matrix_type& J, local_index_type& indRow);

	protected:
	//	common functions
		FunctionGroup m_commonFctGroup;

	//	Function mapping for each disc
		std::vector<FunctionIndexMapping> m_vMap;

	//	current elem discs
		const std::vector<IElemDisc*>* m_pvElemDisc;

	//	data imports which are connected to non-zero derivative ip data
		std::vector<IDataImport*> m_vIDataImport;

	//	Function mapping for import and for the connected data
		std::vector<FunctionIndexMapping> m_vMapImp;
		std::vector<FunctionIndexMapping> m_vMapImpConn;

	//	constant data
		std::vector<IIPData*> m_vConstData;

	//	position dependend data
		std::vector<IIPData*> m_vPosData;

	//	all dependent data
		std::vector<IIPData*> m_vDependData;
		std::vector<FunctionIndexMapping> m_vMapDepend;

	//	data linker
		std::vector<IIPData*> m_vLinkerData;
		std::vector<IDependentIPData*> m_vLinkerDepend;
		std::vector<FunctionIndexMapping> m_vMapLinker;

	//	exports
		std::vector<IDataExport*> m_vIDataExport;
		std::vector<FunctionIndexMapping> m_vMapExp;

		bool m_bUseHanging;
};

} // end namespace ug

#include "data_evaluator_impl.h"

#endif /* __H__UG__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__DATA_EVALUATOR__ */
