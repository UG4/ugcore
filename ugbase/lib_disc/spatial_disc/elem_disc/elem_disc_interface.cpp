/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#include "elem_disc_interface.h"
#include "lib_disc/common/groups_util.h"

namespace ug{

template <typename TDomain>
IElemDiscBase<TDomain>::IElemDiscBase(const char* functions, const char* subsets)
	:	m_spApproxSpace(nullptr), m_spFctPattern(nullptr),
	  	m_timePoint(0), m_pLocalVectorTimeSeries(nullptr), m_bStationaryForced(false)
		//,m_id(ROID_UNKNOWN)
{
	if(functions == nullptr) functions = "";
	if(subsets == nullptr) subsets = "";
	set_functions(functions);
	set_subsets(subsets);
// 	set_default_add_fct();  // OBSOLETE: now part of base-type constructor
}


template <typename TDomain>
IElemDiscBase<TDomain>::
IElemDiscBase(const std::vector<std::string>& vFct,
                              const std::vector<std::string>& vSubset)
	: 	m_spApproxSpace(nullptr), m_spFctPattern(nullptr),
		m_timePoint(0), m_pLocalVectorTimeSeries(nullptr), m_bStationaryForced(false)
		//,m_id(ROID_UNKNOWN)
{
	set_functions(vFct);
	set_subsets(vSubset);
// 	set_default_add_fct(); // OBSOLETE: now part of base-type constructor
}


template <typename TDomain>
void IElemDiscBase<TDomain>::
set_approximation_space(SmartPtr<ApproximationSpace<TDomain> > approxSpace)
{
//	check whether the approximation space has already been set
	bool newApproxSpace = (m_spApproxSpace != approxSpace);

//	remember approx space
	m_spApproxSpace = approxSpace;

//	set function pattern
	set_function_pattern(approxSpace->dof_distribution_info());

//	invoke callback
	if(newApproxSpace)
		approximation_space_changed();
}

template <typename TLeaf, typename TDomain>
void IElemAssembleFuncs<TLeaf, TDomain>::clear_add_fct(ReferenceObjectID id)
{
	m_vPrepareTimestepElemFct[id] = nullptr;
	m_vFinishTimestepElemFct[id] = nullptr;

	m_vPrepareElemLoopFct[id] = nullptr;
	m_vPrepareElemFct[id] = nullptr;
	m_vFinishElemLoopFct[id] = nullptr;

	m_vElemJAFct[id] = nullptr;
	m_vElemJMFct[id] = nullptr;

	m_vElemdAFct[id] = nullptr;
	m_vElemdAExplFct[id] = nullptr;
	m_vElemdMFct[id] = nullptr;

	m_vElemRHSFct[id] = nullptr;
}


template <typename TLeaf, typename TDomain>
void IElemEstimatorFuncs<TLeaf, TDomain>::clear_add_fct(ReferenceObjectID id)
{
	m_vPrepareErrEstElemLoopFct[id] = nullptr;
	m_vPrepareErrEstElemFct[id] = nullptr;
	m_vElemComputeErrEstAFct[id] = nullptr;
	m_vElemComputeErrEstMFct[id] = nullptr;
	m_vElemComputeErrEstRhsFct[id] = nullptr;
	m_vFinishErrEstElemLoopFct[id] = nullptr;
}


//template <typename TDomain>
//void IElemDisc<TDomain>::clear_add_fct(ReferenceObjectID id)
//{
	/*m_vPrepareTimestepElemFct[id] = nullptr;
	m_vFinishTimestepElemFct[id] = nullptr;

	m_vPrepareElemLoopFct[id] = nullptr;
	m_vPrepareElemFct[id] = nullptr;
	m_vFinishElemLoopFct[id] = nullptr;

	m_vElemJAFct[id] = nullptr;
	m_vElemJMFct[id] = nullptr;

	m_vElemdAFct[id] = nullptr;
	m_vElemdAExplFct[id] = nullptr;
	m_vElemdMFct[id] = nullptr;

	m_vElemRHSFct[id] = nullptr;

	m_vPrepareErrEstElemLoopFct[id] = nullptr;
	m_vPrepareErrEstElemFct[id] = nullptr;
	m_vElemComputeErrEstAFct[id] = nullptr;
	m_vElemComputeErrEstMFct[id] = nullptr;
	m_vElemComputeErrEstRhsFct[id] = nullptr;
	m_vFinishErrEstElemLoopFct[id] = nullptr;
*/
//	IElemAssembleFuncs<IElemDisc<TDomain>, TDomain>::clear_add_fct(id);
//	IElemEstimatorFuncs<IElemDisc<TDomain>, TDomain>::clear_add_fct(id);
//}


template <typename TLeaf, typename TDomain>
void IElemAssembleFuncs<TLeaf, TDomain>::
clear_add_fct()
{
	for(size_t i = 0; i < NUM_REFERENCE_OBJECTS; ++i)
			clear_add_fct((ReferenceObjectID) i);

	for (size_t i = 0; i < bridge::NUM_ALGEBRA_TYPES; ++i)
	{
			m_vPrepareTimestepFct[i] = nullptr;
			m_vFinishTimestepFct[i] = nullptr;
	}
}

template <typename TLeaf, typename TDomain>
void IElemEstimatorFuncs<TLeaf, TDomain>::
clear_add_fct()
{
	for(size_t i = 0; i < NUM_REFERENCE_OBJECTS; ++i)
			clear_add_fct((ReferenceObjectID) i);
}
/*
template <typename TDomain>
void IElemDisc<TDomain>::clear_add_fct()
{
	for(size_t i = 0; i < NUM_REFERENCE_OBJECTS; ++i)
		clear_add_fct((ReferenceObjectID) i);

	assemble_base_type::clear_add_fct();
}

*/
template <typename TLeaf, typename TDomain>
void IElemAssembleFuncs<TLeaf, TDomain>::
set_default_add_fct()
{
	for(size_t i = 0; i < NUM_REFERENCE_OBJECTS; ++i)
	{
		m_vPrepareTimestepElemFct[i] = &T::prep_timestep_elem;
		m_vFinishTimestepElemFct[i] = &T::fsh_timestep_elem;

		m_vPrepareElemLoopFct[i] = &T::prep_elem_loop;
		m_vPrepareElemFct[i] = &T::prep_elem;
		m_vFinishElemLoopFct[i] = &T::fsh_elem_loop;

		m_vElemJAFct[i] = &T::add_jac_A_elem;
		m_vElemJMFct[i] = &T::add_jac_M_elem;

		m_vElemdAFct[i] = &T::add_def_A_elem;
		m_vElemdAExplFct[i] = &T::add_def_A_expl_elem;
		m_vElemdMFct[i] = &T::add_def_M_elem;

		m_vElemRHSFct[i] = &T::add_rhs_elem;
	}

	for (size_t i = 0; i < bridge::NUM_ALGEBRA_TYPES; ++i)
	{
		m_vPrepareTimestepFct[i] = &T::prep_timestep;
		m_vFinishTimestepFct[i] = &T::fsh_timestep;
	}

}

template <typename TLeaf, typename TDomain>
void IElemEstimatorFuncs<TLeaf, TDomain>::set_default_add_fct()
{
	for(size_t i = 0; i < NUM_REFERENCE_OBJECTS; ++i)
	{
		m_vPrepareErrEstElemLoopFct[i] = &T::prep_err_est_elem_loop;
		m_vPrepareErrEstElemFct[i] = &T::prep_err_est_elem;
		m_vElemComputeErrEstAFct[i] = &T::compute_err_est_A_elem;
		m_vElemComputeErrEstMFct[i] = &T::compute_err_est_M_elem;
		m_vElemComputeErrEstRhsFct[i] = &T::compute_err_est_rhs_elem;
		m_vFinishErrEstElemLoopFct[i] = &T::fsh_err_est_elem_loop;
	}

}

/*
template <typename TDomain>
void IElemDisc<TDomain>::
set_default_add_fct()
{
	IElemAssembleFuncs<IElemDisc<TDomain>, TDomain>::set_default_add_fct();
	IElemEstimatorFuncs<IElemDisc<TDomain>, TDomain>::set_default_add_fct();
}
*/

template <typename TDomain>
void IElemDiscBase<TDomain>::set_functions(const std::string& fctString)
{
	set_functions(TokenizeString(fctString));
}

template <typename TDomain>
void IElemDiscBase<TDomain>::set_functions(const std::vector<std::string>& functions)
{
	m_vFct = functions;

//	remove white space
	for(size_t i = 0; i < m_vFct.size(); ++i)
		RemoveWhitespaceFromString(m_vFct[i]);

//	if no function passed, clear functions
	if(m_vFct.size() == 1 && m_vFct[0].empty()) m_vFct.clear();

//	if functions passed with separator, but not all tokens filled, throw error
	for(size_t i = 0; i < m_vFct.size(); ++i)
	{
		if(m_vFct.empty())
			UG_THROW("Error while setting functions in an ElemDisc: passed "
							"function string lacks a "
							"function specification at position "<<i<<"(of "
							<<m_vFct.size()-1<<")");
	}

	update_function_index_mapping();
}

template <typename TDomain>
void IElemDiscBase<TDomain>::set_subsets(const std::string& ssString)
{
	set_subsets(TokenizeString(ssString));
}

template <typename TDomain>
void IElemDiscBase<TDomain>::set_subsets(const std::vector<std::string>& subsets)
{
	m_vSubset = subsets;

//	remove white space
	for(size_t i = 0; i < m_vSubset.size(); ++i)
		RemoveWhitespaceFromString(m_vSubset[i]);

//	if no subset passed, clear subsets
	if(m_vSubset.size() == 1 && m_vSubset[0].empty()) m_vSubset.clear();

//	if subsets passed with separator, but not all tokens filled, throw error
	for(size_t i = 0; i < m_vSubset.size(); ++i)
	{
		if(m_vSubset.empty())
			UG_THROW("Error while setting subsets in an ElemDisc: passed "
							"subset string lacks a "
							"subset specification at position "<<i<<"(of "
							<<m_vSubset.size()-1<<")");
	}
}

template <typename TDomain>
void IElemDiscBase<TDomain>::set_function_pattern(ConstSmartPtr<FunctionPattern> fctPatt)
{
	m_spFctPattern = fctPatt;
	update_function_index_mapping();
}

template <typename TDomain>
void IElemDiscBase<TDomain>::update_function_index_mapping()
{
//	without fct pattern, cannot create mappings
	if(m_spFctPattern.invalid()) return;

//	create function group of this elem disc
	try{
		m_fctGrp.set_function_pattern(m_spFctPattern);
		m_fctGrp.add(this->symb_fcts());
	}UG_CATCH_THROW("ElemDisc: Cannot find some symbolic Function Name.");

//	create a mapping between all functions and the function group of this
//	element disc.
	try{CreateFunctionIndexMapping(m_fctIndexMap, m_fctGrp, m_spFctPattern);
	}UG_CATCH_THROW("ElemDisc: Cannot create Function Index Mapping.");

//	set function group at imports
	for(size_t i = 0; i < m_vIImport.size(); ++i){
		m_vIImport[i]->set_map(m_fctIndexMap);
	}
}

template <typename TDomain>
void IElemDiscBase<TDomain>::check_setup(bool bNonRegularGrid)
{
//	check that all functions are defined on chosen subsets
	SubsetGroup discSubsetGrp(m_spFctPattern->subset_handler(), m_vSubset);

//	check that all functions are defined on chosen subsets
	for(size_t fct = 0; fct < m_fctGrp.size(); ++fct){
		for(size_t si = 0; si < discSubsetGrp.size(); ++si){
			if(!m_spFctPattern->is_def_in_subset(m_fctGrp[fct], discSubsetGrp[si])){
				UG_LOG("WARNING in ElemDisc: symbolic Function "<< symb_fcts()[fct]
				 << " is not defined on subset "<< symb_subsets()[si]
				 << ". This may be senseful only in particular cases.\n");
			}
		}
	}

//	request assembling for local finite element id
	std::vector<LFEID> vLfeID(m_fctGrp.size());
	for(size_t f = 0; f < vLfeID.size(); ++f)
		vLfeID[f] = m_fctGrp.local_finite_element_id(f);

	prepare_setting(vLfeID, bNonRegularGrid);
}


template <typename TDomain>
void IElemDiscBase<TDomain>::register_import(IDataImport<dim>& Imp)
{
//	check that not already registered
	for(size_t i = 0; i < m_vIImport.size(); ++i)
		if(m_vIImport[i] == &Imp)
			UG_THROW("Trying to register import twice.");

//	add it
	m_vIImport.push_back(&Imp);

	update_function_index_mapping();
}

template <typename TLeaf, typename TDomain>
void IElemEstimatorFuncs<TLeaf, TDomain>::
check_roid(ReferenceObjectID roid, int discType)
{

}


template <typename TLeaf, typename TDomain>
void IElemEstimatorFuncs<TLeaf, TDomain>::
set_roid(ReferenceObjectID roid, int discType)
{
	m_roid = roid;

	if(roid == ROID_UNKNOWN)
	{
		m_roid = ROID_UNKNOWN;
		UG_THROW("ElemDisc: Reference element type has not been set correctly.");
	}
	check_roid(roid, discType);
};


template <typename TLeaf, typename TDomain>
void IElemAssembleFuncs<TLeaf, TDomain>::
check_roid(ReferenceObjectID roid, int discType)
{
	if(m_vPrepareElemLoopFct[roid]==nullptr)
		UG_THROW("ElemDisc: Missing evaluation method 'prepare_elem_loop' for "<<roid<<" (world dim: "<<dim<<")");
	if(m_vPrepareElemFct[roid]==nullptr)
		UG_THROW("ElemDisc: Missing evaluation method 'prepare_elem' for "<<roid<<" (world dim: "<<dim<<")");
	if(m_vFinishElemLoopFct[roid]==nullptr)
		UG_THROW("ElemDisc: Missing evaluation method 'finish_elem_loop' for "<<roid<<" (world dim: "<<dim<<")");

	if(discType & MASS){
		if(m_vElemJMFct[roid]==nullptr)
			UG_THROW("ElemDisc: Missing evaluation method 'add_jac_M_elem' for "<<roid<<" (world dim: "<<dim<<")");
		if(m_vElemdMFct[roid]==nullptr)
			UG_THROW("ElemDisc: Missing evaluation method 'add_def_M_elem' for "<<roid<<" (world dim: "<<dim<<")");
	}
	if(discType & STIFF){
		if(m_vElemJAFct[roid]==nullptr)
			UG_THROW("ElemDisc: Missing evaluation method 'add_jac_A_elem' for "<<roid<<" (world dim: "<<dim<<")");
		if(m_vElemdAFct[roid]==nullptr)
			UG_THROW("ElemDisc: Missing evaluation method 'add_def_A_elem for' "<<roid<<" (world dim: "<<dim<<")");
	}
	if(discType & RHS){
		if(m_vElemRHSFct[roid]==nullptr)
			UG_THROW("ElemDisc: Missing evaluation method 'add_rhs_elem' for "<<roid<<" (world dim: "<<dim<<")");
	}
}



template <typename TLeaf, typename TDomain>
void IElemAssembleFuncs<TLeaf, TDomain>::
set_roid(ReferenceObjectID roid, int discType)
{
	m_roid = roid;

	if(roid == ROID_UNKNOWN)
	{
		m_roid = ROID_UNKNOWN;
		UG_THROW("ElemDisc: Reference element type has not been set correctly.");
	}
	check_roid(roid, discType);
};


template <typename TDomain>
void IElemDiscBase<TDomain>::
set_time_dependent(LocalVectorTimeSeries& locTimeSeries,
                   const std::vector<number>& vScaleMass,
                   const std::vector<number>& vScaleStiff)
{
	m_pLocalVectorTimeSeries = &locTimeSeries;
	m_vScaleMass = vScaleMass;
	m_vScaleStiff = vScaleStiff;
}

template <typename TDomain>
void IElemDiscBase<TDomain>::set_time_independent()
{
	m_pLocalVectorTimeSeries = nullptr;
	m_vScaleMass.clear();
	m_vScaleStiff.clear();
}

////////////////////////////////////////////////////////////////////////////////
//	throw functions
////////////////////////////////////////////////////////////////////////////////

inline void ThrowMissingVirtualMethod(const char* method, const ReferenceObjectID roid){
	UG_THROW("ElemDisc: No override for the essential assembling function "
			 "'"<<method<<"' for " << roid << " implemented!");
}
inline void ThrowMissingVirtualMethod(const char* method){
	UG_THROW("ElemDisc: No override for the essential assembling function "
			 "'"<<method<<"' implemented!");
}

////////////////////////////////////////////////////////////////////////////////
//	assembling functions dispatches
////////////////////////////////////////////////////////////////////////////////

template <typename TLeaf, typename TDomain>
void IElemAssembleFuncs<TLeaf, TDomain>::
do_prep_timestep(number future_time, const number time, VectorProxyBase* u, size_t algebra_id)
{
	//	call assembling routine
	if (this->m_vPrepareTimestepFct[algebra_id] != nullptr)
		(this->*(m_vPrepareTimestepFct[algebra_id]))(future_time, time, u);
}

template <typename TLeaf, typename TDomain>
void IElemAssembleFuncs<TLeaf, TDomain>::
do_prep_timestep_elem(const number time, LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	//	access by map
	u.access_by_map(asLeaf().map());
	if(asLeaf().local_time_series_needed())
		asLeaf().m_pLocalVectorTimeSeries->access_by_map(asLeaf().map());

	//	call assembling routine
	if (this->m_vPrepareTimestepElemFct[m_roid] != nullptr)
		(this->*(m_vPrepareTimestepElemFct[m_roid]))(time, u, elem, vCornerCoords);
}

template <typename TLeaf, typename TDomain>
void IElemAssembleFuncs<TLeaf, TDomain>::
do_prep_elem(LocalVector& u, GridObject* elem, const ReferenceObjectID roid, const MathVector<dim> vCornerCoords[])
{
	//	access by map
	u.access_by_map(asLeaf().map());
	if(asLeaf().local_time_series_needed())
		asLeaf().m_pLocalVectorTimeSeries->access_by_map(asLeaf().map());

	//	call assembling routine
	UG_ASSERT(m_vPrepareElemFct[m_roid]!=nullptr, "ElemDisc method prepare_elem missing.");
	(this->*(m_vPrepareElemFct[m_roid]))(u, elem, roid, vCornerCoords);
}

template <typename TLeaf, typename TDomain>
void IElemAssembleFuncs<TLeaf, TDomain>::
do_fsh_timestep(const number time, VectorProxyBase* u, size_t algebra_id)
{
	//	call assembling routine
	if (this->m_vFinishTimestepFct[algebra_id] != nullptr)
		(this->*(m_vFinishTimestepFct[algebra_id]))(time, u);
}

template <typename TLeaf, typename TDomain>
void IElemAssembleFuncs<TLeaf, TDomain>::
do_fsh_timestep_elem(const number time, LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	//	access by map
	u.access_by_map(asLeaf().map());
	if(asLeaf().local_time_series_needed())
		asLeaf().m_pLocalVectorTimeSeries->access_by_map(asLeaf().map());

	//	call assembling routine
	if (this->m_vFinishTimestepElemFct[m_roid] != nullptr)
		(this->*(m_vFinishTimestepElemFct[m_roid]))(time, u, elem, vCornerCoords);
}

template <typename TLeaf, typename TDomain>
void IElemAssembleFuncs<TLeaf, TDomain>::
do_prep_elem_loop(const ReferenceObjectID roid, const int si)
{
//	set id and disc part (this checks the assemble functions)
	set_roid(roid, si);

//	remove positions in currently registered imports
	for(size_t i = 0; i < asLeaf().m_vIImport.size(); ++i)
		asLeaf().m_vIImport[i]->clear_ips();

//	call prep_elem_loop (this may set ip-series to imports)
	//	call assembling routine
	UG_ASSERT(m_vPrepareElemLoopFct[m_roid]!=nullptr, "ElemDisc method prepare_elem_loop missing.");
	(this->*m_vPrepareElemLoopFct[m_roid])(roid, si);

//	set roid in imports (for evaluation function)
	for(size_t i = 0; i < asLeaf().m_vIImport.size(); ++i)
		asLeaf().m_vIImport[i]->set_roid(roid);
}

template <typename TLeaf, typename TDomain>
void IElemAssembleFuncs<TLeaf, TDomain>::
do_fsh_elem_loop()
{
//	call finish
	UG_ASSERT(m_vFinishElemLoopFct[m_roid]!=nullptr, "ElemDisc method finish_elem_loop missing.");
	(this->*m_vFinishElemLoopFct[m_roid])();

//	remove positions in currently registered imports
	for(size_t i = 0; i < asLeaf().m_vIImport.size(); ++i)
		asLeaf().m_vIImport[i]->clear_ips();
}

template <typename TLeaf, typename TDomain>
void IElemAssembleFuncs<TLeaf, TDomain>::
do_add_jac_A_elem(LocalMatrix& J, LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	//	access by map
	u.access_by_map(asLeaf().map());
	J.access_by_map(asLeaf().map());
	if(asLeaf().local_time_series_needed())
		asLeaf().m_pLocalVectorTimeSeries->access_by_map(asLeaf().map());

	//	call assembling routine
	UG_ASSERT(m_vElemJAFct[m_roid]!=nullptr, "ElemDisc method add_jac_A missing.");
	(this->*m_vElemJAFct[m_roid])(J, u, elem, vCornerCoords);
}

template <typename TLeaf, typename TDomain>
void IElemAssembleFuncs<TLeaf, TDomain>::
do_add_jac_M_elem(LocalMatrix& J, LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	// check if really needed (may occur in cases, when mixing stat and instat)
	if(asLeaf().m_bStationaryForced) return;

	//	access by map
	u.access_by_map(asLeaf().map());
	J.access_by_map(asLeaf().map());
	if(asLeaf().local_time_series_needed())
		asLeaf().m_pLocalVectorTimeSeries->access_by_map(asLeaf().map());

	//	call assembling routine
	UG_ASSERT(m_vElemJMFct[m_roid]!=nullptr, "ElemDisc method add_jac_M missing.");
	(this->*m_vElemJMFct[m_roid])(J, u, elem, vCornerCoords);
}

template <typename TLeaf, typename TDomain>
void IElemAssembleFuncs<TLeaf, TDomain>::
do_add_def_A_elem(LocalVector& d, LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	//	access by map
	u.access_by_map(asLeaf().map());
	d.access_by_map(asLeaf().map());
	if(asLeaf().local_time_series_needed())
		asLeaf().m_pLocalVectorTimeSeries->access_by_map(asLeaf().map());

	//	call assembling routine
	UG_ASSERT(m_vElemdAFct[m_roid]!=nullptr, "ElemDisc method add_def_A missing.");
	(this->*m_vElemdAFct[m_roid])(d, u, elem, vCornerCoords);
}

template <typename TLeaf, typename TDomain>
void IElemAssembleFuncs<TLeaf, TDomain>::
do_add_def_A_expl_elem(LocalVector& d, LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	//	access by map
	u.access_by_map(asLeaf().map());
	d.access_by_map(asLeaf().map());
	if(asLeaf().local_time_series_needed())
		asLeaf().m_pLocalVectorTimeSeries->access_by_map(asLeaf().map());

	//	call assembling routine
	if(this->m_vElemdAExplFct[m_roid] != nullptr)
		(this->*m_vElemdAExplFct[m_roid])(d, u, elem, vCornerCoords);
}

template <typename TLeaf, typename TDomain>
void IElemAssembleFuncs<TLeaf, TDomain>::
do_add_def_M_elem(LocalVector& d, LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	// check if really needed (may occur in cases, when mixing stat and instat)
	if(asLeaf().m_bStationaryForced) return;

	//	access by map
	u.access_by_map(asLeaf().map());
	d.access_by_map(asLeaf().map());
	if(asLeaf().local_time_series_needed())
		asLeaf().m_pLocalVectorTimeSeries->access_by_map(asLeaf().map());

	//	call assembling routine
	UG_ASSERT(m_vElemdMFct[m_roid]!=nullptr, "ElemDisc method add_def_M missing.");
	(this->*m_vElemdMFct[m_roid])(d, u, elem, vCornerCoords);
}

template <typename TLeaf, typename TDomain>
void IElemAssembleFuncs<TLeaf, TDomain>::
do_add_rhs_elem(LocalVector& rhs, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	//	access by map
	rhs.access_by_map(asLeaf().map());
	if(asLeaf().local_time_series_needed())
		asLeaf().m_pLocalVectorTimeSeries->access_by_map(asLeaf().map());

	//	call assembling routine
	UG_ASSERT(m_vElemRHSFct[m_roid]!=nullptr, "ElemDisc method add_rhs missing.");
	(this->*m_vElemRHSFct[m_roid])(rhs, elem, vCornerCoords);
}

template <typename TLeaf, typename TDomain>
void IElemEstimatorFuncs<TLeaf, TDomain>::
do_prep_err_est_elem_loop(const ReferenceObjectID roid, const int si)
{
//	set id and disc part (this checks the assemble functions)
	set_roid(roid, si);

//	remove positions in currently registered imports
	for(std::size_t i = 0; i < asLeaf().num_imports(); ++i)
		asLeaf().get_import(i).clear_ips();

//	call prep_elem_loop (this may set ip-series to imports)
	UG_ASSERT(m_vPrepareErrEstElemLoopFct[m_roid]!=nullptr, "ElemDisc method prep_err_est_elem_loop missing.");
	(this->*m_vPrepareErrEstElemLoopFct[m_roid])(roid, si);

//	set roid in imports (for evaluation function)
	for(std::size_t i = 0; i < asLeaf().num_imports(); ++i)
		asLeaf().get_import(i).set_roid(roid);
}

template <typename TLeaf, typename TDomain>
void IElemEstimatorFuncs<TLeaf, TDomain>::
do_prep_err_est_elem(LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	//	access by map
	u.access_by_map(asLeaf().map());
	if (asLeaf().local_time_series_needed())
		asLeaf().m_pLocalVectorTimeSeries->access_by_map(asLeaf().map());

	//	call assembling routine
	UG_ASSERT(m_vPrepareErrEstElemFct[m_roid]!=nullptr, "ElemDisc method prepare_err_est_elem missing.");
	(this->*(m_vPrepareErrEstElemFct[m_roid]))(u, elem, vCornerCoords);
}

template <typename TLeaf, typename TDomain>
void IElemEstimatorFuncs<TLeaf, TDomain>::
do_compute_err_est_A_elem(LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[], const number& scale)
{
	//	access by map
	u.access_by_map(asLeaf().map());
	if (asLeaf().local_time_series_needed())
		asLeaf().m_pLocalVectorTimeSeries->access_by_map(asLeaf().map());

	//	call assembling routine
	if (this->m_vElemComputeErrEstAFct[m_roid] != nullptr)
		(this->*(m_vElemComputeErrEstAFct[m_roid]))(u, elem, vCornerCoords, scale);
}

template <typename TLeaf, typename TDomain>
void IElemEstimatorFuncs<TLeaf, TDomain>::
do_compute_err_est_M_elem(LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[], const number& scale)
{
	//	access by map
	u.access_by_map(asLeaf().map());
	if(asLeaf().local_time_series_needed())
		asLeaf().m_pLocalVectorTimeSeries->access_by_map(asLeaf().map());

	//	call assembling routine
	if(this->m_vElemComputeErrEstMFct[m_roid] != nullptr)
		(this->*(m_vElemComputeErrEstMFct[m_roid]))(u, elem, vCornerCoords, scale);
}

template <typename TLeaf, typename TDomain>
void IElemEstimatorFuncs<TLeaf, TDomain>::
do_compute_err_est_rhs_elem(GridObject* elem, const MathVector<dim> vCornerCoords[], const number& scale)
{
	if(asLeaf().local_time_series_needed())
		asLeaf().m_pLocalVectorTimeSeries->access_by_map(asLeaf().map());

	//	call assembling routine
	if(this->m_vElemComputeErrEstRhsFct[m_roid] != nullptr)
		(this->*(m_vElemComputeErrEstRhsFct[m_roid]))(elem, vCornerCoords, scale);
}

template <typename TLeaf, typename TDomain>
void IElemEstimatorFuncs<TLeaf, TDomain>::
do_fsh_err_est_elem_loop()
{
//	call finish
	if (this->m_vFinishErrEstElemLoopFct[m_roid] != nullptr)
		(this->*m_vFinishErrEstElemLoopFct[m_roid])();

//	remove positions in currently registered imports
	for(size_t i = 0; i < asLeaf().num_imports(); ++i)
		asLeaf().get_import(i).clear_ips();
}



////////////////////////////////////////////////////////////////////////////////
//	virtual assembling functions default
////////////////////////////////////////////////////////////////////////////////

template <typename TLeaf, typename TDomain>
void IElemAssembleFuncs<TLeaf, TDomain>::
prep_timestep(number future_time, number time, VectorProxyBase* u)
{
	// do nothing
}

template <typename TLeaf, typename TDomain>
void IElemAssembleFuncs<TLeaf, TDomain>::
prep_timestep_elem(const number time, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	// ThrowMissingVirtualMethod("prep_timestep_elem", elem->reference_object_id ());
}

template <typename TLeaf, typename TDomain>
void IElemAssembleFuncs<TLeaf, TDomain>::
prep_elem(const LocalVector& u, GridObject* elem, const ReferenceObjectID roid, const MathVector<dim> vCornerCoords[])
{
	ThrowMissingVirtualMethod("prep_elem", elem->reference_object_id ());
}

template <typename TLeaf, typename TDomain>
void IElemAssembleFuncs<TLeaf, TDomain>::
fsh_timestep(number time, VectorProxyBase* u)
{
	// do nothing
}

template <typename TLeaf, typename TDomain>
void IElemAssembleFuncs<TLeaf, TDomain>::
fsh_timestep_elem(const number time, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	// ThrowMissingVirtualMethod("fsh_timestep_elem", elem->reference_object_id ());
}

template <typename TLeaf, typename TDomain>
void IElemAssembleFuncs<TLeaf, TDomain>::
prep_elem_loop(const ReferenceObjectID roid, const int si)
{
	ThrowMissingVirtualMethod("prep_elem_loop", roid);
}

template <typename TLeaf, typename TDomain>
void IElemAssembleFuncs<TLeaf, TDomain>::
fsh_elem_loop()
{
	ThrowMissingVirtualMethod("fsh_elem_loop");
}

template <typename TLeaf, typename TDomain>
void IElemAssembleFuncs<TLeaf, TDomain>::
add_jac_A_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	ThrowMissingVirtualMethod("add_jac_A_elem", elem->reference_object_id ());
}

template <typename TLeaf, typename TDomain>
void IElemAssembleFuncs<TLeaf, TDomain>::
add_jac_M_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	ThrowMissingVirtualMethod("add_jac_M_elem", elem->reference_object_id ());
}

template <typename TLeaf, typename TDomain>
void IElemAssembleFuncs<TLeaf, TDomain>::
add_def_A_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	ThrowMissingVirtualMethod("add_def_A_elem", elem->reference_object_id ());
}

template <typename TLeaf, typename TDomain>
void IElemAssembleFuncs<TLeaf, TDomain>::
add_def_A_expl_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	// ThrowMissingVirtualMethod("add_def_A_expl_elem", elem->reference_object_id ());
}

template <typename TLeaf, typename TDomain>
void IElemAssembleFuncs<TLeaf, TDomain>::
add_def_M_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	ThrowMissingVirtualMethod("add_def_M_elem", elem->reference_object_id ());
}

template <typename TLeaf, typename TDomain>
void IElemAssembleFuncs<TLeaf, TDomain>::
add_rhs_elem(LocalVector& rhs, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	ThrowMissingVirtualMethod("add_rhs_elem", elem->reference_object_id ());
}

template <typename TLeaf, typename TDomain>
void IElemEstimatorFuncs<TLeaf, TDomain>::
prep_err_est_elem_loop(const ReferenceObjectID roid, const int si)
{
	// ThrowMissingVirtualMethod("prep_err_est_elem_loop", roi);
}

template <typename TLeaf, typename TDomain>
void IElemEstimatorFuncs<TLeaf, TDomain>::
prep_err_est_elem(const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	//ThrowMissingVirtualMethod("prep_elem", elem->reference_object_id ());
}

template <typename TLeaf, typename TDomain>
void IElemEstimatorFuncs<TLeaf, TDomain>::
compute_err_est_A_elem(const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[], const number& scale)
{
	// ThrowMissingVirtualMethod("compute_err_est_elem", elem->reference_object_id ());
}

template <typename TLeaf, typename TDomain>
void IElemEstimatorFuncs<TLeaf, TDomain>::
compute_err_est_M_elem(const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[], const number& scale)
{
	// ThrowMissingVirtualMethod("compute_err_est_elem", elem->reference_object_id ());
}

template <typename TLeaf, typename TDomain>
void IElemEstimatorFuncs<TLeaf, TDomain>::
compute_err_est_rhs_elem(GridObject* elem, const MathVector<dim> vCornerCoords[], const number& scale)
{
	// ThrowMissingVirtualMethod("compute_err_est_elem", elem->reference_object_id ());
}

template <typename TLeaf, typename TDomain>
void IElemEstimatorFuncs<TLeaf, TDomain>::
fsh_err_est_elem_loop()
{
	// ThrowMissingVirtualMethod("fsh_err_est_elem_loop");
}


////////////////////////////////////////////////////////////////////////////////
//	explicit template instantiations
////////////////////////////////////////////////////////////////////////////////

#ifdef UG_DIM_1
template class IElemDiscBase<Domain1d>;
template class IElemAssembleFuncs<IElemDisc<Domain1d>, Domain1d>;
template class IElemEstimatorFuncs<IElemDisc<Domain1d>, Domain1d>;
template class IElemError<Domain1d>;
template class IElemDisc<Domain1d>;
#endif
#ifdef UG_DIM_2
template class IElemDiscBase<Domain2d>;
template class IElemAssembleFuncs<IElemDisc<Domain2d>, Domain2d>;
template class IElemEstimatorFuncs<IElemDisc<Domain2d>, Domain2d>;
template class IElemError<Domain2d>;
template class IElemDisc<Domain2d>;
#endif
#ifdef UG_DIM_3
template class IElemDiscBase<Domain3d>;
template class IElemAssembleFuncs<IElemDisc<Domain3d>, Domain3d>;
template class IElemEstimatorFuncs<IElemDisc<Domain3d>, Domain3d>;
template class IElemError<Domain3d>;
template class IElemDisc<Domain3d>;
#endif

} // end namespace ug

