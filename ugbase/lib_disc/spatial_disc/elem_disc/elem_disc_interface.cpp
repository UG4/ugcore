/*
 * elem_disc_interface.cpp
 *
 *  Created on: 03.07.2011
 *      Author: andreasvogel
 */

#include "elem_disc_interface.h"
#include "lib_disc/common/groups_util.h"

namespace ug{

template <typename TDomain>
IElemDisc<TDomain>::IElemDisc(const char* functions, const char* subsets)
	: 	m_spApproxSpace(NULL), m_spFctPattern(0),
	  	m_timePoint(0), m_pLocalVectorTimeSeries(NULL), m_bStationaryForced(false),
	  	m_bFastAssembleEnabled(false), m_id(ROID_UNKNOWN)
{
	if(functions == NULL) functions = "";
	if(subsets == NULL) subsets = "";
	set_functions(functions);
	set_subsets(subsets);
	clear_add_fct();
}

template <typename TDomain>
IElemDisc<TDomain>::IElemDisc(const std::vector<std::string>& vFct,
                              const std::vector<std::string>& vSubset)
	: 	m_spApproxSpace(NULL), m_spFctPattern(0),
		m_timePoint(0), m_pLocalVectorTimeSeries(NULL), m_bStationaryForced(false),
		m_bFastAssembleEnabled(false), m_id(ROID_UNKNOWN)
{
	set_functions(vFct);
	set_subsets(vSubset);
	clear_add_fct();
}

template <typename TDomain>
void IElemDisc<TDomain>::
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

template <typename TDomain>
void IElemDisc<TDomain>::clear_add_fct()
{
	for(size_t i = 0; i < NUM_REFERENCE_OBJECTS; ++i)
	{
		m_vPrepareTimestepElemFct[i] = NULL;
		m_vFinishTimestepElemFct[i] = NULL;

		m_vPrepareElemLoopFct[i] = NULL;
		m_vPrepareElemFct[i] = NULL;
		m_vFinishElemLoopFct[i] = NULL;

		m_vElemJAFct[i] = NULL;
		m_vElemJMFct[i] = NULL;

		m_vElemdAFct[i] = NULL;
		m_vElemdAExplFct[i] = NULL;
		m_vElemdMFct[i] = NULL;

		m_vElemRHSFct[i] = NULL;
	}
}


template <typename TDomain>
void IElemDisc<TDomain>::set_functions(const std::string& fctString)
{
	set_functions(TokenizeString(fctString));
}

template <typename TDomain>
void IElemDisc<TDomain>::set_functions(const std::vector<std::string>& functions)
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
void IElemDisc<TDomain>::set_subsets(const std::string& ssString)
{
	set_subsets(TokenizeString(ssString));
}

template <typename TDomain>
void IElemDisc<TDomain>::set_subsets(const std::vector<std::string>& subsets)
{
	m_vSubset = subsets;

//	remove white space
	for(size_t i = 0; i < m_vSubset.size(); ++i)
		RemoveWhitespaceFromString(m_vSubset[i]);

//	if no subset passed, clear subsets
	if(m_vFct.size() == 1 && m_vFct[0].empty()) m_vFct.clear();

//	if subsets passed with separator, but not all tokens filled, throw error
	for(size_t i = 0; i < m_vFct.size(); ++i)
	{
		if(m_vFct.empty())
			UG_THROW("Error while setting subsets in an ElemDisc: passed "
							"subset string lacks a "
							"subset specification at position "<<i<<"(of "
							<<m_vFct.size()-1<<")");
	}
}

template <typename TDomain>
void IElemDisc<TDomain>::set_function_pattern(ConstSmartPtr<FunctionPattern> fctPatt)
{
	m_spFctPattern = fctPatt;
	update_function_index_mapping();
}

template <typename TDomain>
void IElemDisc<TDomain>::update_function_index_mapping()
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
void IElemDisc<TDomain>::check_setup(bool bNonRegularGrid)
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
void IElemDisc<TDomain>::register_import(IDataImport<dim>& Imp)
{
//	check that not already registered
	for(size_t i = 0; i < m_vIImport.size(); ++i)
		if(m_vIImport[i] == &Imp)
			UG_THROW("Trying to register import twice.");

//	add it
	m_vIImport.push_back(&Imp);

	update_function_index_mapping();
}

template <typename TDomain>
void IElemDisc<TDomain>::set_roid(ReferenceObjectID roid, int discType)
{
	m_id = roid;

	if(roid == ROID_UNKNOWN)
	{
		m_id = ROID_UNKNOWN;
		UG_THROW("ElemDisc: Reference element type has not been set correctly.");
	}

	if(fast_add_elem_enabled()){
		if(m_vPrepareElemLoopFct[m_id]==NULL)
			UG_THROW("ElemDisc: Missing evaluation method 'prepare_elem_loop' for "<<roid<<"(world dim: "<<dim<<")");
		if(m_vPrepareElemFct[m_id]==NULL)
			UG_THROW("ElemDisc: Missing evaluation method 'prepare_elem' for "<<roid<<"(world dim: "<<dim<<")");
		if(m_vFinishElemLoopFct[m_id]==NULL)
			UG_THROW("ElemDisc: Missing evaluation method 'finish_elem_loop' for "<<roid<<"(world dim: "<<dim<<")");

		if(discType & MASS){
			if(m_vElemJMFct[m_id]==NULL)
				UG_THROW("ElemDisc: Missing evaluation method 'add_jac_M_elem' for "<<roid<<"(world dim: "<<dim<<")");
			if(m_vElemdMFct[m_id]==NULL)
				UG_THROW("ElemDisc: Missing evaluation method 'add_def_M_elem' for "<<roid<<"(world dim: "<<dim<<")");
		}
		if(discType & STIFF){
			if(m_vElemJAFct[m_id]==NULL)
				UG_THROW("ElemDisc: Missing evaluation method 'add_jac_A_elem' for "<<roid<<"(world dim: "<<dim<<")");
			if(m_vElemdAFct[m_id]==NULL)
				UG_THROW("ElemDisc: Missing evaluation method 'add_def_A_elem for' "<<roid<<"(world dim: "<<dim<<")");
		}
		if(discType & RHS){
			if(m_vElemRHSFct[m_id]==NULL)
				UG_THROW("ElemDisc: Missing evaluation method 'add_rhs_elem' for "<<roid<<"(world dim: "<<dim<<")");
		}
	}
};

template <typename TDomain>
void IElemDisc<TDomain>::
set_time_dependent(LocalVectorTimeSeries& locTimeSeries,
                   const std::vector<number>& vScaleMass,
                   const std::vector<number>& vScaleStiff)
{
	m_pLocalVectorTimeSeries = &locTimeSeries;
	m_vScaleMass = vScaleMass;
	m_vScaleStiff = vScaleStiff;
}

template <typename TDomain>
void IElemDisc<TDomain>::set_time_independent()
{
	m_pLocalVectorTimeSeries = NULL;
	m_vScaleMass.clear();
	m_vScaleStiff.clear();
}

////////////////////////////////////////////////////////////////////////////////
//	assembling functions dispatches
////////////////////////////////////////////////////////////////////////////////

template <typename TDomain>
void IElemDisc<TDomain>::
do_prep_timestep_elem(const number time, LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	//	access by map
	u.access_by_map(map());
	if(local_time_series_needed())
		m_pLocalVectorTimeSeries->access_by_map(map());

	//	call assembling routine
	if(fast_add_elem_enabled()){
		if (this->m_vPrepareTimestepElemFct[m_id] != NULL)
			(this->*(m_vPrepareTimestepElemFct[m_id]))(time, u, elem, vCornerCoords);
	} else {
		prep_timestep_elem(time, u, elem, vCornerCoords);
	}
}

template <typename TDomain>
void IElemDisc<TDomain>::
do_prep_elem(LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	//	access by map
	u.access_by_map(map());
	if(local_time_series_needed())
		m_pLocalVectorTimeSeries->access_by_map(map());

	//	call assembling routine
	if(fast_add_elem_enabled()){
		UG_ASSERT(m_vPrepareElemFct[m_id]!=NULL, "Fast-Assemble Method missing.");
		(this->*(m_vPrepareElemFct[m_id]))(u, elem, vCornerCoords);
	} else {
		prep_elem(u, elem, vCornerCoords);
	}
}

template <typename TDomain>
void IElemDisc<TDomain>::
do_fsh_timestep_elem(const number time, LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	//	access by map
	u.access_by_map(map());
	if(local_time_series_needed())
		m_pLocalVectorTimeSeries->access_by_map(map());

	//	call assembling routine
	if(fast_add_elem_enabled()){
		if (this->m_vFinishTimestepElemFct[m_id] != NULL)
			(this->*(m_vFinishTimestepElemFct[m_id]))(time, u, elem, vCornerCoords);
	} else {
		fsh_timestep_elem(time, u, elem, vCornerCoords);
	}
}

template <typename TDomain>
void IElemDisc<TDomain>::
do_prep_elem_loop(const ReferenceObjectID roid, const int si)
{
//	set id and disc part (this checks and inits fast-assemble functions)
	this->set_roid(roid, si);

//	remove positions in currently registered imports
	for(size_t i = 0; i < m_vIImport.size(); ++i)
		m_vIImport[i]->clear_ips();

//	call prep_elem_loop (this may set ip-series to imports)
	//	call assembling routine
	if(fast_add_elem_enabled()){
		(this->*m_vPrepareElemLoopFct[m_id])(roid, si);
	} else {
		prep_elem_loop(roid, si);
	}

//	set roid in imports (for evaluation function)
	for(size_t i = 0; i < m_vIImport.size(); ++i)
		m_vIImport[i]->set_roid(roid);
}

template <typename TDomain>
void IElemDisc<TDomain>::
do_fsh_elem_loop()
{
	UG_ASSERT(m_vFinishElemLoopFct[m_id]!=NULL, "Fast-Assemble Method missing.");

//	call finish
	if(fast_add_elem_enabled()){
		(this->*m_vFinishElemLoopFct[m_id])();
	} else {
		fsh_elem_loop();
	}

//	remove positions in currently registered imports
	for(size_t i = 0; i < m_vIImport.size(); ++i)
		m_vIImport[i]->clear_ips();
}

template <typename TDomain>
void IElemDisc<TDomain>::
do_add_jac_A_elem(LocalMatrix& J, LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	//	access by map
	u.access_by_map(map());
	J.access_by_map(map());
	if(local_time_series_needed())
		m_pLocalVectorTimeSeries->access_by_map(map());

	//	call assembling routine
	if(fast_add_elem_enabled()){
		UG_ASSERT(m_vElemJAFct[m_id]!=NULL, "Fast-Assemble Method missing.");
		(this->*m_vElemJAFct[m_id])(J, u, elem, vCornerCoords);
	} else {
		add_jac_A_elem(J, u, elem, vCornerCoords);
	}
}

template <typename TDomain>
void IElemDisc<TDomain>::
do_add_jac_M_elem(LocalMatrix& J, LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	// check if really needed (may occur in cases, when mixing stat and instat)
	if(m_bStationaryForced) return;

	//	access by map
	u.access_by_map(map());
	J.access_by_map(map());
	if(local_time_series_needed())
		m_pLocalVectorTimeSeries->access_by_map(map());

	//	call assembling routine
	if(fast_add_elem_enabled()){
		UG_ASSERT(m_vElemJMFct[m_id]!=NULL, "Fast-Assemble Method missing.");
		(this->*m_vElemJMFct[m_id])(J, u, elem, vCornerCoords);
	} else {
		add_jac_M_elem(J, u, elem, vCornerCoords);
	}
}

template <typename TDomain>
void IElemDisc<TDomain>::
do_add_def_A_elem(LocalVector& d, LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	//	access by map
	u.access_by_map(map());
	d.access_by_map(map());
	if(local_time_series_needed())
		m_pLocalVectorTimeSeries->access_by_map(map());

	//	call assembling routine
	if(fast_add_elem_enabled()){
		UG_ASSERT(m_vElemdAFct[m_id]!=NULL, "Fast-Assemble Method missing.");
		(this->*m_vElemdAFct[m_id])(d, u, elem, vCornerCoords);
	} else {
		add_def_A_elem(d, u, elem, vCornerCoords);
	}
}

template <typename TDomain>
void IElemDisc<TDomain>::
do_add_def_A_expl_elem(LocalVector& d, LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	//	access by map
	u.access_by_map(map());
	d.access_by_map(map());
	if(local_time_series_needed())
		m_pLocalVectorTimeSeries->access_by_map(map());

	//	call assembling routine
	if(fast_add_elem_enabled()){
		if(this->m_vElemdAExplFct[m_id] != NULL)
			(this->*m_vElemdAExplFct[m_id])(d, u, elem, vCornerCoords);
	} else {
		add_def_A_expl_elem(d, u, elem, vCornerCoords);
	}
}

template <typename TDomain>
void IElemDisc<TDomain>::
do_add_def_M_elem(LocalVector& d, LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	// check if really needed (may occur in cases, when mixing stat and instat)
	if(m_bStationaryForced) return;

	//	access by map
	u.access_by_map(map());
	d.access_by_map(map());
	if(local_time_series_needed())
		m_pLocalVectorTimeSeries->access_by_map(map());

	//	call assembling routine
	if(fast_add_elem_enabled()){
		UG_ASSERT(m_vElemdMFct[m_id]!=NULL, "Fast-Assemble Method missing.");
		(this->*m_vElemdMFct[m_id])(d, u, elem, vCornerCoords);
	} else {
		add_def_M_elem(d, u, elem, vCornerCoords);
	}
}

template <typename TDomain>
void IElemDisc<TDomain>::
do_add_rhs_elem(LocalVector& rhs, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	//	access by map
	rhs.access_by_map(map());
	if(local_time_series_needed())
		m_pLocalVectorTimeSeries->access_by_map(map());

	//	call assembling routine
	if(fast_add_elem_enabled()){
		UG_ASSERT(m_vElemRHSFct[m_id]!=NULL, "Fast-Assemble Method missing.");
		(this->*m_vElemRHSFct[m_id])(rhs, elem, vCornerCoords);
	} else {
		add_rhs_elem(rhs, elem, vCornerCoords);
	}
}

template <typename TDomain>
void IElemDisc<TDomain>::
do_prep_err_est_elem_loop(const ReferenceObjectID roid, const int si)
{
//	set id and disc part (this checks and inits fast-assemble functions)
	this->set_roid(roid, si);

//	remove positions in currently registered imports
	for(size_t i = 0; i < m_vIImport.size(); ++i)
		m_vIImport[i]->clear_ips();

//	call prep_elem_loop (this may set ip-series to imports)
	//	call assembling routine
	if(fast_add_elem_enabled()){
		if (this->m_vPrepareErrEstElemLoopFct[m_id] != NULL)
			(this->*m_vPrepareErrEstElemLoopFct[m_id])(roid, si);
	} else {
		prep_err_est_elem_loop(roid, si);
	}

//	set roid in imports (for evaluation function)
	for(size_t i = 0; i < m_vIImport.size(); ++i)
		m_vIImport[i]->set_roid(roid);
}

template <typename TDomain>
void IElemDisc<TDomain>::
do_compute_err_est_elem(LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	//	access by map
	u.access_by_map(map());
	if(local_time_series_needed())
		m_pLocalVectorTimeSeries->access_by_map(map());

	//	call assembling routine
	if(fast_add_elem_enabled()){
		if(this->m_vElemComputeErrEstFct[m_id] != NULL)
			(this->*(m_vElemComputeErrEstFct[m_id]))(u, elem, vCornerCoords);
	} else {
		compute_err_est_elem(u, elem, vCornerCoords);
	}
}

template <typename TDomain>
number IElemDisc<TDomain>::
do_get_err_est_elem(LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	//	access by map
	u.access_by_map(map());
	if(local_time_series_needed())
		m_pLocalVectorTimeSeries->access_by_map(map());

	//	call assembling routine
	if(fast_add_elem_enabled()){
		if(this->m_vElemGetErrEstFct[m_id] != NULL)
			return (this->*(m_vElemGetErrEstFct[m_id]))(u, elem, vCornerCoords);
		return 0;
	} else {
		return get_err_est_elem(u, elem, vCornerCoords);
	}
}

template <typename TDomain>
void IElemDisc<TDomain>::
do_fsh_err_est_elem_loop()
{
//	call finish
	if(fast_add_elem_enabled()){
		if (this->m_vFinishErrEstElemLoopFct[m_id] != NULL)
			(this->*m_vFinishErrEstElemLoopFct[m_id])();
	} else {
		fsh_err_est_elem_loop();
	}

//	remove positions in currently registered imports
	for(size_t i = 0; i < m_vIImport.size(); ++i)
		m_vIImport[i]->clear_ips();
}



////////////////////////////////////////////////////////////////////////////////
//	virtual assembling functions default
////////////////////////////////////////////////////////////////////////////////

static void ThrowMissingVirtualMethod(const char* method){
	UG_THROW("ElemDisc: using virtual assembling functions, but no override "
			 "for '"<<method<<"' implemented (although called).");
}

template <typename TDomain>
void IElemDisc<TDomain>::
prep_timestep_elem(const number time, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	ThrowMissingVirtualMethod("prep_timestep_elem");
}

template <typename TDomain>
void IElemDisc<TDomain>::
prep_elem(const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	ThrowMissingVirtualMethod("prep_elem");
}

template <typename TDomain>
void IElemDisc<TDomain>::
fsh_timestep_elem(const number time, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	ThrowMissingVirtualMethod("fsh_timestep_elem");
}

template <typename TDomain>
void IElemDisc<TDomain>::
prep_elem_loop(const ReferenceObjectID roid, const int si)
{
	ThrowMissingVirtualMethod("prep_elem_loop");
}

template <typename TDomain>
void IElemDisc<TDomain>::
fsh_elem_loop()
{
	ThrowMissingVirtualMethod("fsh_elem_loop");
}

template <typename TDomain>
void IElemDisc<TDomain>::
add_jac_A_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	ThrowMissingVirtualMethod("add_jac_A_elem");
}

template <typename TDomain>
void IElemDisc<TDomain>::
add_jac_M_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	ThrowMissingVirtualMethod("add_jac_M_elem");
}

template <typename TDomain>
void IElemDisc<TDomain>::
add_def_A_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	ThrowMissingVirtualMethod("add_def_A_elem");
}

template <typename TDomain>
void IElemDisc<TDomain>::
add_def_A_expl_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	ThrowMissingVirtualMethod("add_def_A_expl_elem");
}

template <typename TDomain>
void IElemDisc<TDomain>::
add_def_M_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	ThrowMissingVirtualMethod("add_def_M_elem");
}

template <typename TDomain>
void IElemDisc<TDomain>::
add_rhs_elem(LocalVector& rhs, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	ThrowMissingVirtualMethod("add_rhs_elem");
}

template <typename TDomain>
void IElemDisc<TDomain>::
prep_err_est_elem_loop(const ReferenceObjectID roid, const int si)
{
	ThrowMissingVirtualMethod("prep_err_est_elem_loop");
}

template <typename TDomain>
void IElemDisc<TDomain>::
compute_err_est_elem(const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	ThrowMissingVirtualMethod("compute_err_est_elem");
}

template <typename TDomain>
number IElemDisc<TDomain>::
get_err_est_elem(const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	ThrowMissingVirtualMethod("get_err_est_elem");
	return 0;
}

template <typename TDomain>
void IElemDisc<TDomain>::
fsh_err_est_elem_loop()
{
	ThrowMissingVirtualMethod("fsh_err_est_elem_loop");
}


////////////////////////////////////////////////////////////////////////////////
//	explicit template instantiations
////////////////////////////////////////////////////////////////////////////////

#ifdef UG_DIM_1
template class IElemDisc<Domain1d>;
#endif
#ifdef UG_DIM_2
template class IElemDisc<Domain2d>;
#endif
#ifdef UG_DIM_3
template class IElemDisc<Domain3d>;
#endif

} // end namespace ug

