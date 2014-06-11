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
	  	m_bDoErrEst(false),
	  	m_id(ROID_UNKNOWN)
{
	if(functions == NULL) functions = "";
	if(subsets == NULL) subsets = "";
	set_functions(functions);
	set_subsets(subsets);
	set_default_add_fct();
}

template <typename TDomain>
IElemDisc<TDomain>::IElemDisc(const std::vector<std::string>& vFct,
                              const std::vector<std::string>& vSubset)
	: 	m_spApproxSpace(NULL), m_spFctPattern(0),
		m_timePoint(0), m_pLocalVectorTimeSeries(NULL), m_bStationaryForced(false),
	  	m_bDoErrEst(false),
		m_id(ROID_UNKNOWN)
{
	set_functions(vFct);
	set_subsets(vSubset);
	set_default_add_fct();
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
void IElemDisc<TDomain>::clear_add_fct(ReferenceObjectID id)
{
	m_vPrepareTimestepElemFct[id] = NULL;
	m_vFinishTimestepElemFct[id] = NULL;

	m_vPrepareElemLoopFct[id] = NULL;
	m_vPrepareElemFct[id] = NULL;
	m_vFinishElemLoopFct[id] = NULL;

	m_vElemJAFct[id] = NULL;
	m_vElemJMFct[id] = NULL;

	m_vElemdAFct[id] = NULL;
	m_vElemdAExplFct[id] = NULL;
	m_vElemdMFct[id] = NULL;

	m_vElemRHSFct[id] = NULL;

	m_vPrepareErrEstElemLoopFct[id] = NULL;
	m_vPrepareErrEstElemFct[id] = NULL;
	m_vElemComputeErrEstAFct[id] = NULL;
	m_vElemComputeErrEstMFct[id] = NULL;
	m_vElemComputeErrEstRhsFct[id] = NULL;
	m_vFinishErrEstElemLoopFct[id] = NULL;
}


template <typename TDomain>
void IElemDisc<TDomain>::clear_add_fct()
{
	for(size_t i = 0; i < NUM_REFERENCE_OBJECTS; ++i)
		clear_add_fct((ReferenceObjectID) i);
}


template <typename TDomain>
void IElemDisc<TDomain>::set_default_add_fct()
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

		m_vPrepareErrEstElemLoopFct[i] = &T::prep_err_est_elem_loop;
		m_vPrepareErrEstElemFct[i] = &T::prep_err_est_elem;
		m_vElemComputeErrEstAFct[i] = &T::compute_err_est_A_elem;
		m_vElemComputeErrEstMFct[i] = &T::compute_err_est_M_elem;
		m_vElemComputeErrEstRhsFct[i] = &T::compute_err_est_rhs_elem;
		m_vFinishErrEstElemLoopFct[i] = &T::fsh_err_est_elem_loop;
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
	if (this->m_vPrepareTimestepElemFct[m_id] != NULL)
		(this->*(m_vPrepareTimestepElemFct[m_id]))(time, u, elem, vCornerCoords);
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
	UG_ASSERT(m_vPrepareElemFct[m_id]!=NULL, "ElemDisc method prepare_elem missing.");
	(this->*(m_vPrepareElemFct[m_id]))(u, elem, vCornerCoords);
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
	if (this->m_vFinishTimestepElemFct[m_id] != NULL)
		(this->*(m_vFinishTimestepElemFct[m_id]))(time, u, elem, vCornerCoords);
}

template <typename TDomain>
void IElemDisc<TDomain>::
do_prep_elem_loop(const ReferenceObjectID roid, const int si)
{
//	set id and disc part (this checks the assemble functions)
	this->set_roid(roid, si);

//	remove positions in currently registered imports
	for(size_t i = 0; i < m_vIImport.size(); ++i)
		m_vIImport[i]->clear_ips();

//	call prep_elem_loop (this may set ip-series to imports)
	//	call assembling routine
	UG_ASSERT(m_vPrepareElemLoopFct[m_id]!=NULL, "ElemDisc method prepare_elem_loop missing.");
	(this->*m_vPrepareElemLoopFct[m_id])(roid, si);

//	set roid in imports (for evaluation function)
	for(size_t i = 0; i < m_vIImport.size(); ++i)
		m_vIImport[i]->set_roid(roid);
}

template <typename TDomain>
void IElemDisc<TDomain>::
do_fsh_elem_loop()
{
//	call finish
	UG_ASSERT(m_vFinishElemLoopFct[m_id]!=NULL, "ElemDisc method finish_elem_loop missing.");
	(this->*m_vFinishElemLoopFct[m_id])();

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
	UG_ASSERT(m_vElemJAFct[m_id]!=NULL, "ElemDisc method add_jac_A missing.");
	(this->*m_vElemJAFct[m_id])(J, u, elem, vCornerCoords);
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
	UG_ASSERT(m_vElemJMFct[m_id]!=NULL, "ElemDisc method add_jac_M missing.");
	(this->*m_vElemJMFct[m_id])(J, u, elem, vCornerCoords);
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
	UG_ASSERT(m_vElemdAFct[m_id]!=NULL, "ElemDisc method add_def_A missing.");
	(this->*m_vElemdAFct[m_id])(d, u, elem, vCornerCoords);
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
	if(this->m_vElemdAExplFct[m_id] != NULL)
		(this->*m_vElemdAExplFct[m_id])(d, u, elem, vCornerCoords);
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
	UG_ASSERT(m_vElemdMFct[m_id]!=NULL, "ElemDisc method add_def_M missing.");
	(this->*m_vElemdMFct[m_id])(d, u, elem, vCornerCoords);
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
	UG_ASSERT(m_vElemRHSFct[m_id]!=NULL, "ElemDisc method add_rhs missing.");
	(this->*m_vElemRHSFct[m_id])(rhs, elem, vCornerCoords);
}

template <typename TDomain>
void IElemDisc<TDomain>::
do_prep_err_est_elem_loop(const ReferenceObjectID roid, const int si)
{
// do nothing if error estimation is turned off
	if (!err_est_enabled()) return;

//	set id and disc part (this checks the assemble functions)
	this->set_roid(roid, si);

//	remove positions in currently registered imports
	for(std::size_t i = 0; i < m_vIImport.size(); ++i)
		m_vIImport[i]->clear_ips();

//	call prep_elem_loop (this may set ip-series to imports)
	UG_ASSERT(m_vPrepareErrEstElemLoopFct[m_id]!=NULL, "ElemDisc method prep_err_est_elem_loop missing.");
	(this->*m_vPrepareErrEstElemLoopFct[m_id])(roid, si);

//	set roid in imports (for evaluation function)
	for(std::size_t i = 0; i < m_vIImport.size(); ++i)
		m_vIImport[i]->set_roid(roid);
}

template <typename TDomain>
void IElemDisc<TDomain>::
do_prep_err_est_elem(LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	// do nothing if error estimation is turned off
	if (!err_est_enabled()) return;

	//	access by map
	u.access_by_map(map());
	if (local_time_series_needed())
		m_pLocalVectorTimeSeries->access_by_map(map());

	//	call assembling routine
	UG_ASSERT(m_vPrepareErrEstElemFct[m_id]!=NULL, "ElemDisc method prepare_err_est_elem missing.");
	(this->*(m_vPrepareErrEstElemFct[m_id]))(u, elem, vCornerCoords);
}

template <typename TDomain>
void IElemDisc<TDomain>::
do_compute_err_est_A_elem(LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[], const number& scale)
{
	// do nothing if error estimation is turned off
	if (!err_est_enabled()) return;

	//	access by map
	u.access_by_map(map());
	if (local_time_series_needed())
		m_pLocalVectorTimeSeries->access_by_map(map());

	//	call assembling routine
	if (this->m_vElemComputeErrEstAFct[m_id] != NULL)
		(this->*(m_vElemComputeErrEstAFct[m_id]))(u, elem, vCornerCoords, scale);
}

template <typename TDomain>
void IElemDisc<TDomain>::
do_compute_err_est_M_elem(LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[], const number& scale)
{
	// do nothing if error estimation is turned off
	if (!err_est_enabled()) return;

	//	access by map
	u.access_by_map(map());
	if(local_time_series_needed())
		m_pLocalVectorTimeSeries->access_by_map(map());

	//	call assembling routine
	if(this->m_vElemComputeErrEstMFct[m_id] != NULL)
		(this->*(m_vElemComputeErrEstMFct[m_id]))(u, elem, vCornerCoords, scale);
}

template <typename TDomain>
void IElemDisc<TDomain>::
do_compute_err_est_rhs_elem(GridObject* elem, const MathVector<dim> vCornerCoords[], const number& scale)
{
	// do nothing if error estimation is turned off
	if (!err_est_enabled()) return;

	if(local_time_series_needed())
		m_pLocalVectorTimeSeries->access_by_map(map());

	//	call assembling routine
	if(this->m_vElemComputeErrEstRhsFct[m_id] != NULL)
		(this->*(m_vElemComputeErrEstRhsFct[m_id]))(elem, vCornerCoords, scale);
}

template <typename TDomain>
void IElemDisc<TDomain>::
do_fsh_err_est_elem_loop()
{
// do nothing if error estimation is turned off
	if (!err_est_enabled()) return;

//	call finish
	if (this->m_vFinishErrEstElemLoopFct[m_id] != NULL)
		(this->*m_vFinishErrEstElemLoopFct[m_id])();

//	remove positions in currently registered imports
	for(size_t i = 0; i < m_vIImport.size(); ++i)
		m_vIImport[i]->clear_ips();
}



////////////////////////////////////////////////////////////////////////////////
//	virtual assembling functions default
////////////////////////////////////////////////////////////////////////////////

inline void ThrowMissingVirtualMethod(const char* method, const ReferenceObjectID roid){
	UG_THROW("ElemDisc: No override for the essential assembling function "
			 "'"<<method<<"' for " << roid << " implemented!");
}
inline void ThrowMissingVirtualMethod(const char* method){
	UG_THROW("ElemDisc: No override for the essential assembling function "
			 "'"<<method<<"' implemented!");
}

template <typename TDomain>
void IElemDisc<TDomain>::
prep_timestep_elem(const number time, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	// ThrowMissingVirtualMethod("prep_timestep_elem", elem->reference_object_id ());
}

template <typename TDomain>
void IElemDisc<TDomain>::
prep_elem(const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	ThrowMissingVirtualMethod("prep_elem", elem->reference_object_id ());
}

template <typename TDomain>
void IElemDisc<TDomain>::
fsh_timestep_elem(const number time, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	// ThrowMissingVirtualMethod("fsh_timestep_elem", elem->reference_object_id ());
}

template <typename TDomain>
void IElemDisc<TDomain>::
prep_elem_loop(const ReferenceObjectID roid, const int si)
{
	ThrowMissingVirtualMethod("prep_elem_loop", roid);
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
	ThrowMissingVirtualMethod("add_jac_A_elem", elem->reference_object_id ());
}

template <typename TDomain>
void IElemDisc<TDomain>::
add_jac_M_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	ThrowMissingVirtualMethod("add_jac_M_elem", elem->reference_object_id ());
}

template <typename TDomain>
void IElemDisc<TDomain>::
add_def_A_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	ThrowMissingVirtualMethod("add_def_A_elem", elem->reference_object_id ());
}

template <typename TDomain>
void IElemDisc<TDomain>::
add_def_A_expl_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	// ThrowMissingVirtualMethod("add_def_A_expl_elem", elem->reference_object_id ());
}

template <typename TDomain>
void IElemDisc<TDomain>::
add_def_M_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	ThrowMissingVirtualMethod("add_def_M_elem", elem->reference_object_id ());
}

template <typename TDomain>
void IElemDisc<TDomain>::
add_rhs_elem(LocalVector& rhs, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	ThrowMissingVirtualMethod("add_rhs_elem", elem->reference_object_id ());
}

template <typename TDomain>
void IElemDisc<TDomain>::
prep_err_est_elem_loop(const ReferenceObjectID roid, const int si)
{
	// ThrowMissingVirtualMethod("prep_err_est_elem_loop", roi);
}

template <typename TDomain>
void IElemDisc<TDomain>::
prep_err_est_elem(const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	//ThrowMissingVirtualMethod("prep_elem", elem->reference_object_id ());
}

template <typename TDomain>
void IElemDisc<TDomain>::
compute_err_est_A_elem(const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[], const number& scale)
{
	// ThrowMissingVirtualMethod("compute_err_est_elem", elem->reference_object_id ());
}

template <typename TDomain>
void IElemDisc<TDomain>::
compute_err_est_M_elem(const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[], const number& scale)
{
	// ThrowMissingVirtualMethod("compute_err_est_elem", elem->reference_object_id ());
}

template <typename TDomain>
void IElemDisc<TDomain>::
compute_err_est_rhs_elem(GridObject* elem, const MathVector<dim> vCornerCoords[], const number& scale)
{
	// ThrowMissingVirtualMethod("compute_err_est_elem", elem->reference_object_id ());
}

template <typename TDomain>
void IElemDisc<TDomain>::
fsh_err_est_elem_loop()
{
	// ThrowMissingVirtualMethod("fsh_err_est_elem_loop");
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

