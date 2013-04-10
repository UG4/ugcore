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
	: 	m_pFctPattern(0), m_timePoint(0),
	  	m_pLocalVectorTimeSeries(NULL), m_bStationaryForced(false),
	  	m_bFastAssembleEnabled(false), m_id(ROID_UNKNOWN), m_spApproxSpace(NULL)
{
	set_functions(functions);
	set_subsets(subsets);
	clear_add_fct();
}

template <typename TDomain>
IElemDisc<TDomain>::IElemDisc(const std::vector<std::string>& vFct,
                              const std::vector<std::string>& vSubset)
	: 	m_pFctPattern(0), m_timePoint(0),
		m_pLocalVectorTimeSeries(NULL), m_bStationaryForced(false),
		m_bFastAssembleEnabled(false), m_id(ROID_UNKNOWN), m_spApproxSpace(NULL)
{
	set_functions(vFct);
	set_subsets(vSubset);
	clear_add_fct();
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
		m_vElemdAFct_explicit[i] = NULL;
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
void IElemDisc<TDomain>::set_function_pattern(const FunctionPattern& fctPatt)
{
	m_pFctPattern = &fctPatt;
	update_function_index_mapping();
}

template <typename TDomain>
void IElemDisc<TDomain>::update_function_index_mapping()
{
//	without fct pattern, cannot create mappings
	if(!m_pFctPattern) return;

//	create function group of this elem disc
	try{
		m_fctGrp.set_function_pattern(*m_pFctPattern);
		m_fctGrp.add(this->symb_fcts());
	}UG_CATCH_THROW("ElemDisc: Cannot find some symbolic Function Name.");

//	create a mapping between all functions and the function group of this
//	element disc.
	try{CreateFunctionIndexMapping(m_fctIndexMap, m_fctGrp, *m_pFctPattern);
	}UG_CATCH_THROW("ElemDisc: Cannot create Function Index Mapping.");

//	set function group at imports
	for(size_t i = 0; i < m_vIImport.size(); ++i){
		m_vIImport[i]->set_map(m_fctIndexMap);
	}

//	set function group at exports
	for(size_t i = 0; i < m_vIExport.size(); ++i){
		m_vIExport[i]->set_function_group(m_fctGrp);
		m_vIExport[i]->set_map(m_fctIndexMap);
	}
}

template <typename TDomain>
void IElemDisc<TDomain>::check_setup()
{
//	check that all functions are defined on chosen subsets
	SubsetGroup discSubsetGrp(m_pFctPattern->subset_handler(), m_vSubset);

//	check that all functions are defined on chosen subsets
	for(size_t fct = 0; fct < m_fctGrp.size(); ++fct){
		for(size_t si = 0; si < discSubsetGrp.size(); ++si){
			if(!m_pFctPattern->is_def_in_subset(m_fctGrp[fct], discSubsetGrp[si])){
				UG_LOG("WARNING in ElemDisc: symbolic Function "<< symb_fcts()[fct]
				 << " is not defined on subset "<< symb_subsets()[si]
				 << ". This may be senseful only in particular cases.\n");
			}
		}
	}

//	check correct number of functions
	if(m_fctGrp.size() != this->num_fct()){
		std::stringstream ss;
		ss << "ElemDisc requires "<< this->num_fct()<<" symbolic "
				"Function Name, but "<< m_fctGrp.size()<<" Functions "
				" specified: ";
		for(size_t f=0; f < symb_fcts().size(); ++f){
			if(f > 0) ss << ", "; ss << symb_fcts()[f];
		}
		UG_THROW(ss.str());
	}

//	request assembling for local finite element id
	std::vector<LFEID> vLfeID(m_fctGrp.size());
	for(size_t f = 0; f < vLfeID.size(); ++f)
		vLfeID[f] = m_fctGrp.local_finite_element_id(f);
	if(!(this->request_finite_element_id(vLfeID)))
	{
		std::stringstream ss;
		ss << "Elem Disc can not assemble the specified local finite element space set:";
		for(size_t f=0; f < symb_fcts().size(); ++f)
		{
			ss << "  Fct "<<f<<": '"<< symb_fcts()[f];
			ss << "' using "<< vLfeID[f];
		}
		UG_THROW(ss.str());
	}
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
void IElemDisc<TDomain>::register_export(SmartPtr<ICplUserData<dim> > Exp)
{
//	check that not already registered
	for(size_t i = 0; i < m_vIExport.size(); ++i)
		if(m_vIExport[i] == Exp)
			UG_THROW("Trying to register export twice.");

//	add it
	m_vIExport.push_back(Exp);

	update_function_index_mapping();
}

template <typename TDomain>
void IElemDisc<TDomain>::set_roid(ReferenceObjectID roid, int discType)
{
	m_id = roid;

	if(roid == ROID_UNKNOWN)
	{
		m_id = ROID_UNKNOWN;
		UG_THROW("Cannot assemble for RefID: "<<roid<<".");
	}

	if(m_vPrepareElemLoopFct[m_id]==NULL)
		UG_THROW("ElemDisc: Missing evaluation method 'prepare_elem_loop' for "<<roid);
	if(m_vPrepareElemFct[m_id]==NULL)
		UG_THROW("ElemDisc: Missing evaluation method 'prepare_elem' for "<<roid);
	if(m_vFinishElemLoopFct[m_id]==NULL)
		UG_THROW("ElemDisc: Missing evaluation method 'finish_elem_loop' for "<<roid);

	if(discType & MASS){
		if(m_vElemJMFct[m_id]==NULL)
			UG_THROW("ElemDisc: Missing evaluation method 'add_jac_M_elem' for "<<roid);
		if(m_vElemdMFct[m_id]==NULL)
			UG_THROW("ElemDisc: Missing evaluation method 'add_def_M_elem' for "<<roid);
	}
	if(discType & STIFF){
		if(m_vElemJAFct[m_id]==NULL)
			UG_THROW("ElemDisc: Missing evaluation method 'add_jac_A_elem' for "<<roid);
		if(m_vElemdAFct[m_id]==NULL)
			UG_THROW("ElemDisc: Missing evaluation method 'add_def_A_elem for' "<<roid);
	}
	if(discType & RHS){
		if(m_vElemRHSFct[m_id]==NULL)
			UG_THROW("ElemDisc: Missing evaluation method 'add_rhs_elem' for "<<roid);
	}
};

template <typename TDomain>
void IElemDisc<TDomain>::set_time_dependent(const LocalVectorTimeSeries& locTimeSeries,
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

