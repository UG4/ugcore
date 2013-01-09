/*
 * neumann_boundary_common.h
 *
 *  Created on: 14.10.2010
 *      Author: andreasvogel
 */

#include "neumann_boundary.h"
#include "common/util/provider.h"
#include "lib_disc/common/groups_util.h"
#include "lib_disc/spatial_disc/user_data/const_user_data.h"

#ifdef UG_FOR_LUA
#include "bindings/lua/lua_user_data.h"
#endif

namespace ug{

////////////////////////////////////////////////////////////////////////////////
// User Data
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
void NeumannBoundary<TDomain>::
extract_data(Data& userData, FunctionGroup& commonFctGrp, std::string& fctNames)
{
//	create Function Group and Subset Group
	FunctionGroup functionGroup;

//	convert strings
	try{
		userData.ssGrp = this->approx_space()->subset_grp_by_name(userData.ssNames.c_str());
	}UG_CATCH_THROW("NeumannBoundary:extract_data':"
					" Subsets '"<<userData.ssNames<<"' not"
					" all contained in ApproximationSpace.");

	try{
		functionGroup = this->approx_space()->fct_grp_by_name(userData.fctName.c_str());
	}UG_CATCH_THROW("NeumannBoundary:extract_data':"
					" Functions '"<<userData.fctName<<"' not"
					" all contained in ApproximationSpace.");

//	check that only one function given
	if(functionGroup.size() != 1)
		UG_THROW("NeumannBoundary:extract_data: Only one function allowed"
						" per neumann value, but passed: " << userData.fctName);

//	get function
	const size_t fct = functionGroup[0];

// 	check if function exist
	if(fct >= this->function_pattern().num_fct())
		UG_THROW("NeumannBoundary:extract_data: Function "<< fct <<
					   " does not exist in pattern.");

//	add to common fct group if not already contained
	if(!commonFctGrp.contains(fct))
	{
		commonFctGrp.add(fct);

	//	build string of functions
		if(!fctNames.empty()) fctNames.append(",");
		fctNames.append(userData.fctName.c_str());
	}

//	set local fct id
	userData.locFct = commonFctGrp.local_index(fct);

//	check subsets and add referenze to data to each segment
	const ISubsetHandler& rSH = *this->function_pattern().subset_handler();
	for(size_t s = 0; s < userData.ssGrp.size(); ++s)
	{
	//	get subset index
		const int si = userData.ssGrp[s];

	// 	check that function is defined for segment
		if(!this->function_pattern().is_def_in_subset(fct, si))
			UG_THROW("NeumannBoundary:extract_data: Function "<<fct<<
						   " not defined on subset "<<si<<".");

	//	check that subsetIndex is valid
		if(si < 0 || si >= rSH.num_subsets())
			UG_THROW("NeumannBoundary:extract_data: Invalid subset "
					"Index " << si <<
					". (Valid is 0, .. , " << rSH.num_subsets() <<").");
	}
}

template<typename TDomain>
void NeumannBoundary<TDomain>::
extract_data()
{
//	a common function group
	FunctionGroup commonFctGrp(this->function_pattern());

//	string of functions
	std::string fctNames;

	for(size_t i = 0; i < m_vNumberData.size(); ++i)
		extract_data(m_vNumberData[i], commonFctGrp, fctNames);
	for(size_t i = 0; i < m_vBNDNumberData.size(); ++i)
		extract_data(m_vBNDNumberData[i], commonFctGrp, fctNames);
	for(size_t i = 0; i < m_vVectorData.size(); ++i)
		extract_data(m_vVectorData[i], commonFctGrp, fctNames);

//	set name of function
	this->set_functions(fctNames.c_str());
}

template<typename TDomain>
void NeumannBoundary<TDomain>::
add(SmartPtr<UserData<number, dim> > data, const char* function, const char* subsets)
{
	m_vNumberData.push_back(NumberData(data, function, subsets));
	if(this->fct_pattern_set()) extract_data();
}

template<typename TDomain>
void NeumannBoundary<TDomain>::
add(SmartPtr<UserData<number, dim, bool> > user, const char* function, const char* subsets)
{
	m_vBNDNumberData.push_back(BNDNumberData(user, function, subsets));
	if(this->fct_pattern_set()) extract_data();
}

template<typename TDomain>
void NeumannBoundary<TDomain>::
add(SmartPtr<UserData<MathVector<dim>, dim> > user, const char* function, const char* subsets)
{
	m_vVectorData.push_back(VectorData(user, function, subsets));
	if(this->fct_pattern_set()) extract_data();
}

template<typename TDomain>
void NeumannBoundary<TDomain>::
add(number val, const char* function, const char* subsets)
{
	SmartPtr<UserData<number, dim> > sp = CreateSmartPtr(new ConstUserNumber<dim>(val));
	add(sp, function, subsets);
}

template<typename TDomain>
void NeumannBoundary<TDomain>::
add(const std::vector<number>& val, const char* function, const char* subsets)
{
	SmartPtr<UserData<MathVector<dim>, dim> > sp = CreateSmartPtr(new ConstUserVector<dim>(val));
	add(sp, function, subsets);
}

#ifdef UG_FOR_LUA
template <typename TDomain>
void NeumannBoundary<TDomain>::
add(const char* name, const char* function, const char* subsets)
{
	if(LuaUserData<number, dim>::check_callback_returns(name)){
		SmartPtr<UserData<number, dim> > sp =
							LuaUserDataFactory<number, dim>::create(name);
		add(sp, function, subsets);
		return;
	}
	if(LuaUserData<number, dim, bool>::check_callback_returns(name)){
		SmartPtr<UserData<number, dim, bool> > sp =
				LuaUserDataFactory<number, dim, bool>::create(name);
		add(sp, function, subsets);
		return;
	}
	if(LuaUserData<MathVector<dim>, dim>::check_callback_returns(name)){
		SmartPtr<UserData<MathVector<dim>, dim> > sp =
				LuaUserDataFactory<MathVector<dim>, dim>::create(name);
		add(sp, function, subsets);
		return;
	}

//	no match found
	if(!CheckLuaCallbackName(name))
		UG_THROW("NeumannBoundary::add: Lua-Callback with name '"<<name<<
		               "' does not exist.");

//	name exists but wrong signature
	UG_THROW("NeumannBoundary::add: Cannot find matching callback "
					"signature. Use one of:\n"
					"a) Number - Callback\n"
					<< (LuaUserData<number, dim>::signature()) << "\n" <<
					"b) Conditional Number - Callback\n"
					<< (LuaUserData<number, dim, bool>::signature()) << "\n" <<
					"c) "<<dim<<"d Vector - Callback\n"
					<< (LuaUserData<MathVector<dim>, dim>::signature()));
}
#endif

////////////////////////////////////////////////////////////////////////////////
// Number Data
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void NeumannBoundary<TDomain>::NumberData::
extract_bip_fv1(const TFVGeom& geo)
{
	typedef typename TFVGeom::BF BF;
	vLocIP.clear();
	vGloIP.clear();
	for(size_t s = 0; s < this->ssGrp.size(); s++)
	{
		const int si = this->ssGrp[s];
		const std::vector<BF>& vBF = geo.bf(si);
		for(size_t i = 0; i < vBF.size(); ++i)
		{
			const BF& bf = vBF[i];
			vLocIP.push_back(bf.local_ip());
			vGloIP.push_back(bf.global_ip());
		}
	}

	import.set_local_ips(&vLocIP[0], vLocIP.size());
	import.set_global_ips(&vGloIP[0], vGloIP.size());
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void NeumannBoundary<TDomain>::NumberData::
extract_bip_fvho(const TFVGeom& geo)
{
	typedef typename TFVGeom::BF BF;
	vLocIP.clear();
	vGloIP.clear();
	for(size_t s = 0; s < this->ssGrp.size(); s++)
	{
		const int si = this->ssGrp[s];
		const std::vector<BF>& vBF = geo.bf(si);
		for(size_t i = 0; i < vBF.size(); ++i)
		{
			const BF& bf = vBF[i];
			for(size_t ip = 0; ip < bf.num_ip(); ++ip){
				vLocIP.push_back(bf.local_ip(ip));
				vGloIP.push_back(bf.global_ip(ip));
			}
		}
	}

	import.set_local_ips(&vLocIP[0], vLocIP.size());
	import.set_global_ips(&vGloIP[0], vGloIP.size());
}

////////////////////////////////////////////////////////////////////////////////
//	Constructor
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
NeumannBoundary<TDomain>::NeumannBoundary(const char* subsets)
 :IDomainElemDisc<TDomain>("", subsets)
{
//	set defaults
	m_order = 1;
	m_discScheme = "fv1";

//	update assemble functions
	set_ass_funcs();
}

template<typename TDomain>
NeumannBoundary<TDomain>::NeumannBoundary(const std::vector<std::string>& vSubset)
 :IDomainElemDisc<TDomain>(std::vector<std::string>(), vSubset)
{
//	set defaults
	m_order = 1;
	m_discScheme = "fv1";

//	update assemble functions
	set_ass_funcs();
}


///	type of trial space for each function used
template<typename TDomain>
bool NeumannBoundary<TDomain>::request_finite_element_id(const std::vector<LFEID>& vLfeID)
{
//	check that Lagrange space
	if(vLfeID[0].type() != LFEID::LAGRANGE)
	{
		UG_LOG("ERROR in 'NeumannBoundary::request_finite_element_id':"
			" Lagrange trial space needed.\n");
		return false;
	}

//	for fv1 only 1st order
	if(m_discScheme == "fv1" && vLfeID[0].order() != 1)
	{
		UG_LOG("ERROR in 'NeumannBoundary::request_finite_element_id':"
				" FV1 Scheme only implemented for 1st order.\n");
		return false;
	}

//	check that not ADAPTIVE
	if(vLfeID[0].order() < 1)
	{
		UG_LOG("ERROR in 'ConvectionDiffusion::request_finite_element_id':"
				" Adaptive or invalid order not implemented.\n");
		return false;
	}

//	remember lfeID;
	m_lfeID = vLfeID[0];

//	set order
	m_order = vLfeID[0].order();

//	update assemble functions
	set_ass_funcs();

//	is supported
	return true;
}

template<typename TDomain>
void NeumannBoundary<TDomain>::set_disc_scheme(const char* c_scheme)
{
//	convert to string
	std::string scheme = c_scheme;

//	check
	if(scheme != std::string("fv1") &&
	   scheme != std::string("fv"))
	{
		UG_THROW("NeumannBoundary: Only 'fv', 'fv1' supported.");
	}

//	remember
	m_discScheme = scheme;

//	update assemble functions
	set_ass_funcs();
}

template<typename TDomain>
void NeumannBoundary<TDomain>::set_ass_funcs()
{
//	switch, which assemble functions to use; both supported.
	if(m_discScheme == "fv1") register_all_fv1_funcs(false);
	else if(m_discScheme == "fv") register_all_fvho_funcs(m_order);
	else UG_THROW("NeumannBoundary: Disc Scheme '"<<m_discScheme<<"' not recognized.");
}

///	switches between non-regular and regular grids
template<typename TDomain>
bool
NeumannBoundary<TDomain>::
request_non_regular_grid(bool bNonRegular)
{
//	switch, which assemble functions to use.
	if(bNonRegular)
	{
		UG_LOG("ERROR in 'FVNeumannBoundaryElemDisc::request_non_regular_grid':"
				" Non-regular grid not implemented.\n");
		return false;
	}

//	this disc supports regular grids
	return true;
}

} // namespace ug

