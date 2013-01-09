/*
 * neumann_boundary_fv1.cpp
 *
 *  Created on: 14.10.2010
 *      Author: andreasvogel
 */

#include "neumann_boundary.h"
#include "common/util/provider.h"
#include "lib_disc/common/groups_util.h"
#include "lib_disc/spatial_disc/disc_util/finite_volume_geometry.h"
#include "lib_disc/spatial_disc/user_data/const_user_data.h"

#ifdef UG_FOR_LUA
#include "bindings/lua/lua_user_data.h"
#endif

namespace ug{

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


template<typename TDomain>
template<typename TElem, template <class Elem, int  Dim> class TFVGeom>
void NeumannBoundary<TDomain>::
prep_elem_loop_fv1()
{
//	register subsetIndex at Geometry
	static TFVGeom<TElem, dim>& geo = Provider<TFVGeom<TElem,dim> >::get();

//	request subset indices as boundary subset. This will force the
//	creation of boundary subsets when calling geo.update

	for(size_t i = 0; i < m_vNumberData.size(); ++i){
		for(size_t s = 0; s < m_vNumberData[i].ssGrp.size(); ++s){
			const int si = m_vNumberData[i].ssGrp[s];
			geo.add_boundary_subset(si);
		}
	}
	for(size_t i = 0; i < m_vBNDNumberData.size(); ++i){
		for(size_t s = 0; s < m_vBNDNumberData[i].ssGrp.size(); ++s){
			const int si = m_vBNDNumberData[i].ssGrp[s];
			geo.add_boundary_subset(si);
		}
	}
	for(size_t i = 0; i < m_vVectorData.size(); ++i){
		for(size_t s = 0; s < m_vVectorData[i].ssGrp.size(); ++s){
			const int si = m_vVectorData[i].ssGrp[s];
			geo.add_boundary_subset(si);
		}
	}

//	clear imports, since we will set them afterwards
	this->clear_imports();

	typedef typename NeumannBoundary<TDomain>::NumberData T;
	ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;

//	set lin defect fct for imports
	for(size_t data = 0; data < m_vNumberData.size(); ++data)
	{
		m_vNumberData[data].import.set_fct(id,
		                                   &m_vNumberData[data],
		                                   &NumberData::template lin_def_fv1<TElem, TFVGeom>);

		this->register_import(m_vNumberData[data].import);
		m_vNumberData[data].import.set_rhs_part();
	}
}

template<typename TDomain>
template<typename TElem, template <class Elem, int  Dim> class TFVGeom>
void NeumannBoundary<TDomain>::
finish_elem_loop_fv1()
{
//	remove subsetIndex from Geometry
	static TFVGeom<TElem, dim>& geo = Provider<TFVGeom<TElem,dim> >::get();


//	unrequest subset indices as boundary subset. This will force the
//	creation of boundary subsets when calling geo.update

	for(size_t i = 0; i < m_vNumberData.size(); ++i){
		for(size_t s = 0; s < m_vNumberData[i].ssGrp.size(); ++s){
			const int si = m_vNumberData[i].ssGrp[s];
			geo.remove_boundary_subset(si);
		}
	}

	for(size_t i = 0; i < m_vBNDNumberData.size(); ++i){
		for(size_t s = 0; s < m_vBNDNumberData[i].ssGrp.size(); ++s){
			const int si = m_vBNDNumberData[i].ssGrp[s];
			geo.remove_boundary_subset(si);
		}
	}

	for(size_t i = 0; i < m_vVectorData.size(); ++i){
		for(size_t s = 0; s < m_vVectorData[i].ssGrp.size(); ++s){
			const int si = m_vVectorData[i].ssGrp[s];
			geo.remove_boundary_subset(si);
		}
	}
}

template<typename TDomain>
template<typename TElem, template <class Elem, int  Dim> class TFVGeom>
void NeumannBoundary<TDomain>::
prep_elem_fv1(TElem* elem, const LocalVector& u)
{
//  update Geometry for this element
	static TFVGeom<TElem, dim>& geo = Provider<TFVGeom<TElem,dim> >::get();
	if(!geo.update(elem,
	               this->template element_corners<TElem>(elem),
	               &(this->subset_handler())))
		UG_THROW("NeumannBoundary::prep_elem_fv1: "
						"Cannot update Finite Volume Geometry.");

	for(size_t i = 0; i < m_vNumberData.size(); ++i)
		m_vNumberData[i].template extract_bip<TElem, TFVGeom>(geo);
}

template<typename TDomain>
template<typename TElem, template <class Elem, int  Dim> class TFVGeom>
void NeumannBoundary<TDomain>::
add_rhs_elem_fv1(LocalVector& d)
{
	const static TFVGeom<TElem, dim>& geo = Provider<TFVGeom<TElem,dim> >::get();
	typedef typename TFVGeom<TElem, dim>::BF BF;

//	Number Data
	for(size_t data = 0; data < m_vNumberData.size(); ++data){
		size_t ip = 0;
		for(size_t s = 0; s < m_vNumberData[data].ssGrp.size(); ++s){
			const int si = m_vNumberData[data].ssGrp[s];
			const std::vector<BF>& vBF = geo.bf(si);

			for(size_t i = 0; i < vBF.size(); ++i, ++ip){
				const int co = vBF[i].node_id();
				d(m_vNumberData[data].locFct, co) -= m_vNumberData[data].import[ip]
				                                    * vBF[i].volume();
			}
		}
	}

//	conditional Number Data
	for(size_t data = 0; data < m_vBNDNumberData.size(); ++data){
		for(size_t s = 0; s < m_vBNDNumberData[data].ssGrp.size(); ++s)	{
			const int si = m_vBNDNumberData[data].ssGrp[s];
			const std::vector<BF>& vBF = geo.bf(si);

			for(size_t i = 0; i < vBF.size(); ++i){
				number val = 0.0;
				if(!(*m_vBNDNumberData[data].functor)(val, vBF[i].global_ip(), this->time(), si))
					continue;

				const int co = vBF[i].node_id();
				d(m_vBNDNumberData[data].locFct, co) -= val * vBF[i].volume();
			}
		}
	}

// 	vector data
	for(size_t data = 0; data < m_vVectorData.size(); ++data){
		for(size_t s = 0; s < m_vVectorData[data].ssGrp.size(); ++s){
			const int si = m_vVectorData[data].ssGrp[s];
			const std::vector<BF>& vBF = geo.bf(si);

			for(size_t i = 0; i < vBF.size(); ++i){
				MathVector<dim> val;
				(*m_vVectorData[data].functor)(val, vBF[i].global_ip(), this->time(), si);

				const int co = vBF[i].node_id();
				d(m_vVectorData[data].locFct, co) -= VecDot(val, vBF[i].normal());
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////////////
// Number Data
////////////////////////////////////////////////////////////////////////////////


template<typename TDomain>
template<typename TElem, template <class Elem, int  Dim> class TFVGeom>
void NeumannBoundary<TDomain>::NumberData::
extract_bip(const TFVGeom<TElem,dim>& geo)
{
	typedef typename TFVGeom<TElem, dim>::BF BF;
	vLocIP.clear();
	vGloIP.clear();
	for(size_t s = 0; s < this->ssGrp.size(); s++)
	{
		const int si = this->ssGrp[s];
		const std::vector<BF>& vBF = geo.bf(si);
		for(size_t i = 0; i < vBF.size(); ++i)
		{
			vLocIP.push_back(vBF[i].local_ip());
			vGloIP.push_back(vBF[i].global_ip());
		}
	}

	import.set_local_ips(&vLocIP[0], vLocIP.size());
	import.set_global_ips(&vGloIP[0], vGloIP.size());
}

template<typename TDomain>
template<typename TElem, template <class Elem, int  Dim> class TFVGeom>
void NeumannBoundary<TDomain>::NumberData::
lin_def_fv1(const LocalVector& u,
            std::vector<std::vector<number> > vvvLinDef[],
            const size_t nip)
{
//  get finite volume geometry
	const static TFVGeom<TElem, dim>& geo = Provider<TFVGeom<TElem,dim> >::get();
	typedef typename TFVGeom<TElem, dim>::BF BF;

	size_t ip = 0;
	for(size_t s = 0; s < this->ssGrp.size(); ++s)
	{
		const int si = this->ssGrp[s];
		const std::vector<BF>& vBF = geo.bf(si);
		for(size_t i = 0; i < vBF.size(); ++i){
			const int co = vBF[i].node_id();
			vvvLinDef[ip][this->locFct][co] -= vBF[i].volume();
		}
	}
}


////////////////////////////////////////////////////////////////////////////////
//	Constructor
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
NeumannBoundary<TDomain>::NeumannBoundary(const char* subsets)
 :IDomainElemDisc<TDomain>("", subsets)
{
	register_all_fv1_funcs(false);
}


///	type of trial space for each function used
template<typename TDomain>
bool
NeumannBoundary<TDomain>::
request_finite_element_id(const std::vector<LFEID>& vLfeID)
{
//	check that Lagrange 1st order
	for(size_t i = 0; i < vLfeID.size(); ++i)
		if(vLfeID[i] != LFEID(LFEID::LAGRANGE, 1)) return false;
	return true;
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

////////////////////////////////////////////////////////////////////////////////
//	register assemble functions
////////////////////////////////////////////////////////////////////////////////

// register for 1D
template<typename TDomain>
void
NeumannBoundary<TDomain>::
register_all_fv1_funcs(bool bHang)
{
//	get all grid element types in this dimension and below
	typedef typename domain_traits<dim>::DimElemList ElemList;

//	switch assemble functions
	if(!bHang) boost::mpl::for_each<ElemList>( RegisterFV1<FV1Geometry>(this) );
	else throw (UGError("Not implemented."));
}

template<typename TDomain>
template<typename TElem, template <class Elem, int WorldDim> class TFVGeom>
void
NeumannBoundary<TDomain>::
register_fv1_func()
{
	ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;
	typedef this_type T;

	this->enable_fast_ass_elem(true);
	this->set_prep_elem_loop_fct(id, &T::template prep_elem_loop_fv1<TElem, TFVGeom>);
	this->set_prep_elem_fct(	 id, &T::template prep_elem_fv1<TElem, TFVGeom>);
	this->set_fsh_elem_loop_fct( id, &T::template finish_elem_loop_fv1<TElem, TFVGeom>);
	this->set_ass_JA_elem_fct(	 id, &T::template add_JA_elem_fv1<TElem, TFVGeom>);
	this->set_ass_JM_elem_fct(	 id, &T::template add_JM_elem_fv1<TElem, TFVGeom>);
	this->set_ass_dA_elem_fct(	 id, &T::template add_dA_elem_fv1<TElem, TFVGeom>);
	this->set_ass_dM_elem_fct(	 id, &T::template add_dM_elem_fv1<TElem, TFVGeom>);
	this->set_ass_rhs_elem_fct(	 id, &T::template add_rhs_elem_fv1<TElem, TFVGeom>);
}

////////////////////////////////////////////////////////////////////////////////
//	explicit template instantiations
////////////////////////////////////////////////////////////////////////////////

template class NeumannBoundary<Domain1d>;
template class NeumannBoundary<Domain2d>;
template class NeumannBoundary<Domain3d>;

} // namespace ug

