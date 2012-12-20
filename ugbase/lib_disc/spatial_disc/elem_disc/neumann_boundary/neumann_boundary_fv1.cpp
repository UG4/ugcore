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
template <typename TUserData>
void NeumannBoundary<TDomain>::
extract_data(std::map<int, std::vector<TUserData*> >& mvUserDataBndSegment,
             std::vector<TUserData>& vUserData,
             FunctionGroup& commonFctGrp, std::string& fctNames)
{
//	clear all extracted data
	mvUserDataBndSegment.clear();

//	loop all scheduled data
	for(size_t i = 0; i < vUserData.size(); ++i)
	{
	//	create Function Group and Subset Group
		FunctionGroup functionGroup;

	//	convert strings
		try{
			vUserData[i].ssGrp = this->approx_space()->subset_grp_by_name(vUserData[i].ssNames.c_str());
		}UG_CATCH_THROW("NeumannBoundary:extract_data':"
						" Subsets '"<<vUserData[i].ssNames<<"' not"
						" all contained in ApproximationSpace.");

		try{
			functionGroup = this->approx_space()->fct_grp_by_name(vUserData[i].fctName.c_str());
		}UG_CATCH_THROW("NeumannBoundary:extract_data':"
						" Functions '"<<vUserData[i].fctName<<"' not"
						" all contained in ApproximationSpace.");

	//	check that only one function given
		if(functionGroup.size() != 1)
			UG_THROW("NeumannBoundary:extract_data: Only one function allowed"
							" per neumann value, but passed: " << vUserData[i].fctName);

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
			fctNames.append(vUserData[i].fctName.c_str());
		}

	//	set local fct id
		vUserData[i].locFct = commonFctGrp.local_index(fct);

	//	check subsets and add referenze to data to each segment
		const ISubsetHandler& rSH = *this->function_pattern().subset_handler();
		for(size_t si = 0; si < vUserData[i].ssGrp.size(); ++si)
		{
		//	get subset index
			const int subsetIndex = vUserData[i].ssGrp[si];

		// 	check that function is defined for segment
			if(!this->function_pattern().is_def_in_subset(fct, subsetIndex))
				UG_THROW("NeumannBoundary:extract_data: Function "<<fct<<
				               " not defined on subset "<<subsetIndex<<".");

		//	check that subsetIndex is valid
			if(subsetIndex < 0 || subsetIndex >= rSH.num_subsets())
				UG_THROW("NeumannBoundary:extract_data: Invalid subset "
						"Index " << subsetIndex <<
						". (Valid is 0, .. , " << rSH.num_subsets() <<").");

		//	add to segment
			mvUserDataBndSegment[subsetIndex].push_back(&vUserData[i]);
		}
	}
}

template<typename TDomain>
void NeumannBoundary<TDomain>::
extract_data()
{
//	a common function group
	FunctionGroup commonFctGrp;
	commonFctGrp.set_function_pattern(this->function_pattern());

//	string of functions
	std::string fctNames;

	extract_data(m_mNumberBndSegment, m_vNumberData, commonFctGrp, fctNames);
	extract_data(m_mBNDNumberBndSegment, m_vBNDNumberData, commonFctGrp, fctNames);
	extract_data(m_mVectorBndSegment, m_vVectorData, commonFctGrp, fctNames);

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
prepare_element_loop()
{
//	register subsetIndex at Geometry
	static TFVGeom<TElem, dim>& geo = Provider<TFVGeom<TElem,dim> >::get();

	for(size_t i = 0; i < m_vNumberData.size(); ++i)
	{
		for(size_t s = 0; s < m_vNumberData[i].ssGrp.size(); ++s)
		{
		//	get subset index
			const int bndSubset = m_vNumberData[i].ssGrp[s];

		//	request this subset index as boundary subset. This will force the
		//	creation of boundary subsets when calling geo.update
			geo.add_boundary_subset(bndSubset);
		}
	}

	{
	typename std::map<int, std::vector<BNDNumberData*> >::const_iterator subsetIter;
	for(subsetIter = m_mBNDNumberBndSegment.begin();
			subsetIter != m_mBNDNumberBndSegment.end(); ++subsetIter)
	{
	//	get subset index
		const int bndSubset = (*subsetIter).first;

	//	request this subset index as boundary subset. This will force the
	//	creation of boundary subsets when calling geo.update
		geo.add_boundary_subset(bndSubset);
	}
	}

	{
	typename std::map<int, std::vector<VectorData*> >::const_iterator subsetIter;
	for(subsetIter = m_mVectorBndSegment.begin();
			subsetIter != m_mVectorBndSegment.end(); ++subsetIter)
	{
	//	get subset index
		const int bndSubset = (*subsetIter).first;

	//	request this subset index as boundary subset. This will force the
	//	creation of boundary subsets when calling geo.update
		geo.add_boundary_subset(bndSubset);
	}
	}

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
finish_element_loop()
{
//	remove subsetIndex from Geometry
	static TFVGeom<TElem, dim>& geo = Provider<TFVGeom<TElem,dim> >::get();

	for(size_t i = 0; i < m_vNumberData.size(); ++i)
	{
		for(size_t s = 0; s < m_vNumberData[i].ssGrp.size(); ++s)
		{
		//	get subset index
			const int bndSubset = m_vNumberData[i].ssGrp[s];

		//	request this subset index as boundary subset. This will force the
		//	creation of boundary subsets when calling geo.update
			geo.remove_boundary_subset(bndSubset);
		}
	}

	{
	typename std::map<int, std::vector<BNDNumberData*> >::const_iterator subsetIter;
	for(subsetIter = m_mBNDNumberBndSegment.begin();
			subsetIter != m_mBNDNumberBndSegment.end(); ++subsetIter)
	{
	//	get subset index
		const int bndSubset = (*subsetIter).first;

	//	remove requested bnd subset
		geo.remove_boundary_subset(bndSubset);
	}
	}

	{
	typename std::map<int, std::vector<VectorData*> >::const_iterator subsetIter;
	for(subsetIter = m_mVectorBndSegment.begin();
			subsetIter != m_mVectorBndSegment.end(); ++subsetIter)
	{
	//	get subset index
		const int bndSubset = (*subsetIter).first;

	//	remove requested bnd subset
		geo.remove_boundary_subset(bndSubset);
	}
	}
}

template<typename TDomain>
template<typename TElem, template <class Elem, int  Dim> class TFVGeom>
void NeumannBoundary<TDomain>::
prepare_element(TElem* elem, const LocalVector& u)
{
//	get corners
	m_vCornerCoords = this->template element_corners<TElem>(elem);

//  update Geometry for this element
	static TFVGeom<TElem, dim>& geo = Provider<TFVGeom<TElem,dim> >::get();
	if(!geo.update(elem, &m_vCornerCoords[0], &(this->subset_handler())))
		UG_THROW("NeumannBoundary::prepare_element: "
						"Cannot update Finite Volume Geometry.");

	for(size_t i = 0; i < m_vNumberData.size(); ++i)
		m_vNumberData[i].template extract_bip<TElem, TFVGeom>(geo);
}

template<typename TDomain>
template<typename TElem, template <class Elem, int  Dim> class TFVGeom>
void NeumannBoundary<TDomain>::
ass_rhs_elem(LocalVector& d)
{
// 	get finite volume geometry
	const static TFVGeom<TElem, dim>& geo = Provider<TFVGeom<TElem,dim> >::get();
	typedef typename TFVGeom<TElem, dim>::BF BF;

	for(size_t data = 0; data < m_vNumberData.size(); ++data)
	{
		size_t ip = 0;
		for(size_t s = 0; s < m_vNumberData[data].ssGrp.size(); ++s)
		{
		//	get subset index
			const int bndSubset = m_vNumberData[data].ssGrp[s];

			if(geo.num_bf(bndSubset) == 0) continue;

			const std::vector<BF>& vBF = geo.bf(bndSubset);

			for(size_t i = 0; i < vBF.size(); ++i, ++ip)
			{
			// 	get associated node
				const int co = vBF[i].node_id();

			// 	add to local matrix
				d(m_vNumberData[data].locFct, co) -= m_vNumberData[data].import[ip]
				                                    * vBF[i].volume();
			}
		}
	}

// 	loop registered boundary segments
	if(!m_vBNDNumberData.empty())
	{
		typename std::map<int, std::vector<BNDNumberData*> >::const_iterator subsetIter;
		for(subsetIter = m_mBNDNumberBndSegment.begin(); subsetIter != m_mBNDNumberBndSegment.end(); ++subsetIter)
		{
		//	get subset index corresponding to boundary
			const int bndSubset = (*subsetIter).first;

		//	check for boundary faces
			if(geo.num_bf(bndSubset) == 0) continue;

		//	get evaluation function vector
			const std::vector<BNDNumberData*>& vSegmentData = (*subsetIter).second;

		//	get boundary faces
			const std::vector<BF>& vBF = geo.bf(bndSubset);

		// 	loop Boundary Faces
			for(size_t i = 0; i < vBF.size(); ++i)
			{
			// get current BF
				const BF& bf = vBF[i];

			//	loop functions, where neumann bnd is set for this bndSubset
				for(size_t fct = 0; fct < vSegmentData.size(); ++fct)
				{
				// 	get neumann value
					number val = 0.0;
					if(!(*vSegmentData[fct]->functor)(val, bf.global_ip(), this->time(), bndSubset))
						continue;

				// 	get associated node
					const int co = bf.node_id();

				// 	add to local matrix
					d(vSegmentData[fct]->locFct, co) -= val * bf.volume();
				}
			}
		}
	}

// 	loop registered boundary segments
	if(!m_vVectorData.empty())
	{
		typename std::map<int, std::vector<VectorData*> >::const_iterator subsetIter;
		for(subsetIter = m_mVectorBndSegment.begin(); subsetIter != m_mVectorBndSegment.end(); ++subsetIter)
		{
		//	get subset index corresponding to boundary
			const int bndSubset = (*subsetIter).first;

		//	check for boundary faces
			if(geo.num_bf(bndSubset) == 0) continue;

		//	get evaluation function vector
			const std::vector<VectorData*>& vSegmentData = (*subsetIter).second;

		//	get boundary faces
			const std::vector<BF>& vBF = geo.bf(bndSubset);

		// 	loop Boundary Faces
			for(size_t i = 0; i < vBF.size(); ++i)
			{
			// get current BF
				const BF& bf = vBF[i];

			//	loop functions, where neumann bnd is set for this bndSubset
				for(size_t fct = 0; fct < vSegmentData.size(); ++fct)
				{
				// 	get neumann value
					MathVector<dim> val;
					(*vSegmentData[fct]->functor)(val, bf.global_ip(), this->time(), bndSubset);

				// 	get associated node
					const int co = bf.node_id();

				// 	add to local matrix
					d(vSegmentData[fct]->locFct, co) -= VecDot(val, bf.normal());
				}
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
	for(size_t s = 0; s < ssGrp.size(); s++)
	{
		const int bndSubset = ssGrp[s];
		if(geo.num_bf(bndSubset) == 0) continue;
		const std::vector<BF>& vBF = geo.bf(bndSubset);
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
	for(size_t s = 0; s < ssGrp.size(); ++s)
	{
		const int bndSubset = ssGrp[s];
		if(geo.num_bf(bndSubset) == 0) continue;
		const std::vector<BF>& vBF = geo.bf(bndSubset);
		for(size_t i = 0; i < vBF.size(); ++i)
		{
		// 	get associated node
			const int co = vBF[i].node_id();

		// 	set lin defect
			vvvLinDef[ip][locFct][co] -= vBF[i].volume();
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
	m_mBNDNumberBndSegment.clear();
	m_mVectorBndSegment.clear();
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
	this->set_prep_elem_loop_fct(id, &T::template prepare_element_loop<TElem, TFVGeom>);
	this->set_prep_elem_fct(	 id, &T::template prepare_element<TElem, TFVGeom>);
	this->set_fsh_elem_loop_fct( id, &T::template finish_element_loop<TElem, TFVGeom>);
	this->set_ass_JA_elem_fct(		 id, &T::template ass_JA_elem<TElem, TFVGeom>);
	this->set_ass_JM_elem_fct(		 id, &T::template ass_JM_elem<TElem, TFVGeom>);
	this->set_ass_dA_elem_fct(		 id, &T::template ass_dA_elem<TElem, TFVGeom>);
	this->set_ass_dM_elem_fct(		 id, &T::template ass_dM_elem<TElem, TFVGeom>);
	this->set_ass_rhs_elem_fct(	 id, &T::template ass_rhs_elem<TElem, TFVGeom>);
}

////////////////////////////////////////////////////////////////////////////////
//	explicit template instantiations
////////////////////////////////////////////////////////////////////////////////

template class NeumannBoundary<Domain1d>;
template class NeumannBoundary<Domain2d>;
template class NeumannBoundary<Domain3d>;

} // namespace ug

