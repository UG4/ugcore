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

namespace ug{

template<typename TDomain>
template <typename TUserData>
void FV1NeumannBoundaryElemDisc<TDomain>::
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
		if(functionGroup.num_fct() != 1)
			UG_THROW_FATAL("NeumannBoundary:extract_data: Only one function allowed"
							" per neumann value, but passed: " << vUserData[i].fctName);

	//	get function
		const size_t fct = functionGroup[0];

	// 	check if function exist
		if(fct >= this->function_pattern().num_fct())
			UG_THROW_FATAL("NeumannBoundary:extract_data: Function "<< fct <<
			               " does not exist in pattern.");

	//	add to common fct group if not already contained
		if(!commonFctGrp.contains(fct)) commonFctGrp.add(fct);

	//	build string of functions
		if(!fctNames.empty()) fctNames.append(",");
		fctNames.append(vUserData[i].fctName.c_str());

	//	set local fct id
		vUserData[i].locFct = commonFctGrp.local_index(fct);

	//	check subsets and add referenze to data to each segment
		const ISubsetHandler& rSH = *this->function_pattern().subset_handler();
		for(size_t si = 0; si < vUserData[i].ssGrp.num_subsets(); ++si)
		{
		//	get subset index
			const int subsetIndex = vUserData[i].ssGrp[si];

		// 	check that function is defined for segment
			if(!this->function_pattern().is_def_in_subset(fct, subsetIndex))
				UG_THROW_FATAL("NeumannBoundary:extract_data: Function "<<fct<<
				               " not defined on subset "<<subsetIndex<<".");

		//	check that subsetIndex is valid
			if(subsetIndex < 0 || subsetIndex >= rSH.num_subsets())
				UG_THROW_FATAL("NeumannBoundary:extract_data: Invalid subset "
						"Index " << subsetIndex <<
						". (Valid is 0, .. , " << rSH.num_subsets() <<").");

		//	add to segment
			mvUserDataBndSegment[subsetIndex].push_back(&vUserData[i]);
		}
	}
}

template<typename TDomain>
void FV1NeumannBoundaryElemDisc<TDomain>::
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
void FV1NeumannBoundaryElemDisc<TDomain>::
add(SmartPtr<IPData<number, dim> > data, const char* function, const char* subsets)
{
	m_vNumberData.push_back(NumberData(data, function, subsets));

	if(this->fct_pattern_set()) extract_data();
}

template<typename TDomain>
void FV1NeumannBoundaryElemDisc<TDomain>::
add(BNDNumberFunctor& user, const char* function, const char* subsets)
{
	m_vBNDNumberData.push_back(BNDNumberData(user, function, subsets));

	if(this->fct_pattern_set()) extract_data();
}

template<typename TDomain>
void FV1NeumannBoundaryElemDisc<TDomain>::
add(VectorFunctor& user, const char* function, const char* subsets)
{
	m_vVectorData.push_back(VectorData(user, function, subsets));

	if(this->fct_pattern_set()) extract_data();
}

template<typename TDomain>
template<typename TElem, template <class Elem, int  Dim> class TFVGeom>
inline
bool
FV1NeumannBoundaryElemDisc<TDomain>::
prepare_element_loop()
{
//	register subsetIndex at Geometry
	static TFVGeom<TElem, dim>& geo = Provider<TFVGeom<TElem,dim> >::get();

	for(size_t i = 0; i < m_vNumberData.size(); ++i)
	{
		for(size_t s = 0; s < m_vNumberData[i].ssGrp.num_subsets(); ++s)
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

	typedef typename FV1NeumannBoundaryElemDisc<TDomain>::NumberData T;
	ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;

//	set lin defect fct for imports
	for(size_t data = 0; data < m_vNumberData.size(); ++data)
	{
		m_vNumberData[data].import.set_fct(id,
		                                   &m_vNumberData[data],
		                                   &NumberData::template lin_def_fv1<TElem, TFVGeom>);

		this->register_import(m_vNumberData[data].import);
		m_vNumberData[data].import.set_rhs_part(true);
	}

//	we're done
	return true;
}

template<typename TDomain>
template<typename TElem, template <class Elem, int  Dim> class TFVGeom>
bool
FV1NeumannBoundaryElemDisc<TDomain>::
finish_element_loop()
{
//	remove subsetIndex from Geometry
	static TFVGeom<TElem, dim>& geo = Provider<TFVGeom<TElem,dim> >::get();

	for(size_t i = 0; i < m_vNumberData.size(); ++i)
	{
		for(size_t s = 0; s < m_vNumberData[i].ssGrp.num_subsets(); ++s)
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

//	we're done
	return true;
}

template<typename TDomain>

template<typename TElem, template <class Elem, int  Dim> class TFVGeom>
inline
bool
FV1NeumannBoundaryElemDisc<TDomain>::
prepare_element(TElem* elem, const LocalVector& u)
{
//	get corners
	m_vCornerCoords = this->template element_corners<TElem>(elem);

//  update Geometry for this element
	static TFVGeom<TElem, dim>& geo = Provider<TFVGeom<TElem,dim> >::get();
	if(!geo.update(elem, &m_vCornerCoords[0], &(this->subset_handler())))
	{
		UG_LOG("ERROR in 'FVNeumannBoundaryElemDisc::prepare_element': "
				"Cannot update Finite Volume Geometry.\n");
		return false;
	}

	for(size_t i = 0; i < m_vNumberData.size(); ++i)
		m_vNumberData[i].template extract_bip<TElem, TFVGeom>(geo);

//	we're done
	return true;
}

template<typename TDomain>
template<typename TElem, template <class Elem, int  Dim> class TFVGeom>
inline
bool
FV1NeumannBoundaryElemDisc<TDomain>::
assemble_JA(LocalMatrix& J, const LocalVector& u)
{
	// we're done
	return true;
}


template<typename TDomain>
template<typename TElem, template <class Elem, int  Dim> class TFVGeom>
inline
bool
FV1NeumannBoundaryElemDisc<TDomain>::
assemble_JM(LocalMatrix& J, const LocalVector& u)
{
	// we're done
	return true;
}


template<typename TDomain>
template<typename TElem, template <class Elem, int  Dim> class TFVGeom>
inline
bool
FV1NeumannBoundaryElemDisc<TDomain>::
assemble_A(LocalVector& d, const LocalVector& u)
{
	// we're done
	return true;
}


template<typename TDomain>
template<typename TElem, template <class Elem, int  Dim> class TFVGeom>
inline
bool
FV1NeumannBoundaryElemDisc<TDomain>::
assemble_M(LocalVector& d, const LocalVector& u)
{
	// we're done
	return true;
}


template<typename TDomain>
template<typename TElem, template <class Elem, int  Dim> class TFVGeom>
inline
bool
FV1NeumannBoundaryElemDisc<TDomain>::
assemble_f(LocalVector& d)
{
// 	get finite volume geometry
	const static TFVGeom<TElem, dim>& geo = Provider<TFVGeom<TElem,dim> >::get();
	typedef typename TFVGeom<TElem, dim>::BF BF;

	size_t ip = 0;
	for(size_t data = 0; data < m_vNumberData.size(); ++data)
	{
		for(size_t s = 0; s < m_vNumberData[data].ssGrp.num_subsets(); ++s)
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
					vSegmentData[fct]->functor(val, bf.global_ip(), time());

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
					vSegmentData[fct]->functor(val, bf.global_ip(), time());

				// 	get associated node
					const int co = bf.node_id();

				// 	add to local matrix
					d(vSegmentData[fct]->locFct, co) -= VecDot(val, bf.normal());
				}
			}
		}
	}

// 	we're done
	return true;
}

////////////////////////////////////////////////////////////////////////////////
//	Constructor
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
FV1NeumannBoundaryElemDisc<TDomain>::FV1NeumannBoundaryElemDisc(const char* subsets)
 :IDomainElemDisc<TDomain>("", subsets)
{
	m_mBNDNumberBndSegment.clear();
	m_mVectorBndSegment.clear();
	register_all_fv1_funcs(false);
}


///	type of trial space for each function used
template<typename TDomain>
bool
FV1NeumannBoundaryElemDisc<TDomain>::
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
FV1NeumannBoundaryElemDisc<TDomain>::
treat_non_regular_grid(bool bNonRegular)
{
//	switch, which assemble functions to use.
	if(bNonRegular)
	{
		UG_LOG("ERROR in 'FVNeumannBoundaryElemDisc::treat_non_regular_grid':"
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
FV1NeumannBoundaryElemDisc<TDomain>::
register_all_fv1_funcs(bool bHang)
{
//	get all grid element types in this dimension and below
	typedef typename domain_traits<dim>::DimElemList ElemList;

//	switch assemble functions
	if(!bHang) boost::mpl::for_each<ElemList>( RegisterFV1<FV1Geometry>(this) );
	else throw (UGFatalError("Not implemented."));
}

template<typename TDomain>
template<typename TElem, template <class Elem, int WorldDim> class TFVGeom>
void
FV1NeumannBoundaryElemDisc<TDomain>::
register_fv1_func()
{
	ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;
	typedef this_type T;

	this->set_prep_elem_loop_fct(id, &T::template prepare_element_loop<TElem, TFVGeom>);
	this->set_prep_elem_fct(	 id, &T::template prepare_element<TElem, TFVGeom>);
	this->set_fsh_elem_loop_fct( id, &T::template finish_element_loop<TElem, TFVGeom>);
	this->set_ass_JA_elem_fct(		 id, &T::template assemble_JA<TElem, TFVGeom>);
	this->set_ass_JM_elem_fct(		 id, &T::template assemble_JM<TElem, TFVGeom>);
	this->set_ass_dA_elem_fct(		 id, &T::template assemble_A<TElem, TFVGeom>);
	this->set_ass_dM_elem_fct(		 id, &T::template assemble_M<TElem, TFVGeom>);
	this->set_ass_rhs_elem_fct(	 id, &T::template assemble_f<TElem, TFVGeom>);
}

////////////////////////////////////////////////////////////////////////////////
//	explicit template instantiations
////////////////////////////////////////////////////////////////////////////////

template class FV1NeumannBoundaryElemDisc<Domain1d>;
template class FV1NeumannBoundaryElemDisc<Domain2d>;
template class FV1NeumannBoundaryElemDisc<Domain3d>;

} // namespace ug

