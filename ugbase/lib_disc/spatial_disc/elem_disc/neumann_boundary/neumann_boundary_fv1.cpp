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
template <typename TUserData, typename TScheduledUserData>
bool FV1NeumannBoundaryElemDisc<TDomain>::
extract_scheduled_data(std::map<int, std::vector<TUserData> >& mvUserDataBndSegment,
                       const std::vector<TScheduledUserData>& vScheduledUserData,
                       FunctionGroup& commonFctGrp, std::string& fctNames)
{
//	clear all extracted data
	mvUserDataBndSegment.clear();

//	loop all scheduled data
	for(size_t i = 0; i < vScheduledUserData.size(); ++i)
	{
	//	create Function Group and Subset Group
		FunctionGroup functionGroup;
		SubsetGroup subsetGroup;

	//	convert strings
		if(!ConvertStringToSubsetGroup(subsetGroup, this->get_fct_pattern(),
		                               vScheduledUserData[i].ssName.c_str()))
		{
			UG_LOG("ERROR in 'FV1NeumannBoundaryElemDisc:extract_scheduled_data':"
					" Subsets '"<<vScheduledUserData[i].ssName<<"' not"
					" all contained in ApproximationSpace.\n");
		}

		if(!ConvertStringToFunctionGroup(functionGroup, this->get_fct_pattern(),
		                                 vScheduledUserData[i].fctName.c_str()))
		{
			UG_LOG("ERROR in 'FV1NeumannBoundaryElemDisc:extract_scheduled_data':"
					" Functions '"<<vScheduledUserData[i].fctName<<"' not"
					" all contained in ApproximationSpace.\n");
			return false;
		}

	//	check that only one function given
		if(functionGroup.num_fct() != 1)
		{
			UG_LOG("ERROR in 'FV1NeumannBoundaryElemDisc:extract_scheduled_data':"
					" Only one function allowed for a neumann value.\n");
			return false;
		}

	//	get function
		const size_t fct = functionGroup[0];

	// 	check if function exist
		if(fct >= this->get_fct_pattern().num_fct())
		{
			UG_LOG("ERROR in 'LagrangeDirichletBoundary:extract_scheduled_data':"
					" Function "<< fct << " does not exist in pattern.\n");
			return false;
		}

	//	add to common fct group if not already contained
		if(commonFctGrp.contains(fct))
		{
			UG_LOG("ERROR in 'LagrangeDirichletBoundary:extract_scheduled_data':"
				" More than one neumann value specified for function "<<
				vScheduledUserData[i].fctName.c_str()<<".\n");
			return false;
		}
		else
		{
			commonFctGrp.add(fct);
		}

	//	build string of functions
		if(!fctNames.empty()) fctNames.append(",");
		fctNames.append(vScheduledUserData[i].fctName.c_str());

	//	get subsethandler
		const ISubsetHandler* pSH = this->get_fct_pattern().get_subset_handler();

	// 	loop subsets
		for(size_t si = 0; si < subsetGroup.num_subsets(); ++si)
		{
		//	get subset index
			const int subsetIndex = subsetGroup[si];

		// 	check that function is defined for segment
			if(!this->get_fct_pattern().is_def_in_subset(fct, subsetIndex))
			{
				UG_LOG("ERROR in 'LagrangeDirichletBoundary:extract_scheduled_data':"
					" Function "<<fct<<" not defined on subset "<<subsetIndex<<".\n");
				return false;
			}

		//	check that subsetIndex is valid
			if(subsetIndex < 0 || subsetIndex >= pSH->num_subsets())
			{
				UG_LOG("ERROR in 'FVNeumannBoundaryElemDisc:extract_scheduled_data':"
						" Invalid subset Index " << subsetIndex <<
						". (Valid is 0, .. , " << pSH->num_subsets() <<").\n");
				return false;
			}

		//	get Boundary segment from map
			std::vector<TUserData>& vSegmentFunction = mvUserDataBndSegment[subsetIndex];

		//	remember functor and function
			vSegmentFunction.push_back(TUserData(commonFctGrp.num_fct()-1,
			                                     vScheduledUserData[i].functor));
		}
	}

	return true;
}

template<typename TDomain>
bool FV1NeumannBoundaryElemDisc<TDomain>::
extract_scheduled_data()
{
//	a common function group
	FunctionGroup commonFctGrp;
	commonFctGrp.set_function_pattern(this->get_fct_pattern());

//	string of functions
	std::string fctNames;

	bool bRet = true;
	bRet &= extract_scheduled_data(m_mBNDNumberBndSegment, m_vScheduledBNDNumberData,
	                               commonFctGrp, fctNames);
	bRet &= extract_scheduled_data(m_mVectorBndSegment, m_vScheduledVectorData,
	                               commonFctGrp, fctNames);

//	set name of function
	this->set_functions(fctNames.c_str());

//	done
	return bRet;
}

template<typename TDomain>
void FV1NeumannBoundaryElemDisc<TDomain>::
add(BNDNumberFunctor& user, const char* function, const char* subsets)
{
	m_vScheduledBNDNumberData.push_back(ScheduledBNDNumberData(user, function, subsets));

	if(this->fct_pattern_set()) extract_scheduled_data();
}

template<typename TDomain>
void FV1NeumannBoundaryElemDisc<TDomain>::
add(VectorFunctor& user, const char* function, const char* subsets)
{
	m_vScheduledVectorData.push_back(ScheduledVectorData(user, function, subsets));

	if(this->fct_pattern_set()) extract_scheduled_data();
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

	{
	typename std::map<int, std::vector<BNDNumberData> >::const_iterator subsetIter;
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
	typename std::map<int, std::vector<VectorData> >::const_iterator subsetIter;
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

	{
	typename std::map<int, std::vector<BNDNumberData> >::const_iterator subsetIter;
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
	typename std::map<int, std::vector<VectorData> >::const_iterator subsetIter;
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
prepare_element(TElem* elem, const local_vector_type& u)
{
//	get corners
	m_vCornerCoords = this->template get_element_corners<TElem>(elem);

//  update Geometry for this element
	static TFVGeom<TElem, dim>& geo = Provider<TFVGeom<TElem,dim> >::get();
	if(!geo.update(elem, &m_vCornerCoords[0], &(this->get_subset_handler())))
	{
		UG_LOG("ERROR in 'FVNeumannBoundaryElemDisc::prepare_element': "
				"Cannot update Finite Volume Geometry.\n");
		return false;
	}

//	we're done
	return true;
}

template<typename TDomain>
template<typename TElem, template <class Elem, int  Dim> class TFVGeom>
inline
bool
FV1NeumannBoundaryElemDisc<TDomain>::
assemble_JA(local_matrix_type& J, const local_vector_type& u)
{
	// we're done
	return true;
}


template<typename TDomain>
template<typename TElem, template <class Elem, int  Dim> class TFVGeom>
inline
bool
FV1NeumannBoundaryElemDisc<TDomain>::
assemble_JM(local_matrix_type& J, const local_vector_type& u)
{
	// we're done
	return true;
}


template<typename TDomain>
template<typename TElem, template <class Elem, int  Dim> class TFVGeom>
inline
bool
FV1NeumannBoundaryElemDisc<TDomain>::
assemble_A(local_vector_type& d, const local_vector_type& u)
{
	// we're done
	return true;
}


template<typename TDomain>
template<typename TElem, template <class Elem, int  Dim> class TFVGeom>
inline
bool
FV1NeumannBoundaryElemDisc<TDomain>::
assemble_M(local_vector_type& d, const local_vector_type& u)
{
	// we're done
	return true;
}


template<typename TDomain>
template<typename TElem, template <class Elem, int  Dim> class TFVGeom>
inline
bool
FV1NeumannBoundaryElemDisc<TDomain>::
assemble_f(local_vector_type& d)
{
// 	get finite volume geometry
	const static TFVGeom<TElem, dim>& geo = Provider<TFVGeom<TElem,dim> >::get();

// 	loop registered boundary segments
	{
	typename std::map<int, std::vector<BNDNumberData> >::const_iterator subsetIter;
	for(subsetIter = m_mBNDNumberBndSegment.begin(); subsetIter != m_mBNDNumberBndSegment.end(); ++subsetIter)
	{
	//	get subset index corresponding to boundary
		const int bndSubset = (*subsetIter).first;

	//	get evaluation function vector
		const std::vector<BNDNumberData>& vSegmentFunction = (*subsetIter).second;

	// 	loop Boundary Faces
		for(size_t i = 0; i < geo.num_bf(bndSubset); ++i)
		{
		// get current BF
			const typename TFVGeom<TElem, dim>::BF& bf = geo.bf(bndSubset, i);

		//	loop functions, where neumann bnd is set for this bndSubset
			for(size_t fct = 0; fct < vSegmentFunction.size(); ++fct)
			{
			// 	get neumann value
				number val = 0.0;
				vSegmentFunction[fct].functor(val, bf.global_ip(), time());

			// 	get associated node
				const int co = bf.node_id();

			// 	add to local matrix
				d(vSegmentFunction[fct].loc_fct, co) -= val * bf.volume();
			}
		}
	}
	}

// 	loop registered boundary segments
	{
	typename std::map<int, std::vector<VectorData> >::const_iterator subsetIter;
	for(subsetIter = m_mVectorBndSegment.begin(); subsetIter != m_mVectorBndSegment.end(); ++subsetIter)
	{
	//	get subset index corresponding to boundary
		const int bndSubset = (*subsetIter).first;

	//	get evaluation function vector
		const std::vector<VectorData>& vSegmentFunction = (*subsetIter).second;

	// 	loop Boundary Faces
		for(size_t i = 0; i < geo.num_bf(bndSubset); ++i)
		{
		// get current BF
			const typename TFVGeom<TElem, dim>::BF& bf = geo.bf(bndSubset, i);

		//	loop functions, where neumann bnd is set for this bndSubset
			for(size_t fct = 0; fct < vSegmentFunction.size(); ++fct)
			{
			// 	get neumann value
				MathVector<dim> val;
				vSegmentFunction[fct].functor(val, bf.global_ip(), time());

			// 	get associated node
				const int co = bf.node_id();

			// 	add to local matrix
				d(vSegmentFunction[fct].loc_fct, co) -= VecDot(val, bf.normal());
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
FV1NeumannBoundaryElemDisc<TDomain>::FV1NeumannBoundaryElemDisc()
{
	m_mBNDNumberBndSegment.clear();
	register_all_fv1_funcs(false);
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
	typedef typename domain_traits<dim>::AllElemList ElemList;

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

	set_prep_elem_loop_fct(id, &T::template prepare_element_loop<TElem, TFVGeom>);
	set_prep_elem_fct(	 id, &T::template prepare_element<TElem, TFVGeom>);
	set_fsh_elem_loop_fct( id, &T::template finish_element_loop<TElem, TFVGeom>);
	set_ass_JA_elem_fct(		 id, &T::template assemble_JA<TElem, TFVGeom>);
	set_ass_JM_elem_fct(		 id, &T::template assemble_JM<TElem, TFVGeom>);
	set_ass_dA_elem_fct(		 id, &T::template assemble_A<TElem, TFVGeom>);
	set_ass_dM_elem_fct(		 id, &T::template assemble_M<TElem, TFVGeom>);
	set_ass_rhs_elem_fct(	 id, &T::template assemble_f<TElem, TFVGeom>);
}

////////////////////////////////////////////////////////////////////////////////
//	explicit template instantiations
////////////////////////////////////////////////////////////////////////////////

template class FV1NeumannBoundaryElemDisc<Domain1d>;
template class FV1NeumannBoundaryElemDisc<Domain2d>;
template class FV1NeumannBoundaryElemDisc<Domain3d>;

} // namespace ug

