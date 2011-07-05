/*
 * neumann_boundary_fv1.cpp
 *
 *  Created on: 14.10.2010
 *      Author: andreasvogel
 */

#include "neumann_boundary.h"
#include "lib_discretization/spatial_discretization/disc_util/geometry_provider.h"
#include "lib_discretization/spatial_discretization/disc_util/finite_volume_geometry.h"

namespace ug{


template<typename TDomain>
bool
FV1NeumannBoundaryElemDisc<TDomain>::
add_boundary_value(IBoundaryData<number, dim>& user, const char* function, const char* subsets)
{
//	create Function Group and Subset Group
	FunctionGroup functionGroup;
	SubsetGroup subsetGroup;

//	convert strings
	if(!ConvertStringToSubsetGroup(subsetGroup, this->get_fct_pattern(), subsets))
	{
		UG_LOG("ERROR while parsing Subsets.\n");
		return false;
	}
	if(!ConvertStringToFunctionGroup(functionGroup, this->get_fct_pattern(), function))
	{
		UG_LOG("ERROR while parsing Functions.\n");
		return false;
	}

//	only one function allowed
	if(functionGroup.num_fct() != 1)
	{
		UG_LOG("ERROR in 'FVNeumannBoundaryElemDisc:add_boundary_value':"
				" Exactly one function needed, but given '"<<function<<"' as functions.\n");
		return false;
	}

//	set name of function
	m_numFct = 1;
	this->set_functions(function);

//	forward request
	return add_boundary_value(user, functionGroup.unique_id(0), subsetGroup);
}

template<typename TDomain>
bool
FV1NeumannBoundaryElemDisc<TDomain>::
add_boundary_value(IBoundaryData<number, dim>& user, size_t fct, SubsetGroup bndSubsetGroup)
{
//	get subsethandler
	const ISubsetHandler* pSH = this->get_fct_pattern().get_subset_handler();

// 	check if function exist
	if(fct >= this->get_fct_pattern().num_fct())
	{
		UG_LOG("ERROR in 'FVNeumannBoundaryElemDisc:add_boundary_value': Function "
				<< fct << " does not exist in pattern.\n");
		return false;
	}

	//\TODO: THIS IS NOW BROKEN, HANDLE WHEN POSSIBLE  !!!!!!
//	get position of function in function group
	const size_t index = 0;

// 	loop subsets
	for(size_t si = 0; si < bndSubsetGroup.num_subsets(); ++si)
	{
	//	get subset index
		const int subsetIndex = bndSubsetGroup[si];

	//	check that subsetIndex is valid
		if(subsetIndex < 0 || subsetIndex >= pSH->num_subsets())
		{
			UG_LOG("ERROR in 'FVNeumannBoundaryElemDisc:add_boundary_value':"
					" Invalid subset Index " << subsetIndex <<
					". (Valid is 0, .. , " << pSH->num_subsets() <<").\n");
			return false;
		}

	//	get Boundary segment from map
		std::vector<UserDataFunction>& vSegmentFunction = m_mBoundarySegment[subsetIndex];

	//	remember functor and function
		vSegmentFunction.push_back(UserDataFunction(index, user.get_functor()));
	}

//	we're done
	return true;
}


template<typename TDomain>
template<typename TElem, template <class Elem, int  Dim> class TFVGeom>
inline
bool
FV1NeumannBoundaryElemDisc<TDomain>::
prepare_element_loop()
{
//	register subsetIndex at Geometry
	static TFVGeom<TElem, dim>& geo = GeomProvider::get<TFVGeom<TElem,dim> >();

	typename std::map<int, std::vector<UserDataFunction> >::const_iterator subsetIter;
	for(subsetIter = m_mBoundarySegment.begin();
			subsetIter != m_mBoundarySegment.end(); ++subsetIter)
	{
		const int bndSubset = (*subsetIter).first;

		geo.add_boundary_subset(bndSubset);
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
	static TFVGeom<TElem, dim>& geo = GeomProvider::get<TFVGeom<TElem,dim> >();

	typename std::map<int, std::vector<UserDataFunction> >::const_iterator subsetIter;
	for(subsetIter = m_mBoundarySegment.begin();
			subsetIter != m_mBoundarySegment.end(); ++subsetIter)
	{
		const int bndSubset = (*subsetIter).first;

		geo.remove_boundary_subset(bndSubset);
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
	static TFVGeom<TElem, dim>& geo = GeomProvider::get<TFVGeom<TElem,dim> >();
	if(!geo.update(elem, this->get_subset_handler(), &m_vCornerCoords[0]))
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
	// get finite volume geometry
	const static TFVGeom<TElem, dim>& geo = GeomProvider::get<TFVGeom<TElem,dim> >();

	// loop registered boundary segments
	typename std::map<int, std::vector<UserDataFunction> >::const_iterator subsetIter;
	for(subsetIter = m_mBoundarySegment.begin(); subsetIter != m_mBoundarySegment.end(); ++subsetIter)
	{
		const int bndSubset = (*subsetIter).first;
		const std::vector<UserDataFunction>& vSegmentFunction = (*subsetIter).second;

		// loop Boundary Faces
		for(size_t i = 0; i < geo.num_bf(bndSubset); ++i)
		{
		// get current BF
			const typename TFVGeom<TElem, dim>::BF& bf = geo.bf(bndSubset, i);

		//	loop functions, where neumann bnd is set for this bndSubset
			for(size_t fct = 0; fct < vSegmentFunction.size(); ++fct)
			{
				// first value
				number val = 0.0;
				vSegmentFunction[fct].functor(val, bf.global_ip(), time());

				// get associated node
				const int co = bf.node_id();

				// Add to local matrix
				d(vSegmentFunction[fct].loc_fct, co) -= val * bf.volume();
			}
		}
	}

	// we're done
	return true;
}

////////////////////////////////////////////////////////////////////////////////
//	Constructor
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
FV1NeumannBoundaryElemDisc<TDomain>::FV1NeumannBoundaryElemDisc()
: m_numFct(0)
{
	m_mBoundarySegment.clear();
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
	typedef typename GridElemTypes<dim>::AllElemList ElemList;

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

	reg_prepare_elem_loop_fct(id, &T::template prepare_element_loop<TElem, TFVGeom>);
	reg_prepare_elem_fct(	 id, &T::template prepare_element<TElem, TFVGeom>);
	reg_finish_elem_loop_fct( id, &T::template finish_element_loop<TElem, TFVGeom>);
	reg_ass_JA_elem_fct(		 id, &T::template assemble_JA<TElem, TFVGeom>);
	reg_ass_JM_elem_fct(		 id, &T::template assemble_JM<TElem, TFVGeom>);
	reg_ass_dA_elem_fct(		 id, &T::template assemble_A<TElem, TFVGeom>);
	reg_ass_dM_elem_fct(		 id, &T::template assemble_M<TElem, TFVGeom>);
	reg_ass_rhs_elem_fct(	 id, &T::template assemble_f<TElem, TFVGeom>);
}

////////////////////////////////////////////////////////////////////////////////
//	explicit template instantiations
////////////////////////////////////////////////////////////////////////////////

template class FV1NeumannBoundaryElemDisc<Domain1d>;
template class FV1NeumannBoundaryElemDisc<Domain2d>;
template class FV1NeumannBoundaryElemDisc<Domain3d>;

} // namespace ug

