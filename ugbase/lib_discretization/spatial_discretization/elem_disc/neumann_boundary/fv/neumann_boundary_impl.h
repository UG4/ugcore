/*
 * neumann_boundary_impl.h
 *
 *  Created on: 14.10.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__NEUMANN_BOUNDARY__FV1__NEUMANN_BOUNDARY_IMPL__
#define __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__NEUMANN_BOUNDARY__FV1__NEUMANN_BOUNDARY_IMPL__

#include "neumann_boundary.h"

namespace ug{


template<typename TDomain, typename TAlgebra>
bool
FV1NeumannBoundaryElemDisc<TDomain, TAlgebra>::
add_boundary_value(IBoundaryData<number, dim>& user, const char* function, const char* subsets)
{
//	check that function pattern exists
	if(this->m_pPattern == NULL)
	{
		UG_LOG("ERROR in 'FVNeumannBoundaryElemDisc:add_boundary_value':"
				" Function Pattern not set.\n");
		return false;
	}

//	create Function Group and Subset Group
	FunctionGroup functionGroup;
	SubsetGroup subsetGroup;

//	convert strings
	if(!ConvertStringToSubsetGroup(subsetGroup, *this->m_pPattern, subsets))
	{
		UG_LOG("ERROR while parsing Subsets.\n");
		return false;
	}
	if(!ConvertStringToFunctionGroup(functionGroup, *this->m_pPattern, function))
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

//	forward request
	return add_boundary_value(user, functionGroup.unique_id(0), subsetGroup);
}

template<typename TDomain, typename TAlgebra>
bool
FV1NeumannBoundaryElemDisc<TDomain, TAlgebra>::
add_boundary_value(IBoundaryData<number, dim>& user, size_t fct, SubsetGroup bndSubsetGroup)
{
//	check that function pattern exists
	if(this->m_pPattern == NULL)
	{
		UG_LOG("ERROR in 'FVNeumannBoundaryElemDisc:add_boundary_value': "
				"Function Pattern not set.\n");
		return false;
	}

//	get subsethandler
	const ISubsetHandler* pSH = this->m_pPattern->get_subset_handler();

// 	check if function exist
	if(fct >= this->m_pPattern->num_fct())
	{
		UG_LOG("ERROR in 'FVNeumannBoundaryElemDisc:add_boundary_value': Function "
				<< fct << " does not exist in pattern.\n");
		return false;
	}

//	add function to function group if needed
//	this will force, that the local vector contains the function
//	note: 	The local indices for allready contained functions are not changed.
//			Therefore an update in the UserDataFunctions is not necessary
	if(!this->m_FunctionGroup.contains(fct))
	{
		this->m_FunctionGroup.add(fct);
		m_numFct = this->m_FunctionGroup.num_fct();
	}

//	get position of function in function group
	const size_t index = this->m_FunctionGroup.local_index(fct);

// 	check that function is defined on inner subset
	if(!this->m_SubsetGroup.empty())
	{
		for(size_t si = 0; si < this->m_SubsetGroup.num_subsets(); ++si)
		{
			const int subsetIndex = this->m_SubsetGroup[si];
			if(!this->m_pPattern->is_def_in_subset(fct, subsetIndex))
			{
				UG_LOG("ERROR in 'FVNeumannBoundaryElemDisc:add_boundary_value':"
						" Function " << fct << " not defined in subset "
						<< subsetIndex << ".\n");
				return false;
			}
		}
	}

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


template<typename TDomain, typename TAlgebra>
template<typename TElem, template <class Elem, int  Dim> class TFVGeom>
inline
bool
FV1NeumannBoundaryElemDisc<TDomain, TAlgebra>::
prepare_element_loop()
{
//	register subsetIndex at Geometry
	TFVGeom<TElem, dim>& geo = FVGeometryProvider::get_geom<TFVGeom, TElem,dim>();

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

template<typename TDomain, typename TAlgebra>
template<typename TElem, template <class Elem, int  Dim> class TFVGeom>
bool
FV1NeumannBoundaryElemDisc<TDomain, TAlgebra>::
finish_element_loop()
{
//	remove subsetIndex from Geometry
	TFVGeom<TElem, dim>& geo = FVGeometryProvider::get_geom<TFVGeom, TElem,dim>();

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

template<typename TDomain, typename TAlgebra>

template<typename TElem, template <class Elem, int  Dim> class TFVGeom>
inline
bool
FV1NeumannBoundaryElemDisc<TDomain, TAlgebra>::
prepare_element(TElem* elem, const local_vector_type& u, const local_index_type& glob_ind)
{
//	get corners
	m_vCornerCoords = this->template get_element_corners<TElem>(elem);

//  update Geometry for this element
	TFVGeom<TElem, dim>& geo = FVGeometryProvider::get_geom<TFVGeom, TElem,dim>();
	if(!geo.update(elem, this->get_subset_handler(), &m_vCornerCoords[0]))
	{
		UG_LOG("ERROR in 'FVNeumannBoundaryElemDisc::prepare_element': "
				"Cannot update Finite Volume Geometry.\n");
		return false;
	}

//	we're done
	return true;
}

template<typename TDomain, typename TAlgebra>
template<typename TElem, template <class Elem, int  Dim> class TFVGeom>
inline
bool
FV1NeumannBoundaryElemDisc<TDomain, TAlgebra>::
assemble_JA(local_matrix_type& J, const local_vector_type& u)
{
	// we're done
	return true;
}


template<typename TDomain, typename TAlgebra>
template<typename TElem, template <class Elem, int  Dim> class TFVGeom>
inline
bool
FV1NeumannBoundaryElemDisc<TDomain, TAlgebra>::
assemble_JM(local_matrix_type& J, const local_vector_type& u)
{
	// we're done
	return true;
}


template<typename TDomain, typename TAlgebra>
template<typename TElem, template <class Elem, int  Dim> class TFVGeom>
inline
bool
FV1NeumannBoundaryElemDisc<TDomain, TAlgebra>::
assemble_A(local_vector_type& d, const local_vector_type& u)
{
	// we're done
	return true;
}


template<typename TDomain, typename TAlgebra>
template<typename TElem, template <class Elem, int  Dim> class TFVGeom>
inline
bool
FV1NeumannBoundaryElemDisc<TDomain, TAlgebra>::
assemble_M(local_vector_type& d, const local_vector_type& u)
{
	// we're done
	return true;
}


template<typename TDomain, typename TAlgebra>
template<typename TElem, template <class Elem, int  Dim> class TFVGeom>
inline
bool
FV1NeumannBoundaryElemDisc<TDomain, TAlgebra>::
assemble_f(local_vector_type& d)
{
	// get finite volume geometry
	TFVGeom<TElem, dim>& geo = FVGeometryProvider::get_geom<TFVGeom, TElem, dim>();

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

} // namespace ug


#endif /*__H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__NEUMANN_BOUNDARY__FV1__NEUMANN_BOUNDARY_IMPL__*/
