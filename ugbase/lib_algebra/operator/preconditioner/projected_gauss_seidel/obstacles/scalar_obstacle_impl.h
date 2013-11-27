/*
 * scalar_obstacle_impl.h
 *
 *  Created on: 25.11.2013
 *      Author: raphaelprohl
 */

#ifndef __H__UG__LIB_ALGEBRA__OPERATOR__PRECONDITIONER__PROJECTED_GAUSS_SEIDEL__SCALAR_OBSTACLE_IMPL__
#define __H__UG__LIB_ALGEBRA__OPERATOR__PRECONDITIONER__PROJECTED_GAUSS_SEIDEL__SCALAR_OBSTACLE_IMPL__

#include "scalar_obstacle.h"

namespace ug{

template <typename TDomain, typename TAlgebra>
template <typename TElem, typename TIterator>
void
ScalarObstacle<TDomain, TAlgebra>::
obstacle_indices_on_subset(TIterator iterBegin,
		TIterator iterEnd,
		function_type& u)
{
	// 	check if at least an element exists, else return
	if(iterBegin == iterEnd) return;

	// 	local indices
	LocalIndices ind;

	// 	Loop over all elements on subset
	for(TIterator iter = iterBegin; iter != iterEnd; ++iter)
	{
		TElem* elem = *iter;

		// 	get global indices
		u.indices(elem, ind);
	}
}




template <typename TDomain, typename TAlgebra>
void
ScalarObstacle<TDomain, TAlgebra>::
init(function_type& u)
{
	//	create Subset Group
	try{
		m_ssGrp = m_spApproxSpace->subset_grp_by_name(m_ssName.c_str());
	}UG_CATCH_THROW(" Subsets '"<<m_ssName<<"' not"
					" all contained in ApproximationSpace.");

	//	loop all subsets contained in the Subset Group
	//for (vector<int>::iterator activeSI = m_vActiveSubsets.begin();
	//			activeSI != m_vActiveSubsets.end(); ++activeSI)
	for (size_t i = 0; i < m_ssGrp.size(); i++)
	{
		int si = m_ssGrp[i];
		const int subsetDim = DimensionOfSubset(*m_spDD->subset_handler(), si);

		//	store the indices in m_vObsIndices
		switch(subsetDim)
		{
		case 0:
			break;
		case 1:
			obstacle_indices_on_subset<Edge>(m_spDD->template begin<Edge>(si),
						m_spDD->template end<Edge>(si), u);
			break;
		case 2:
			obstacle_indices_on_subset<Triangle>(m_spDD->template begin<Triangle>(si),
						m_spDD->template end<Triangle>(si), u);
			obstacle_indices_on_subset<Quadrilateral>(m_spDD->template begin<Quadrilateral>(si),
						m_spDD->template end<Quadrilateral>(si), u);
			break;
		case 3:
			obstacle_indices_on_subset<Tetrahedron>(m_spDD->template begin<Tetrahedron>(si),
						m_spDD->template end<Tetrahedron>(si), u);
			obstacle_indices_on_subset<Pyramid>(m_spDD->template begin<Pyramid>(si),
						m_spDD->template end<Pyramid>(si), u);
			obstacle_indices_on_subset<Prism>(m_spDD->template begin<Prism>(si),
						m_spDD->template end<Prism>(si), u);
			obstacle_indices_on_subset<Hexahedron>(m_spDD->template begin<Hexahedron>(si),
						m_spDD->template end<Hexahedron>(si), u);
			break;
		default:
			UG_THROW("ScalarObstacle::init:"
				"SubsetDimension "<< subsetDim <<" (subset="<< si <<") not supported.");
		}
	}
}

template <typename TDomain, typename TAlgebra>
void
ScalarObstacle<TDomain, TAlgebra>::
correction_for_lower_obs(vector_type& c, vector_type& lastSol, const size_t index, const value_type& tmpSol)
{
	//	TODO: in order to handle such cases, in which the obstacle is only set on some boundary-subset
	//	enter here something like
	//	if(obsSubset) do //oder if (indexIsInObsSubset)
	//	{	if (tmpSol - lowerObs) < 0.0 -> projection
	//		else lastSol = tmpSol;
	//	}
	//	else lastSol = tmpSol;


	//	get index-th lower obstacle value
	const value_type& lowerObsVal = (*m_spVecOfLowObsValues)[index];

	if(GetSize(tmpSol) != GetSize(lowerObsVal))
		UG_THROW("size of tmpSol and size of upperObsVal need to be the same");

	for(size_t j = 0; j < GetSize(tmpSol); j++)
	{
		if ( (BlockRef(tmpSol, j) - BlockRef(lowerObsVal, j)) < 0.0)
		{
			//	u_{s-1/2} < lowerObsValue (:the lower constraint is not fulfilled)

			//	adjust correction c := u_s - u_{s-1} = m_obsVal - u_{s-1}
			BlockRef(c[index], j) = BlockRef(lowerObsVal, j) - BlockRef(lastSol[index], j);

			//	set new solution u_s to the obstacle value
			//	and store the current index in a vector for further treatment
			BlockRef(lastSol[index], j) = BlockRef(lowerObsVal, j);
			(*m_spLowerActiveInd).push_back(MultiIndex<2>(index, j) );
		}
		else
		{
			//	the 'tmpSol' is valid with respect to the lower constraints
			BlockRef(lastSol[index], j) = BlockRef(tmpSol, j);
		}
	}
}

template <typename TDomain, typename TAlgebra>
void
ScalarObstacle<TDomain, TAlgebra>::
correction_for_upper_obs(vector_type& c, vector_type& lastSol, const size_t index, const value_type& tmpSol)
{
	//	get index-th upper obstacle value
	const value_type& upperObsVal = (*m_spVecOfUpObsValues)[index];

	if(GetSize(tmpSol) != GetSize(upperObsVal))
		UG_THROW("size of tmpSol and size of upperObsVal need to be the same");

	for(size_t j = 0; j < GetSize(tmpSol); j++)
	{
		if ( (BlockRef(tmpSol, j) - BlockRef(upperObsVal, j)) > 0.0)
		{
			//	u_{s-1/2} > upperObsValue (:the upper constraint is not fulfilled)

			//	adjust correction c := u_s - u_{s-1} = m_obsVal - u_{s-1}
			BlockRef(c[index], j) = BlockRef(upperObsVal, j) - BlockRef(lastSol[index], j);

			//	set new solution u_s to the obstacle value
			//	and store the current index in a vector for further treatment
			BlockRef(lastSol[index], j) = BlockRef(upperObsVal, j);
			(*m_spUpperActiveInd).push_back(MultiIndex<2>(index, j) );
		}
		else
		{
			//	the 'tmpSol' is valid with respect to the upper constraints
			BlockRef(lastSol[index], j) = BlockRef(tmpSol, j);
		}
	}
}

template <typename TDomain, typename TAlgebra>
void
ScalarObstacle<TDomain, TAlgebra>::
correction_for_lower_and_upper_obs(vector_type& c, vector_type& lastSol, const size_t index, const value_type& tmpSol)
{
	//	get index-th lower obstacle value
	const value_type& upperObsVal = (*m_spVecOfUpObsValues)[index];
	const value_type& lowerObsVal = (*m_spVecOfLowObsValues)[index];

	if(GetSize(tmpSol) != GetSize(upperObsVal))
		UG_THROW("size of tmpSol and size of upperObsVal need to be the same");
	if(GetSize(tmpSol) != GetSize(lowerObsVal))
		UG_THROW("size of tmpSol and size of upperObsVal need to be the same");

	for(size_t j = 0; j < GetSize(tmpSol); j++)
	{
		if ( (BlockRef(tmpSol, j) - BlockRef(upperObsVal, j)) > 0.0)
		{
			//	u_{s-1/2} > upperObsValue (:the upper constraint is not fulfilled)

			//	adjust correction c := u_s - u_{s-1} = m_obsVal - u_{s-1}
			BlockRef(c[index], j)  = BlockRef(upperObsVal, j) - BlockRef(lastSol[index], j);

			//	set new solution u_s to the obstacle value
			//	and store the current index in a vector for further treatment
			BlockRef(lastSol[index], j) = BlockRef(upperObsVal, j);
			(*m_spUpperActiveInd).push_back(MultiIndex<2>(index, j) );
		}
		else
		{
			if ( (BlockRef(tmpSol, j) - BlockRef(lowerObsVal, j)) < 0.0)
			{
				//	u_{s-1/2} < lowerObsValue (:the lower constraint is not fulfilled)

				//	adjust correction c := u_s - u_{s-1} = m_obsVal - u_{s-1}
				BlockRef(c[index], j) = BlockRef(lowerObsVal, j) - BlockRef(lastSol[index], j);

				//	set new solution u_s to the obstacle value
				//	and store the current index in a vector for further treatment
				BlockRef(lastSol[index], j) = BlockRef(lowerObsVal, j);
				(*m_spLowerActiveInd).push_back(MultiIndex<2>(index, j) );
			}
			else
			{
				//	the 'tmpSol' is valid with respect to all constraints
				BlockRef(lastSol[index], j) = BlockRef(tmpSol, j);
			}
		}
	}
}

} // end namespace ug

#endif /* __H__UG__LIB_ALGEBRA__OPERATOR__PRECONDITIONER__PROJECTED_GAUSS_SEIDEL__SCALAR_OBSTACLE_IMPL__ */
