/*
 * obstacle_constraint_interface_impl.h
 *
 *  Created on: 26.11.2013
 *      Author: raphaelprohl
 */

#ifndef __H__UG__LIB_ALGEBRA__OPERATOR__PRECONDITIONER__PROJECTED_GAUSS_SEIDEL__OBSTACLE_CONSTRAINT_INTERFACE_IMPL__
#define __H__UG__LIB_ALGEBRA__OPERATOR__PRECONDITIONER__PROJECTED_GAUSS_SEIDEL__OBSTACLE_CONSTRAINT_INTERFACE_IMPL__

#include "obstacle_constraint_interface.h"

namespace ug{

template <typename TAlgebra>
template <typename TElem>
void
IObstacleConstraint<TAlgebra>::
obstacle_indices_on_subset(const size_t si)
{
	typedef typename DoFDistribution::traits<TElem>::const_iterator iter_type;
	iter_type iter = m_spDD->begin<TElem>(si);
	iter_type iterEnd = m_spDD->end<TElem>(si);

//	loop all elements of type
	for( ; iter != iterEnd; ++iter)
	{
		//	get elem
		TElem* elem = *iter;

		m_spDD->inner_algebra_indices(elem, m_vIndicesOfObsSubsets, false);
		//m_spDD->inner_dof_indices(elem, fct, m_vIndicesOfObsSubset, false);
	}
}

/*template <typename TDomain, typename TAlgebra>
template <typename TElem, typename TIterator>
void
IObstacleConstraint<TDomain,TAlgebra>::
obstacle_indices_on_subset(TIterator iterBegin,
		TIterator iterEnd,
		const function_type& u)
{
	// 	local indices
	LocalIndices ind;

	// 	Loop over all elements on subset
	for(TIterator iter = iterBegin; iter != iterEnd; ++iter)
	{
		TElem* elem = *iter;

		// 	get global indices
		u.indices(elem, ind);

		UG_LOG("Elem: " << elem << "\n");

	//	TODO: maybe this is not the best order to push the obstacle-indices
	//	into the m_vIndicesOfObsSubset-vector, since the obsIndices are called
	//	wrt algebra indices later on!

		for (size_t fct = 0; fct < ind.num_fct(); fct++)
		{
			for (size_t dof = 0; dof < ind.num_dof(fct); dof++)
			{
				size_t index = ind.index(fct, dof);
				size_t comp = ind.comp(fct, dof);

				//	create vector of DoFindices in obstacle-subset.
				//	Only those pairs should be attached which are not
				//	already a member of the 'm_vIndicesOfObsSubset'-vector
				bool bDoFAlreadyInVec = false;

				UG_LOG("#DoFIndicesOfObsSubset global: " << m_vIndicesOfObsSubset.size() << "\n");

				for (vector<DoFIndex>::iterator itObsInd = m_vIndicesOfObsSubset.begin();
						itObsInd < m_vIndicesOfObsSubset.end(); itObsInd++)
				{
					UG_LOG("(*itObsInd)[0]: " << (*itObsInd)[0] << "\n");
					UG_LOG("(*itObsInd)[1]: " << (*itObsInd)[1] << "\n");

					if ((*itObsInd)[0] == index
							&& (*itObsInd)[1] == comp)
						bDoFAlreadyInVec = true;

					if (!bDoFAlreadyInVec)
					{
						UG_LOG("in !bDoFAlreadyInVec \n");

						if((*itObsInd)[0] > index)
						{
							UG_LOG("current index "<<index<<" is smaller \n");
							//	current index 'index' is smaller than the current index
							//	in m_vIndicesOfObsSubset, the iterator is pointing to
							m_vIndicesOfObsSubset.insert(itObsInd, DoFIndex(index, comp));
							bDoFAlreadyInVec = true;
						}
						else
						{
							if((*itObsInd)[0] == index)
							{
								if((*itObsInd)[1] > comp)
								{
									//	current comp 'comp' is smaller than the current comp
									//	in m_vIndicesOfObsSubset, the iterator is pointing to
									m_vIndicesOfObsSubset.insert(itObsInd, DoFIndex(index, comp));
									bDoFAlreadyInVec = true;
								}
							}
						}
					} //end if (!bAlreadyStoredDoFIndex)
				} //end for(vector)

				if (!bDoFAlreadyInVec)
				{
					//	insert (index,comp) at the end of the vector
					m_vIndicesOfObsSubset.push_back(DoFIndex(index, comp));
				}

			} //end(dof)
		} //end(fct)
	}//end(elem)

	UG_LOG("#DoFIndicesOfObsSubset global: " << m_vIndicesOfObsSubset.size() << "\n");

	for (vector<DoFIndex>::iterator itObsInd = m_vIndicesOfObsSubset.begin();
							itObsInd < m_vIndicesOfObsSubset.end(); ++itObsInd)
	{
		UG_LOG("globObsIndex: "<<(*itObsInd)[0]<<"\n");
		UG_LOG("globObsComp: "<<(*itObsInd)[1]<<"\n");
		UG_LOG("\n");
		UG_LOG("\n");
	}
}*/

template <typename TAlgebra>
void
IObstacleConstraint<TAlgebra>::init(const vector_type& u)
{
	//	TODO: the init of the vector of obstacle values should not depend on the solution vector u
	//	which is usually defined on the whole domain!

// 	init vector of obstacle values and init values with zero
	if (!m_bLowerObs)
		m_spVecOfLowObsValues = u.clone_without_values();
	if (!m_bUpperObs)
		m_spVecOfUpObsValues = u.clone_without_values();

//	check, that lower obstacle is <= upper obstacle (for all indices)
	if (m_bLowerObs && m_bUpperObs)
	{
		if ((*m_spVecOfLowObsValues).size() != (*m_spVecOfUpObsValues).size())
			UG_THROW("In IObstacleConstraint::init(u) :Vector of lower obstacle values [size= "
					<<(*m_spVecOfLowObsValues).size()<<"] and "
					" Vector of upper obstacle values [size= "
					<<(*m_spVecOfUpObsValues).size()<<"] sizes have to match!");

		for(size_t i = 0; i < (*m_spVecOfLowObsValues).size(); i++)
		{
			const value_type& lowerObsVal = (*m_spVecOfLowObsValues)[i];
			const value_type& upperObsVal = (*m_spVecOfUpObsValues)[i];
			for(size_t j = 0; j < GetSize(lowerObsVal); j++)
			{
				if (BlockRef(lowerObsVal, j) - BlockRef(upperObsVal, j) > 0.0)
					UG_THROW("In IObstacleConstraint::init(u) " <<i<<"-th index and "<<j<<"-th"
						" component of vector of lower obstacle [value= "<<lowerObsVal<<"] needs "
						"to be lower equal the "<<i<<"-th value of vector of upper obstacle "
						"[value= "<<upperObsVal<<"]!");
			}
		}
	}

	//	init pointer to vector of active indices
		m_spLowerActiveInd = &m_vActiveIndicesLow;
		m_spUpperActiveInd = &m_vActiveIndicesUp;


	//	create Subset Group
	try{
		m_ssGrp = m_spDD->subset_grp_by_name(m_ssName.c_str());
	}UG_CATCH_THROW(" Subsets '"<<m_ssName<<"' not"
					" all contained in DofDistribution.");

	//	reset vector of indices in obstacle-subsets
	m_vIndicesOfObsSubsets.resize(0); //m_spDD->num_indices());

	//	loop all subsets contained in the Subset Group
	//for (vector<int>::iterator activeSI = m_vActiveSubsets.begin();
	//			activeSI != m_vActiveSubsets.end(); ++activeSI)
	for (size_t i = 0; i < m_ssGrp.size(); i++)
	{
		int si = m_ssGrp[i];

		if(m_spDD->max_dofs(VERTEX)) {
			obstacle_indices_on_subset<VertexBase>(si);
			UG_LOG("VERTEX-BASE \n");
		}
		if(m_spDD->max_dofs(EDGE)) 	 {
			obstacle_indices_on_subset<EdgeBase>(si);
			UG_LOG("EDGE-BASE \n");
		}
		if(m_spDD->max_dofs(FACE))	 {
			obstacle_indices_on_subset<Face>(si);
			UG_LOG("FACE \n");
		}
		if(m_spDD->max_dofs(VOLUME)) {
			obstacle_indices_on_subset<Volume>(si);
			UG_LOG("VOLUME \n");
		}

		//	sort m_vIndicesOfObsSubsets
		sort(m_vIndicesOfObsSubsets.begin(), m_vIndicesOfObsSubsets.end());

		for (vector<size_t>::iterator itObsInd = m_vIndicesOfObsSubsets.begin();
			itObsInd < m_vIndicesOfObsSubsets.end(); ++itObsInd)
		{
			UG_LOG("algIndex: "<<(*itObsInd)<<"\n");
			UG_LOG("\n");
			UG_LOG("\n");
		}
	}

	/*const int subsetDim = DimensionOfSubset(*m_spDD->subset_handler(), si);

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
	}*/

}

} // end namespace ug

#endif /* __H__UG__LIB_ALGEBRA__OPERATOR__PRECONDITIONER__PROJECTED_GAUSS_SEIDEL__OBSTACLE_CONSTRAINT_INTERFACE_IMPL__ */
