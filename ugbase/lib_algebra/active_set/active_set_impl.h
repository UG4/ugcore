/*
 * active_set_impl.h
 *
 *  Created on: 15.02.2013
 *      Author: raphaelprohl
 */

#ifndef ACTIVE_SET_IMPL_H_
#define ACTIVE_SET_IMPL_H_

#include "active_set.h"

namespace ug {

template <typename TAlgebra>
ActiveSet<TAlgebra>::ActiveSet() : m_bCons(false)
{};


template <typename TAlgebra>
void ActiveSet<TAlgebra>::prepare(vector_type& u)
{
	m_vActiveSet.resize(0); m_vActiveSetOld.resize(0);
}

template <typename TAlgebra>
bool ActiveSet<TAlgebra>::check_dist_to_obs(vector_type& u)
{
	//	STILL IN PROGRESS
	value_type dist;
	//	get number of unknowns per value_type
	//	(e.g. if CPU == 3 -> nrFcts = 3!)
	size_t nrFcts = GetSize(dist);

	bool geometry_cut_by_cons = false;

	for(size_t i = 0; i < u.size(); i++)
	{
		dist = u[i] - m_ConsVec[i];
		//TODO: anstatt u muss hier die geometrische Info einflie§en!
		for (size_t fct = 0; fct < nrFcts; fct++)
		{
			if (BlockRef(dist,fct) < 0.0)
			{
				geometry_cut_by_cons = true;
				break;
			}
		}

		if (geometry_cut_by_cons)
			break;

	}

	return geometry_cut_by_cons;
}

template <typename TAlgebra>
bool ActiveSet<TAlgebra>::active_index(vector_type& u,
		vector_type& lambda)
{
	if(!m_bCons)
		UG_THROW("No constraint set in ActiveSet \n");

	//	remember old ActiveSet for convergence check
	//	TODO: avoid this vector copy; m_vActiveSetOld really necessary?
	m_vActiveSetOld = m_vActiveSet;
	m_vActiveSet.resize(0);

	if (u.size() != lambda.size())
		UG_THROW("Temporarily u and lambda need to be "
				"of same size in ActiveSet:active_index \n");

	value_type complementaryVal;

	//	get number of unknowns per value_type
	//	(e.g. if CPU == 3 -> nrFcts = 3!)
	size_t nrFcts = GetSize(complementaryVal);

	bool one_fct_is_active = false;

	for(size_t i = 0; i < u.size(); i++)
	{
		//	note: complementaryVal, lambda[i], etc. are blocks here
		complementaryVal = lambda[i] + u[i] - m_ConsVec[i];

		for (size_t fct = 0; fct < nrFcts; fct++)
		{
			if (BlockRef(complementaryVal,fct) <= 0.0)
			{
				//	multiindex (i,fct) is inactive!
				//	temporarily this is only valid
				//	for a constraint of type u <= *m_spConsVec
				BlockRef(lambda[i],fct) = 0.0;
			}
			else
			{
				one_fct_is_active = true;

				//	mark MultiIndex-pair (i,fct) as active
				MultiIndex<2> activeMultiIndex(i,fct);

				//	create list of active MultiIndex-pairs
				m_vActiveSet.push_back(activeMultiIndex);

				//	this corresponds to adjust_solution
				//	in the context of Dirichlet-nodes:
				BlockRef(u[i],fct) = BlockRef(m_ConsVec[i],fct);
			}
		}
	}

	return one_fct_is_active;
}

template <typename TAlgebra>
void ActiveSet<TAlgebra>::comp_lambda(vector_type& lambda,
		const matrix_type& mat,
		const vector_type& u,
		const vector_type& rhs)
{
	if (u.size() != lambda.size())
		UG_THROW("Temporarily u and lambda need to be "
				"of same size in ActiveSet:comp_lambda \n");

	vector_type mat_u; // mat_u2;
	mat_u.resize(u.size());
	//mat_u2.resize(u.size());

	// 	only if some indices are active we need to compute contact forces
	if(m_vActiveSet.size() != 0.0)
	{
		//	we only need *it-th row of mat -> number mat_u = 0.0;
		#ifdef UG_PARALLEL
			MatMultDirect(mat_u, 1.0, mat, u);
		#else
			MatMult(mat_u, 1.0, mat, u);
		#endif

		//	loop MultiIndex-pairs in activeSet-vector
		for (vector<MultiIndex<2> >::iterator it = m_vActiveSet.begin();
				it < m_vActiveSet.end(); ++it)
		{
			/*#ifdef UG_PARALLEL
				MatMultDirect(u2[*it], 1.0, mat(*it), u[*it]);
			#else
				MatMult(mat_u2[*it], 1.0, mat.row_index(*it), u[*it]);
			#endif*/

			//	compute contact forces (lambda) for active multiIndices

			//	get active (DoF,fct)-pairs out of m_vvActiveSet
			MultiIndex<2> activeMultiIndex = *it;

			size_t dof = activeMultiIndex[0];
			size_t fct = activeMultiIndex[1];

			//	lambda = rhs - Mat * u;
			BlockRef(lambda[dof],fct) = BlockRef(rhs[dof],fct) - BlockRef(mat_u[dof],fct);
		}

		UG_LOG("new lambda-values computed \n");
	}
	else{
		UG_LOG("no active index in comp_lambda \n");
	}

	/*for(size_t i = 0; i < u.size(); i++){
		rhs[i] -= lambda[i];
	}
	UG_LOG("rhs updated \n");*/
}

template <typename TAlgebra>
bool ActiveSet<TAlgebra>::check_conv(const vector_type& u, const size_t step)
{
	//	ensure that at least one activeSet-iteration is performed
	if (step <= 1)
		return false;

	//	NOW TWO CHECKS WILL BE PERFORMED TO ENSURE CONVERGENCE:
	//	1. 	Is the constraint violated by any multiIndex?
	//	2. 	Did some multiIndices change from 'active' to 'inactive' or vice versa
	//		in the last iteration-step?

	UG_LOG(m_vActiveSet.size() << " indices are active in step " << step << " ! \n");

	value_type penetration;
	//	get number of unknowns per value_type
	//	(e.g. if CPU == 3 -> nrFcts = 3!)
	size_t nrFcts = GetSize(penetration);

	//	check if constraint is satisfied for all multiIndices
	for(size_t i = 0; i < u.size(); i++)
	{
		penetration = u[i] - m_ConsVec[i];

		for (size_t fct = 0; fct < nrFcts; fct++){
			if (BlockRef(penetration,fct) > 0.0)
				return false;
		}
	}

	//	check if activeSet has changed
	if (m_vActiveSet.size() == m_vActiveSetOld.size())
	{
		UG_LOG("Old and new active Set have the same number of members \n");

		vector<MultiIndex<2> >::iterator it = m_vActiveSet.begin();

		for (vector<MultiIndex<2> >::iterator itOld = m_vActiveSetOld.begin();
				itOld < m_vActiveSetOld.end(); ++itOld)
		{
			MultiIndex<2> multiIndexOld = *itOld;
			MultiIndex<2> multiIndex = *it;

			if ((multiIndex[0] != multiIndexOld[0])
					|| (multiIndex[1] != multiIndexOld[1]))
				return false;

			++it;

		} // itOld

		//	activeSet remains unchanged & constraint is fulfilled for all indices
		return true;
	}

	return false;
}

template <typename TAlgebra>
void ActiveSet<TAlgebra>::createVecOfPointers()
{
	m_vActiveSetSP.resize(m_vActiveSet.size());

	vector<MultiIndex<2> >::iterator it = m_vActiveSet.begin();

	for (vector<SmartPtr<MultiIndex<2> > >::iterator itSP = m_vActiveSetSP.begin();
				itSP < m_vActiveSetSP.end(); ++itSP)
	{
		*itSP = &*it;
		++it;
	}

}

}; // namespace ug

#endif /* ACTIVE_SET_IMPL_H_ */
