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

	//	TODO: note: these are blocks here!!!
	vector_type complementaryVal;
	complementaryVal.resize(0);
	complementaryVal.resize(1);

	bool index_is_active = false;

	for(size_t i = 0; i < u.size(); i++)
	{
		//	TODO: consider blocks here!
		complementaryVal[0] = lambda[i] + u[i] - m_ConsVec[i];
		if (complementaryVal[0] <= 0.0)
		{
			//	index i is inactive!
			//	temporarily this is only valid
			//	for a constraint of type u <= *m_spConsVec
			//UG_LOG("Index " << i << " is inactive! \n");
			lambda[i] = 0.0;
		}
		else{
			//UG_LOG("Index " << i << " is active! \n");
			u[i] = m_ConsVec[i];

			//	create list of active indices
			m_vActiveSet.push_back(i);
			index_is_active = true;
		}
	}

	return index_is_active;
}

template <typename TAlgebra>
void ActiveSet<TAlgebra>::comp_lambda(matrix_type& mat,
		vector_type& u,
		vector_type& rhs,
		vector_type& lambda)
{
	if (u.size() != lambda.size())
		UG_THROW("Temporarily u and lambda need to be "
				"of same size in ActiveSet:comp_lambda \n");

	//	fŸr die inaktive Menge -> lambda = 0 setzen oder die untere Operation
	//	nur fŸr aktive Indizes durchfŸhren!!!
	// TODO: next call is only valid for ParallelVector/-Matrix!
	//MatMultAddDirect(lambda, 1.0, rhs, -1.0, mat, u);

	vector_type mat_u;
	mat_u.resize(u.size());

	// 	only if some indices are active we need to compute contact forces
	if(m_vActiveSet.size() != 0.0)
	{
		// TODO: next call is only valid for ParallelVector/-Matrix!
		//	we only need *it-th row of mat -> number mat_u = 0.0;
		#ifdef UG_PARALLEL
			MatMultDirect(mat_u, 1.0, mat, u);
		#else
			for(size_t i = 0; i < u.size(); i++)
				mat_u[i] = 0.0;
			UG_LOG("mat_u set to zero!!! \n");
		#endif


		//	loop indices in activeSet-vector
		for (vector<size_t>::iterator it = m_vActiveSet.begin();
				it < m_vActiveSet.end(); ++it)
		{
			lambda[*it] = rhs[*it] - mat_u[*it];
			//UG_LOG("lambda[ " << *it << " ]" << lambda[*it] << " \n");
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
bool ActiveSet<TAlgebra>::check_conv(vector_type& u, size_t step)
{
	//	ensure that at least one activeSet-iteration is performed
	if (step <= 1)
		return false;

	//	check if constraint is satisfied for all indices
	vector_type penetration;
	penetration.resize(0);
	penetration.resize(1);

	for(size_t i = 0; i < u.size(); i++)
	{
		penetration[0] = u[i] - m_ConsVec[i];
		if (penetration[0] > 0.0) return false;
	}

	if (m_vActiveSet.size() == m_vActiveSetOld.size())
	{
		size_t ind = 0;
		for (vector<size_t>::iterator it = m_vActiveSet.begin();
				it < m_vActiveSet.end(); ++it)
		{
			if (*it != m_vActiveSetOld[ind]) return false;
			/*UG_LOG("index " << m_vActiveSetOld[ind] <<
					" remains unchanged in activeSet \n");*/
			ind++;
		}
		//	activeSet remains unchanged & constraint is fulfilled for all indices
		return true;
	}

	return false;
}

}; // namespace ug

#endif /* ACTIVE_SET_IMPL_H_ */
