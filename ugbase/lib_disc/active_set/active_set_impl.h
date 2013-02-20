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

template <typename TDomain, typename TAlgebra>
ActiveSet<TDomain, TAlgebra>::ActiveSet() : m_bCons(false)
{};

template <typename TDomain, typename TAlgebra>
void ActiveSet<TDomain, TAlgebra>::prepare(vector_type& u)
{
	m_vActiveSet.resize(0); m_vActiveSetOld.resize(0); m_vInactiveSet.resize(0);
	//m_vActiveSet.resize(u.size());

	if(!m_bCons)
		UG_THROW("No constraint set in ActiveSet \n");
}

template <typename TDomain, typename TAlgebra>
bool ActiveSet<TDomain, TAlgebra>::active_index(vector_type& u,
		vector_type& lambda)
{
	UG_LOG("ActiveSet:apply called \n");

	//	TODO: avoid this vector copy; m_vActiveSetOld really necessary?
	m_vActiveSetOld = m_vActiveSet;
	m_vActiveSet.resize(0); m_vInactiveSet.resize(0);

	if (u.size() != lambda.size())
		UG_THROW("Temporarily u and lambda need to be "
				"of same size in ActiveSet:is_an_index_active \n");

	//	TODO: note: these are blocks here!!!
	vector_type complementaryVal;
	complementaryVal.resize(0);
	complementaryVal.resize(1);

	bool index_is_active = false;

	for(size_t i = 0; i < u.size(); i++)
	{
		complementaryVal[0] = lambda[i] + u[i] - m_ConsVec[i];
		if (complementaryVal[0] <= 0.0)
		{
			// index i is inactive!
			// this is only valid for a constraint of type u <= *m_spConsVec
			UG_LOG("Index " << i << " is inactive! \n");
			lambda[i] = 0.0;
			m_vInactiveSet.push_back(i);
		}
		else{
			UG_LOG("Index " << i << " is active! \n");
			u[i] = m_ConsVec[i];
			//	fŸr diesen Index eine Dirichlet-Reihe in Matrix setzen!
			//	create list of Dirichlet-indizes!
			//	Problem: was passiert mit Dirichlet-Indizes auf Boundary?
			//lambda[i] = rhs[i] - u[i]; // or A * u instead of u?
			m_vActiveSet.push_back(i); //liste immer wieder neu aufbauen!
			index_is_active = true;
		}
	}

	return index_is_active;
}

template <typename TDomain, typename TAlgebra>
void ActiveSet<TDomain, TAlgebra>::comp_lambda(matrix_type& mat, vector_type& u,
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
	for (vector<size_t>::iterator it = m_vActiveSet.begin(); it < m_vActiveSet.end(); ++it)
	{
		UG_LOG("iterator: " << *it << " \n");

		/*
		for(size_t j = 0; j < u.size(); j++)
		{
			mat_u += mat[*it][j] * u[j];
		}*/

		lambda[*it] = rhs[*it] - mat_u[*it];
		UG_LOG("lambda[ " << *it << " ]" << lambda[*it] << " \n");

	}

	UG_LOG("new lambda computed \n");

	/*for(size_t i = 0; i < u.size(); i++){
		rhs[i] -= lambda[i];
	}

	UG_LOG("rhs updated \n");*/
}

template <typename TDomain, typename TAlgebra>
bool ActiveSet<TDomain, TAlgebra>::check_conv()
{
	if (m_vActiveSet.size() == m_vActiveSetOld.size())
	{
		size_t ind = 0;
		for (vector<size_t>::iterator it = m_vActiveSet.begin(); it < m_vActiveSet.end(); ++it)
		{
			if (*it != m_vActiveSetOld[ind]) return false;
			ind++;
		}
		//	activeSet not changed // +check, ob Nebenbedingung fŸr alle DoFs erfŸllt ist?!
		return true;
	}

	return false;
}

}; // namespace ug

#endif /* ACTIVE_SET_IMPL_H_ */
