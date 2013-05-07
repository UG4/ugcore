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
void ActiveSet<TDomain, TAlgebra>::prepare(vector_type& u)
{
	m_vActiveSet.resize(0); m_vActiveSetOld.resize(0);
}

template <typename TDomain, typename TAlgebra>
bool ActiveSet<TDomain, TAlgebra>::check_dist_to_obs(vector_type& u)
{
	//	STILL IN PROGRESS: u sollte hier reference-position + Startlšsung sein!
	value_type dist;

	bool geometry_cut_by_cons = false;

	for(size_t i = 0; i < u.size(); i++)
	{
		UG_LOG("u(" << i << "):" << u[i] << "\n");
		UG_LOG("m_ConsVec(" << i << "):" << m_ConsVec[i] << "\n");
		dist = m_ConsVec[i] - u[i];
		UG_LOG("dist:" << dist << "\n");
		//TODO: anstatt u muss hier die geometrische Info einflie§en!
		for (size_t fct = 0; fct < m_nrFcts; fct++)
		{
			if (BlockRef(dist,fct) < 0.0) // i.e.: u < m_ConsVec
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

template <typename TDomain, typename TAlgebra>
bool ActiveSet<TDomain, TAlgebra>::active_index(vector_type& u,
		vector_type& contactForce)
{
	if(!m_bCons)
		UG_THROW("No constraint set in ActiveSet \n");

	//	remember old ActiveSet for convergence check
	//	TODO: avoid this vector copy; m_vActiveSetOld really necessary?
	m_vActiveSetOld = m_vActiveSet;
	m_vActiveSet.resize(0);

	if (u.size() != contactForce.size())
		UG_THROW("Temporarily u and contactForce need to be "
				"of same size in ActiveSet:active_index \n");

	value_type complementaryVal;

	bool one_fct_is_active = false;

	for(size_t i = 0; i < u.size(); i++)
	{
		//	note: complementaryVal, contactForce[i], etc. are blocks here
		complementaryVal = contactForce[i] + u[i] - m_ConsVec[i];

		for (size_t fct = 0; fct < m_nrFcts; fct++)
		{
			if (BlockRef(complementaryVal,fct) <= 0.0)
			{
				//	multiindex (i,fct) is inactive!
				//	temporarily this is only valid
				//	for a constraint of type u <= *m_spConsVec
				BlockRef(contactForce[i],fct) = 0.0;
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

//	computes the contact forces for a given contact discretization
template <typename TDomain, typename TAlgebra>
void ActiveSet<TDomain, TAlgebra>::contactForces(vector_type& contactforce,
		const vector_type& u)
{
	if(m_vActiveSet.size() != 0.0)
	{
		//	check that contact disc is set
		if (m_spContactDisc.invalid())
			UG_THROW("No contact discretization set in "
						"ActiveSet:contactForces \n");

		//	TODO: restrict contactforce-gridfunction.size() to possible contact subset
		if (u.size() != contactforce.size())
			UG_THROW("Temporarily u and contactForce need to be "
					"of same size in ActiveSet:contactForces \n");

		m_spContactDisc->contactForces(contactforce, u, m_vActiveSet);
	}
}


//	computes the contact forces via the residuum: contactforce = rhs - mat * u;
template <typename TDomain, typename TAlgebra>
void ActiveSet<TDomain, TAlgebra>::contactForcesRes(vector_type& contactforce,
		const matrix_type& mat,
		const vector_type& u,
		const vector_type& rhs)
{
	// 	only if some indices are active we need to compute contact forces
	if(m_vActiveSet.size() != 0.0)
	{
		if (u.size() != contactforce.size())
			UG_THROW("Temporarily u and contactForce need to be "
					"of same size in ActiveSet:contactForcesRes \n");

		vector_type mat_u;
		mat_u.resize(u.size());

		#ifdef UG_PARALLEL
			MatMultDirect(mat_u, 1.0, mat, u);
		#else
			MatMult(mat_u, 1.0, mat, u);
		#endif

		//	loop MultiIndex-pairs in activeSet-vector
		for (vector<MultiIndex<2> >::iterator it = m_vActiveSet.begin();
				it < m_vActiveSet.end(); ++it)
		{
			//	compute contact forces for active multiIndices

			//	get active (DoF,fct)-pairs out of m_vActiveSet
			MultiIndex<2> activeMultiIndex = *it;

			size_t dof = activeMultiIndex[0];
			size_t fct = activeMultiIndex[1];

			//	contactForce = rhs - Mat * u;
			BlockRef(contactforce[dof],fct) = BlockRef(rhs[dof],fct) - BlockRef(mat_u[dof],fct);
		}

		UG_LOG("new contactforce-values computed \n");
	}
	else{
		UG_LOG("no active index in contactForcesRes \n");
	}

	/*for(size_t i = 0; i < u.size(); i++){
		rhs[i] -= contactforce[i];
	}
	UG_LOG("rhs updated \n");*/
}

template <typename TDomain, typename TAlgebra>
bool ActiveSet<TDomain, TAlgebra>::check_conv(const vector_type& u, const size_t step)
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

	//	check if constraint is satisfied for all multiIndices
	for(size_t i = 0; i < u.size(); i++)
	{
		penetration = u[i] - m_ConsVec[i];

		for (size_t fct = 0; fct < m_nrFcts; fct++){
			if (BlockRef(penetration,fct) > 0.0) //	i.e.: u > m_ConsVec
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

template <typename TDomain, typename TAlgebra>
void ActiveSet<TDomain, TAlgebra>::createVecOfPointers()
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
