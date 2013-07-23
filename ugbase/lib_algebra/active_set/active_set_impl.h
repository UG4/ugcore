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
	m_vActiveSetGlob.resize(0); m_vActiveSetGlobOld.resize(0);
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
		UG_LOG("m_spConsGF(" << i << "):" << (*m_spConsGF)[i] << "\n");
		dist = (*m_spConsGF)[i] - u[i];
		UG_LOG("dist:" << dist << "\n");
		//TODO: anstatt u muss hier die geometrische Info einflie§en!
		for (size_t fct = 0; fct < m_nrFcts; fct++)
		{
			if (BlockRef(dist,fct) < 0.0) // i.e.: u < m_spConsGF
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
template <typename TElem, typename TIterator>
void ActiveSet<TDomain, TAlgebra>::ActiveIndexElem(TIterator iterBegin,
		TIterator iterEnd,
		const function_type& u,
		function_type& contactForce)
{
// 	check if at least an element exists, else return
	if(iterBegin == iterEnd) return;

	int countElem = 0;
// 	local indices and local algebra
	LocalIndices indU, indCF, indCons;
	LocalVector locU, locCF, locCons;

	number complementaryVal;

	//	TODO: eventuell ist es sinnvoll auch die aktiven Elemente zu speichern,
	//	damit folgende Elem-loops (z.B. zur Berechnung der contactForces) verkleinert werden kšnnen!

// 	Loop over all elements on active subsets
	for(TIterator iter = iterBegin; iter != iterEnd; ++iter)
	{
		countElem++;

	// 	get global indices
		u.indices(*iter, indU); contactForce.indices(*iter, indCF); (*m_spConsGF).indices(*iter, indCons);

	// 	adapt local algebra
		locU.resize(indU); locCF.resize(indCF); locCons.resize(indCons);

	// 	read local values of u and contactForce
		GetLocalVector(locU, u); GetLocalVector(locCF, contactForce); GetLocalVector(locCons, *m_spConsGF);

	//	loop over DoFs in element and store all activeMultiIndex-pairs in vector
		size_t nFctElem = indU.num_fct();
		UG_LOG("#FctElem: " << nFctElem << "\n");
		for(size_t fct = 0; fct < nFctElem; ++fct)
		{
			size_t nDoFsPerFctElem = indU.num_dof(fct);
			UG_LOG("#nDoFsPerFctElem: " << nDoFsPerFctElem << "\n");

			for(size_t dof = 0; dof < nDoFsPerFctElem; ++dof)
			{
				UG_LOG("#locU( " << fct << "," << dof << "): " << locU(fct, dof) << "\n");
				UG_LOG("#locCF( " << fct << "," << dof << "): " << locCF(fct, dof) << "\n");
				UG_LOG("#locCons( " << fct << "," << dof << "): " << locCons(fct, dof) << "\n");

				//	TODO: complementaryVal only holds for displacements locU along
				//	a POSITIVE normal and contactForces along a NEGATIVE normal!
				complementaryVal = - 1.0 * locCF(fct ,dof) + locU(fct, dof) - locCons(fct, dof);
				UG_LOG("complementaryVal: " << complementaryVal << "\n");

				if (complementaryVal <= -1e-10)
				{
					//	TODO: ist das notwendig hier?
					//	multiindex (i,fct) is inactive!
					//	temporarily this is only valid
					//	for a constraint of type u <= *m_spConsVec
					locCF(fct,dof) = 0.0;
				}
				else{
					//	mark MultiIndex-pair (dof,fct) as active
					size_t globIndex = indU.index(fct, dof);
					size_t globComp = indU.comp(fct, dof);
					MultiIndex<2> activeMultiIndex(globIndex, globComp);

					bool bAlreadyActive = false;

					//	create list of active global MultiIndex-pairs. Only those pairs should be attached
					//	which are not already a member of the activeSet
					for (vector<MultiIndex<2> >::iterator itSet = m_vActiveSetGlob.begin();
							itSet < m_vActiveSetGlob.end(); ++itSet)
					{
						MultiIndex<2> multiIndexSet = *itSet;
						if ((multiIndexSet[0] == activeMultiIndex[0])
								&& (multiIndexSet[1] == activeMultiIndex[1]))
							bAlreadyActive = true;
					}

					if (!bAlreadyActive)
						m_vActiveSetGlob.push_back(activeMultiIndex);
				}
			} // end(dof)
		} // end(fct)

		// 	send local to global contactForce
		AddLocalVector(contactForce, locCF);

		UG_LOG("#activeDoFFctPairs global: " << m_vActiveSetGlob.size() << "\n");
	}
	UG_LOG("#elems: " << countElem << "\n");

}

template <typename TDomain, typename TAlgebra>
bool ActiveSet<TDomain, TAlgebra>::active_index(function_type& u,
		function_type& contactForce)
{
	if(!m_bCons)
		UG_THROW("No constraint set in ActiveSet \n");

	//	remember old ActiveSet for convergence check
	//	TODO: avoid this vector copy; m_vActiveSetGlobOld really necessary?
	m_vActiveSetGlobOld = m_vActiveSetGlob;
	m_vActiveSetGlob.resize(0);

	SmartPtr<DoFDistribution> dd = contactForce.dof_distribution();
	UG_LOG("#subsets: " << dd->num_subsets() << "\n");

	//	1.) get all subsets on which contactForce is defined!
	//	-> store them in vSubsetsOfContactForces
	m_vSubsetsOfContact.resize(0);
	for (int si = 0; si < dd->num_subsets(); si++){
		for (size_t fct = 0; fct < contactForce.num_fct(); fct++)
			if (contactForce.is_def_in_subset(fct,si))
			{
				m_vSubsetsOfContact.push_back(si);
				//	'break' is necessary to ensure that 'si' is
				//	only added once when several fcts of
				//	'contactForce' are defined in subset 'si'!
				break;
			}
	}

	UG_LOG("#sizeOfvSubsetsOfContactForces: " << m_vSubsetsOfContact.size() << "\n");

	//	case: contactForce and ConsVec are only defined on some subsets of the original domain
	//	2.) loop over all elements of the possible contact subsets
	for (std::vector<int>::iterator siContact = m_vSubsetsOfContact.begin();
			siContact != m_vSubsetsOfContact.end(); ++siContact)
	{
		UG_LOG("siContact: " << *siContact << "\n");
		const int subsetDim = DimensionOfSubset(*dd->subset_handler(), *siContact);
		UG_LOG("subsetDim: " << subsetDim << "\n");

		//	3.) get localU out of u for each element and
		//	4.) store the active (local) dof,fct-pair in m_vActiveSet
		switch(subsetDim)
		{
		case 0:
			break;
		case 1:
			ActiveIndexElem<Edge>
				(dd->template begin<Edge>(*siContact), dd->template end<Edge>(*siContact), u, contactForce);
			break;
		case 2:
			ActiveIndexElem<Triangle>
				(dd->template begin<Triangle>(*siContact), dd->template end<Triangle>(*siContact), u, contactForce);
			ActiveIndexElem<Quadrilateral>
				(dd->template begin<Quadrilateral>(*siContact), dd->template end<Quadrilateral>(*siContact), u, contactForce);
			break;
		case 3:
			ActiveIndexElem<Tetrahedron>
				(dd->template begin<Tetrahedron>(*siContact), dd->template end<Tetrahedron>(*siContact), u, contactForce);
			ActiveIndexElem<Pyramid>
				(dd->template begin<Pyramid>(*siContact), dd->template end<Pyramid>(*siContact), u, contactForce);
			ActiveIndexElem<Prism>
				(dd->template begin<Prism>(*siContact), dd->template end<Prism>(*siContact), u, contactForce);
			ActiveIndexElem<Hexahedron>
				(dd->template begin<Hexahedron>(*siContact), dd->template end<Hexahedron>(*siContact), u, contactForce);
			break;
		default:
			UG_THROW("ActiveSet::active_index:"
				"SubsetDimension "<<subsetDim<<" (subset="<<*siContact<<") not supported.");
		}
	}

	if ((*m_spConsGF).size() != contactForce.size())
		UG_THROW("ConstraintGridFunction and contactForce need to be "
				"of same size in ActiveSet:active_index \n");

	/*if (u.size() != contactForce.size())
		UG_THROW("Temporarily u and contactForce need to be "
				"of same size in ActiveSet:active_index \n");

	value_type complementaryVal;

	for(size_t i = 0; i < u.size(); i++)
	{
		//	note: complementaryVal, contactForce[i], etc. are blocks here
		//	TODO: complementaryVal only holds for displacements u along
		//	a POSITIVE normal and contactForces along a NEGATIVE normal!
		complementaryVal = - 1.0 * contactForce[i] + u[i] - (*m_spConsGF)[i];
		//complementaryVal = contactForce[i] + u[i] - *m_spConsGF[i];

		UG_LOG("for DoF i: " << i << "contactForce: "
				<< contactForce[i] << "u: " << u[i] << "Obs: " << (*m_spConsGF)[i] << "\n");
		UG_LOG("complementaryVal: " << complementaryVal << "\n");

		for (size_t fct = 0; fct < m_nrFcts; fct++)
		{
			if (BlockRef(complementaryVal, fct) <= -1e-10)
			{
				//	multiindex (i,fct) is inactive!
				//	temporarily this is only valid
				//	for a constraint of type u <= *m_spConsVec
				BlockRef(contactForce[i], fct) = 0.0;
			}
			else
			{
				//	mark MultiIndex-pair (i,fct) as active
				MultiIndex<2> activeMultiIndex(i, fct);

				//	create list of active MultiIndex-pairs
				m_vActiveSetGlob.push_back(activeMultiIndex);

				//UG_LOG("active (dof,fct) in ActiveSet::active_index: ("
				//		<< i << "," << fct << ")\n");

				//	this corresponds to adjust_solution
				//	in the context of Dirichlet-nodes:
				//	by uncommenting the next line, u is allowed to penetrate the obstacle
				//	for computing the contactForces later on!
				//BlockRef(u[i],fct) = BlockRef(*m_spConsGF[i],fct);
			}
		}
	}*/

	if (m_vActiveSetGlob.size() > 0) return true;
	else return false;
}

//	computes the contact forces for a given contact discretization
template <typename TDomain, typename TAlgebra>
void ActiveSet<TDomain, TAlgebra>::contactForces(function_type& contactforce,
		function_type& rhs, const function_type& u)
{
	//	check that contact disc is set
	if (m_spContactDisc.invalid())
		UG_THROW("No contact discretization set in "
					"ActiveSet:contactForces \n");

	if(m_vActiveSetGlob.size() != 0.0)
	{
		UG_LOG("activeDoFs in ActiveSet:contactForces " << m_vActiveSetGlob.size() << "\n");

		if ((*m_spConsGF).size() != contactforce.size())
			UG_THROW("ConstraintGridFunction and contactForce need to be "
					"of same size in ActiveSet:contactForces \n");

		m_spContactDisc->contactForces(contactforce, u, m_vActiveSetGlob, m_vSubsetsOfContact);

		/*for (vector<MultiIndex<2> >::iterator it = m_vActiveSetGlob.begin();
						it < m_vActiveSetGlob.end(); ++it)
		{
			//	update rhs with contact forces for active multiIndices

			//	get active (DoF,fct)-pairs out of m_vActiveSetGlob
			MultiIndex<2> activeMultiIndex = *it;

			size_t dof = activeMultiIndex[0];
			size_t fct = activeMultiIndex[1];

			//	TODO: pass the contactForces to ass_rhs instead of adding it here!
			//	temporarily this is only valid for initial rhs = 0!
			//	rhs = rhs + contactForce;
			//BlockRef(rhs[dof],fct) = BlockRef(rhs[dof],fct) + BlockRef(contactforce[dof],fct);
			BlockRef(rhs[dof],fct) = BlockRef(contactforce[dof],fct);
		}*/
	}
}


//	computes the contact forces via the residuum: contactforce = rhs - mat * u;
template <typename TDomain, typename TAlgebra>
void ActiveSet<TDomain, TAlgebra>::contactForcesRes(vector_type& contactforce,
		const matrix_type& mat,
		const vector_type& u,
		vector_type& rhs)
{
	// 	only if some indices are active we need to compute contact forces
	if(m_vActiveSetGlob.size() != 0.0)
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
		for (vector<MultiIndex<2> >::iterator it = m_vActiveSetGlob.begin();
				it < m_vActiveSetGlob.end(); ++it)
		{
			//	compute contact forces for active multiIndices

			//	get active (DoF,fct)-pairs out of m_vActiveSetGlob
			MultiIndex<2> activeMultiIndex = *it;

			size_t dof = activeMultiIndex[0];
			size_t fct = activeMultiIndex[1];

			//	contactForce = rhs - Mat * u;
			BlockRef(contactforce[dof],fct) = BlockRef(rhs[dof],fct) - BlockRef(mat_u[dof],fct);
			//BlockRef(rhs[dof],fct) = BlockRef(rhs[dof],fct) + BlockRef(contactforce[dof],fct);
		}

		UG_LOG("new contactforce-values computed \n");
		//UG_LOG("rhs updated \n");
	}
	else{
		UG_LOG("no active index in contactForcesRes \n");
	}

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

	UG_LOG(m_vActiveSetGlob.size() << " indices are active at the begin "
			"of step " << step << " ! \n");

	value_type gap;

	//	check if constraint is satisfied for all multiIndices
	for(size_t i = 0; i < u.size(); i++)
	{
		//	TODO: the following gap-computation holds only for displacements along
		//	a POSITIVE normal and contactForces along a NEGATIVE normal!
		gap = (*m_spConsGF)[i] - u[i];

		UG_LOG("for DoF i: " << i <<  " u: " << u[i] << "Obs: " << (*m_spConsGF)[i] << "\n");
		UG_LOG("gap: " << gap << "\n");

		for (size_t fct = 0; fct < m_nrFcts; fct++){
			if (BlockRef(gap,fct) < -1e-10) //	i.e.: m_spConsGF < u
				return false;
		}
	}

	//	check if activeSet has changed
	if (m_vActiveSetGlob.size() == m_vActiveSetGlobOld.size())
	{
		UG_LOG("Old and new active Set have the same number of members \n");

		vector<MultiIndex<2> >::iterator it = m_vActiveSetGlob.begin();

		for (vector<MultiIndex<2> >::iterator itOld = m_vActiveSetGlobOld.begin();
				itOld < m_vActiveSetGlobOld.end(); ++itOld)
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
	m_vActiveSetGlobSP.resize(m_vActiveSetGlob.size());

	vector<MultiIndex<2> >::iterator it = m_vActiveSetGlob.begin();

	for (vector<SmartPtr<MultiIndex<2> > >::iterator itSP = m_vActiveSetGlobSP.begin();
				itSP < m_vActiveSetGlobSP.end(); ++itSP)
	{
		*itSP = &*it;
		++it;
	}

}

}; // namespace ug

#endif /* ACTIVE_SET_IMPL_H_ */
