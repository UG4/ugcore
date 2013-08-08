/*
 * active_set_impl.h
 *
 *  Created on: 15.02.2013
 *      Author: raphaelprohl
 */

#ifndef ACTIVE_SET_IMPL_H_
#define ACTIVE_SET_IMPL_H_

#include "active_set.h"
#include "lib_disc/common/geometry_util.h"

namespace ug {

template <int dim> struct face_type_traits
{
    typedef void face_type0;
	typedef void face_type1;
};

template <> struct face_type_traits<1>
{
    typedef ReferenceVertex face_type0;
	typedef ReferenceVertex face_type1;
};

template <> struct face_type_traits<2>
{
    typedef ReferenceEdge face_type0;
	typedef ReferenceEdge face_type1;
};

template <> struct face_type_traits<3>
{
    typedef ReferenceTriangle face_type0;
	typedef ReferenceQuadrilateral face_type1;
};

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

//	determines a vector of global (dof,fct)-pairs, which are active (wrt. an obstacle-constraint)
template <typename TDomain, typename TAlgebra>
template <typename TElem, typename TIterator>
void ActiveSet<TDomain, TAlgebra>::ActiveIndexElem(TIterator iterBegin,
		TIterator iterEnd,
		function_type& u,
		function_type& rhs,
		function_type& contactForce)
{
	// 	check if at least an element exists, else return
	if(iterBegin == iterEnd) return;

	int countElem = 0;

	static const int dim = function_type::dim;
	typedef typename function_type::domain_type domain_type;
	typedef typename domain_type::grid_type grid_type;
	typedef typename function_type::element_type element_type;
	typedef typename element_type::side side_type;

	typename grid_type::template traits<side_type>::secure_container sides;
	typename grid_type::template traits<element_type>::secure_container assoElements;

	typedef typename face_type_traits<dim>::face_type0 face_type0;
	typedef typename face_type_traits<dim>::face_type1 face_type1;

	//	get position accessor
	typename domain_type::position_accessor_type& aaPos
			= u.domain()->position_accessor();

	//ConstSmartPtr<TDomain> dom = u.domain();

	//	storage for corner coordinates
	MathVector<dim> normal;
	std::vector<MathVector<dim> > vCorner;
	MathVector<dim> sideCoPos[dim+1];
	//std::vector<MathVector<dim> > sideCorners;
	//MathVector<dim> vCornerCoords[TElem::NUM_VERTICES];

	// 	local indices and local algebra
	LocalIndices indU, indCF, indCons, indRhs;
	LocalVector locU, locCF, locCons, locRhs;

	number complementaryVal, complementaryVal_n_scaled;

	//	TODO: eventuell ist es sinnvoll auch die aktiven Elemente zu speichern,
	//	damit folgende Elem-loops (z.B. zur Berechnung der contactForces) verkleinert werden kšnnen!

	// 	Loop over all elements on active subsets
	for(TIterator iter = iterBegin; iter != iterEnd; ++iter)
	{
		TElem* elem = *iter;

	//  get sides of element
		//grid_type& grid = *u.domain()->grid();
		//grid.get_associated(elements, elem);

		/*grid.associated_elements_sorted(sides, elem);
		for (size_t s = 0; s < sides.size(); s++)
		{
			grid.associated_elements(assoElements, sides[s]);
			if (assoElements.size() > 1)
				UG_THROW("Multiple associated elements found in ActiveSet::ActiveIndexElem! \n");

			element_type* assoE = assoElements[0];

			MathVector<dim> normal;
			CollectCornerCoordinates(sideCorners, *sides[s], aaPos);
			for (size_t i = 0; i < sideCorners.size(); i++)
				for (int d = 0; d < dim; d++)
					sideCoPos[i][d] = sideCorners[i][d];

			//	faces have dim corners in 1d, 2d
			//	in 3d they have dim corners (triangle) or dim+1 corners (quadrilateral)
			//	TODO: ensure that we use the outer normal here!
			if ((int)sideCorners.size() == dim)
				ElementNormal<face_type0, dim>(normal, sideCoPos);
			else
				ElementNormal<face_type1, dim>(normal, sideCoPos);

			for (size_t i = 0; i < sideCorners.size(); i++)
				for (int d = 0; d < dim; d++)
					UG_LOG("sideCoPos(" << i << "," << d << "): " << sideCoPos[i][d] << "\n");
			UG_LOG("normal: " << normal << "\n");
		}*/

	//	reference object type
		ReferenceObjectID roid = elem->reference_object_id();

		const DimReferenceElement<dim-1>& rRefElem
				= ReferenceElementProvider::get<dim-1>(roid);

		//UG_LOG("NORMAL-COMPUTATION VIA VCORNER \n");

	//	get corners of element
		CollectCornerCoordinates(vCorner, *elem, aaPos);

	//	here the ordering of the corners in the reference element is exploited
	//	in order to compute the outer normal later on
		for (int i = 0; i < (int)vCorner.size(); ++i)
			sideCoPos[i] = vCorner[rRefElem.id(dim-1, 0, 0, i)];

		if ((int)vCorner.size() == dim)
		{
			ElementNormal<face_type0, dim>(normal, sideCoPos);
			UG_LOG("face_type0 \n");
		}
		else
		{
			ElementNormal<face_type1, dim>(normal, sideCoPos);
			UG_LOG("face_type1 \n");
		}

		for (int i = 0; i < (int)vCorner.size(); ++i)
			UG_LOG("sideCoPos: " << sideCoPos[i] << "\n");
		UG_LOG("normal: " << normal << "\n");

		/*UG_LOG("NORMAL-COMPUTATION VIA VCORNERCOORDS \n");
	//	get corner coordinates
		FillCornerCoordinates(vCornerCoords, *elem, *dom);

		//	call 'ElementNormal' depending on elem-type
		//	faces have dim corners in 1d, 2d
		//	in 3d they have dim corners (triangle) or dim+1 corners (quadrilateral)
		//	TODO: ensure that we use the outer normal here!
		//if ((int)(*vCornerCoords).size() == dim)
		if (sizeof(vCornerCoords)/sizeof(vCornerCoords[0]) == dim)
		{
			ElementNormal<face_type0, dim>(normal, vCornerCoords);
			UG_LOG("face_type0 \n");
		}
		else
		{
			ElementNormal<face_type1, dim>(normal, vCornerCoords);
			UG_LOG("face_type1 \n");
		}

		for (int i = 0; i <= (int)(*vCornerCoords).size(); ++i)
			UG_LOG("vCornerCoord: " << vCornerCoords[i] << "\n");
		UG_LOG("normal: " << normal << "\n");*/

		countElem++;

	// 	get global indices
		u.indices(*iter, indU); contactForce.indices(*iter, indCF);
		(*m_spConsGF).indices(*iter, indCons); rhs.indices(*iter, indRhs);

	// 	adapt local algebra
		locU.resize(indU); locCF.resize(indCF);
		locCons.resize(indCons); locRhs.resize(indRhs);

	//	reset contribution of this element
		locCF = 0.0; locRhs = 0.0;

	// 	read local values of u and contactForce
		GetLocalVector(locU, u); GetLocalVector(locCF, contactForce);
		GetLocalVector(locCons, *m_spConsGF); GetLocalVector(locRhs, rhs);

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

				MathVector<dim> locUDof;
				for(int i = 0; i < dim; ++i)
					locUDof[i] = locU(i, dof);

				number locU_n = VecDot(locUDof, normal);
				number NormNormal = VecLength(normal);
				UG_LOG("NormNormal:" << NormNormal << "\n");
				number locU_n_scaled = locU_n / NormNormal;
				if (locU_n != 0.0)
				{
					UG_LOG("locU_n:" << locU_n << "\n");
					UG_LOG("locU_n_scaled:" << locU_n_scaled << "\n");
				}

				//	-1.0 * locCF corresponds to the lagrange multiplier lambda, c.f. Hintermueller/Ito/Kunisch(2003)
				complementaryVal = -1.0 * locCF(fct ,dof) + locU(fct, dof) - locCons(fct, dof);
				UG_LOG("complementaryVal: " << complementaryVal << "\n");
				complementaryVal_n_scaled = -1.0 * locCF(fct ,dof) + locU_n_scaled - locCons(fct, dof);
				UG_LOG("complementaryVal_n_scaled: " << complementaryVal_n_scaled << "\n");

				if (complementaryVal_n_scaled <= 1e-10)
				{
					//	TODO: ist das notwendig hier?
					//	multiindex (dof,fct0) is inactive!
					locCF(fct,dof) = 0.0;
				}
				else
				{
					//	mark MultiIndex-pair (dof,fct) as active
					const size_t globIndex = indU.index(fct, dof);
					const size_t globComp = indU.comp(fct, dof);
					MultiIndex<2> activeMultiIndex(globIndex, globComp);

					bool bAlreadyActive = false;

					//	create list of active global MultiIndex-pairs. Only those pairs should be attached
					//	which are not already a member of the 'activeSetGlob'-vector
					for (vector<MultiIndex<2> >::iterator itSet = m_vActiveSetGlob.begin();
							itSet < m_vActiveSetGlob.end(); ++itSet)
					{
						MultiIndex<2> multiIndexSet = *itSet;
						if ((multiIndexSet[0] == activeMultiIndex[0])
								&& (multiIndexSet[1] == activeMultiIndex[1]))
							bAlreadyActive = true;
					}

					if (!bAlreadyActive)
					{
						m_vActiveSetGlob.push_back(activeMultiIndex);
						BlockRef(rhs[globIndex], globComp) = locCons.value(fct,dof);
					}
				}
			} // end(dof)
		} // end(fct)

		// 	send local to global contactForce
		//AddLocalVector(contactForce, locCF);

		UG_LOG("#activeDoFFctPairs global: " << m_vActiveSetGlob.size() << "\n");
	} // end(elem)
	UG_LOG("#elems: " << countElem << "\n");

}

template <typename TDomain, typename TAlgebra>
bool ActiveSet<TDomain, TAlgebra>::active_index(function_type& u,
		function_type& rhs, function_type& contactForce, function_type& gap)
{
	//	note: first the active-index search is restricted to those cases,
	//	in which the constraint is of the form u * n <= consGF

	if(!m_bCons)
		UG_THROW("No constraint set in ActiveSet \n");

	//	remember old ActiveSet for convergence check
	//	TODO: avoid this vector copy; m_vActiveSetGlobOld really necessary?
	m_vActiveSetGlobOld = m_vActiveSetGlob;
	m_vActiveSetGlob.resize(0);

	SmartPtr<DoFDistribution> dd = gap.dof_distribution();
	UG_LOG("#subsets: " << dd->num_subsets() << "\n");

	//	1.) get all subsets on which the gap-gridfunction is defined!
	//	-> store them in vSubsetsOfContactForces
	m_vSubsetsOfContact.resize(0);
	//TODO: it is only necessary to loop over all BoundarySubsets!
	for (int si = 0; si < dd->num_subsets(); si++){
		for (size_t fct = 0; fct < gap.num_fct(); fct++)
			if (gap.is_def_in_subset(fct, si))
			{
				m_vSubsetsOfContact.push_back(si);
				//	'break' is necessary to ensure that 'si' is
				//	only added once when several fcts of
				//	'gap' are defined in subset 'si'!
				break;
			}
	}

	UG_LOG("#sizeOfvSubsetsOfGap: " << m_vSubsetsOfContact.size() << "\n");

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
				(dd->template begin<Edge>(*siContact), dd->template end<Edge>(*siContact), u, rhs, contactForce);
			break;
		case 2:
			ActiveIndexElem<Triangle>
				(dd->template begin<Triangle>(*siContact), dd->template end<Triangle>(*siContact), u, rhs, contactForce);
			ActiveIndexElem<Quadrilateral>
				(dd->template begin<Quadrilateral>(*siContact), dd->template end<Quadrilateral>(*siContact), u, rhs, contactForce);
			break;
		default:
			UG_THROW("ActiveSet::active_index:"
				"SubsetDimension "<<subsetDim<<" (subset="<<*siContact<<") not supported.");
		}
	}

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
		m_spContactDisc->contactForces(contactforce, u, m_vActiveSetGlob, m_vSubsetsOfContact);
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
template <typename TElem, typename TIterator>
bool ActiveSet<TDomain, TAlgebra>::ActiveSetConvCheckElem(TIterator iterBegin,
		TIterator iterEnd, function_type& u)
{
	static const int dim = function_type::dim;
	typedef typename function_type::domain_type domain_type;
	typedef typename domain_type::grid_type grid_type;
	typedef typename function_type::element_type element_type;
	typedef typename element_type::side side_type;

	typename grid_type::template traits<side_type>::secure_container sides;
	typename grid_type::template traits<element_type>::secure_container assoElements;

	typedef typename face_type_traits<dim>::face_type0 face_type0;
	typedef typename face_type_traits<dim>::face_type1 face_type1;

	//	get position accessor
	typename domain_type::position_accessor_type& aaPos
			= u.domain()->position_accessor();

	//	storage for corner coordinates
	MathVector<dim> normal;
	std::vector<MathVector<dim> > vCorner;
	MathVector<dim> sideCoPos[dim+1];

	// 	local indices and local algebra
	LocalIndices indU, indCons;
	LocalVector locU, locCons;

	// 	Loop over all elements on active subsets
	for(TIterator iter = iterBegin; iter != iterEnd; ++iter)
	{
		TElem* elem = *iter;

	//	reference object type
		ReferenceObjectID roid = elem->reference_object_id();

		const DimReferenceElement<dim-1>& rRefElem
				= ReferenceElementProvider::get<dim-1>(roid);

	//	get corners of element
		CollectCornerCoordinates(vCorner, *elem, aaPos);

	//	here the ordering of the corners in the reference element is exploited
	//	in order to compute the outer normal later on
		for (int i = 0; i < (int)vCorner.size(); ++i)
			sideCoPos[i] = vCorner[rRefElem.id(dim-1, 0, 0, i)];

		if ((int)vCorner.size() == dim)
			ElementNormal<face_type0, dim>(normal, sideCoPos);
		else
			ElementNormal<face_type1, dim>(normal, sideCoPos);

	// 	get global indices
		u.indices(*iter, indU); (*m_spConsGF).indices(*iter, indCons);

	// 	adapt local algebra
		locU.resize(indU); locCons.resize(indCons);

	// 	read local values of u and contactForce
		GetLocalVector(locU, u); GetLocalVector(locCons, *m_spConsGF);

	//	loop over DoFs in element and store all activeMultiIndex-pairs in vector
		size_t nFctElem = indU.num_fct();

		number gap_value, locU_n, NormNormal, locU_n_scaled;

		for(size_t fct = 0; fct < nFctElem; ++fct)
		{
			size_t nDoFsPerFctElem = indU.num_dof(fct);
			for(size_t dof = 0; dof < nDoFsPerFctElem; ++dof)
			{
				MathVector<dim> locUDof;
				for(int i = 0; i < dim; ++i)
					locUDof[i] = locU(i, dof);

				locU_n = VecDot(locUDof, normal);
				NormNormal = VecLength(normal);
				locU_n_scaled = locU_n / NormNormal;

				gap_value = locCons(fct, dof) - locU_n_scaled;

				if (gap_value < -1e-10){
					//	i.e.: m_spConsGF < u
					//	constraint is violated
					return false;
				}
			} // end(dof)
		} // end(fct)
	} // end(elem)

	return true;
}

template <typename TDomain, typename TAlgebra>
bool ActiveSet<TDomain, TAlgebra>::check_conv(function_type& u, const size_t step)
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

	//	check if constraint is violated
	SmartPtr<function_type> spConsGF = m_spConsGF.cast_const();
	SmartPtr<DoFDistribution> dd = (*spConsGF).dof_distribution();

	bool bConstraintViolated = false;
	for (std::vector<int>::iterator siContact = m_vSubsetsOfContact.begin();
			siContact != m_vSubsetsOfContact.end(); ++siContact)
	{
		const int subsetDim = DimensionOfSubset(*dd->subset_handler(), *siContact);
		switch(subsetDim)
		{
		case 0:
			break;
		case 1:
			if (!ActiveSetConvCheckElem<Edge>
				(dd->template begin<Edge>(*siContact), dd->template end<Edge>(*siContact), u))
			{bConstraintViolated = true;}

			break;
		case 2:
			if (!ActiveSetConvCheckElem<Triangle>
				(dd->template begin<Triangle>(*siContact), dd->template end<Triangle>(*siContact), u))
			{bConstraintViolated = true;}

			if (!ActiveSetConvCheckElem<Quadrilateral>
				(dd->template begin<Quadrilateral>(*siContact), dd->template end<Quadrilateral>(*siContact), u))
			{bConstraintViolated = true;}

			break;
		default:
			UG_THROW("ActiveSet::check_conv:"
				"SubsetDimension "<<subsetDim<<" (subset="<<*siContact<<") not supported.");
		}

		if (bConstraintViolated)
			return false;
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
