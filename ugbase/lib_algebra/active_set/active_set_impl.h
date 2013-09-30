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
#include "lib_disc/spatial_disc/disc_util/fe_geom.h"

namespace ug {

template <int dim> struct face_type_traits
{
    typedef void face_type0;
	typedef void face_type1;
	typedef void DimFEGeo;
};

template <> struct face_type_traits<1>
{
    typedef ReferenceVertex face_type0;
	typedef ReferenceVertex face_type1;
	typedef DimFEGeometry<1, 1> DimFEGeo;
};

template <> struct face_type_traits<2>
{
    typedef ReferenceEdge face_type0;
	typedef ReferenceEdge face_type1;
	typedef DimFEGeometry<2, 1> DimFEGeo;
};

template <> struct face_type_traits<3>
{
    typedef ReferenceTriangle face_type0;
	typedef ReferenceQuadrilateral face_type1;
	typedef DimFEGeometry<3, 2> DimFEGeo;
};

template <typename TDomain, typename TAlgebra>
void ActiveSet<TDomain, TAlgebra>::prepare(function_type& u)
{
	m_vActiveSetGlob.resize(0); m_vActiveSetGlobOld.resize(0);
	m_spDD = u.dof_distribution();
	m_spDom = u.domain();
}

/*template <typename TDomain, typename TAlgebra>
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
}*/

//	builds up a vector of global (dof,fct)-pairs, which are active (wrt. an obstacle-constraint)
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

	//	storage for corner coordinates
	MathVector<dim> normal;
	vector<MathVector<dim> > vCorner;
	MathVector<dim> sideCoPos[dim+1];

	// 	local indices and local algebra
	LocalIndices ind, indCons;
	LocalVector locU, locCF, locCons;

	//number complementaryVal;
	number complementaryVal_n_scaled;

	//	TODO: eventuell ist es sinnvoll auch die aktiven Elemente zu speichern,
	//	damit folgende Elem-loops (z.B. zur Berechnung der contactForces) verkleinert werden kšnnen!

	// 	Loop over all elements on active subsets
	for(TIterator iter = iterBegin; iter != iterEnd; ++iter)
	{
		TElem* sideElem = *iter;

	//	reference object type
		ReferenceObjectID roid = sideElem->reference_object_id();

		const DimReferenceElement<dim-1>& rRefElem
				= ReferenceElementProvider::get<dim-1>(roid);

	//	get corners of element
		CollectCornerCoordinates(vCorner, *sideElem, aaPos);

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

		/*for (int i = 0; i < (int)vCorner.size(); ++i)
			UG_LOG("sideCoPos: " << sideCoPos[i] << "\n");
		UG_LOG("normal: " << normal << "\n");*/

		countElem++;

	// 	get global indices
		u.indices(*iter, ind); (*m_spConsGF).indices(*iter, indCons);

	// 	adapt local algebra
		locU.resize(ind); locCF.resize(ind); locCons.resize(indCons);

	// 	read local values of u and contactForce
		GetLocalVector(locU, u); GetLocalVector(locCF, contactForce);
		GetLocalVector(locCons, *m_spConsGF);

	//	loop over DoFs in element and store all activeMultiIndex-pairs in vector
		size_t nFctElem = ind.num_fct();
		number NormNormal = VecLength(normal);

		for(size_t fct = 0; fct < nFctElem; ++fct)
		{
			size_t nDoFsPerFctElem = ind.num_dof(fct);
			for(size_t dof = 0; dof < nDoFsPerFctElem; ++dof)
			{
				MathVector<dim> locUDof;
				for(int i = 0; i < dim; ++i)
					locUDof[i] = locU(i, dof);

				//	locCF corresponds to the lagrange multiplier lambda,
				//	c.f. Hintermueller/Ito/Kunisch(2003)
				number locU_n_scaled = VecDot(locUDof, normal) / NormNormal;
				complementaryVal_n_scaled = locCF(fct ,dof) + locU_n_scaled - locCons(fct, dof);

				if (complementaryVal_n_scaled <= 1e-10){
					locCF(fct,dof) = 0.0;
				}
				else
				{
					//	create list of active global MultiIndex-pairs. Only those pairs should be attached
					//	which are not already a member of the 'activeSetGlob'-vector
					bool bAlreadyActive = false;
					for (vector<DoFIndex>::iterator itActiveInd = m_vActiveSetGlob.begin();
							itActiveInd < m_vActiveSetGlob.end(); ++itActiveInd)
					{
						if ((*itActiveInd)[0] == ind.index(fct, dof)
								&& (*itActiveInd)[1] == ind.comp(fct, dof))
							bAlreadyActive = true;
					}

					if (!bAlreadyActive)
					{
						DoFIndex newActiveIndex(ind.index(fct, dof), ind.comp(fct, dof));
						m_vActiveSetGlob.push_back(newActiveIndex);

						BlockRef(rhs[ind.index(fct, dof)], ind.comp(fct, dof)) = locCons.value(fct,dof);
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

//	builds up a vector of global (dof,fct)-pairs, which are active (wrt. an obstacle-constraint)
template <typename TDomain, typename TAlgebra>
bool ActiveSet<TDomain, TAlgebra>::active_index(function_type& u,
		function_type& rhs, function_type& contactForce, function_type& gap)
{
	//	note: the active-index search is restricted
	//	to constraint of the form u * n <= consGF

	if(!m_bCons)
		UG_THROW("No constraint set in ActiveSet \n");

	if (u.num_indices() != rhs.num_indices() || u.num_indices() != contactForce.num_indices() )
		UG_THROW("GridFunctions u, rhs and contactForce need to be defined "
				"for the same domain and of the same size in 'ActiveSet:active_index' \n");

	//	remember old ActiveSet for convergence check
	//	TODO: avoid this vector copy; m_vActiveSetGlobOld really necessary?
	m_vActiveSetGlobOld = m_vActiveSetGlob;
	m_vActiveSetGlob.resize(0);

	//	1.) get all subsets on which the 'gap'-gridfunction is defined!
	//	-> store them in vSubsetsOfContactForces
	m_vSubsetsOfContact.resize(0);
	//TODO: it is only necessary to loop over all BoundarySubsets!
	for (int si = 0; si < m_spDD->num_subsets(); si++){
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

	if (m_vSubsetsOfContact.size() == 0)
		UG_LOG("No subsets chosen as possible contact subsets. \n");

	UG_LOG("#sizeOfvSubsetsOfGap: " << m_vSubsetsOfContact.size() << "\n");

	//	2.) loop over all elements of the possible contact subsets
	for (vector<int>::iterator siContact = m_vSubsetsOfContact.begin();
			siContact != m_vSubsetsOfContact.end(); ++siContact)
	{
		//UG_LOG("siContact: " << *siContact << "\n");
		const int subsetDim = DimensionOfSubset(*m_spDD->subset_handler(), *siContact);
		//UG_LOG("subsetDim: " << subsetDim << "\n");

		//	3.) get localU out of u for each element and
		//	4.) store the active (global) dof,fct-pair in m_vActiveSetGlob
		switch(subsetDim)
		{
		case 0:
			break;
		case 1:
			ActiveIndexElem<Edge>
				(m_spDD->template begin<Edge>(*siContact),
						m_spDD->template end<Edge>(*siContact), u, rhs, contactForce);
			break;
		case 2:
			ActiveIndexElem<Triangle>
				(m_spDD->template begin<Triangle>(*siContact),
						m_spDD->template end<Triangle>(*siContact), u, rhs, contactForce);
			ActiveIndexElem<Quadrilateral>
				(m_spDD->template begin<Quadrilateral>(*siContact),
						m_spDD->template end<Quadrilateral>(*siContact), u, rhs, contactForce);
			break;
		default:
			UG_THROW("ActiveSet::active_index:"
				"SubsetDimension "<< subsetDim <<" (subset="<< *siContact <<") not supported.");
		}
	}

	if (m_vActiveSetGlob.size() > 0) return true;
	else return false;
}

//	sets a Dirichlet row for all active Indices
template <typename TDomain, typename TAlgebra>
void ActiveSet<TDomain, TAlgebra>::
adjust_matrix(matrix_type& mat, vector<SmartPtr<DoFIndex> > vActiveIndices)
{
	for (vector<SmartPtr<DoFIndex> >::iterator itActiveInd = vActiveIndices.begin();
			itActiveInd < vActiveIndices.end(); ++itActiveInd)
	{
		SetDirichletRow(mat, **itActiveInd);
	}
}

//	computes the contact forces for a given contact discretization
template <typename TDomain, typename TAlgebra>
void ActiveSet<TDomain, TAlgebra>::contactForces(function_type& contactforce,
		const function_type& u)
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

template <typename TDomain, typename TAlgebra>
template <typename TElem, typename TIterator>
void ActiveSet<TDomain, TAlgebra>::ass_lagrangeMatIElem(TIterator iterBegin,
		TIterator iterEnd, matrix_type& lagrangeMatI)
{
	typedef typename face_type_traits<dim>::face_type0 face_type0;
	typedef typename face_type_traits<dim>::face_type1 face_type1;
	typedef typename face_type_traits<dim>::DimFEGeo sideGeo;
	typename TDomain::position_accessor_type& aaPos = m_spDom->position_accessor();

	//	some storage
	MathVector<dim> sideCoPos[dim+1];
	vector<MathVector<dim> > vCorner;

	// 	local indices and local algebra
	LocalIndices ind; LocalMatrix locLagrangeMatI;

	// 	Loop over all elements on active subsets
	for(TIterator iter = iterBegin; iter != iterEnd; ++iter)
	{
		TElem* elem = *iter;

		// 	get global indices
		m_spDD->indices(elem, ind);

		locLagrangeMatI.resize(ind, ind);
		locLagrangeMatI = 0.0;

		//	reference object type and geometry
		ReferenceObjectID sideRoid = elem->reference_object_id();
		sideGeo geo(sideRoid, 3, LFEID(LFEID::LAGRANGE, dim, 1));

		//	prepare geometry for type and order
		try{
			geo.update_local(sideRoid, LFEID(LFEID::LAGRANGE, dim, 1), 1);
		}
		UG_CATCH_THROW("ActiveSet::ass_lagrangeMatIElem:"
				" Cannot update local values of finite element geometry.");

		//	get corners of element
		CollectCornerCoordinates(vCorner, *elem, aaPos);

		const DimReferenceElement<dim-1>& rRefElem
				= ReferenceElementProvider::get<dim-1>(sideRoid);

		for (int i = 0; i < (int)vCorner.size(); ++i)
			sideCoPos[i] = vCorner[rRefElem.id(dim-1, 0, 0, i)];

		MathVector<dim> normal;
		if ((int)vCorner.size() == dim)
			ElementNormal<face_type0, dim>(normal, sideCoPos);
		else
			ElementNormal<face_type1, dim>(normal, sideCoPos);
		UG_LOG("normal in ass_lagrangeMatIElem: " << normal << "\n");
		//number normNormal = VecLength(normal);

		try{
			geo.update(elem, sideCoPos, LFEID(LFEID::LAGRANGE, dim, 1), 1);
		}
		UG_CATCH_THROW("ActiveSet::ass_lagrangeMatIElem:"
						" Cannot update finite element geometry.");

		for (size_t ip = 0; ip < geo.num_ip(); ++ip)
			for(size_t sh = 0; sh < geo.num_sh(); ++sh)
				for (size_t i = 0; i < (size_t) dim; ++i)
					locLagrangeMatI(i, sh, i, sh) += 1.0/(geo.weight(ip) * geo.shape(ip, sh));
			// * normal[i] / normNormal;

		AddLocalMatrixToGlobal(lagrangeMatI, locLagrangeMatI);
	}
}

template <typename TDomain, typename TAlgebra>
void ActiveSet<TDomain, TAlgebra>::ass_lagrangeMatI(matrix_type& lagrangeMatI)
{
	for (vector<int>::iterator siContact = m_vSubsetsOfContact.begin();
				siContact != m_vSubsetsOfContact.end(); ++siContact)
	{
		UG_LOG("siContact: " << *siContact << "\n");
		const int subsetDim = DimensionOfSubset(*m_spDD->subset_handler(), *siContact);
		UG_LOG("subsetDim: " << subsetDim << "\n");

		//	3.) get localU out of u for each element and
		//	4.) store the active (local) dof,fct-pair in m_vActiveSet
		switch(subsetDim)
		{
		case 0:
			break;
		case 1:
			ass_lagrangeMatIElem<Edge>
				(m_spDD->template begin<Edge>(*siContact), m_spDD->template end<Edge>(*siContact), lagrangeMatI);
			break;
		case 2:
			ass_lagrangeMatIElem<Triangle>
				(m_spDD->template begin<Triangle>(*siContact), m_spDD->template end<Triangle>(*siContact), lagrangeMatI);
			ass_lagrangeMatIElem<Quadrilateral>
				(m_spDD->template begin<Quadrilateral>(*siContact), m_spDD->template end<Quadrilateral>(*siContact), lagrangeMatI);
			break;
		default:
			UG_THROW("ActiveSet::ass_lagrangeMat:"
				"SubsetDimension "<<subsetDim<<" (subset="<<*siContact<<") not supported.");
		}
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

		/*matrix_type lagrangeMatI;
		//SmartPtr<AssembledLinearOperator<algebra_type> > spLagrangeMatI;
		ass_lagrangeMatI(lagrangeMatI);*/

		SmartPtr<vector_type> spMat_u = u.clone_without_values();
		(*spMat_u).resize(u.size());

		#ifdef UG_PARALLEL
			MatMultDirect(*spMat_u, 1.0, mat, u);
			spMat_u->set_storage_type(u.get_storage_mask());
		#else
			MatMult(*spMat_u, 1.0, mat, u);
		#endif

		//SmartPtr<vector_type> spRes = u.clone_without_values();
		//(*spRes).resize(u.size());

		//	loop MultiIndex-pairs in activeSet-vector
		for (vector<DoFIndex>::iterator itActiveInd = m_vActiveSetGlob.begin();
				itActiveInd < m_vActiveSetGlob.end(); ++itActiveInd)
		{
			//	compute contact forces for active multiIndices
			//	contactForce = rhs - Mat * u;
			DoFRef(contactforce, *itActiveInd) = DoFRef(rhs, *itActiveInd) - DoFRef((*spMat_u), *itActiveInd);
		}

		/*#ifdef UG_PARALLEL
			MatMultDirect(contactforce, 1.0, lagrangeMatI, *spRes);
		#else
			MatMult(contactforce, 1.0, lagrangeMatI, *spRes);
		#endif*/

		UG_LOG("new contactforce-values computed \n");
		//UG_LOG("rhs updated \n");
	}
	else{
		UG_LOG("no active index in contactForcesRes \n");
	}
}


template <typename TDomain, typename TAlgebra>
template <typename TElem, typename TIterator>
bool ActiveSet<TDomain, TAlgebra>::ConvCheckElem(TIterator iterBegin,
		TIterator iterEnd, function_type& u, const function_type& lambda)
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
	vector<MathVector<dim> > vCorner;
	MathVector<dim> sideCoPos[dim+1], normal;

	// 	local indices and local algebra
	LocalIndices ind, indCons;
	LocalVector locU, locCons, locLambda;

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
		u.indices(*iter, ind); (*m_spConsGF).indices(*iter, indCons);

	// 	adapt local algebra
		locU.resize(ind); locCons.resize(indCons); locLambda.resize(ind);

	// 	read local values of u and contactForce
		GetLocalVector(locU, u); GetLocalVector(locCons, *m_spConsGF);
		GetLocalVector(locLambda, lambda);

		size_t nFctElem = ind.num_fct();
		number gap_value, locU_n_scaled, kktcond;
		number NormNormal = VecLength(normal);

		for(size_t fct = 0; fct < nFctElem; ++fct)
		{
			size_t nDoFsPerFctElem = ind.num_dof(fct);
			for(size_t dof = 0; dof < nDoFsPerFctElem; ++dof)
			{
				MathVector<dim> locUDof;
				for(int i = 0; i < dim; ++i)
					locUDof[i] = locU(i, dof);
				locU_n_scaled = VecDot(locUDof, normal) / NormNormal;

				gap_value =	locU_n_scaled - locCons(fct, dof);
				if (gap_value > 1e-06){
					//	i.e.: m_spConsGF < u
					//	constraint is violated
					return false;
				}

				kktcond = gap_value * locLambda(fct, dof);
				if (kktcond > 1e-06 || kktcond < -1e-06)
					return false;

			} // end(dof)
		} // end(fct)
	} // end(elem)

	return true;
}

template <typename TDomain, typename TAlgebra>
bool ActiveSet<TDomain, TAlgebra>::check_conv(function_type& u, const function_type& lambda,
		const size_t step)
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

	bool bConstraintViolated = false;
	for (vector<int>::iterator siContact = m_vSubsetsOfContact.begin();
			siContact != m_vSubsetsOfContact.end(); ++siContact)
	{
		const int subsetDim = DimensionOfSubset(*m_spDD->subset_handler(), *siContact);
		switch(subsetDim)
		{
		case 0:
			break;
		case 1:
			if (!ConvCheckElem<Edge>
				(m_spDD->template begin<Edge>(*siContact), m_spDD->template end<Edge>(*siContact),
						u, lambda))
			{bConstraintViolated = true;}

			break;
		case 2:
			if (!ConvCheckElem<Triangle>
				(m_spDD->template begin<Triangle>(*siContact),
						m_spDD->template end<Triangle>(*siContact), u, lambda))
			{bConstraintViolated = true;}

			if (!ConvCheckElem<Quadrilateral>
				(m_spDD->template begin<Quadrilateral>(*siContact),
						m_spDD->template end<Quadrilateral>(*siContact), u, lambda))
			{bConstraintViolated = true;}

			break;
		default:
			UG_THROW("ActiveSet::check_conv:"
				"SubsetDimension "<< subsetDim <<" (subset="<< *siContact <<") not supported.");
		}

		if (bConstraintViolated)
			return false;
	}

	//	check if activeSet has changed
	if (m_vActiveSetGlob.size() == m_vActiveSetGlobOld.size())
	{
		UG_LOG("Old and new active Set have the same number of members \n");

		vector<DoFIndex>::iterator itActiveInd = m_vActiveSetGlob.begin();

		for (vector<DoFIndex>::iterator itActiveIndOld = m_vActiveSetGlobOld.begin();
				itActiveIndOld < m_vActiveSetGlobOld.end(); ++itActiveIndOld)
		{
			if (*itActiveInd != *itActiveIndOld)
				return false;

			++itActiveInd;
		} // itActiveIndOld

		//	activeSet remains unchanged & constraint is fulfilled for all indices
		return true;
	}
	else{
		return false;
	}
}


template <typename TDomain, typename TAlgebra>
bool ActiveSet<TDomain, TAlgebra>::checkInequ(const matrix_type& mat,
				const vector_type& u,
				const vector_type& lambda,
				const vector_type& rhs)
{
	if (u.size() != lambda.size())
				UG_THROW("Temporarily u and lambda need to be "
						"of same size in ActiveSet:checkInequ \n");

	SmartPtr<vector_type> spMat_u = u.clone_without_values();
	SmartPtr<vector_type> spRes = u.clone_without_values();
	(*spMat_u).resize(u.size()); (*spRes).resize(u.size());

	#ifdef UG_PARALLEL
		MatMultDirect(*spMat_u, 1.0, mat, u);
		spMat_u->set_storage_type(u.get_storage_mask());
		spRes->set_storage_type(u.get_storage_mask());
	#else
		MatMult(*spMat_u, 1.0, mat, u);
	#endif


	for (size_t i = 0; i < u.size(); ++i)
	{
		//if (lambda[i] < -1e-06)
		//	return false;

		(*spRes)[i] = (*spMat_u)[i] - rhs[i]; //+ lambda[i] - rhs[i]; //
		UG_LOG("lambda["<< i << "]: " << lambda[i] <<", res[" << i << "]: " << (*spRes)[i] << "\n");
		//if ((*spRes)[i] > 1e-06)
		//	return false;
	}

	return true;
}

template <typename TDomain, typename TAlgebra>
void ActiveSet<TDomain, TAlgebra>::createVecOfPointers()
{
	m_vActiveSetGlobSP.resize(m_vActiveSetGlob.size());

	vector<DoFIndex>::iterator it = m_vActiveSetGlob.begin();

	for (vector<SmartPtr<DoFIndex> >::iterator itSP = m_vActiveSetGlobSP.begin();
				itSP < m_vActiveSetGlobSP.end(); ++itSP)
	{
		*itSP = &*it;
		++it;
	}

}

}; // namespace ug

#endif /* ACTIVE_SET_IMPL_H_ */
