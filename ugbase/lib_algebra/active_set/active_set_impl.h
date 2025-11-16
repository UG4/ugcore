/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
 * Author: Raphael Prohl
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#ifndef ACTIVE_SET_IMPL_H_
#define ACTIVE_SET_IMPL_H_

#include "active_set.h"
#include "lib_disc/common/geometry_util.h"
#include "lib_disc/spatial_disc/disc_util/fe_geom.h"

using namespace std;

namespace ug {

template <int dim> struct face_type_traits
{
	using face_type0 = void;
	using face_type1 = void;
	using DimFEGeo = void;
};

template <> struct face_type_traits<1>
{
	using face_type0 = ReferenceVertex;
	using face_type1 = ReferenceVertex;
	using DimFEGeo = DimFEGeometry<1, 1>;
};

template <> struct face_type_traits<2>
{
	using face_type0 = ReferenceEdge;
	using face_type1 = ReferenceEdge;
	using DimFEGeo = DimFEGeometry<2, 1>;
};

template <> struct face_type_traits<3>
{
	using face_type0 = ReferenceTriangle;
	using face_type1 = ReferenceQuadrilateral;
	using DimFEGeo = DimFEGeometry<3, 2>;
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
	//	STILL IN PROGRESS: u sollte hier reference-position + Startl�sung sein!
	value_type dist;

	bool geometry_cut_by_cons = false;

	for(size_t i = 0; i < u.size(); i++)
	{
		UG_LOG("u(" << i << "):" << u[i] << "\n");
		UG_LOG("m_spObs(" << i << "):" << (*m_spObs)[i] << "\n");
		dist = (*m_spObs)[i] - u[i];
		UG_LOG("dist:" << dist << "\n");
		//TODO: anstatt u muss hier die geometrische Info einflie�en!
		for (size_t fct = 0; fct < m_nrFcts; fct++)
		{
			if (BlockRef(dist,fct) < 0.0) // i.e.: u < m_spObs
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
void ActiveSet<TDomain, TAlgebra>::active_index_elem(TIterator iterBegin,
		TIterator iterEnd,
		function_type& u,
		function_type& rhs,
		function_type& lagrangeMult)
{
	// 	check if at least an element exists, else return
	if(iterBegin == iterEnd) return;

	int elemCounter = 0;

	static constexpr int dim = function_type::dim;
	using domain_type = typename function_type::domain_type;
	using face_type0 = typename face_type_traits<dim>::face_type0;
	using face_type1 = typename face_type_traits<dim>::face_type1;

	//	get position accessor
	typename domain_type::position_accessor_type& aaPos
			= u.domain()->position_accessor();

	//	storage for corner coordinates
	vector<MathVector<dim> > vCorner;
	vector<MathVector<dim> > vSideCoPos;

	MathVector<dim> normal;

	// 	local indices and local algebra
	LocalIndices ind, indObs;
	LocalVector locU, locLM, locObs;
	number complementaryVal;

	//	TODO: it could be more efficient to store the active elements
	//	and its active local indices as well

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
		int nCorner = (int)vCorner.size();
		for (int i = 0; i < nCorner; ++i)
			vSideCoPos.push_back(vCorner[rRefElem.id(dim-1, 0, 0, i)]);

		if (nCorner == dim)
		{
			ElementNormal<face_type0, dim>(normal, &vSideCoPos[0]);
			UG_LOG("face_type0 \n");
		}
		else
		{
			ElementNormal<face_type1, dim>(normal, &vSideCoPos[0]);
			UG_LOG("face_type1 \n");
		}

		/*for (int i = 0; i < (int)vCorner.size(); ++i)
			UG_LOG("sideCoPos: " << sideCoPos[i] << "\n");
		UG_LOG("normal: " << normal << "\n");*/

		elemCounter++;

	// 	get global indices
		u.indices(*iter, ind); (*m_spObs).indices(*iter, indObs);

	// 	adapt local algebra
		locU.resize(ind); locLM.resize(ind); locObs.resize(indObs);

	// 	read local values of u and lagrangeMult
		GetLocalVector(locU, u); GetLocalVector(locLM, lagrangeMult);
		GetLocalVector(locObs, *m_spObs);

	//	loop over DoFs in element and store all activeMultiIndex-pairs in vector
		size_t nrFctElem = ind.num_fct();
		number normOfNormal = VecLength(normal);
		MathVector<dim> locUDof;

		for(size_t fct = 0; fct < nrFctElem; ++fct)
		{
			size_t nrDoFsPerFctElem = ind.num_dof(fct);
			for(size_t dof = 0; dof < nrDoFsPerFctElem; ++dof)
			{
				for(int i = 0; i < dim; ++i)
					locUDof[i] = locU(i, dof);

				//	locLM corresponds to the lagrange multiplier lambda,
				//	c.f. Hintermueller/Ito/Kunisch(2003)
				number locUNormal = VecDot(locUDof, normal) / normOfNormal;
				complementaryVal = locLM(fct ,dof) + locUNormal - locObs(fct, dof);

				if (complementaryVal <= 1e-10){
					locLM(fct,dof) = 0.0;
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

						BlockRef(rhs[ind.index(fct, dof)], ind.comp(fct, dof)) = locObs.value(fct,dof);
					}
				}
			} // end(dof)
		} // end(fct)

		// 	send local to global lagrangeMult
		//AddLocalVector(lagrangeMult, locLM);

		UG_LOG("#activeDoFFctPairs global: " << m_vActiveSetGlob.size() << "\n");
	} // end(elem)
	UG_LOG("#elems: " << elemCounter << "\n");

}

//	builds up a vector of global (dof,fct)-pairs, which are active (wrt. an obstacle-constraint)
//	and sets dirichlet values in rhs for active indices
template <typename TDomain, typename TAlgebra>
bool ActiveSet<TDomain, TAlgebra>::active_index(function_type& u,
		function_type& rhs, function_type& lagrangeMult, function_type& gap)
{
	//	note: the active-index search is restricted
	//	to constraint of the form u * n <= consGF
	if(!m_bObs)
		UG_THROW("No constraint set in ActiveSet \n");

	if (u.num_indices() != rhs.num_indices() || u.num_indices() != lagrangeMult.num_indices() )
		UG_THROW("GridFunctions u, rhs and lagrangeMult need to be defined "
				"for the same domain and of the same size in 'ActiveSet:active_index' \n");

	//	remember old ActiveSet for convergence check
	//	TODO: avoid this vector copy; m_vActiveSetGlobOld really necessary?
	m_vActiveSetGlobOld = m_vActiveSetGlob;
	m_vActiveSetGlob.resize(0);

	//	1.) get all subsets on which the 'gap'-gridfunction is defined!
	//	-> store them in vSubsetsOflagrangeMults
	m_vActiveSubsets.resize(0);
	//TODO: it is only necessary to loop over all BoundarySubsets!
	for (int si = 0; si < m_spDD->num_subsets(); si++){
		for (size_t fct = 0; fct < gap.num_fct(); fct++)
			if (gap.is_def_in_subset(fct, si))
			{
				m_vActiveSubsets.push_back(si);
				//	'break' is necessary to ensure that 'si' is
				//	only added once when several fcts of
				//	'gap' are defined in subset 'si'!
				break;
			}
	}

	if (m_vActiveSubsets.size() == 0)
		UG_LOG("No subsets chosen as possible active subsets. \n");

	UG_LOG("#sizeOfvActiveSubsets: " << m_vActiveSubsets.size() << "\n");

	//	2.) loop over all elements of the possible active subsets
	for (vector<int>::iterator activeSI = m_vActiveSubsets.begin();
			activeSI != m_vActiveSubsets.end(); ++activeSI)
	{
		//UG_LOG("activeSI: " << *activeSI << "\n");
		const int subsetDim = DimensionOfSubset(*m_spDD->subset_handler(), *activeSI);
		//UG_LOG("subsetDim: " << subsetDim << "\n");

		//	3.) get localU out of u for each element and
		//	4.) store the active global (dof,fct)-pair in m_vActiveSetGlob
		switch(subsetDim)
		{
		case 0:
			break;
		case 1:
			active_index_elem<RegularEdge>
				(m_spDD->template begin<RegularEdge>(*activeSI),
						m_spDD->template end<RegularEdge>(*activeSI), u, rhs, lagrangeMult);
			break;
		case 2:
			active_index_elem<Triangle>
				(m_spDD->template begin<Triangle>(*activeSI),
						m_spDD->template end<Triangle>(*activeSI), u, rhs, lagrangeMult);
			active_index_elem<Quadrilateral>
				(m_spDD->template begin<Quadrilateral>(*activeSI),
						m_spDD->template end<Quadrilateral>(*activeSI), u, rhs, lagrangeMult);
			break;
		default:
			UG_THROW("ActiveSet::active_index:"
				"SubsetDimension "<< subsetDim <<" (subset="<< *activeSI <<") not supported.");
		}
	}

	if (m_vActiveSetGlob.size() > 0) return true;
	else return false;
}

//	sets a Dirichlet row for all active Indices
template <typename TDomain, typename TAlgebra>
void ActiveSet<TDomain, TAlgebra>::
set_dirichlet_rows(matrix_type& mat)
{
	for (vector<DoFIndex>::iterator itActiveInd = m_vActiveSetGlob.begin();
			itActiveInd < m_vActiveSetGlob.end(); ++itActiveInd){
		SetDirichletRow(mat, *itActiveInd);
	}
}

//	computes the lagrange multiplier for a given discretization
template <typename TDomain, typename TAlgebra>
void ActiveSet<TDomain, TAlgebra>::lagrange_multiplier(function_type& lagrangeMult,
		const function_type& u)
{
	//	check that lagrange multiplier disc is set
	if (m_spLagMultDisc.invalid())
		UG_THROW("No discretization set to compute the lagrange multiplier (in "
					"ActiveSet:lagrange_multiplier) \n");

	if(m_vActiveSetGlob.size() != 0.0)
	{
		UG_LOG("activeDoFs in ActiveSet:lagrange_multiplier " << m_vActiveSetGlob.size() << "\n");
		//TODO: pass localActiveElemAndIndex here!!!
		m_spLagMultDisc->lagrange_multiplier(lagrangeMult, u, m_vActiveSetGlob, m_vActiveSubsets);
	}
}

template <typename TDomain, typename TAlgebra>
template <typename TElem, typename TIterator>
void ActiveSet<TDomain, TAlgebra>::lagrange_mat_inv_elem(TIterator iterBegin,
		TIterator iterEnd, matrix_type& lagrangeMatInv)
{
	using face_type0 = typename face_type_traits<dim>::face_type0;
	using face_type1 = typename face_type_traits<dim>::face_type1;
	using sideGeo = typename face_type_traits<dim>::DimFEGeo;
	typename TDomain::position_accessor_type& aaPos = m_spDom->position_accessor();

	//	some storage
	MathVector<dim> normal;
	vector<MathVector<dim> > vCorner;
	vector<MathVector<dim> > vSideCoPos;

	// 	local indices and local algebra
	LocalIndices ind; LocalMatrix loclagrangeMatInv;

	// 	Loop over all elements on active subsets
	for(TIterator iter = iterBegin; iter != iterEnd; ++iter)
	{
		TElem* elem = *iter;

		// 	get global indices
		m_spDD->indices(elem, ind);

		loclagrangeMatInv.resize(ind, ind);
		loclagrangeMatInv = 0.0;

		//	reference object type and geometry
		ReferenceObjectID sideRoid = elem->reference_object_id();
		sideGeo geo(sideRoid, 3, LFEID(LFEID::LAGRANGE, dim, 1));

		//	prepare geometry for type and order
		try{
			geo.update_local(sideRoid, LFEID(LFEID::LAGRANGE, dim, 1), 1);
		}
		UG_CATCH_THROW("ActiveSet::lagrange_mat_inv_elem:"
				" Cannot update local values of finite element geometry.");

		//	get corners of element
		CollectCornerCoordinates(vCorner, *elem, aaPos);

		const DimReferenceElement<dim-1>& rRefElem
				= ReferenceElementProvider::get<dim-1>(sideRoid);

		int nCorner = (int)vCorner.size();
		for (int i = 0; i < nCorner; ++i)
			vSideCoPos.push_back(vCorner[rRefElem.id(dim-1, 0, 0, i)]);

		if (nCorner == dim)
			ElementNormal<face_type0, dim>(normal, &vSideCoPos[0]);
		else
			ElementNormal<face_type1, dim>(normal, &vSideCoPos[0]);

		UG_LOG("normal in lagrange_mat_inv_elem: " << normal << "\n");
		//number normOfNormal = VecLength(normal);

		try{
			geo.update(elem, &vSideCoPos[0], LFEID(LFEID::LAGRANGE, dim, 1), 1);
		}
		UG_CATCH_THROW("ActiveSet::lagrange_mat_inv_elem:"
						" Cannot update finite element geometry.");

		for (size_t ip = 0; ip < geo.num_ip(); ++ip)
			for(size_t sh = 0; sh < geo.num_sh(); ++sh)
				for (size_t i = 0; i < (size_t) dim; ++i)
					loclagrangeMatInv(i, sh, i, sh) += 1.0/(geo.weight(ip) * geo.shape(ip, sh));
			// * normal[i] / normOfNormal;

		AddLocalMatrixToGlobal(lagrangeMatInv, loclagrangeMatInv);
	}
}

template <typename TDomain, typename TAlgebra>
void ActiveSet<TDomain, TAlgebra>::lagrange_mat_inv(matrix_type& lagrangeMatInv)
{
	for (vector<int>::iterator activeSI = m_vActiveSubsets.begin();
				activeSI != m_vActiveSubsets.end(); ++activeSI)
	{
		UG_LOG("activeSI: " << *activeSI << "\n");
		const int subsetDim = DimensionOfSubset(*m_spDD->subset_handler(), *activeSI);
		UG_LOG("subsetDim: " << subsetDim << "\n");

		switch(subsetDim)
		{
		case 0:
			break;
		case 1:
			lagrange_mat_inv_elem<RegularEdge>
				(m_spDD->template begin<RegularEdge>(*activeSI), m_spDD->template end<RegularEdge>(*activeSI), lagrangeMatInv);
			break;
		case 2:
			lagrange_mat_inv_elem<Triangle>
				(m_spDD->template begin<Triangle>(*activeSI), m_spDD->template end<Triangle>(*activeSI), lagrangeMatInv);
			lagrange_mat_inv_elem<Quadrilateral>
				(m_spDD->template begin<Quadrilateral>(*activeSI), m_spDD->template end<Quadrilateral>(*activeSI), lagrangeMatInv);
			break;
		default:
			UG_THROW("ActiveSet::lagrange_mat_inv:"
				"SubsetDimension "<<subsetDim<<" (subset="<<*activeSI<<") not supported.");
		}
	}
}

template <typename TDomain, typename TAlgebra>
void ActiveSet<TDomain, TAlgebra>::residual_lagrange_mult(vector_type& lagMult,
		const matrix_type& mat,
		const vector_type& u,
		vector_type& rhs)
{
	// 	only if some indices are active the lagrange multiplier is computed
	if(m_vActiveSetGlob.size() != 0.0)
	{
		if (u.size() != lagMult.size())
			UG_THROW("Temporarily u and lagMult need to be "
					"of same size in ActiveSet:residual_lagrange_mult \n");

		/*matrix_type lagrangeMatInv;
		//SmartPtr<AssembledLinearOperator<algebra_type> > splagrangeMatInv;
		ass_lagrangeMatInv(lagrangeMatInv);*/

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
			//	compute lagrange multiplier for active multiIndices
			//	lagMult = rhs - Mat * u;
			DoFRef(lagMult, *itActiveInd) = DoFRef(rhs, *itActiveInd) - DoFRef((*spMat_u), *itActiveInd);
		}

		/*#ifdef UG_PARALLEL
			MatMultDirect(lagMult, 1.0, lagrangeMatInv, *spRes);
		#else
			MatMult(lagMult, 1.0, lagrangeMatInv, *spRes);
		#endif*/

		UG_LOG("new lagMult-values computed \n");
		//UG_LOG("rhs updated \n");
	}
	else{
		UG_LOG("no active index in residual_lagrange_mult \n");
	}
}


template <typename TDomain, typename TAlgebra>
template <typename TElem, typename TIterator>
bool ActiveSet<TDomain, TAlgebra>::check_conv_elem(TIterator iterBegin,
		TIterator iterEnd, function_type& u, const function_type& lambda)
{
	static constexpr int dim = function_type::dim;
	using domain_type = typename function_type::domain_type;
	using face_type0 = typename face_type_traits<dim>::face_type0;
	using face_type1 = typename face_type_traits<dim>::face_type1;

	//	get position accessor
	typename domain_type::position_accessor_type& aaPos
			= u.domain()->position_accessor();

	//	storage for corner coordinates
	vector<MathVector<dim> > vCorner;
	vector<MathVector<dim> > vSideCoPos;
	MathVector<dim> normal;

	// 	local indices and local algebra
	LocalIndices ind, indObs;
	LocalVector locU, locObs, locLambda;

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
		int nCorner = (int)vCorner.size();
		for (int i = 0; i < nCorner; ++i)
			vSideCoPos.push_back(vCorner[rRefElem.id(dim-1, 0, 0, i)]);

		if (nCorner == dim)
			ElementNormal<face_type0, dim>(normal, &vSideCoPos[0]);
		else
			ElementNormal<face_type1, dim>(normal, &vSideCoPos[0]);

	// 	get global indices
		u.indices(*iter, ind); (*m_spObs).indices(*iter, indObs);

	// 	adapt local algebra
		locU.resize(ind); locObs.resize(indObs); locLambda.resize(ind);

	// 	read local values of u and lagrangeMult
		GetLocalVector(locU, u); GetLocalVector(locObs, *m_spObs);
		GetLocalVector(locLambda, lambda);

		size_t nrFctElem = ind.num_fct();
		number gapValue, locUNormal, kktcond;
		number normOfNormal = VecLength(normal);
		MathVector<dim> locUDof;

		for(size_t fct = 0; fct < nrFctElem; ++fct)
		{
			size_t nrDoFsPerFctElem = ind.num_dof(fct);
			for(size_t dof = 0; dof < nrDoFsPerFctElem; ++dof)
			{
				for(int i = 0; i < dim; ++i)
					locUDof[i] = locU(i, dof);
				locUNormal = VecDot(locUDof, normal) / normOfNormal;

				gapValue =	locUNormal - locObs(fct, dof);
				if (gapValue > 1e-06){
					//	i.e.: m_spObs < u
					//	constraint is violated
					return false;
				}

				kktcond = gapValue * locLambda(fct, dof);
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
	//	1. 	Did some multiIndices change from 'active' to 'inactive' or vice versa
	//		in the last iteration-step?
	//	2. 	Is the constraint violated for any DoFIndex?

	UG_LOG(m_vActiveSetGlob.size() << " indices are active at the begin "
			"of step " << step << " ! \n");

	//	check if activeSet has changed
	if (m_vActiveSetGlob == m_vActiveSetGlobOld)
	{
		//	check if constraint is fulfilled
		bool bConstraintViolated = false;
		for (vector<int>::iterator activeSI = m_vActiveSubsets.begin();
				activeSI != m_vActiveSubsets.end(); ++activeSI)
		{
			const int subsetDim = DimensionOfSubset(*m_spDD->subset_handler(), *activeSI);
			switch(subsetDim)
			{
			case 0:
				break;
			case 1:
				if (!check_conv_elem<RegularEdge>
					(m_spDD->template begin<RegularEdge>(*activeSI), m_spDD->template end<RegularEdge>(*activeSI),
							u, lambda))
				{bConstraintViolated = true;}

				break;
			case 2:
				if (!check_conv_elem<Triangle>
					(m_spDD->template begin<Triangle>(*activeSI),
							m_spDD->template end<Triangle>(*activeSI), u, lambda))
				{bConstraintViolated = true;}

				if (!check_conv_elem<Quadrilateral>
					(m_spDD->template begin<Quadrilateral>(*activeSI),
							m_spDD->template end<Quadrilateral>(*activeSI), u, lambda))
				{bConstraintViolated = true;}

				break;
			default:
				UG_THROW("ActiveSet::check_conv:"
					"SubsetDimension "<< subsetDim <<" (subset="<< *activeSI <<") not supported.");
			}

			if (bConstraintViolated)
				return false;
		}

		//	activeSet remains unchanged & constraint is fulfilled for all indices
		return true;
	}
	else{
		return false;
	}
}


template <typename TDomain, typename TAlgebra>
bool ActiveSet<TDomain, TAlgebra>::check_inequ(const matrix_type& mat,
				const vector_type& u,
				const vector_type& lambda,
				const vector_type& rhs)
{
	if (u.size() != lambda.size())
		UG_THROW("Temporarily u and lambda need to be "
				"of same size in ActiveSet:check_inequ \n");

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

}; // namespace ug

#endif /* ACTIVE_SET_IMPL_H_ */
