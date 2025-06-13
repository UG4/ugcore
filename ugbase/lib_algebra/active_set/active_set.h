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

#ifndef ACTIVE_SET_H_
#define ACTIVE_SET_H_

#include "lib_algebra/active_set/lagrange_multiplier_disc_interface.h"
#include "lib_disc/function_spaces/grid_function.h"

using namespace std;

namespace ug {

/// Active Set method
/**
 *	The active Set method is a well-known method in constrained optimization theory.
 *	A general formulation for these problems reads
 *  \f{eqnarray*}{
 *		min_{x \in \mathbb{R}^n} f(x)
 *  \f}
 *	s.t.
 *  \f{eqnarray*}{
 *			c_i(x) = 0, \qquad i \in E \\
 *			c_i(x) \ge 0 \qquad i \in I,
 *  \f}
 *
 *	where \f$ f \f$, \f$ c_i \f$ are smooth, real-valued functions on a subset of \f$ \mathbb{R}^n \f$.
 *	\f$ I \f$ (set of inequality constraints) and \f$ E \f$ (set of equality constraints)
 *	are two finite sets of indices. \f$ f \f$ is called objective function.
 *
 *	The active Set \f$ A \f$ is defined as:
 *  \f{eqnarray*}{
 *		A(x) := E \cup \{ i \in I | c_i(x) = 0 \},
 *  \f}
 *	i.e. for \f$ i \in I \f$ the inequality constraint is said to be active, if \f$ c_i(x) = 0 \f$.
 *	Otherwise (\f$ c_i(x) > 0 \f$) it is called inactive.
 *
 *	A common approach to treat the inequality constraints is its reformulation as
 *	equations by using so called complementarity functions, see e.g. C.Hager und B. I. Wohlmuth:
 *	"Hindernis- und Kontaktprobleme" for a simple introduction into this topic.
 *	By means of complementarity functions, constraints of the form
 *  \f{eqnarray*}{
 *		a \ge 0, \, b \ge 0, \, a b = 0
 *  \f}
 *	with \f$ a, b \in \mathbb{R}^n \f$ can be reformulated as
 *  \f{eqnarray*}{
 *		C(a,b) = 0,
 *  \f}
 *	with \f$ C: \mathbb{R}^n x \mathbb{R}^n \to \mathbb{R}^n \f$ being an appropriate complementarity function.
 *	The value of \f$ C \f$ indicates, whether the index is active or inactive (see method 'active_index').
 *	Thus, it determines in every step the set of active indices, for which the original system
 *	of equations needs to be adapted. The influence of these active indices on the original system
 *	of equations ( \f$ K u = f \f$, with \f$ K \f$: system-matrix; \f$ u, f \f$ vectors) can be modelled by means of a
 *	lagrange multiplier '\f$ \lambda \f$'. \f$ \lambda \f$ can either be defined by the residual ( \f$ \lambda := f - K u \f$,
 *	see method 'residual_lagrange_mult') or by a problem-dependent computation (see method 'lagrange_multiplier').
 *
 *	In every Active Set step a linear or linearized system is solved. The algorithm stops when the active
 *	and inactive Set remains unchanged (see method 'check_conv').
 *
 *
 * References:
 * <ul>
 * <li> J. Nocedal and S. J. Wright. Numerical optimization.(2000)
 * </ul>
 *
 *  \tparam 	TDomain			Domain type
 *  \tparam 	TAlgebra		Algebra type
 */
template <typename TDomain, typename TAlgebra>
class ActiveSet
{
	public:
	///	Type of algebra
		typedef TAlgebra algebra_type;

	///	Type of algebra matrix
		typedef typename algebra_type::matrix_type matrix_type;

	///	Type of algebra vector
		typedef typename algebra_type::vector_type vector_type;

	///	Type of algebra value
		typedef typename vector_type::value_type value_type;

	///	Type of grid function
		typedef GridFunction<TDomain, TAlgebra> function_type;

	///	base element type of associated domain
		typedef typename domain_traits<TDomain::dim>::grid_base_object TBaseElem;

	///	domain dimension
		static const int dim = TDomain::dim;

	public:
	///	constructor
		ActiveSet() : m_bObs(false), m_spLagMultDisc(nullptr) {
			//	specifies the number of fcts
			//value_type u_val;
			//m_nrFcts = GetSize(u_val);  //ToDo: This field is only used in check_dist_to_obs which is commented out. Remove it completely?
		};

	///	sets obstacle gridfunction, which limits the solution u
		void set_obstacle(ConstSmartPtr<function_type> obs) {
			m_spObs = obs; m_bObs = true;
			//	if 'obs'-gridfunction is defined on a subset,
			//	which is not a boundary-subset -> UG_LOG
		}

	///	sets a discretization in order to compute the lagrange multiplier
		void set_lagrange_multiplier_disc(
				SmartPtr<ILagrangeMultiplierDisc<TDomain, function_type> > lagMultDisc)
		{m_spLagMultDisc = lagMultDisc;};

		void prepare(function_type& u);

	///	checks the distance to the prescribed obstacle/constraint
		//bool check_dist_to_obs(vector_type& u);

		template <typename TElem, typename TIterator>
		void active_index_elem(TIterator iterBegin,
				TIterator iterEnd, function_type& u,
				function_type& rhs, function_type& lagrangeMult);

	///	determines the active indices, stores them in a vector and sets dirichlet values in rhs for active indices
		bool active_index(function_type& u, function_type& rhs, function_type& lagrangeMult,
				function_type& gap);

		void set_dirichlet_rows(matrix_type& mat);

	///	computes the lagrange multiplier for a given disc
		void lagrange_multiplier(function_type& lagrangeMult, const function_type& u);

	///	computes the lagrange multiplier by means of
	/// the residuum (lagMult = rhs - mat * u)
		void residual_lagrange_mult(vector_type& lagMult, const matrix_type& mat,
				const vector_type& u, vector_type& rhs);

		template <typename TElem, typename TIterator>
		bool check_conv_elem(TIterator iterBegin,
				TIterator iterEnd, function_type& u, const function_type& lambda);

	///	checks if all constraints are fulfilled & the activeSet remained unchanged
		bool check_conv(function_type& u, const function_type& lambda, const size_t step);

	///	checks if all inequalities are fulfilled
		bool check_inequ(const matrix_type& mat, const vector_type& u,
					const vector_type& lagrangeMult, const vector_type& rhs);

		template <typename TElem, typename TIterator>
		void lagrange_mat_inv_elem(TIterator iterBegin,
				TIterator iterEnd, matrix_type& lagrangeMatInv);

		void lagrange_mat_inv(matrix_type& lagrangeMatInv);

	private:
		///	pointer to the DofDistribution on the whole domain
		SmartPtr<DoFDistribution> m_spDD;
		SmartPtr<TDomain> m_spDom;

		///	number of functions
		//size_t m_nrFcts; //ToDo: This field is only used in check_dist_to_obs which is commented out. Remove it completely?

		/// pointer to a gridfunction describing an obstacle-constraint
		ConstSmartPtr<function_type> m_spObs;
		bool m_bObs;

		///	pointer to a lagrangeMultiplier-Disc
		SmartPtr<ILagrangeMultiplierDisc<TDomain, function_type> > m_spLagMultDisc;

		///	vector of possible active subsets
		vector<int> m_vActiveSubsets;

		/*template <typename TElem>
		struct activeElemAndLocInd{
			TElem* pElem; 						// pointer to active elem
			vector<vector<size_t> > vlocInd; 	// vector of local active indices
		};*/

		///	vector of the current active set of global DoFIndices
		vector<DoFIndex> m_vActiveSetGlob;
		///	vector remembering the active set of global DoFIndices
		vector<DoFIndex> m_vActiveSetGlobOld;
};

} // namespace ug

#include "active_set_impl.h"

#endif /* ACTIVE_SET_H_ */
