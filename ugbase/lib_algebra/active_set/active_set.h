/*
 * active_set.h
 *
 *  Created on: 15.02.2013
 *      Author: raphaelprohl
 */

#ifndef ACTIVE_SET_H_
#define ACTIVE_SET_H_

#include "lib_algebra/active_set/lagrange_multiplier_disc_interface.h"
#include "lib_disc/function_spaces/grid_function.h"

using namespace std;

namespace ug {

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
		typedef typename domain_traits<TDomain::dim>::geometric_base_object TBaseElem;

	///	domain dimension
		static const int dim = TDomain::dim;

	public:
	///	constructor
		ActiveSet() : m_bObs(false), m_spLagMultDisc(NULL) {
			//	specifies the number of fcts
			value_type u_val;
			m_nrFcts = GetSize(u_val);
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

		void set_dirichlet_rows(matrix_type& mat, vector<SmartPtr<DoFIndex> > vActiveIndices);

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

	///	method used for lua-call in order to pass the ActiveSet to assemble-funcs
		vector<SmartPtr<DoFIndex> > active_dof_indices()
		{
			create_vec_of_pointers();
			return m_vActiveSetGlobSP;
		};

	private:
	///	creates a list of pointers to the active Indices for lua-registry
		void create_vec_of_pointers();

	private:
		///	pointer to the DofDistribution on the whole domain
		SmartPtr<DoFDistribution> m_spDD;
		SmartPtr<TDomain> m_spDom;

		///	number of functions
		size_t m_nrFcts;

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
		///	vector of pointers to active set needed for lua-call
		vector<SmartPtr<DoFIndex> > m_vActiveSetGlobSP;
};

} // namespace ug

#include "active_set_impl.h"

#endif /* ACTIVE_SET_H_ */
