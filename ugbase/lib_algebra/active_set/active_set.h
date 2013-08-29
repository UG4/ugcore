/*
 * active_set.h
 *
 *  Created on: 15.02.2013
 *      Author: raphaelprohl
 */

#ifndef ACTIVE_SET_H_
#define ACTIVE_SET_H_

#include "lib_disc/spatial_disc/elem_disc/contact_boundary/contact_interface.h"
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
		ActiveSet() : m_bCons(false), m_spContactDisc(NULL) {
			//	specifies the number of fcts
			value_type u_val;
			m_nrFcts = GetSize(u_val);
		};

	///	sets constraint/obstacle
		void set_constraint(ConstSmartPtr<function_type> cons) {
			m_spConsGF = cons; m_bCons = true;
			//	if cons-GF is defined on a subset,
			//	which is not a boundary-subset -> UG_LOG
		}

	///	sets a contact discretization in order to use contact-disc dependent funcs
		void set_contactDisc(SmartPtr<IContactDisc<TDomain, function_type> > contact)
		{m_spContactDisc = contact;};

		void prepare(function_type& u);

	///	checks the distance to the prescribed obstacle/constraint
		bool check_dist_to_obs(vector_type& u);

		template <typename TElem, typename TIterator>
		void ActiveIndexElem(TIterator iterBegin,
				TIterator iterEnd, function_type& u,
				function_type& rhs, function_type& contactForce);

	///	determines the active indices
		bool active_index(function_type& u, function_type& rhs, function_type& contactForce,
				function_type& gap);

	///	computes the contact forces for a given contact disc
		void contactForces(function_type& contactForce, const function_type& u);

	///	computes the contact forces via the residuum
		void contactForcesRes(vector_type& contactForce, const matrix_type& mat,
				const vector_type& u, vector_type& rhs);

		template <typename TElem, typename TIterator>
		bool ConvCheckElem(TIterator iterBegin,
				TIterator iterEnd, function_type& u, const function_type& lambda);

	///	checks if all constraints are fulfilled & the activeSet remained unchanged
		bool check_conv(function_type& u, const function_type& lambda, const size_t step);

		bool checkEqu(const matrix_type& mat,
						const vector_type& u,
						const vector_type& contactforce,
						const vector_type& rhs);

		template <typename TElem, typename TIterator>
		void ass_lagrangeMatIElem(TIterator iterBegin,
				TIterator iterEnd, matrix_type& lagrangeMatI);

		void ass_lagrangeMatI(matrix_type& lagrangeMatI);

	///	method used for lua-call in order to pass the ActiveSet to assemble-funcs
		vector<SmartPtr<MultiIndex<2> > >  activeMultiIndices()
		{
			createVecOfPointers();
			return m_vActiveSetGlobSP;
		};

	private:
	///	creates a list of pointers to the active Indices for lua-registry
		void createVecOfPointers();

	private:
		///	point to the DofDistribution on the whole domain
		SmartPtr<DoFDistribution> m_spDD;
		SmartPtr<TDomain> m_spDom;

		///	#fcts for value_type
		size_t m_nrFcts;

		/// pointer to a gridfunction describing a constraint
		ConstSmartPtr<function_type> m_spConsGF;
		bool m_bCons;

		///	pointer to a contact-Disc
		SmartPtr<IContactDisc<TDomain, function_type> > m_spContactDisc;

		///	vector of possible contact subsets
		vector<int> m_vSubsetsOfContact;

		///	vector of the current active set of MultiIndices (DoF,Fct)
		vector<MultiIndex<2> > m_vActiveSetGlob;
		///	vector remembering the active set of MultiIndices (DoF,Fct)
		vector<MultiIndex<2> > m_vActiveSetGlobOld;
		///	vector of pointers to active set needed for lua-call
		vector<SmartPtr<MultiIndex<2> > > m_vActiveSetGlobSP;
};

} // namespace ug

#include "active_set_impl.h"

#endif /* ACTIVE_SET_H_ */
