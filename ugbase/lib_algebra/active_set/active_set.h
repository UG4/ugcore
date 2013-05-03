/*
 * active_set.h
 *
 *  Created on: 15.02.2013
 *      Author: raphaelprohl
 */

#ifndef ACTIVE_SET_H_
#define ACTIVE_SET_H_

using namespace std;

namespace ug {

template <typename TAlgebra>
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

	public:
	///	constructor
		ActiveSet() : m_bCons(false) {
			//	specifies the number of fcts
			value_type u_val;
			m_nrFcts = GetSize(u_val);
		};

	///	sets constraint/obstacle
		void set_constraint(vector_type& cons) {m_ConsVec = cons; m_bCons = true;}

		void prepare(vector_type& u);

	///	checks the distance to the prescribed obstacle/constraint
		bool check_dist_to_obs(vector_type& u);

	///	determines the active indices
		bool active_index(vector_type& u,
				vector_type& lambda);

	///	computes the contact forces lambda
		void contactForces(vector_type& lambda, const matrix_type& mat,
				const vector_type& u, const vector_type& rhs);

	///	checks if all constraints are fulfilled & the activeSet remained unchanged
		bool check_conv(const vector_type& u, const size_t step);

	///	creates a list of pointers to the active Indices for lua-registry
		void createVecOfPointers();

	///	method used for lua-call in order to pass the ActiveSet to assemble-funcs
		vector<SmartPtr<MultiIndex<2> > >  activeMultiIndices()
		{
			createVecOfPointers();
			return m_vActiveSetSP;
		};

	private:
		///	#fcts for value_type
		size_t m_nrFcts;

		/// vector describing a constraint
		vector_type m_ConsVec;
		bool m_bCons;

		///	vector of the current active set of MultiIndices (DoF,Fct)
		vector<MultiIndex<2> > m_vActiveSet;
		///	vector remembering the active set of MultiIndices (DoF,Fct)
		vector<MultiIndex<2> > m_vActiveSetOld;
		///	vector of pointers to active set needed for lua-call
		vector<SmartPtr<MultiIndex<2> > > m_vActiveSetSP;
};

} // namespace ug

#include "active_set_impl.h"

#endif /* ACTIVE_SET_H_ */
