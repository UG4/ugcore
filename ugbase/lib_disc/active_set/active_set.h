/*
 * active_set.h
 *
 *  Created on: 15.02.2013
 *      Author: raphaelprohl
 */

#ifndef ACTIVE_SET_H_
#define ACTIVE_SET_H_

namespace ug {

//	TODO: maybe the dependency on TDomain could be removed.
//	In that case move the class in lib_algebra!
template <typename TDomain, typename TAlgebra>
class ActiveSet
{
	public:
	///	Type of domain
		typedef TDomain domain_type;

	///	world Dimension
		static const int dim = domain_type::dim;

	///	Type of position coordinates (e.g. position_type)
		typedef typename domain_type::position_type position_type;

	///	Type of algebra
		typedef TAlgebra algebra_type;

	///	Type of algebra matrix
		typedef typename algebra_type::matrix_type matrix_type;

	///	Type of algebra vector
		typedef typename algebra_type::vector_type vector_type;

	public:
	///	constructor
		ActiveSet();

		void set_constraint(vector_type& cons) {
			//	note: temporarily only constraints
			//	which do not differ for the different fcts are valid here!!!
			m_ConsVec = cons; m_bCons = true;
		}

		void prepare(vector_type& u);

		bool active_index(vector_type& u,
				vector_type& lambda);

		void comp_lambda(matrix_type& mat, vector_type& u, vector_type& rhs, vector_type& lambda);

		vector<size_t> get_activeSet() { return m_vActiveSet;};

		//SmartPtr<vector_type> get_constraint() { return m_spConsVec;};

		bool check_conv();

	private:
		//	smart pointer to a vector describing a constraint
		//SmartPtr<vector_type> m_spConsVec;
		vector_type m_ConsVec;
		bool m_bCons;

		//	vector remembering the active set of DoFs
		vector<size_t> m_vActiveSet;
		vector<size_t> m_vInactiveSet;
		vector<size_t> m_vActiveSetOld;
};

} // namespace ug

#include "active_set_impl.h";

#endif /* ACTIVE_SET_H_ */
