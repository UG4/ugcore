/*
 * p1_continuity_constraints.h
 *
 *  Created on: 01.03.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__CONSTRAINTS__CONTINUITY_CONSTRAINTS__P1_CONTINUITY_CONSTRAINTS__
#define __H__UG__LIB_DISC__SPATIAL_DISC__CONSTRAINTS__CONTINUITY_CONSTRAINTS__P1_CONTINUITY_CONSTRAINTS__

#include "lib_disc/assemble_interface.h"
#include "lib_grid/algorithms/geom_obj_util/vertex_util.h"
#include "lib_disc/spatial_disc/constraints/constraint_base.h"

namespace ug {

template <typename TDomain, typename TAlgebra>
class SymP1Constraints
	: public ConstraintBase<TDomain, TAlgebra,
	  	  	  	  	  	  	  SymP1Constraints<TDomain, TAlgebra> >
{
	public:
	// 	Algebra type
		typedef TAlgebra algebra_type;

	// 	Type of algebra matrix
		typedef typename algebra_type::matrix_type matrix_type;

	// 	Type of algebra vector
		typedef typename algebra_type::vector_type vector_type;

	public:
		virtual int type() {return CT_CONSTRAINTS;}

		template <typename TDD>
		void adjust_defect(vector_type& d, const vector_type& u,
		                   ConstSmartPtr<TDD> dd, number time = 0.0)
		{
			UG_THROW_FATAL("not implemented.");
		}

		template <typename TDD>
		void adjust_rhs(vector_type& rhs, const vector_type& u,
		                ConstSmartPtr<TDD> dd, number time = 0.0)
		{
			UG_THROW_FATAL("not implemented.");
		}

		template <typename TDD>
		void adjust_jacobian(matrix_type& J, const vector_type& u,
		                     ConstSmartPtr<TDD> dd, number time = 0.0)
		{
		//  \todo: Implement correctly
		//	dummy for rhs
			vector_type rhsDummy; rhsDummy.resize(u.size());

			adjust_linear(J, rhsDummy, dd, time);
		}

		template <typename TDD>
		void adjust_linear(matrix_type& mat, vector_type& rhs,
		                   ConstSmartPtr<TDD> dd, number time);

		template <typename TDD>
		void adjust_solution(vector_type& u, ConstSmartPtr<TDD> dd,
		                     number time);

	protected:
		void SplitAddRow(matrix_type& A	,
		                 std::vector<size_t> & constrainedIndex,
		                 std::vector<std::vector<size_t> >& vConstrainingIndices);

		void SetInterpolation(matrix_type& A,
		                      std::vector<size_t> & constrainedIndex,
		                      std::vector<std::vector<size_t> >& vConstrainingIndices);

		void HandleRhs(vector_type& rhs,
		               std::vector<size_t> & constrainedIndex,
		               std::vector<std::vector<size_t> >& vConstrainingIndices);

		void InterpolateValues(vector_type& u,
		                       std::vector<size_t> & constrainedIndex,
		                       std::vector<std::vector<size_t> >& vConstrainingIndices);
};



template <typename TDomain, typename TAlgebra>
class OneSideP1Constraints
	: public ConstraintBase<TDomain, TAlgebra,
	  	  	  	  	  	  OneSideP1Constraints<TDomain, TAlgebra> >
{
	public:
	// 	Algebra type
		typedef TAlgebra algebra_type;

	// 	Type of algebra matrix
		typedef typename algebra_type::matrix_type matrix_type;

	// 	Type of algebra vector
		typedef typename algebra_type::vector_type vector_type;

	public:
		virtual int type() {return CT_CONSTRAINTS;}

		template <typename TDD>
		void adjust_jacobian(matrix_type& J,
		                     const vector_type& u,
		                     ConstSmartPtr<TDD> dd,
		                     number time = 0.0)
		{
		//  \todo: Implement correctly
		//	dummy for rhs
			vector_type rhsDummy; rhsDummy.resize(u.size());

			adjust_linear(J, rhsDummy, dd, time);
		}

		template <typename TDD>
		void adjust_defect(vector_type& d, const vector_type& u,
		                   ConstSmartPtr<TDD> dd, number time = 0.0)
		{
			UG_THROW_FATAL("not implemented.");
		}

		template <typename TDD>
		void adjust_rhs(vector_type& rhs, const vector_type& u,
		                ConstSmartPtr<TDD> dd, number time = 0.0)
		{
			UG_THROW_FATAL("not implemented.");
		}

		template <typename TDD>
		void adjust_solution(vector_type& u, ConstSmartPtr<TDD> dd,
		                     number time = 0.0)
		{
			UG_THROW_FATAL("not implemented.");
		}

		template <typename TDD>
		void adjust_linear(matrix_type& mat, vector_type& rhs,
		                   ConstSmartPtr<TDD> dd, number time);

	protected:
		void SplitAddRow(matrix_type& A,
		                 std::vector<size_t> & constrainedIndex,
		                 std::vector<std::vector<size_t> >& vConstrainingIndices);

		void SetInterpolation(matrix_type& A,
		                      std::vector<size_t> & constrainedIndex,
		                      std::vector<std::vector<size_t> >& vConstrainingIndices);

		void HandleRhs(vector_type& rhs,
		               std::vector<size_t> & constrainedIndex,
		               std::vector<std::vector<size_t> >& vConstrainingIndices);
};

}; // namespace ug

#include "p1_continuity_constraints_impl.h"

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__CONSTRAINTS__CONTINUITY_CONSTRAINTS__P1_CONTINUITY_CONSTRAINTS__ */
