/*
 * amg_solver.h
 *
 *  Created on: 16.06.2010
 *      Author: mrupp
 */

#ifndef __H__LIB_DISCRETIZATION__LAPACK_LU__LAPACK_LU__
#define __H__LIB_DISCRETIZATION__LAPACK_LU__LAPACK_LU__

namespace ug{

#define AMG_MAX_LEVELS 32

template <typename TDiscreteFunction>
class LapackLUOperator : public ILinearizedOperatorInverse<TDiscreteFunction, TDiscreteFunction>
{

public:
	// domain function type
	typedef TDiscreteFunction domain_function_type;
	// codomain function type
	typedef TDiscreteFunction codomain_function_type;

	typedef typename TDiscreteFunction::algebra_type algebra_type;
	typedef typename algebra_type::matrix_type Matrix_type;


public:
	LapackLUOperator() : m_lapacklu()
	{};

	virtual bool init(ILinearizedOperator<TDiscreteFunction, TDiscreteFunction>& Op)
	{
		AssembledLinearizedOperator<TDiscreteFunction>* A
			= dynamic_cast<AssembledLinearizedOperator<TDiscreteFunction>*>(&Op);
		UG_ASSERT(A != NULL, 	"Operator used does not use a matrix."
								" Currently only matrix based operators can be inverted.\n");

		// remember operator
		m_pOperator = A;
		m_pMatrix = &m_pOperator->get_matrix();

		m_lapacklu.init(*m_pMatrix);

		return true;
	}

		// prepare Operator
	virtual bool prepare(domain_function_type& u, domain_function_type& d, codomain_function_type& c)
	{
		typename domain_function_type::vector_type& d_vec = d.get_vector();
		typename codomain_function_type::vector_type& c_vec = c.get_vector();

		UG_ASSERT(d_vec.size() == m_pMatrix->num_rows(),	"Vector and Row sizes have to match!");
		UG_ASSERT(c_vec.size() == m_pMatrix->num_cols(), "Vector and Column sizes have to match!");

		m_lapacklu.prepare(d_vec, c_vec);
		return true;
	}

	// compute new correction c = B*d
	//    AND
	// update defect: d := d - A*c
	virtual bool apply(domain_function_type& d, codomain_function_type& c)
	{
		typename domain_function_type::vector_type& d_vec = d.get_vector();
		typename codomain_function_type::vector_type& c_vec = c.get_vector();

		UG_ASSERT(d_vec.size() == m_pMatrix->num_rows(),	"Vector and Row sizes have to match!");
		UG_ASSERT(c_vec.size() == m_pMatrix->num_cols(), "Vector and Column sizes have to match!");
		UG_ASSERT(d_vec.size() == c_vec.size(), "Vector sizes have to match!");

		m_lapacklu.apply(d_vec, c_vec);

		return true;
	}



	// destructor
	virtual ~LapackLUOperator() {};

protected:
	AssembledLinearizedOperator<TDiscreteFunction>* m_pOperator;

	LapackLU m_lapacklu;
	SparseMatrix<double>* m_pMatrix;
	bool m_bOpChanged;
};

}

#endif /* __H__LIB_DISCRETIZATION__LAPACK_LU__LAPACK_LU__ */
