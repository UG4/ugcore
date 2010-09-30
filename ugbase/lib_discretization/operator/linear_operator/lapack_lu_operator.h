/*
 * amg_solver.h
 *
 *  Created on: 16.06.2010
 *      Author: mrupp
 */

#ifndef __H__LIB_DISCRETIZATION__LAPACK_LU__LAPACK_LU__
#define __H__LIB_DISCRETIZATION__LAPACK_LU__LAPACK_LU__

namespace ug{

template <typename TFunction>
class LapackLUOperator : public ILinearizedOperatorInverse<TFunction, TFunction>
{

public:
	// domain function type
	typedef TFunction domain_function_type;
	// codomain function type
	typedef TFunction codomain_function_type;

	typedef typename TFunction::algebra_type algebra_type;
	typedef typename algebra_type::matrix_type Matrix_type;


public:
	LapackLUOperator() : m_lapacklu()
	{};

	virtual bool init(ILinearizedOperator<TFunction, TFunction>& Op)
	{
		AssembledLinearizedOperator<TFunction>* A
			= dynamic_cast<AssembledLinearizedOperator<TFunction>*>(&Op);
		UG_ASSERT(A != NULL, 	"Operator used does not use a matrix."
								" Currently only matrix based operators can be inverted.\n");

		// remember operator
		m_pOperator = A;
		m_pMatrix = &m_pOperator->get_matrix();

		m_lapacklu.init(*m_pMatrix);

		return true;
	}

		// prepare Operator
	virtual bool prepare(TFunction& dOut, TFunction& uIn, TFunction& cIn)
	{
		typename domain_function_type::vector_type& d_vec = dOut.get_vector();
		typename codomain_function_type::vector_type& c_vec = cIn.get_vector();

 		UG_ASSERT(d_vec.size() == m_pMatrix->num_rows(),	"Vector and Row sizes have to match!");
		UG_ASSERT(c_vec.size() == m_pMatrix->num_cols(), "Vector and Column sizes have to match!");

		// TODO: This must be inverted
		m_lapacklu.prepare(d_vec, c_vec);
		return true;
	}

	// compute new correction c = B*d
	//    AND
	// update defect: d := d - A*c
	virtual bool apply(TFunction& cOut, TFunction& dIn)
	{
		typename domain_function_type::vector_type& d_vec = dIn.get_vector();
		typename codomain_function_type::vector_type& c_vec = cOut.get_vector();

		UG_ASSERT(d_vec.size() == m_pMatrix->num_rows(),	"Vector and Row sizes have to match!");
		UG_ASSERT(c_vec.size() == m_pMatrix->num_cols(), "Vector and Column sizes have to match!");
		UG_ASSERT(d_vec.size() == c_vec.size(), "Vector sizes have to match!");

		// TODO: This must be inverted
		m_lapacklu.apply(d_vec, c_vec);
		
		return true;
	}



	// destructor
	virtual ~LapackLUOperator() {};

protected:
	AssembledLinearizedOperator<TFunction>* m_pOperator;

	LapackLU m_lapacklu;
	typename algebra_type::matrix_type *m_pMatrix;
	bool m_bOpChanged;
};

}

#endif /* __H__LIB_DISCRETIZATION__LAPACK_LU__LAPACK_LU__ */
