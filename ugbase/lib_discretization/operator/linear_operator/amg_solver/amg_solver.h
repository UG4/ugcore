/**
 * \file amg_rs_prolongation.h
 *
 * \author Martin Rupp
 *
 * \date 16.06.2010
 *
 * Goethe-Center for Scientific Computing 2009-2010.
 */


#ifndef __H__LIB_DISCRETIZATION__AMG_SOLVER__AMG_SOLVER__
#define __H__LIB_DISCRETIZATION__AMG_SOLVER__AMG_SOLVER__

#include "amg.h"

namespace ug{


template <typename TDiscreteFunction>
class AMGSolver : 	public ILinearizedIteratorOperator<TDiscreteFunction, TDiscreteFunction>,
					public amg<SparseMatrix<double>, Vector<double> > //amg<TDiscreteFunction::algebra_type::Matrix_type, TDiscreteFunction::algebra_type::Vector_type>
{

public:
	// domain function type
	typedef TDiscreteFunction domain_function_type;
	// codomain function type
	typedef TDiscreteFunction codomain_function_type;

	typedef typename TDiscreteFunction::algebra_type algebra_type;
	typedef typename algebra_type::matrix_type Matrix_type;

private:
	typedef ILinearizedOperatorInverse<domain_function_type, domain_function_type> base_solver_type;
	typedef ILinearizedIteratorOperator<domain_function_type, domain_function_type> smoother_type;

public:
	AMGSolver() : m_amg()
	{
	};

	bool init(ILinearizedOperator<domain_function_type,codomain_function_type>& Op)
	{
		AssembledLinearizedOperator<TDiscreteFunction>* A
			= dynamic_cast<AssembledLinearizedOperator<TDiscreteFunction>*>(&Op);
		UG_ASSERT(A != NULL, 	"Operator used does not use a matrix."
								" Currently only matrix based operators can be inverted.\n");

		// remember operator
		m_pOperator = A;
		m_pMatrix = &m_pOperator->get_matrix();

		m_amg.init(*m_pMatrix);

		m_bOpChanged = true;

		return true;
	}

	void set_debug_positions(const MathVector<2> *positions2d, size_t size)
	{
		m_amg.set_debug_positions(positions2d, size);
	}

	void set_debug_positions(const MathVector<3> *positions3d, size_t size)
	{
		m_amg.set_debug_positions(positions3d, size);
	}


	// prepare Operator
	virtual bool prepare(domain_function_type& u, domain_function_type& d, codomain_function_type& c)
	{
		return true;
	}

	// compute new correction c = B*d
	//    AND
	// update defect: d := d - A*c
	virtual bool apply(domain_function_type& d, codomain_function_type& c, bool updateDefect)
	{
		typename domain_function_type::vector_type& d_vec = d.get_vector();
		typename codomain_function_type::vector_type& c_vec = c.get_vector();

		UG_ASSERT(d_vec.size() == m_pMatrix->num_rows(),	"Vector and Row sizes have to match!");
		UG_ASSERT(c_vec.size() == m_pMatrix->num_cols(), "Vector and Column sizes have to match!");
		UG_ASSERT(d_vec.size() == c_vec.size(), "Vector sizes have to match!");

		m_amg.get_correction_and_update_defect(d_vec, c_vec);

		return true;
	}

	// destructor
	virtual ~AMGSolver() {};

	virtual ILinearizedIteratorOperator<TDiscreteFunction,TDiscreteFunction>* clone()
	{
		return new AMGSolver;
	}

protected:
	AssembledLinearizedOperator<TDiscreteFunction>* m_pOperator;

	SparseMatrix<double>* m_pMatrix;
	bool m_bOpChanged;
	amg<SparseMatrix<double>, Vector<double> > m_amg;
	int m_nu1, m_nu2;
};


} // namespace ug


#endif /* __H__LIB_DISCRETIZATION__AMG_SOLVER__AMG_SOLVER__ */
