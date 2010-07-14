/*
 * jacobi.h
 *
 *  Created on: 04.07.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__JACOBI__
#define __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__JACOBI__

#include "lib_algebra/lib_algebra.h"
#include "lib_discretization/operator/operator_interface.h"

namespace ug{

/* This Operator type behaves different on application. It not only computes v = L*u, but also changes u. */
/* It is used in iterative schemes. */
template <typename TDiscreteFunction>
class AssembledJacobiOperator : public ILinearizedIteratorOperator<TDiscreteFunction, TDiscreteFunction>
{
	public:
		// export types:

		// domain function type
		typedef TDiscreteFunction domain_function_type;

		// codomain function type
		typedef TDiscreteFunction codomain_function_type;

	private:
		typedef typename TDiscreteFunction::algebra_type algebra_type;

	public:
		AssembledJacobiOperator(number damp) : m_damp(damp), m_bOpChanged(true)
		{};

		bool init(ILinearizedOperator<domain_function_type,codomain_function_type>& Op)
		{
			AssembledLinearizedOperator<TDiscreteFunction>* A
				= dynamic_cast<AssembledLinearizedOperator<TDiscreteFunction>*>(&Op);
			UG_ASSERT(A != NULL, 	"Operator used does not use a matrix."
									" Currently only matrix based operators can be inverted by this Jacobi.\n");

			// remember operator
			m_pOperator = A;
			m_pMatrix = &m_pOperator->get_matrix();

			m_bOpChanged = true;

			return true;
		}

		// prepare Operator
		virtual bool prepare(domain_function_type& u, domain_function_type& d, codomain_function_type& c)
		{
#ifdef UG_PARALLEL
			if(m_bOpChanged)
			{
				// create help vector to apply diagonal
				size_t size = m_pMatrix->num_rows();
				if(size != m_pMatrix->num_cols())
				{
					UG_LOG("Square Matrix needed for Jacobi Iteration.\n");
					return false;
				}

				if(m_diagInv.size() == size) m_diagInv.set(0.0);
				else {
					m_diagInv.destroy();
					m_diagInv.create(size);
				}

				m_diagInv.set_slave_layout(d.get_slave_layout());
				m_diagInv.set_master_layout(d.get_master_layout());

				typename algebra_type::matrix_type::local_matrix_type locMat(1, 1);
				typename algebra_type::matrix_type::local_index_type locInd(1);
				typename domain_function_type::vector_type::local_vector_type locVec(1);

				// collect diagonal
				for(size_t i = 0; i < m_diagInv.size(); ++i){
					locInd[0][0] = i;
					m_pMatrix->get(locMat, locInd, locInd);

					locVec[0] = locMat(0, 0);
					m_diagInv.set(locVec, locInd);
				}

				//	make diagonal consistent
				m_diagInv.set_storage_type(PST_ADDITIVE);
				m_diagInv.change_storage_type(PST_CONSISTENT);

				// invert diagonal and multiply by damping
				for(size_t i = 0; i < m_diagInv.size(); ++i){
					locInd[0][0] = i;
					m_diagInv.get(locVec, locInd);

					locVec[0] = m_damp / locVec[0];
					m_diagInv.set(locVec, locInd);
				}
				m_bOpChanged = false;
			}
#endif

			// TODO: Do we assume, that m_Op has been prepared?
			//m_Op->prepare(c,d);
			return true;
		}

		// compute new correction c = B*d
		//    AND
		// update defect: d := d - A*c
		virtual bool apply(domain_function_type& d, codomain_function_type& c)
		{
			if(!d.has_storage_type(PST_ADDITIVE)) return false;

			typename domain_function_type::vector_type& d_vec = d.get_vector();
			typename codomain_function_type::vector_type& c_vec = c.get_vector();

			UG_ASSERT(d_vec.size() == m_pMatrix->num_rows(),	"Vector and Row sizes have to match!");
			UG_ASSERT(c_vec.size() == m_pMatrix->num_cols(), "Vector and Column sizes have to match!");
			UG_ASSERT(d_vec.size() == c_vec.size(), "Vector sizes have to match!");

#ifdef UG_PARALLEL
			typename algebra_type::matrix_type::local_matrix_type locMat(1, 1);
			typename algebra_type::matrix_type::local_index_type locInd(1);
			typename domain_function_type::vector_type::local_vector_type locVec(1);

			// multiply defect with diagonal, c = damp * D^{-1} * d
			for(size_t i = 0; i < m_diagInv.size(); ++i){
				locInd[0][0] = i;
				m_diagInv.get(locVec, locInd);
				number a = locVec[0];

				d_vec.get(locVec, locInd);
				locVec[0] *= a;

				c_vec.set(locVec, locInd);
			}

			// make correction consistent
			c.set_storage_type(PST_ADDITIVE);
			if(c.change_storage_type(PST_CONSISTENT) != true)
				return false;
#else
			// apply iterator: c = B*d (damp is not used)
			diag_step(*m_pMatrix, c_vec, d_vec, m_damp);

			// damp correction
			c *= m_damp;
#endif
			// update defect
			// TODO: Check that matrix has correct type (additive)
			m_pMatrix->matmul_minus(d_vec, c_vec);

			// defect is now no longer unique (maybe we should handle this in matrix multiplication)
			d.set_storage_type(PST_ADDITIVE);
			return true;
		}

		// clone
		virtual ILinearizedIteratorOperator<TDiscreteFunction,TDiscreteFunction>* clone()
		{
			AssembledJacobiOperator<TDiscreteFunction>* clone = new AssembledJacobiOperator<TDiscreteFunction>(m_damp);

			return dynamic_cast<ILinearizedIteratorOperator<TDiscreteFunction,TDiscreteFunction>* >(clone);
		}


		// destructor
		virtual ~AssembledJacobiOperator() {};

	protected:
		AssembledLinearizedOperator<TDiscreteFunction>* m_pOperator;

		typename algebra_type::matrix_type* m_pMatrix;

		typename domain_function_type::vector_type m_diagInv;

		number m_damp;

		bool m_bOpChanged;
};


} // end namespace ug

#endif
