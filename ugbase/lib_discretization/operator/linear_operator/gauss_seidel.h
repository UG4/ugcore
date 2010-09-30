/*
 * gauss_seidel.h
 *
 *  Created on: 14.07.2010
 *      Author: Martin Rupp
 */

#ifndef __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__GAUSS_SEIDEL__
#define __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__GAUSS_SEIDEL__

#include "lib_algebra/lib_algebra.h"
#include "lib_discretization/operator/operator_interface.h"

namespace ug{

/* This Operator type behaves different on application. It not only computes v = L*u, but also changes u. */
/* It is used in iterative schemes. */
template <typename TDiscreteFunction>
class AssembledGSOperator : public ILinearizedIteratorOperator<TDiscreteFunction, TDiscreteFunction>
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
		AssembledGSOperator() : m_bOpChanged(true)
		{};

		bool init(ILinearizedOperator<TDiscreteFunction,TDiscreteFunction>& Op)
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
		virtual bool prepare(TDiscreteFunction& cOut, TDiscreteFunction& uIn, TDiscreteFunction& dIn)
		{

			// TODO: Do we assume, that m_Op has been prepared?
			//m_Op->prepare(c,d);
			return true;
		}

		// compute new correction c = B*d
		//    AND
		// update defect: d := d - A*c
		virtual bool apply(TDiscreteFunction& cOut, TDiscreteFunction& dInOut, bool updateDefect)
		{
#ifdef UG_PARALLEL
			if(!dInOut.has_storage_type(PST_ADDITIVE)) return false;
#endif

			typename domain_function_type::vector_type& d_vec = dInOut.get_vector();
			typename codomain_function_type::vector_type& c_vec = cOut.get_vector();

			UG_ASSERT(d_vec.size() == m_pMatrix->num_rows(),	"Vector and Row sizes have to match!");
			UG_ASSERT(c_vec.size() == m_pMatrix->num_cols(), "Vector and Column sizes have to match!");
			UG_ASSERT(d_vec.size() == c_vec.size(), "Vector sizes have to match!");


			// apply iterator: c = B*d (damp is not used)
			gs_step_LL(*m_pMatrix, c_vec, d_vec);

			// update defect
			// TODO: Check that matrix has correct type (additive)
			if(updateDefect) m_pMatrix->matmul_minus(d_vec, c_vec);

			// defect is now no longer unique (maybe we should handle this in matrix multiplication)
#ifdef UG_PARALLEL
			dInOut.set_storage_type(PST_ADDITIVE);
#endif
			return true;
		}

		// clone
		virtual ILinearizedIteratorOperator<TDiscreteFunction,TDiscreteFunction>* clone()
		{
			AssembledGSOperator<TDiscreteFunction>* clone = new AssembledGSOperator<TDiscreteFunction>();
			return dynamic_cast<ILinearizedIteratorOperator<TDiscreteFunction,TDiscreteFunction>* >(clone);
		}

		// destructor
		virtual ~AssembledGSOperator() {};

	protected:
		AssembledLinearizedOperator<TDiscreteFunction>* m_pOperator;

		typename algebra_type::matrix_type* m_pMatrix;

		bool m_bOpChanged;
};


template <typename TDiscreteFunction>
class AssembledBGSOperator : public ILinearizedIteratorOperator<TDiscreteFunction, TDiscreteFunction>
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
		AssembledBGSOperator() : m_bOpChanged(true)
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

			// TODO: Do we assume, that m_Op has been prepared?
			//m_Op->prepare(c,d);
			return true;
		}

		// compute new correction c = B*d
		//    AND
		// update defect: d := d - A*c
		virtual bool apply(domain_function_type& d, codomain_function_type& c, bool updateDefect)
		{
#ifdef UG_PARALLEL
			if(!d.has_storage_type(PST_ADDITIVE)) return false;
#endif

			typename domain_function_type::vector_type& d_vec = d.get_vector();
			typename codomain_function_type::vector_type& c_vec = c.get_vector();

			UG_ASSERT(d_vec.size() == m_pMatrix->num_rows(),	"Vector and Row sizes have to match!");
			UG_ASSERT(c_vec.size() == m_pMatrix->num_cols(), "Vector and Column sizes have to match!");
			UG_ASSERT(d_vec.size() == c_vec.size(), "Vector sizes have to match!");


			// apply iterator: c = B*d (damp is not used)
			gs_step_UR(*m_pMatrix, c_vec, d_vec);

			// update defect
			// TODO: Check that matrix has correct type (additive)
			if(updateDefect) m_pMatrix->matmul_minus(d_vec, c_vec);

			// defect is now no longer unique (maybe we should handle this in matrix multiplication)
#ifdef UG_PARALLEL
			d.set_storage_type(PST_ADDITIVE);
#endif
			return true;
		}

		// clone
		virtual ILinearizedIteratorOperator<TDiscreteFunction,TDiscreteFunction>* clone()
		{
			AssembledBGSOperator<TDiscreteFunction>* clone = new AssembledBGSOperator<TDiscreteFunction>();
			return dynamic_cast<ILinearizedIteratorOperator<TDiscreteFunction,TDiscreteFunction>* >(clone);
		}

		// destructor
		virtual ~AssembledBGSOperator() {};

	protected:
		AssembledLinearizedOperator<TDiscreteFunction>* m_pOperator;

		typename algebra_type::matrix_type* m_pMatrix;

		bool m_bOpChanged;
};



template <typename TDiscreteFunction>
class AssembledSGSOperator : public ILinearizedIteratorOperator<TDiscreteFunction, TDiscreteFunction>
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
		AssembledSGSOperator() : m_bOpChanged(true)
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

			// TODO: Do we assume, that m_Op has been prepared?
			//m_Op->prepare(c,d);
			return true;
		}

		// compute new correction c = B*d
		//    AND
		// update defect: d := d - A*c
		virtual bool apply(domain_function_type& d, codomain_function_type& c, bool updateDefect)
		{
#ifdef UG_PARALLEL
			if(!d.has_storage_type(PST_ADDITIVE)) return false;
#endif
			typename domain_function_type::vector_type& d_vec = d.get_vector();
			typename codomain_function_type::vector_type& c_vec = c.get_vector();

			UG_ASSERT(d_vec.size() == m_pMatrix->num_rows(),	"Vector and Row sizes have to match!");
			UG_ASSERT(c_vec.size() == m_pMatrix->num_cols(), "Vector and Column sizes have to match!");
			UG_ASSERT(d_vec.size() == c_vec.size(), "Vector sizes have to match!");


			// apply iterator: c = B*d (damp is not used)
			sgs_step(*m_pMatrix, c_vec, d_vec);

			// update defect
			// TODO: Check that matrix has correct type (additive)
			if(updateDefect) m_pMatrix->matmul_minus(d_vec, c_vec);

			// defect is now no longer unique (maybe we should handle this in matrix multiplication)
#ifdef UG_PARALLEL
			d.set_storage_type(PST_ADDITIVE);
#endif
			return true;
		}

		// clone
		virtual ILinearizedIteratorOperator<TDiscreteFunction,TDiscreteFunction>* clone()
		{
			AssembledSGSOperator<TDiscreteFunction>* clone = new AssembledSGSOperator<TDiscreteFunction>();
			return dynamic_cast<ILinearizedIteratorOperator<TDiscreteFunction,TDiscreteFunction>* >(clone);
		}

		// destructor
		virtual ~AssembledSGSOperator() {};

	protected:
		AssembledLinearizedOperator<TDiscreteFunction>* m_pOperator;

		typename algebra_type::matrix_type* m_pMatrix;

		bool m_bOpChanged;
};

} // end namespace ug

#endif // __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__GAUSS_SEIDEL__
