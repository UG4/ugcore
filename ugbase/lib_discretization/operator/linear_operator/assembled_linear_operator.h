
#ifndef __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__ASSEMBLED_LINEAR_OPERATOR__
#define __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__ASSEMBLED_LINEAR_OPERATOR__

#include "lib_algebra/lib_algebra.h"

namespace ug{

template <typename TDiscreteFunction>
class AssembledDiscreteLinearizedOperator : public IDiscreteLinearizedOperator<TDiscreteFunction, TDiscreteFunction>
{
	public:
		// export types:

		// domain function type
		typedef TDiscreteFunction domain_function_type;

		// codomain function type
		typedef TDiscreteFunction codomain_function_type;

		// type of algebra
		typedef typename TDiscreteFunction::algebra_type algebra_type;

	public:
		AssembledDiscreteLinearizedOperator(IAssemble<algebra_type, domain_function_type>& ass) :
			m_ass(ass)
		{};

		virtual bool init()
		{
			return true;
		}

		// prepare the operator for application (e.g. compute an intern Matrix J(u))
		virtual bool prepare(domain_function_type& u, domain_function_type& c, codomain_function_type& d)
		{
			const typename domain_function_type::vector_type& c_vec = c.get_vector();
			typename codomain_function_type::vector_type& d_vec = d.get_vector();

			UG_DLOG(LIB_ALG_LINEAR_OPERATOR, 2, "Creating Matrix of size (" << d_vec.size() <<", " << c_vec.size() <<").\n");
			if(m_J.row_size() == d_vec.size() && m_J.col_size() == c_vec.size())
			{
				UG_DLOG(LIB_ALG_LINEAR_OPERATOR, 3, " ---- 'prepare': Reset matrix to 0.0.\n");
				if(m_J.set(0.0) != true) return false;
			}
			else
			{
				UG_DLOG(LIB_ALG_LINEAR_OPERATOR, 3, " ---- 'prepare': Calling destroy matrix.\n");
				if(m_J.destroy() != true) return false;
				UG_DLOG(LIB_ALG_LINEAR_OPERATOR, 3, " ---- 'prepare': Calling allocate matrix.\n");
				if(m_J.create(d_vec.size(), c_vec.size()) != true) return false;
			}

			UG_DLOG(LIB_ALG_LINEAR_OPERATOR, 2, " ---- 'prepare': Calling 'm_ass.assemble_jacobian'.\n");
			if(m_ass.assemble_jacobian(m_J, u) != IAssemble_OK) return false;

			UG_DLOG(LIB_DISC_ASSEMBLE, 10, "Matrix:\n" << m_J << "\n");

			return true;
		}

		// compute d = J(u)*c (here, J(u) is a Matrix)
		virtual bool apply(domain_function_type& c, codomain_function_type& d)
		{
			const typename domain_function_type::vector_type& x = c.get_vector();
			typename codomain_function_type::vector_type& b = d.get_vector();

			UG_ASSERT(x.size() == m_J.row_size(), "Row size '" << m_J.row_size() << "' of Matrix J and size '" << x.size() << "' of Vector x do not match. Cannot calculate L*x.");
			UG_ASSERT(b.size() == m_J.col_size(), "Column size '" << m_J.row_size() << "' of Matrix J and size  '" << b.size() << "' of Vector b do not match. Cannot calculate b := L*x.");

			return m_J.apply(b,x);
		}

		// d := d - J(u)*c
		virtual bool apply_sub(domain_function_type& c, codomain_function_type& d)
		{
			const typename domain_function_type::vector_type& x = c.get_vector();
			typename codomain_function_type::vector_type& b = d.get_vector();

			UG_ASSERT(x.size() == m_J.row_size(), "Row size '" << m_J.row_size() << "' of Matrix J and size '" << x.size() << "' of Vector x do not match. Cannot calculate L*x.");
			UG_ASSERT(b.size() == m_J.col_size(), "Column size '" << m_J.row_size() << "' of Matrix J and size  '" << b.size() << "' of Vector b do not match. Cannot calculate b := b - L*x.");

			return m_J.matmul_minus(b,x);
		}

		virtual typename algebra_type::matrix_type& get_matrix()
		{
			return m_J;
		}

		// destructor
		virtual ~AssembledDiscreteLinearizedOperator()
		{
			m_J.destroy();
		};

	protected:
		// matrix type used
		typedef typename algebra_type::matrix_type matrix_type;

	protected:
		// assembling procedure
		IAssemble<algebra_type, domain_function_type>& m_ass;

		// matrix storage
		matrix_type m_J;
};


template <typename TDiscreteFunction>
class AssembledDiscreteLinearOperator : public AssembledDiscreteLinearizedOperator<TDiscreteFunction>
{
	public:
		// export types:

		// domain function type
		typedef TDiscreteFunction domain_function_type;

		// codomain function type
		typedef TDiscreteFunction codomain_function_type;

		// type of algebra
		typedef typename TDiscreteFunction::algebra_type algebra_type;

	public:
		AssembledDiscreteLinearOperator(IAssemble<algebra_type, domain_function_type>& ass, bool assemble_rhs = false) :
			AssembledDiscreteLinearizedOperator<TDiscreteFunction>(ass),
			m_assemble_rhs(assemble_rhs), m_ass(ass)
		{};

		virtual bool init()
		{
			return true;
		}

		// prepare the operator for application (e.g. compute an intern Matrix L)
		virtual bool prepare(domain_function_type& u, codomain_function_type& f)
		{
			const typename domain_function_type::vector_type& x = u.get_vector();
			typename codomain_function_type::vector_type& b = f.get_vector();

			UG_DLOG(LIB_ALG_LINEAR_OPERATOR, 2, "Creating Matrix of size (" << b.size() <<", " << x.size() <<").\n");
			if(m_Matrix.row_size() == b.size() && m_Matrix.col_size() == x.size())
			{
				UG_DLOG(LIB_ALG_LINEAR_OPERATOR, 3, " ---- 'prepare': Reset matrix to 0.0.\n");
				if(m_Matrix.set(0.0) != true) return false;
			}
			else
			{
				UG_DLOG(LIB_ALG_LINEAR_OPERATOR, 3, " ---- 'prepare': Calling destroy matrix.\n");
				if(m_Matrix.destroy() != true) return false;
				UG_DLOG(LIB_ALG_LINEAR_OPERATOR, 3, " ---- 'prepare': Calling allocate matrix.\n");
				if(m_Matrix.create(b.size(), x.size()) != true) return false;
			}

			if(m_assemble_rhs)
			{
				UG_DLOG(LIB_ALG_LINEAR_OPERATOR, 2, " ---- 'prepare': Calling 'm_ass.assemble_linear'.\n");
				if(m_ass.assemble_linear(m_Matrix, b, u) != IAssemble_OK) return false;
			}
			else
			{
				UG_DLOG(LIB_ALG_LINEAR_OPERATOR, 2, " ---- 'prepare': Calling 'm_ass.assemble_jacobian'.\n");
				if(m_ass.assemble_jacobian(m_Matrix, u) != IAssemble_OK) return false;
			}

			if(m_ass.assemble_solution(u) != IAssemble_OK) return false;

			UG_DLOG(LIB_DISC_ASSEMBLE, 10, "Matrix:\n" << m_Matrix << "\n");

			return true;
		}

		// compute f = L*u (here, L is a Matrix)
		virtual bool apply(domain_function_type& u, codomain_function_type& f)
		{
			const typename domain_function_type::vector_type& x = u.get_vector();
			typename codomain_function_type::vector_type& b = f.get_vector();

			UG_ASSERT(x.size() == m_Matrix.row_size(), "Row size '" << m_Matrix.row_size() << "' of Matrix L and size '" << x.size() << "' of Vector x do not match. Cannot calculate L*x.");
			UG_ASSERT(b.size() == m_Matrix.col_size(), "Column size '" << m_Matrix.row_size() << "' of Matrix L and size  '" << b.size() << "' of Vector b do not match. Cannot calculate b := L*x.");

			return m_Matrix.apply(b,x);
		}

		// f := f - L*u
		virtual bool apply_sub(domain_function_type& u, codomain_function_type& f)
		{
			const typename domain_function_type::vector_type& x = u.get_vector();
			typename codomain_function_type::vector_type& b = f.get_vector();

			UG_ASSERT(x.size() == m_Matrix.row_size(), "Row size '" << m_Matrix.row_size() << "' of Matrix L and size '" << x.size() << "' of Vector x do not match. Cannot calculate L*x.");
			UG_ASSERT(b.size() == m_Matrix.col_size(), "Column size '" << m_Matrix.row_size() << "' of Matrix L and size  '" << b.size() << "' of Vector b do not match. Cannot calculate b := b - L*x.");

			return m_Matrix.matmul_minus(b,x);
		}

		virtual typename algebra_type::matrix_type& get_matrix()
		{
			return m_Matrix;
		}

		// destructor
		virtual ~AssembledDiscreteLinearOperator()
		{
			m_Matrix.destroy();
		};

	protected:
		// matrix type used
		typedef typename algebra_type::matrix_type matrix_type;

	protected:
		// choose, whether rhs should be assembled as well.
		bool m_assemble_rhs;

		// assembling procedure
		IAssemble<algebra_type, domain_function_type>& m_ass;

		// matrix storage
		matrix_type m_Matrix;
};


/* This Operator type behaves different on application. It not only computes v = L*u, but also changes u. */
/* It is used in iterative schemes. */
template <typename TDiscreteFunction>
class AssembledJacobiOperator : public IDiscreteLinearizedIteratorOperator<TDiscreteFunction, TDiscreteFunction>
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
		AssembledJacobiOperator(number damp) : m_damp(damp)
		{};

		bool init(IDiscreteLinearizedOperator<domain_function_type,codomain_function_type>& Op)
		{
			AssembledDiscreteLinearizedOperator<TDiscreteFunction>* A = dynamic_cast<AssembledDiscreteLinearizedOperator<TDiscreteFunction>*>(&Op);
			UG_ASSERT(A != NULL, "Operator used does not use a matrix. Currently only matrix based operators can be inverted by this Jacobi.\n");

			m_Op = A;
			m_matrix = &m_Op->get_matrix();
			return true;
		}

		// prepare Operator
		virtual bool prepare(domain_function_type& u, domain_function_type& d, codomain_function_type& c)
		{
			UG_ASSERT(d.get_vector().size() == m_matrix->row_size(), "Row size '" << m_matrix->row_size() << "' of Matrix B and size '" << d.get_vector().size() << "' of Vector d do not match. Cannot calculate B*d.");
			UG_ASSERT(c.get_vector().size() == m_matrix->col_size(), "Column size '" << m_matrix->row_size() << "' of Matrix B and size  '" << c.get_vector().size() << "' of Vector c do not match. Cannot calculate c := B*d.");

			// TODO: Do we assume, that m_Op has been prepared?
			//m_Op->prepare(c,d);
			return true;
		}

		// compute new correction c = B*d
		//    AND
		// update defect: d := d - A*c
		virtual bool apply(domain_function_type& d, codomain_function_type& c)
		{
			typename domain_function_type::vector_type& d_vec = d.get_vector();
			typename codomain_function_type::vector_type& c_vec = c.get_vector();

			UG_ASSERT(d_vec.size() == m_matrix->row_size(), "Row size '" << m_matrix->row_size() << "' of Matrix B and size '" << d_vec.size() << "' of Vector d do not match. Cannot calculate B*d.");
			UG_ASSERT(c_vec.size() == m_matrix->col_size(), "Column size '" << m_matrix->row_size() << "' of Matrix B and size  '" << c_vec.size() << "' of Vector c do not match. Cannot calculate c := B*d.");

			return diag_step(*m_matrix, c_vec, d_vec, m_damp);
		}

		// destructor
		virtual ~AssembledJacobiOperator() {};

	protected:
		AssembledDiscreteLinearizedOperator<TDiscreteFunction>* m_Op;

		typename algebra_type::matrix_type* m_matrix;

		number m_damp;
};



} // namespace ug

#endif /* __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__ASSEMBLED_LINEAR_OPERATOR__ */
