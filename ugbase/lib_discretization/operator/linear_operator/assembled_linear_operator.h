


#ifndef __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__ASSEMBLED_LINEAR_OPERATOR__
#define __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__ASSEMBLED_LINEAR_OPERATOR__

#include "lib_algebra/lib_algebra.h"

namespace ug{

template <typename TDiscreteFunction>
class AssembledLinearizedOperator : public ILinearizedOperator<TDiscreteFunction, TDiscreteFunction>
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
		AssembledLinearizedOperator(IAssemble<domain_function_type, algebra_type>& ass) :
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

			if(m_J.row_size() == d_vec.size() && m_J.col_size() == c_vec.size())
			{
				if(m_J.set(0.0) != true) return false;
			}
			else
			{
				if(m_J.destroy() != true) return false;
				if(m_J.create(d_vec.size(), c_vec.size()) != true) return false;
			}

			if(m_ass.assemble_jacobian(m_J, u) != IAssemble_OK) return false;

			return true;
		}

		// compute d = J(u)*c (here, J(u) is a Matrix)
		virtual bool apply(domain_function_type& c, codomain_function_type& d)
		{
			const typename domain_function_type::vector_type& x = c.get_vector();
			typename codomain_function_type::vector_type& b = d.get_vector();

			UG_ASSERT(x.size() == m_J.row_size(), "Row size '" << m_J.row_size() << "' of Matrix J and size '" << x.size() << "' of Vector x do not match. Cannot calculate L*x.");
			UG_ASSERT(b.size() == m_J.col_size(), "Column size '" << m_J.row_size() << "' of Matrix J and size  '" << b.size() << "' of Vector b do not match. Cannot calculate b := L*x.");

			if(!c.has_storage_type(PST_CONSISTENT)) return false;
			d.set_storage_type(PST_ADDITIVE);
			return m_J.apply(b,x);
		}

		// d := d - J(u)*c
		virtual bool apply_sub(domain_function_type& c, codomain_function_type& d)
		{
			const typename domain_function_type::vector_type& x = c.get_vector();
			typename codomain_function_type::vector_type& b = d.get_vector();

			UG_ASSERT(x.size() == m_J.row_size(), "Row size '" << m_J.row_size() << "' of Matrix J and size '" << x.size() << "' of Vector x do not match. Cannot calculate L*x.");
			UG_ASSERT(b.size() == m_J.col_size(), "Column size '" << m_J.row_size() << "' of Matrix J and size  '" << b.size() << "' of Vector b do not match. Cannot calculate b := b - L*x.");

			// TODO: Check that matrix has correct type (additive)
			if(!c.has_storage_type(PST_CONSISTENT)) return false;
			if(!d.has_storage_type(PST_ADDITIVE)) return false;
			d.set_storage_type(PST_ADDITIVE);
			return m_J.matmul_minus(b,x);
		}

		virtual typename algebra_type::matrix_type& get_matrix()
		{
			return m_J;
		}

		// destructor
		virtual ~AssembledLinearizedOperator()
		{
			m_J.destroy();
		};

	protected:
		// matrix type used
		typedef typename algebra_type::matrix_type matrix_type;

	protected:
		// assembling procedure
		IAssemble<domain_function_type, algebra_type>& m_ass;

		// matrix storage
		matrix_type m_J;
};


template <typename TDiscreteFunction>
class AssembledLinearOperator : public ILinearOperator<TDiscreteFunction, TDiscreteFunction>, public AssembledLinearizedOperator<TDiscreteFunction>
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
		AssembledLinearOperator(IAssemble<domain_function_type, algebra_type>& ass, bool assemble_rhs = false) :
			AssembledLinearizedOperator<TDiscreteFunction>(ass),
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

			if(m_Matrix.row_size() == b.size() && m_Matrix.col_size() == x.size())
			{
				if(m_Matrix.set(0.0) != true) return false;
			}
			else
			{
				if(m_Matrix.destroy() != true) return false;
				if(m_Matrix.create(b.size(), x.size()) != true) return false;
			}

			if(m_assemble_rhs)
			{
				b.set(0.0);
				if(m_ass.assemble_linear(m_Matrix, b, u) != IAssemble_OK) return false;
				f.set_storage_type(PST_ADDITIVE);
			}
			else
			{
				if(m_ass.assemble_jacobian(m_Matrix, u) != IAssemble_OK) return false;
			}

			if(m_ass.assemble_solution(u) != IAssemble_OK) return false;

			return true;
		}

		// compute f = L*u (here, L is a Matrix)
		virtual bool apply(domain_function_type& u, codomain_function_type& f)
		{
			const typename domain_function_type::vector_type& x = u.get_vector();
			typename codomain_function_type::vector_type& b = f.get_vector();

			UG_ASSERT(x.size() == m_Matrix.row_size(), "Row size '" << m_Matrix.row_size() << "' of Matrix L and size '" << x.size() << "' of Vector x do not match. Cannot calculate L*x.");
			UG_ASSERT(b.size() == m_Matrix.col_size(), "Column size '" << m_Matrix.row_size() << "' of Matrix L and size  '" << b.size() << "' of Vector b do not match. Cannot calculate b := L*x.");

			// TODO: Check that matrix has correct type (additive)
			if(!u.has_storage_type(PST_CONSISTENT)) return false;
			f.set_storage_type(PST_ADDITIVE);
			return m_Matrix.apply(b,x);
		}

		// f := f - L*u
		virtual bool apply_sub(domain_function_type& u, codomain_function_type& f)
		{
			const typename domain_function_type::vector_type& x = u.get_vector();
			typename codomain_function_type::vector_type& b = f.get_vector();

			UG_ASSERT(x.size() == m_Matrix.row_size(), "Row size '" << m_Matrix.row_size() << "' of Matrix L and size '" << x.size() << "' of Vector x do not match. Cannot calculate L*x.");
			UG_ASSERT(b.size() == m_Matrix.col_size(), "Column size '" << m_Matrix.row_size() << "' of Matrix L and size  '" << b.size() << "' of Vector b do not match. Cannot calculate b := b - L*x.");

			// TODO: Check that matrix has correct type (additive)
			if(!u.has_storage_type(PST_CONSISTENT)) return false;
			if(!f.has_storage_type(PST_ADDITIVE)) return false;
			f.set_storage_type(PST_ADDITIVE);
			return m_Matrix.matmul_minus(b,x);
		}

		virtual typename algebra_type::matrix_type& get_matrix()
		{
			return m_Matrix;
		}

		// destructor
		virtual ~AssembledLinearOperator()
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
		IAssemble<domain_function_type, algebra_type>& m_ass;

		// matrix storage
		matrix_type m_Matrix;
};


} // namespace ug

#endif /* __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__ASSEMBLED_LINEAR_OPERATOR__ */
