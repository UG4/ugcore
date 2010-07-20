


#ifndef __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__ASSEMBLED_LINEAR_OPERATOR__
#define __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__ASSEMBLED_LINEAR_OPERATOR__

#include "lib_algebra/lib_algebra.h"

namespace ug{

template <typename TFunction>
class AssembledLinearizedOperator : public ILinearizedOperator<TFunction, TFunction>
{
	public:
		// domain function type
		typedef TFunction domain_function_type;

		// codomain function type
		typedef TFunction codomain_function_type;

		// type of algebra
		typedef typename TFunction::algebra_type algebra_type;

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

			if(m_J.num_rows() == d_vec.size() && m_J.num_cols() == c_vec.size())
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
			#ifdef UG_PARALLEL
			if(!c.has_storage_type(PST_CONSISTENT))
				{
					UG_LOG("WARNING: In 'AssembledLinearizedOperator::apply':Inadequate storage format of Vector c."
										" Use consistent to avoid internal type conversion.\n");
					if(!c.change_storage_type(PST_CONSISTENT)) return false;
				}
			#endif

			const typename domain_function_type::vector_type& x = c.get_vector();
			typename codomain_function_type::vector_type& b = d.get_vector();

			UG_ASSERT(x.size() == m_J.num_rows(), "Row size '" << m_J.num_rows() << "' of Matrix J and size '" << x.size() << "' of Vector x do not match. Cannot calculate L*x.");
			UG_ASSERT(b.size() == m_J.num_cols(), "Column size '" << m_J.num_rows() << "' of Matrix J and size  '" << b.size() << "' of Vector b do not match. Cannot calculate b := L*x.");

			// set storage type to additiv, since it could been additive unique before
			// TODO: Handle this in matrix multiplication
			d.set_storage_type(PST_ADDITIVE);
			return m_J.apply(b,x);
		}

		// d := d - J(u)*c
		virtual bool apply_sub(domain_function_type& c, codomain_function_type& d)
		{
			#ifdef UG_PARALLEL
			if(!d.has_storage_type(PST_ADDITIVE) || !c.has_storage_type(PST_CONSISTENT))
				{
					UG_LOG("WARNING: In 'AssembledLinearizedOperator::apply_sub':Inadequate storage format of Vectors.\n");
					UG_LOG("                                             	 use: d additive and c consistent to avoid internal type conversion.\n");
					if(!d.change_storage_type(PST_ADDITIVE)) return false;
					if(!c.change_storage_type(PST_CONSISTENT)) return false;
				}
			#endif
			const typename domain_function_type::vector_type& x = c.get_vector();
			typename codomain_function_type::vector_type& b = d.get_vector();

			UG_ASSERT(x.size() == m_J.num_rows(), "Row size '" << m_J.num_rows() << "' of Matrix J and size '" << x.size() << "' of Vector x do not match. Cannot calculate L*x.");
			UG_ASSERT(b.size() == m_J.num_cols(), "Column size '" << m_J.num_rows() << "' of Matrix J and size  '" << b.size() << "' of Vector b do not match. Cannot calculate b := b - L*x.");

			// set storage type to additiv, since it could been additive unique before
			// TODO: Handle this in matrix multiplication
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


template <typename TFunction>
class AssembledLinearOperator : public ILinearOperator<TFunction, TFunction>, public AssembledLinearizedOperator<TFunction>
{
	public:
		// domain function type
		typedef TFunction domain_function_type;

		// codomain function type
		typedef TFunction codomain_function_type;

		// type of algebra
		typedef typename TFunction::algebra_type algebra_type;

	public:
		AssembledLinearOperator(IAssemble<domain_function_type, algebra_type>& ass, bool assemble_rhs = false) :
			AssembledLinearizedOperator<TFunction>(ass),
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

			if(m_Matrix.num_rows() == b.size() && m_Matrix.num_cols() == x.size())
			{
				if(m_Matrix.set(0.0) != true) return false;
			}
			else
			{
				if(m_Matrix.destroy() != true)
					{UG_LOG("Cannot destroy matrix.\n"); return false;}
				if(m_Matrix.create(b.size(), x.size()) != true)
					{UG_LOG("Cannot create matrix of size " << b.size() << "x" << x.size() << ".\n"); return false;}
			}

			if(m_assemble_rhs)
			{
				b.set(0.0);
				if(m_ass.assemble_linear(m_Matrix, b, u) != IAssemble_OK)
					{UG_LOG("Error while assembling Matrix and rhs.\n"); return false;}
				f.set_storage_type(PST_ADDITIVE);
			}
			else
			{
				if(m_ass.assemble_jacobian(m_Matrix, u) != IAssemble_OK)
					{UG_LOG("Error while assembling Matrix.\n"); return false;}
			}

			if(m_ass.assemble_solution(u) != IAssemble_OK)
				{UG_LOG("Error while assembling solution.\n"); return false;}

			return true;
		}

		// compute f = L*u (here, L is a Matrix)
		virtual bool apply(domain_function_type& u, codomain_function_type& f)
		{
			#ifdef UG_PARALLEL
			if(!u.has_storage_type(PST_CONSISTENT))
				{
					UG_LOG("WARNING: In 'AssembledLinearOperator::apply':Inadequate storage format of Vector u."
							" Use consistent to avoid internal type conversion.\n");
					if(!u.change_storage_type(PST_CONSISTENT)) return false;
				}
			#endif

			const typename domain_function_type::vector_type& x = u.get_vector();
			typename codomain_function_type::vector_type& b = f.get_vector();

			UG_ASSERT(x.size() == m_Matrix.num_rows(), "Row size '" << m_Matrix.num_rows() << "' of Matrix L and size '" << x.size() << "' of Vector x do not match. Cannot calculate L*x.");
			UG_ASSERT(b.size() == m_Matrix.num_cols(), "Column size '" << m_Matrix.num_rows() << "' of Matrix L and size  '" << b.size() << "' of Vector b do not match. Cannot calculate b := L*x.");

			// set storage type to additiv, since it could been additive unique before
			// TODO: Handle this in matrix multiplication
			f.set_storage_type(PST_ADDITIVE);
			return m_Matrix.apply(b,x);
		}

		// f := f - L*u
		virtual bool apply_sub(domain_function_type& u, codomain_function_type& f)
		{
			#ifdef UG_PARALLEL
			if(!f.has_storage_type(PST_ADDITIVE) || !u.has_storage_type(PST_CONSISTENT))
				{
					UG_LOG("WARNING: In 'AssembledLinearOperator::apply_sub':Inadequate storage format of Vectors.\n");
					UG_LOG("                                             use: f additive and u consistent to avoid internal type conversion.\n");
					if(!f.change_storage_type(PST_ADDITIVE)) return false;
					if(!u.change_storage_type(PST_CONSISTENT)) return false;
				}
			#endif
			const typename domain_function_type::vector_type& x = u.get_vector();
			typename codomain_function_type::vector_type& b = f.get_vector();

			UG_ASSERT(x.size() == m_Matrix.num_rows(), "Row size '" << m_Matrix.num_rows() << "' of Matrix L and size '" << x.size() << "' of Vector x do not match. Cannot calculate L*x.");
			UG_ASSERT(b.size() == m_Matrix.num_cols(), "Column size '" << m_Matrix.num_rows() << "' of Matrix L and size  '" << b.size() << "' of Vector b do not match. Cannot calculate b := b - L*x.");

			// set storage type to additiv, since it could been additive unique before
			// TODO: Handle this in matrix multiplication
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
