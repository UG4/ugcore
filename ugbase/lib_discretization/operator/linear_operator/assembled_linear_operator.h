


#ifndef __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__ASSEMBLED_LINEAR_OPERATOR__
#define __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__ASSEMBLED_LINEAR_OPERATOR__

#include "lib_algebra/lib_algebra.h"

namespace ug{

template <typename TFunction>
class AssembledLinearizedOperator : virtual public ILinearizedOperator<TFunction, TFunction>
{
	public:
		// domain function type
		typedef TFunction domain_function_type;

		// codomain function type
		typedef TFunction codomain_function_type;

		// type of algebra
		typedef typename TFunction::algebra_type algebra_type;

	public:
		AssembledLinearizedOperator(IAssemble<TFunction, algebra_type>& ass) :
			m_ass(ass)
		{};

		virtual bool init()
		{
			return true;
		}

		// prepare the operator for application (e.g. compute an intern Matrix J(u))
		virtual bool prepare(TFunction& dOut, TFunction& uIn, TFunction& cIn)
		{
			const typename domain_function_type::vector_type& c_vec = cIn.get_vector();
			typename codomain_function_type::vector_type& d_vec = dOut.get_vector();

			if(m_J.num_rows() == d_vec.size() && m_J.num_cols() == c_vec.size())
			{
				if(m_J.set(0.0) != true) return false;
			}
			else
			{
				if(m_J.destroy() != true) return false;
				if(m_J.create(d_vec.size(), c_vec.size()) != true) return false;
			}

			if(m_ass.assemble_jacobian(m_J, uIn) != IAssemble_OK) return false;

			return true;
		}

		// compute d = J(u)*c (here, J(u) is a Matrix)
		virtual bool apply(TFunction& dOut, TFunction& cIn)
		{
			#ifdef UG_PARALLEL
			if(!cIn.has_storage_type(PST_CONSISTENT))
				{
					UG_LOG("WARNING: In 'AssembledLinearizedOperator::apply':Inadequate storage format of Vector c."
										" Use consistent to avoid internal type conversion.\n");
					if(!cIn.change_storage_type(PST_CONSISTENT)) return false;
				}
			#endif

			const typename domain_function_type::vector_type& x = cIn.get_vector();
			typename codomain_function_type::vector_type& b = dOut.get_vector();

			UG_ASSERT(x.size() == m_J.num_rows(), "Row size '" << m_J.num_rows() << "' of Matrix J and size '" << x.size() << "' of Vector x do not match. Cannot calculate L*x.");
			UG_ASSERT(b.size() == m_J.num_cols(), "Column size '" << m_J.num_rows() << "' of Matrix J and size  '" << b.size() << "' of Vector b do not match. Cannot calculate b := L*x.");

			// set storage type to additiv, since it could been additive unique before
			// TODO: Handle this in matrix multiplication
			#ifdef UG_PARALLEL
			dOut.set_storage_type(PST_ADDITIVE);
			#endif

			return m_J.apply(b,x);
		}

		// d := d - J(u)*c
		virtual bool apply_sub(TFunction& dOut, TFunction& cIn)
		{
			#ifdef UG_PARALLEL
			if(!dOut.has_storage_type(PST_ADDITIVE) || !cIn.has_storage_type(PST_CONSISTENT))
				{
					UG_LOG("WARNING: In 'AssembledLinearizedOperator::apply_sub':Inadequate storage format of Vectors.\n");
					UG_LOG("                                             	 use: d additive and c consistent to avoid internal type conversion.\n");
					if(!dOut.change_storage_type(PST_ADDITIVE)) return false;
					if(!cIn.change_storage_type(PST_CONSISTENT)) return false;
				}
			#endif
			const typename domain_function_type::vector_type& x = cIn.get_vector();
			typename codomain_function_type::vector_type& b = dOut.get_vector();

			UG_ASSERT(x.size() == m_J.num_rows(), "Row size '" << m_J.num_rows() << "' of Matrix J and size '" << x.size() << "' of Vector x do not match. Cannot calculate L*x.");
			UG_ASSERT(b.size() == m_J.num_cols(), "Column size '" << m_J.num_rows() << "' of Matrix J and size  '" << b.size() << "' of Vector b do not match. Cannot calculate b := b - L*x.");

			// set storage type to additiv, since it could been additive unique before
			// TODO: Handle this in matrix multiplication
			#ifdef UG_PARALLEL
			dOut.set_storage_type(PST_ADDITIVE);
			#endif

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
class AssembledLinearOperator : public AssembledLinearizedOperator<TFunction>, virtual public ILinearOperator<TFunction, TFunction>
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

		// otherwise not clear: AssembledLinearizedOperator<TFunction>::prepare(u, c, d) or ILinearOperator(u,c,d)
		virtual bool prepare(TFunction& dOut, TFunction& uIn, TFunction& cIn)
		{
			return AssembledLinearizedOperator<TFunction>::prepare(dOut, uIn, cIn);
		}


		// prepare the operator for application (e.g. compute an intern Matrix L)
		virtual bool prepare(TFunction& fOut, TFunction& uIn)
		{
			const typename domain_function_type::vector_type& x = uIn.get_vector();
			typename codomain_function_type::vector_type& b = fOut.get_vector();

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
				if(m_ass.assemble_linear(m_Matrix, b, uIn) != IAssemble_OK)
					{UG_LOG("Error while assembling Matrix and rhs.\n"); return false;}
				#ifdef UG_PARALLEL
				fOut.set_storage_type(PST_ADDITIVE);
				#endif
			}
			else
			{
				if(m_ass.assemble_jacobian(m_Matrix, uIn) != IAssemble_OK)
					{UG_LOG("Error while assembling Matrix.\n"); return false;}
			}

			if(m_ass.assemble_solution(uIn) != IAssemble_OK)
				{UG_LOG("Error while assembling solution.\n"); return false;}

			return true;
		}

		// compute f = L*u (here, L is a Matrix)
		virtual bool apply(TFunction& fOut, TFunction& uIn)
		{
			#ifdef UG_PARALLEL
			if(!uIn.has_storage_type(PST_CONSISTENT))
				{
					UG_LOG("WARNING: In 'AssembledLinearOperator::apply':Inadequate storage format of Vector u."
							" Use consistent to avoid internal type conversion.\n");
					if(!uIn.change_storage_type(PST_CONSISTENT)) return false;
				}
			#endif

			const typename domain_function_type::vector_type& x = uIn.get_vector();
			typename codomain_function_type::vector_type& b = fOut.get_vector();

			UG_ASSERT(x.size() == m_Matrix.num_rows(), "Row size '" << m_Matrix.num_rows() << "' of Matrix L and size '" << x.size() << "' of Vector x do not match. Cannot calculate L*x.");
			UG_ASSERT(b.size() == m_Matrix.num_cols(), "Column size '" << m_Matrix.num_rows() << "' of Matrix L and size  '" << b.size() << "' of Vector b do not match. Cannot calculate b := L*x.");

			// set storage type to additiv, since it could been additive unique before
			// TODO: Handle this in matrix multiplication
			#ifdef UG_PARALLEL
			fOut.set_storage_type(PST_ADDITIVE);
			#endif

			return m_Matrix.apply(b,x);
		}

		// f := f - L*u
		virtual bool apply_sub(TFunction& fOut, TFunction& uIn)
		{
			#ifdef UG_PARALLEL
			if(!fOut.has_storage_type(PST_ADDITIVE) || !uIn.has_storage_type(PST_CONSISTENT))
				{
					UG_LOG("WARNING: In 'AssembledLinearOperator::apply_sub':Inadequate storage format of Vectors.\n");
					UG_LOG("                                             use: f additive and u consistent to avoid internal type conversion.\n");
					if(!fOut.change_storage_type(PST_ADDITIVE)) return false;
					if(!uIn.change_storage_type(PST_CONSISTENT)) return false;
				}
			#endif
			const typename domain_function_type::vector_type& x = uIn.get_vector();
			typename codomain_function_type::vector_type& b = fOut.get_vector();

			UG_ASSERT(x.size() == m_Matrix.num_rows(), "Row size '" << m_Matrix.num_rows() << "' of Matrix L and size '" << x.size() << "' of Vector x do not match. Cannot calculate L*x.");
			UG_ASSERT(b.size() == m_Matrix.num_cols(), "Column size '" << m_Matrix.num_rows() << "' of Matrix L and size  '" << b.size() << "' of Vector b do not match. Cannot calculate b := b - L*x.");

			// set storage type to additiv, since it could been additive unique before
			// TODO: Handle this in matrix multiplication
			#ifdef UG_PARALLEL
			fOut.set_storage_type(PST_ADDITIVE);
			#endif

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
