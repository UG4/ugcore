
#ifndef __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__ASSEMBLED_LINEAR_OPERATOR__
#define __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__ASSEMBLED_LINEAR_OPERATOR__

#include "lib_algebra/lib_algebra.h"
#include "lib_discretization/io/vtkoutput.h"


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

			if(!c.has_storage_type(GFST_CONSISTENT)) return false;
			d.set_storage_type(GFST_ADDITIVE);
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
			if(!c.has_storage_type(GFST_CONSISTENT)) return false;
			if(!d.has_storage_type(GFST_ADDITIVE)) return false;
			d.set_storage_type(GFST_ADDITIVE);
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
				f.set_storage_type(GFST_ADDITIVE);
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

			// TODO: Check that matrix has correct type (additive)
			if(!u.has_storage_type(GFST_CONSISTENT)) return false;
			f.set_storage_type(GFST_ADDITIVE);
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
			if(!u.has_storage_type(GFST_CONSISTENT)) return false;
			if(!f.has_storage_type(GFST_ADDITIVE)) return false;
			f.set_storage_type(GFST_ADDITIVE);
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
		AssembledJacobiOperator(number damp) : m_damp(damp), m_bOpChanged(true)
		{};

		bool init(IDiscreteLinearizedOperator<domain_function_type,codomain_function_type>& Op)
		{
			AssembledDiscreteLinearizedOperator<TDiscreteFunction>* A
				= dynamic_cast<AssembledDiscreteLinearizedOperator<TDiscreteFunction>*>(&Op);
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
				size_t size = m_pMatrix->row_size();
				if(size != m_pMatrix->col_size())
				{
					UG_LOG("Square Matrix needed for Jacobi Iteration.\n");
					return false;
				}

				if(m_diagInv.size() == size) m_diagInv.set(0.0);
				else {
					m_diagInv.destroy();
					m_diagInv.create(size);
				}

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

				typename domain_function_type::dof_manager_type& dofManager = d.get_dof_manager();

				//	make diagonal consistent
				AdditiveToConsistent(	&m_diagInv,
										dofManager.get_master_layout(d.get_level()),
										dofManager.get_slave_layout(d.get_level()));

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
			if(!d.has_storage_type(GFST_ADDITIVE)) return false;

			typename domain_function_type::vector_type& d_vec = d.get_vector();
			typename codomain_function_type::vector_type& c_vec = c.get_vector();

			UG_ASSERT(d_vec.size() == m_pMatrix->row_size(),	"Vector and Row sizes have to match!");
			UG_ASSERT(c_vec.size() == m_pMatrix->col_size(), "Vector and Column sizes have to match!");
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
			c.set_storage_type(GFST_ADDITIVE);
			if(c.change_storage_type(GFST_CONSISTENT) != true)
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
			d.set_storage_type(GFST_ADDITIVE);
			return true;
		}

		// destructor
		virtual ~AssembledJacobiOperator() {};

	protected:
		AssembledDiscreteLinearizedOperator<TDiscreteFunction>* m_pOperator;

		typename algebra_type::matrix_type* m_pMatrix;

		typename domain_function_type::vector_type m_diagInv;

		number m_damp;

		bool m_bOpChanged;
};



} // namespace ug

#endif /* __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__ASSEMBLED_LINEAR_OPERATOR__ */
