


#ifndef __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__ASSEMBLED_LINEAR_OPERATOR__
#define __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__ASSEMBLED_LINEAR_OPERATOR__

#include "lib_algebra/operator/operator_interface.h"

#ifdef UG_PARALLEL
#include "lib_discretization/parallelization/parallelization_util.h"
#endif

namespace ug{

template <typename TDoFDistribution, typename TAlgebra>
class AssembledLinearOperator :
	public virtual IMatrixOperator<	typename TAlgebra::vector_type,
									typename TAlgebra::vector_type,
									typename TAlgebra::matrix_type>
{
	public:
	// 	Type of algebra
		typedef TAlgebra algebra_type;

	//	Type of Vector
		typedef typename TAlgebra::vector_type vector_type;

	//	Type of Vector
		typedef typename TAlgebra::matrix_type matrix_type;

	public:
		AssembledLinearOperator() :
			m_bInit(false), m_bAssembleRhs(false),
			m_pAss(NULL), m_pDoFDistribution(NULL)
			{};

		AssembledLinearOperator(IAssemble<TDoFDistribution, algebra_type>& ass,
		                        bool assemble_rhs = false)
		:	m_bInit(false), m_bAssembleRhs(assemble_rhs),
			m_pAss(&ass), m_pDoFDistribution(NULL)
		{};


		void set_discretization(IAssemble<TDoFDistribution, algebra_type>& ass)
		{
			m_pAss = &ass;
		}

		void export_rhs(bool assemble_rhs) {m_bAssembleRhs = assemble_rhs;}

		bool set_dof_distribution(const IDoFDistribution<TDoFDistribution>& dofDistr)
		{
			m_pDoFDistribution = &dofDistr;
			return true;
		}

		const IDoFDistribution<TDoFDistribution>* get_dof_distribution()
		{
			return m_pDoFDistribution;
		}

		virtual bool init(const vector_type& uIn)
		{
			if(m_pDoFDistribution == NULL)
			{
				UG_LOG("ERROR in AssembledLinearizedOperator::init: "
						"DoF Distribution not set.\n");
				return false;
			}

			if(m_pAss == NULL)
			{
				UG_LOG("ERROR in AssembledLinearizedOperator::init:"
						" Assembling rountine not set.\n");
				return false;
			}

		//	get number of dofs
			const size_t numDoFs = m_pDoFDistribution->num_dofs();

		//	resize matrix and set to zero
			if(m_J.num_rows() == numDoFs && m_J.num_cols() == numDoFs)
			{
				if(!m_J.set(0.0)) return false;
			}
			else
			{
				if(!m_J.destroy()) return false;
				if(!m_J.create(numDoFs, numDoFs)) return false;
			}

		//	assemble matrix (depending on u, i.e. J(u))
			if(m_pAss->assemble_jacobian(m_J, uIn, *m_pDoFDistribution) != IAssemble_OK)
			{
				UG_LOG("ERROR in AssembledLinearizedOperator::init:"
						" Cannot assemble Jacobi matrix.\n");
				return false;
			}

		//	Remember parallel storage type
			#ifdef UG_PARALLEL
				m_J.set_storage_type(PST_ADDITIVE);
				IDoFDistribution<TDoFDistribution>* dist = const_cast<IDoFDistribution<TDoFDistribution>*>(m_pDoFDistribution);
				CopyLayoutsAndCommunicatorIntoMatrix(m_J, *dist);
			#endif

		//	remember that operator is initialized
			m_bInit = true;
			return true;
		}

	//	Initialize the operator
		virtual bool init()
		{
		//	todo: check that assembling is linear

		//	check if DoF Distribution is set
			if(m_pDoFDistribution == NULL)
			{
				UG_LOG("ERROR in AssembledLinearOperator::prepare:"
						" DoF Distribution not set.\n");
				return false;
			}

		//	get number of dofs
			const size_t numDoFs = m_pDoFDistribution->num_dofs();

		//	Resize Matrix and set to zero
			if(m_J.num_rows() == numDoFs && m_J.num_cols() == numDoFs)
			{
				if(!m_J.set(0.0))
				{
					UG_LOG("ERROR in AssembledLinearOperator::init: "
							"Cannot set matrix to zero.\n");
					return false;
				}
			}
			else
			{
				if(!m_J.resize(numDoFs, numDoFs))
				{
					UG_LOG("ERROR in AssembledLinearOperator::prepare:"
							" Cannot resize matrix.\n");
					return false;
				}
			}

		//	Resize rhs
			if(m_rhs.size() != numDoFs)
			{
				if(!m_rhs.resize(numDoFs))
				{
					UG_LOG("ERROR in AssembledLinearOperator::init:"
							" Cannot resize rhs.\n");
					return false;
				}
			}

		//	Compute matrix (and rhs if needed)
			if(m_bAssembleRhs)
			{
				m_rhs.set(0.0);
				if(m_pAss->assemble_linear(m_J, m_rhs, m_rhs, *m_pDoFDistribution) != IAssemble_OK)
				{
					UG_LOG("ERROR in AssembledLinearOperator::init:"
							" Cannot assemble Matrix and Rhs.\n");
					return false;
				}
				#ifdef UG_PARALLEL
				m_rhs.set_storage_type(PST_ADDITIVE);
				#endif
			}
			else
			{
				if(m_pAss->assemble_jacobian(m_J, m_rhs, *m_pDoFDistribution) != IAssemble_OK)
				{
					UG_LOG("ERROR in AssembledLinearOperator::init:"
							" Cannot assemble Matrix.\n");
					return false;
				}
			}

			#ifdef UG_PARALLEL
			m_J.set_storage_type(PST_ADDITIVE);
			IDoFDistribution<TDoFDistribution>* dist =
					const_cast<IDoFDistribution<TDoFDistribution>*>(m_pDoFDistribution);
			CopyLayoutsAndCommunicatorIntoMatrix(m_J, *dist);
			#endif

		//	Remember that operator has been initialized
			m_bInit = true;

			return true;
		}


	//	Compute d = J(u)*c (here, J(u) is a Matrix)
		virtual bool apply(vector_type& dOut, const vector_type& cIn)
		{
			if(!m_bInit)
			{
				UG_LOG("ERROR in AssembledLinearizedOperator::apply: Operator not initialized.\n");
				return false;
			}

			#ifdef UG_PARALLEL
			if(!cIn.has_storage_type(PST_CONSISTENT))
				{
					UG_LOG("WARNING: In 'AssembledLinearizedOperator::apply':Inadequate storage format of Vector c.");
					return false;
				}
			#endif

			UG_ASSERT(cIn.size() == m_J.num_rows(), "Row size '" << m_J.num_rows() << "' of Matrix J and size '"
													<< cIn.size() << "' of Vector x do not match. Cannot calculate L*x.");
			UG_ASSERT(dOut.size() == m_J.num_cols(), "Column size '" << m_J.num_rows() << "' of Matrix J and size  '"
													<< dOut.size() << "' of Vector b do not match. Cannot calculate b := L*x.");

		//	Apply Matrix
			return m_J.apply(dOut, cIn);
		}

	//	Compute d := d - J(u)*c
		virtual bool apply_sub(vector_type& dOut, const vector_type& cIn)
		{
			if(!m_bInit)
			{
				UG_LOG("ERROR in AssembledLinearizedOperator::apply_sub: Operator not initialized.\n");
				return false;
			}

#ifdef UG_PARALLEL
			if(!dOut.has_storage_type(PST_ADDITIVE))
			{
				UG_LOG("ERROR: In 'AssembledLinearizedOperator::apply_sub':Inadequate storage format of Vector d.\n");
				return false;
			}
			if(!cIn.has_storage_type(PST_CONSISTENT))
			{
				UG_LOG("ERROR: In 'AssembledLinearizedOperator::apply_sub':Inadequate storage format of Vector c.\n");
				return false;
			}
#endif

			UG_ASSERT(cIn.size() == m_J.num_rows(), "Row size '" << m_J.num_rows() << "' of Matrix J and size '"
													<< cIn.size() << "' of Vector x do not match. Cannot calculate L*x.");
			UG_ASSERT(dOut.size() == m_J.num_cols(), "Column size '" << m_J.num_rows() << "' of Matrix J and size  '"
													<< dOut.size() << "' of Vector b do not match. Cannot calculate b := b - L*x.");

		//	Apply Matrix
			return m_J.matmul_minus(dOut,cIn);
		}

	//	Export matrix
		virtual matrix_type& get_matrix()
		{
			return m_J;
		}

	//	Export assembled rhs
		const vector_type& get_rhs() const {return m_rhs;}

	//	Set Dirichlet values
		bool set_dirichlet_values(vector_type& u)
		{
			if(!m_bInit)
			{
				UG_LOG("ERROR in AssembledLinearizedOperator::set_dirichlet_values:"
						" Operator not initialized.\n");
				return false;
			}

			if(m_pAss->assemble_solution(u, *m_pDoFDistribution) != IAssemble_OK)
			{
				UG_LOG("ERROR in AssembledLinearOperator::set_dirichlet_values:"
						" Cannot assemble solution.\n");
				return false;
			}
			return true;
		}

	// 	Destructor
		virtual ~AssembledLinearOperator()
		{
			m_J.destroy();
			m_rhs.destroy();
		};

	protected:
		// init flag
		bool m_bInit;

		// assemble rhs flag
		bool m_bAssembleRhs;

		// assembling procedure
		IAssemble<TDoFDistribution, algebra_type>* m_pAss;

		// DoF Distribution used
		const IDoFDistribution<TDoFDistribution>* m_pDoFDistribution;

		// matrix storage
		matrix_type m_J;

		// vector storage
		vector_type m_rhs;
};


} // namespace ug

#endif /* __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__ASSEMBLED_LINEAR_OPERATOR__ */
