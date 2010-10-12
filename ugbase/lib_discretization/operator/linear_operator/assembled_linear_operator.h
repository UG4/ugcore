


#ifndef __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__ASSEMBLED_LINEAR_OPERATOR__
#define __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__ASSEMBLED_LINEAR_OPERATOR__

#include "lib_algebra/operator/operator_interface.h"

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

	//	Type of DoFDistribution
		typedef TDoFDistribution dof_distribution_type;

	public:
		AssembledLinearOperator() :
			m_bInit(false), m_bAssembleRhs(false),
			m_pAss(NULL), m_pDoFDistribution(NULL)
			{};

		AssembledLinearOperator(IAssemble<dof_distribution_type, algebra_type>& ass, bool assemble_rhs = false) :
			m_bInit(false), m_bAssembleRhs(assemble_rhs),
			m_pAss(&ass), m_pDoFDistribution(NULL)
		{};


		void set_discretization(IAssemble<dof_distribution_type, algebra_type>& ass) {m_pAss = &ass;}
		void export_rhs(bool assemble_rhs) {m_bAssembleRhs = assemble_rhs;}

		bool set_dof_distribution(const TDoFDistribution& dofDistr)
		{
			m_pDoFDistribution = &dofDistr;
			return true;
		}

		const TDoFDistribution* get_dof_distribution()
		{
			return m_pDoFDistribution;
		}

		virtual bool init(const vector_type& uIn)
		{
			if(m_pDoFDistribution == NULL)
			{
				UG_LOG("ERROR in AssembledLinearizedOperator::init: DoF Distribution not set.\n");
				return false;
			}

			if(m_pAss == NULL)
			{
				UG_LOG("ERROR in AssembledLinearizedOperator::init: Assembling rountine not set.\n");
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
				UG_LOG("ERROR in AssembledLinearizedOperator::init: Cannot assemble Jacobi matrix.\n");
				return false;
			}

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
				UG_LOG("ERROR in AssembledLinearOperator::prepare: DoF Distribution not set.\n");
				return false;
			}

		//	get number of dofs
			const size_t numDoFs = m_pDoFDistribution->num_dofs();

		//	Resize Matrix and set to zero
			if(m_J.num_rows() == numDoFs && m_J.num_cols() == numDoFs)
			{
				if(!m_J.set(0.0))
				{
					UG_LOG("ERROR in AssembledLinearOperator::init: Cannot set matrix to zero.\n");
					return false;
				}
			}
			else
			{
				if(!m_J.resize(numDoFs, numDoFs))
				{
					UG_LOG("ERROR in AssembledLinearOperator::prepare: Cannot resize matrix.\n");
					return false;
				}
			}

		//	Resize rhs
			if(m_rhs.size() != numDoFs)
			{
				if(!m_rhs.resize(numDoFs))
				{
					UG_LOG("ERROR in AssembledLinearOperator::init: Cannot resize rhs.\n");
					return false;
				}
			}

		//	Compute matrix (and rhs if needed)
			if(m_bAssembleRhs)
			{
				m_rhs.set(0.0);
				if(m_pAss->assemble_linear(m_J, m_rhs, m_rhs, *m_pDoFDistribution) != IAssemble_OK)
					{UG_LOG("Error while assembling Matrix and rhs.\n"); return false;}
				#ifdef UG_PARALLEL
				m_rhs.set_storage_type(PST_ADDITIVE);
				#endif
			}
			else
			{
				if(m_pAss->assemble_jacobian(m_J, m_rhs, *m_pDoFDistribution) != IAssemble_OK)
					{UG_LOG("Error while assembling Matrix.\n"); return false;}
			}

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

			UG_ASSERT(cIn.size() == m_J.num_rows(), "Row size '" << m_J.num_rows() << "' of Matrix J and size '" << cIn.size() << "' of Vector x do not match. Cannot calculate L*x.");
			UG_ASSERT(dOut.size() == m_J.num_cols(), "Column size '" << m_J.num_rows() << "' of Matrix J and size  '" << dOut.size() << "' of Vector b do not match. Cannot calculate b := L*x.");

			// set storage type to additiv, since it could been additive unique before
			// TODO: Handle this in matrix multiplication
			#ifdef UG_PARALLEL
			dOut.set_storage_type(PST_ADDITIVE);
			#endif

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
				UG_LOG("ERROR: In 'LaplackLUSolver::apply':Inadequate storage format of Vector d.\n");
				return false;
			}
			if(!cIn.has_storage_type(PST_CONSISTENT))
			{
				UG_LOG("ERROR: In 'LaplackLUSolver::apply':Inadequate storage format of Vector c.\n");
				return false;
			}
#endif

			UG_ASSERT(cIn.size() == m_J.num_rows(), "Row size '" << m_J.num_rows() << "' of Matrix J and size '" << cIn.size() << "' of Vector x do not match. Cannot calculate L*x.");
			UG_ASSERT(dOut.size() == m_J.num_cols(), "Column size '" << m_J.num_rows() << "' of Matrix J and size  '" << dOut.size() << "' of Vector b do not match. Cannot calculate b := b - L*x.");

			// set storage type to additiv, since it could been additive unique before
			// TODO: Handle this in matrix multiplication
			#ifdef UG_PARALLEL
			dOut.set_storage_type(PST_ADDITIVE);
			#endif

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
				UG_LOG("ERROR in AssembledLinearizedOperator::set_dirichlet_values: Operator not initialized.\n");
				return false;
			}

			if(m_pAss->assemble_solution(u, *m_pDoFDistribution) != IAssemble_OK)
				{UG_LOG("Error while assembling solution.\n"); return false;}
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
		IAssemble<dof_distribution_type, algebra_type>* m_pAss;

		// DoF Distribution used
		const TDoFDistribution* m_pDoFDistribution;

		// matrix storage
		matrix_type m_J;

		// vector storage
		vector_type m_rhs;
};

	/*

template <typename TDoFDistribution, typename TAlgebra>
class AssembledLinearOperator :
	public AssembledLinearizedOperator<TDoFDistribution, TAlgebra>,
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

	//	Type of DoFDistribution
		typedef TDoFDistribution dof_distribution_type;

	public:
		AssembledLinearOperator(IAssemble<dof_distribution_type, algebra_type>& ass, bool assemble_rhs = false) :
			AssembledLinearizedOperator<TDoFDistribution, TAlgebra>(ass),
			m_bInit(false), m_assemble_rhs(assemble_rhs),
			m_pDoFDistribution(NULL), m_ass(ass)
		{};

	//	Set the DoF Distribution (i.e. the trial space)
		bool set_dof_distribution(const TDoFDistribution& dofDistr)
		{
			m_pDoFDistribution = &dofDistr;
			return true;
		}

		virtual bool init(const vector_type& u)
		{
			return init();
		}


	// 	Prepare functions
		virtual bool prepare(vector_type& fOut, vector_type& uIn)
		{
			if(!m_bInit)
			{
				UG_LOG("ERROR in AssembledLinearOperator::prepare: Operator not initialized.\n");
				return false;
			}

			if(m_ass.assemble_solution(uIn, *m_pDoFDistribution) != IAssemble_OK)
				{UG_LOG("Error while assembling solution.\n"); return false;}

			return true;
		}

	// 	Compute f = L*u (here, L is a Matrix)
		virtual bool apply(vector_type& fOut, const vector_type& uIn)
		{
			if(!m_bInit)
			{
				UG_LOG("ERROR in AssembledLinearOperator::apply: Operator not initialized.\n");
				return false;
			}

			#ifdef UG_PARALLEL
			if(!uIn.has_storage_type(PST_CONSISTENT))
			{
				UG_LOG("ERROR: In 'AssembledLinearOperator::apply':Inadequate storage format of Vector u.");
				return false;
			}
			#endif

			UG_ASSERT(uIn.size() == m_Matrix.num_rows(), "Row size '" << m_Matrix.num_rows() << "' of Matrix L and size '" << uIn.size() << "' of Vector x do not match. Cannot calculate L*x.");
			UG_ASSERT(fOut.size() == m_Matrix.num_cols(), "Column size '" << m_Matrix.num_rows() << "' of Matrix L and size  '" << fOut.size() << "' of Vector b do not match. Cannot calculate b := L*x.");

			// set storage type to additiv, since it could been additive unique before
			// TODO: Handle this in matrix multiplication
			#ifdef UG_PARALLEL
			fOut.set_storage_type(PST_ADDITIVE);
			#endif

			return m_Matrix.apply(fOut,uIn);
		}

	// 	Compute f := f - L*u
		virtual bool apply_sub(vector_type& fOut, const vector_type& uIn)
		{
			if(!m_bInit)
			{
				UG_LOG("ERROR in AssembledLinearOperator::apply_sub: Operator not initialized.\n");
				return false;
			}

			#ifdef UG_PARALLEL
			if(!fOut.has_storage_type(PST_ADDITIVE) || !uIn.has_storage_type(PST_CONSISTENT))
			{
				UG_LOG("ERROR: In 'AssembledLinearOperator::apply_sub':Inadequate storage format of Vectors.\n");
				return false;
			}
			#endif

			UG_ASSERT(uIn.size() == m_Matrix.num_rows(), "Row size '" << m_Matrix.num_rows() << "' of Matrix L and size '" << uIn.size() << "' of Vector x do not match. Cannot calculate L*x.");
			UG_ASSERT(fOut.size() == m_Matrix.num_cols(), "Column size '" << m_Matrix.num_rows() << "' of Matrix L and size  '" << fOut.size() << "' of Vector b do not match. Cannot calculate b := b - L*x.");

			// set storage type to additiv, since it could been additive unique before
			// TODO: Handle this in matrix multiplication
			#ifdef UG_PARALLEL
			fOut.set_storage_type(PST_ADDITIVE);
			#endif

			return m_Matrix.matmul_minus(fOut, uIn);
		}

	//	Export matrix
		virtual matrix_type& get_matrix() {return m_Matrix;}

	//	Export assembled rhs
		const vector_type& get_rhs() const {return m_rhs;}

	// 	Destructor
		virtual ~AssembledLinearOperator() {m_Matrix.destroy();};

	protected:
		// init flag
		bool m_bInit;

		// choose, weather rhs should be assembled as well.
		bool m_assemble_rhs;

		// DoF Distribution used
		const TDoFDistribution* m_pDoFDistribution;

		// assembling procedure
		IAssemble<dof_distribution_type, algebra_type>& m_ass;

		// matrix storage
		matrix_type m_Matrix;

		// vector storage
		vector_type m_rhs;
};


*/



} // namespace ug

#endif /* __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__ASSEMBLED_LINEAR_OPERATOR__ */
