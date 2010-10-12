/*
 * operator_iterator_interface.h
 *
 *  Created on: 22.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_ALGEBRA__OPERATOR__OPERATOR_ITERATOR_INTERFACE__
#define __H__LIB_ALGEBRA__OPERATOR__OPERATOR_ITERATOR_INTERFACE__

#include "operator_interface.h"
#ifdef UG_PARALLEL
	#include "lib_algebra/parallelization/parallelization.h"
#endif

namespace ug{

/* This Operator type behaves different on application. It not only computes c = B*d, but also changes d. */
/* It is used in iterative schemes. */
template <typename X, typename Y>
class ILinearIterator
{
	public:
	// 	Domain space
		typedef X domain_function_type;

	// 	Range space
		typedef Y codomain_function_type;

	public:
	// 	Prepare for Operator J(u) and linearization point u (current solution)
		virtual bool init(ILinearOperator<Y,X>& J, const Y& u) = 0;

	//	Prepare for Linear Operartor L
		virtual bool init(ILinearOperator<Y,X>& L) = 0;

	//	Compute new correction c = B*d
		virtual bool apply(Y& c, const X& d) = 0;

	//	Compute new correction c = B*d and return new defect d := d - A*c
		virtual bool apply_update_defect(Y& c, X& d) = 0;

	//	Clone
		virtual ILinearIterator<X,Y>* clone() = 0;

	// 	Destructor
		virtual ~ILinearIterator() {};
};


template <typename TAlgebra>
class IPreconditioner :
	public virtual ILinearIterator<	typename TAlgebra::vector_type,
									typename TAlgebra::vector_type>
{
	public:
	//	Algebra type
		typedef TAlgebra algebra_type;

	//	Vector type
		typedef typename TAlgebra::vector_type vector_type;

	//	Matrix type
		typedef typename TAlgebra::matrix_type matrix_type;

	public:
		IPreconditioner() :
			m_pOperator(NULL), m_pMatrix(NULL), m_bInit(false)
		{};

	protected:
	//	Name of preconditioner
		virtual const char* name() const = 0;

	//	Preprocess routine
		virtual bool preprocess(matrix_type& mat) = 0;

	//	Stepping routine
		virtual bool step(matrix_type& mat, vector_type& c, const vector_type& d)  = 0;

	//	Postprocess routine
		virtual bool postprocess() = 0;

	public:
	//	Implement general interface
		virtual bool init(ILinearOperator<vector_type, vector_type>& J, const vector_type& u)
		{
			IMatrixOperator<vector_type, vector_type, matrix_type>* Op =
						dynamic_cast<IMatrixOperator<vector_type, vector_type, matrix_type>*>(&J);

		//	Check that matrix if of correct type
			if(Op == NULL)
			{
				UG_LOG("ERROR in '" << name() << "::init': Passed Operator is not based on matrix.\n"
						"This Preconditioner can only handle matrix-based operators. Aborting.\n");
				return false;
			}

		//	Pass to own implementation
			return init(*Op);
		}

	//	Implement general interface
		virtual bool init(ILinearOperator<vector_type, vector_type>& L)
		{
			IMatrixOperator<vector_type, vector_type, matrix_type>* Op =
						dynamic_cast<IMatrixOperator<vector_type, vector_type, matrix_type>*>(&L);

		//	Check that matrix if of correct type
			if(Op == NULL)
			{
				UG_LOG("ERROR in '" << name() << "::init': Passed Operator is not based on matrix.\n"
						"This Preconditioner can only handle matrix-based operators. Aborting.\n");
				return false;
			}

		//	Pass to own implementation
			return init(*Op);
		}

	//	Init operator
		bool init(IMatrixOperator<vector_type, vector_type, matrix_type>& Op)
		{
		// 	Remember operator
			m_pOperator = &Op;

		//	Remember matrix
			m_pMatrix = &m_pOperator->get_matrix();

		//	Check that matrix exists
			if(m_pMatrix == NULL)
			{
				UG_LOG("ERROR in '"<< name() << "::init': Matrix not found, though matrix-based operator given.\n");
				return false;
			}

		//	Preprocess
			if(!preprocess(*m_pMatrix))
			{
				UG_LOG("ERROR in '"<< name() << "::init': Preprocess failed.\n");
				return false;
			}

		//	Remember, that operator has been initialized
			m_bInit = true;

			return true;
		}

	// 	Compute new correction c = B*d
		virtual bool apply(vector_type& c, const vector_type& d)
		{
		//	Check that operator is initialized
			if(!m_bInit)
			{
				UG_LOG("ERROR in '" << name() << "::prepare': Iterator not initialized.\n");
				return false;
			}

		//	Check parallel status
#ifdef UG_PARALLEL
			if(!d.has_storage_type(PST_ADDITIVE))
			{
				UG_LOG("ERROR in '" << name() << "::apply': Wrong parallel storage format. Defect must be additive.\n");
				return false;
			}
#endif

		//	Some Assertions
			UG_ASSERT(d.size() == m_pMatrix->num_rows(),"Vector [size= " << d.size() << "] and Row [size= " << m_pMatrix->num_rows() <<"] sizes have to match!");
			UG_ASSERT(c.size() == m_pMatrix->num_cols(),"Vector [size= " << c.size() << "] and Column [size= " << m_pMatrix->num_cols() <<"] sizes have to match!");
			UG_ASSERT(d.size() == c.size(), 			"Vector [d size= " << d.size() << ", c size = " << c.size() << "] sizes have to match!");

		// 	apply iterator: c = B*d (damp is not used)
			if(!step(*m_pMatrix, c, d))
			{
				UG_LOG("ERROR in '" << name() << "::apply': Step Routine failed.\n");
				return false;
			}

#ifdef 	UG_PARALLEL
		//	Correction is always consistent
			if(!c.change_storage_type(PST_CONSISTENT))
			{
				UG_LOG("ERROR in '" << name() << "::apply': Cannot change parallel storage type of correction to consistent.\n");
				return false;
			}
#endif

		//	we're done
			return true;
		}

	// 	Compute new correction c = B*d and update defect
		virtual bool apply_update_defect(vector_type& c, vector_type& d)
		{
		//	compute new correction
			if(!apply(c, d)) return false;

		// 	update defect d := d - A*c
		// 	TODO: Check that matrix has correct type (additive)
		//	TODO: check return value
			m_pMatrix->matmul_minus(d, c);

#ifdef UG_PARALLEL
			// defect is now no longer unique (maybe we should handle this in matrix multiplication)
			d.set_storage_type(PST_ADDITIVE);
#endif
			return true;
		}

		// destructor
		virtual ~IPreconditioner() {};

	protected:
	//	Operator
		IMatrixOperator<vector_type, vector_type, matrix_type>* m_pOperator;

	//	Matrix
		matrix_type* m_pMatrix;

	// 	init flag
		bool m_bInit;
};

} // end namespace ug

#endif /* __H__LIB_ALGEBRA__OPERATOR__OPERATOR_ITERATOR_INTERFACE__ */

