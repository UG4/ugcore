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

///////////////////////////////////////////////////////////////////////////////
// Iterator Operator
///////////////////////////////////////////////////////////////////////////////

/// describes a linear iterator
/**
 * This class is the base class for all linear iterators. Iterators (also called
 * preconditioners) are used in iterative schemes when solving a linear system.
 * Usually, a linear problem like L*u = f is intended to be solved. This is
 * done in an iterative way by performing an iteration of
 *
 * 	start: compute d := f - L*u
 * 	iterate: 	- 	c := B*d 		(compute correction)
 * 				-	u := u + c		(update solution)
 * 				- 	d := d - L*c	(update defect)
 *
 * This iterator class describes the application of B in the scheme above.
 * The application has been split up into two parts.
 *
 * 1. init(L, u) or init(L):
 * 		These methods initialize the iterator and one of these methods has to
 * 		be called before any of the apply methods can be used. Passing the
 * 		linear operator indicates that this operator is used as underlying
 * 		for the iterator.
 *
 * 2. apply or apply_return_defect:
 * 		These methods are used to compute the correction (and to update the
 * 		defect at the same time). Note, that these methods can only be called
 * 		when the iterator has been initialized.
 *
 * This splitting has been made, since initialization may be computationally
 * expansive. Thus, the user of this class has the choice when to call this
 * initialization. E.g. when the operator is applied several times the init of
 * the iterator is only needed once.
 *
 * \tparam	X 	Domain space function
 * \tparam	Y	Range space function
 */
template <typename X, typename Y>
class ILinearIterator
{
	public:
	///	Domain space
		typedef X domain_function_type;

	///	Range space
		typedef Y codomain_function_type;

	public:
	///	returns the name of iterator
	/**
	 * This method returns the name of the iterator operator. This function is
	 * typically needed, when the iterator operator is used inside of another
	 * operator and some debug output should be printed
	 *
	 * \returns 	const char* 	name of inverse operator
	 */
		virtual const char* name() const = 0;

	///	initialize for operator J(u) and linearization point u
	/**
	 * This method passes the linear operator J(u) that should be used as
	 * underlying by this iterator. As second argument the linearization point
	 * is passed. This is needed e.g. for the geometric multigrid method.
	 *
	 * \param[in]	J		linearized operator to use as underlying
	 * \param[in]	u		linearization point
	 * \returns		bool	success flag
	 */
		virtual bool init(ILinearOperator<Y,X>& J, const Y& u) = 0;

	///	initialize for linear operator L
	/**
	 * This method passes the operator L that used as underlying by this
	 * operator. In addition some preparation step can be made.
	 *
	 * \param[in]	L		linear operator to use as underlying
	 * \returns		bool	success flag
	 */
		virtual bool init(ILinearOperator<Y,X>& L) = 0;

	///	compute new correction c = B*d
	/**
	 * This method applies the iterator operator, i.e. c = B*d. The
	 * domain function d remains unchanged.
	 * Note, that this method can always be implemented by creating a copy of
	 * d and calling apply_update_defect with this copy.
	 *
	 * \param[in]	d		defect
	 * \param[out]	c		correction
	 * \returns		bool	success flag
	 */
		virtual bool apply(Y& c, const X& d) = 0;

	///	compute new correction c = B*d and update defect d := d - A*c
	/**
	 * This method applies the inverse operator, i.e. c = B*d. The
	 * domain function d is changed in the way, that the defect d := d - A*c
	 * is returned in the function. This is always useful, when the iterating
	 * algorithm can (or must) update the defect during computation (this is
	 * e.g. the case for the geometric multigrid method).
	 * Note, that this method can always be implemented by calling apply and
	 * then computing d := d - A*c.
	 *
	 * \param[in,out]	d		defect
	 * \param[out]		u		correction
	 * \returns			bool	success flag
	 */
		virtual bool apply_update_defect(Y& c, X& d) = 0;

	///	clone
		virtual ILinearIterator<X,Y>* clone() = 0;

	/// virtual destructor
		virtual ~ILinearIterator() {};
};

///////////////////////////////////////////////////////////////////////////////
// Matrix Iterator Operator (Preconditioner)
///////////////////////////////////////////////////////////////////////////////

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
			if(!m_pMatrix->matmul_minus(d, c))
			{
				UG_LOG("ERROR in 'IPreconditioner::apply_update_defect': Cannot execute matmul_minus to compute d:=d-A*c.\n");
				return false;
			}

		//	we're done
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

