/*
 * preconditioner.h
 *
 *  Created on: 23.04.2013
 *      Author: andreasvogel
 */

#ifndef __H__LIB_ALGEBRA__OPERATOR__INTERFACE__PRECONDITIONER__
#define __H__LIB_ALGEBRA__OPERATOR__INTERFACE__PRECONDITIONER__

#include "linear_iterator.h"
#include "matrix_operator.h"
#include "lib_algebra/operator/damping.h"
#include "lib_algebra/operator/debug_writer.h"

#ifdef UG_PARALLEL
	#include "lib_algebra/parallelization/parallelization.h"
#endif

namespace ug{
///////////////////////////////////////////////////////////////////////////////
// Matrix Iterator Operator (Preconditioner)
///////////////////////////////////////////////////////////////////////////////

/// describes a linear iterator that is based on a matrix operator
/**
 * This class is the base class for all linear iterators, that act on matrix
 * based linear operators. We call the matrix based iterator a 'Preconditioner'.
 * They are used in iterative schemes when solving a linear system.
 * Usually, a linear problem like L*u = f is intended to be solved. This is
 * done in an iterative way by performing an iteration of
 *
 * 	start: compute d := f - L*u
 * 	iterate: 	- 	c := B*d 		(compute correction)
 * 				-	u := u + c		(update solution)
 * 				- 	d := d - L*c	(update defect)
 *
 * This iterator class describes the application of B in the scheme above, where
 * L is internally represented by a matrix.
 *
 * This method derives from the ILinearIterator-interface and thus must
 * implement these to steps:
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
 * In order to facilitate the implementation of derived classes a default
 * implementation of these virtual method is given. Thus, the implementational
 * part for a matrix based iterator is to implement the three functions:
 *
 * 1. preprocess: Initializes the preconditioner for a matrix, that is used as
 * 				  the underlying matrix. Calling these method again with a
 * 				  different matrix results in a different preconditioner.
 *
 * 2. step:	Computes the new correction.
 *
 * 3. postprocess:	Clean up (if needed)
 *
 * \tparam	TAlgebra 	Type of Algebra
 */
template <typename TAlgebra>
class IPreconditioner :
	public virtual ILinearIterator<typename TAlgebra::vector_type>,
	public DebugWritingObject<TAlgebra>
{
	public:
	///	Algebra type
		typedef TAlgebra algebra_type;

	///	Vector type
		typedef typename TAlgebra::vector_type vector_type;

	///	Matrix type
		typedef typename TAlgebra::matrix_type matrix_type;

	///	Matrix Operator type
		typedef MatrixOperator<matrix_type, vector_type> matrix_operator_type;
		using DebugWritingObject<TAlgebra>::set_debug;
	protected:
		using ILinearIterator<vector_type>::damping;

		using DebugWritingObject<TAlgebra>::debug_writer;
		using DebugWritingObject<TAlgebra>::write_debug;

	public:
	///	default constructor
		IPreconditioner() :
			m_spDefectOperator(NULL), m_spApproxOperator(NULL), m_bInit(false), m_bOtherApproxOperator(false)
		{};

	///	constructor setting debug writer
		IPreconditioner(SmartPtr<IDebugWriter<algebra_type> > spDebugWriter) :
			DebugWritingObject<TAlgebra>(spDebugWriter),
			m_spDefectOperator(NULL), m_spApproxOperator(NULL), m_bInit(false), m_bOtherApproxOperator(false)
		{};

	/// clone constructor
		IPreconditioner( IPreconditioner<TAlgebra> *parent ) :
			ILinearIterator<vector_type>(parent),
			DebugWritingObject<TAlgebra>(parent),
			m_spDefectOperator(NULL), m_spApproxOperator(NULL), m_bInit(false), m_bOtherApproxOperator(false)
		{
		}
	protected:
	///	returns the name of iterator
	/**
	 * This method returns the name of the iterator operator. This function is
	 * typically needed, when the iterator operator is used inside of another
	 * operator and some debug output should be printed
	 *
	 * \returns 	const char* 	name of inverse operator
	 */
		virtual const char* name() const = 0;

	///	initializes the preconditioner
	/**
	 * This method is used to initialize the preconditioner. Usually, here
	 * are performed computationally expensive operations, that should only be
	 * computed once for an underlying matrix (e.g. LU factorization),  while
	 * the preconditioner will by applied (using 'step'-method) several times.
	 *
	 * \param[in]	mat			underlying matrix (i.e. L in L*u = f)
	 * \returns		bool		success flag
	 */
		virtual bool preprocess(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp) = 0;

	///	computes a new correction c = B*d
	/**
	 * This method computes a new correction c = B*d. It can only be called,
	 * when the preprocess has been done.
	 *
	 * \param[in]	mat			underlying matrix (i.e. L in L*u = f)
	 * \param[out]	c			correction
	 * \param[in]	d			defect
	 * \returns		bool		success flag
	 */
		virtual bool step(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp, vector_type& c, const vector_type& d)  = 0;

	///	cleans the operator
		virtual bool postprocess() = 0;

	public:
	///	implements the ILinearIterator-interface for matrix based preconditioner
	/**
	 * This method implements the ILinearIterator interface. It check if the
	 * passed linear operator is matrix based (otherwise this preconditioner can
	 * not be used for the linear operator). Then the request is forwarded to
	 * the implementation of matrix based operators.
	 *
	 * \param[in]	J		linear operator
	 * \param[in]	u		linearization point
	 * \returns		bool	success flag
	 */
		virtual bool init(SmartPtr<ILinearOperator<vector_type> > J,
		                  const vector_type& u)
		{
		//	cast to matrix based operator
			SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp =
					J.template cast_dynamic<MatrixOperator<matrix_type, vector_type> >();

		//	Check that matrix if of correct type
			if(pOp.invalid())
				UG_THROW(name() << "::init': Passed Operator is "
						"not based on matrix. This Preconditioner can only "
						"handle matrix-based operators.");

		//	forward request to matrix based implementation
			return init(pOp);
		}

	///	implements the ILinearIterator-interface for matrix based preconditioner
	/**
	 * This method implements the ILinearIterator interface. It check if the
	 * passed linear operator is matrix based (otherwise this preconditioner can
	 * not be used for the linear operator). Then the request is forwarded to
	 * the implementation of matrix based operators.
	 *
	 * \param[in]	L		linear operator
	 * \returns		bool	success flag
	 */
		bool init(SmartPtr<ILinearOperator<vector_type> > L)
		{
		// 	Remember operator
			m_spDefectOperator = L;
			if(m_bOtherApproxOperator) return true;

		//	cast to matrix based operator
			SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp =
					L.template cast_dynamic<MatrixOperator<matrix_type, vector_type> >();

		//	Check that matrix if of correct type
			if(pOp.invalid())
				UG_THROW(name() << "::init': Passed Operator is "
						"not based on matrix. This Preconditioner can only "
						"handle matrix-based operators.");

		//	forward request to matrix based implementation
			return init(pOp);

		}


	///	initializes the preconditioner for a matrix based operator
	/**
	 * This method initializes the preconditioner for matrix based operators.
	 * It performs some default checks and then forwards internally the
	 * initialization to the (virtual) 'preprocess'-method
	 *
	 * \param[in]	Op		matrix based operator
	 * \returns		bool	success flag
	 */
		bool init(SmartPtr<MatrixOperator<matrix_type, vector_type> > Op)
		{
		// 	Remember operator
			m_spApproxOperator = Op;
			m_spDefectOperator = Op;

		//	Check that matrix exists
			if(m_spApproxOperator.invalid())
				UG_THROW(name() << "::init': Passed Operator is invalid.");

		//	Preprocess
			if(!preprocess(m_spApproxOperator))
			{
				UG_LOG("ERROR in '"<<name()<<"::init': Preprocess failed.\n");
				return false;
			}

		//	Remember, that operator has been initialized
			m_bInit = true;

		//	we're done
			return true;
		}

	///	compute new correction c = B*d
	/**
	 * This method implements the virtual method of the ILinearIterator-interface.
	 * Basically, besides some common checks the request is forwarded to the
	 * (virtual) 'step'-method.
	 *
	 * \param[out]	c		correction
	 * \param[in]	d		defect
	 * \returns		bool	success flag
	 */
		virtual bool apply(vector_type& c, const vector_type& d)
		{
		//	Check that operator is initialized
			if(!m_bInit)
			{
				UG_LOG("ERROR in '"<<name()<<"::apply': Iterator not initialized.\n");
				return false;
			}

		//	Check parallel status
			#ifdef UG_PARALLEL
			if(!d.has_storage_type(PST_ADDITIVE))
				UG_THROW(name() << "::apply: Wrong parallel "
				               "storage format. Defect must be additive.");
			#endif

		//	Check sizes
			THROW_IF_NOT_EQUAL_4(c.size(), d.size(),
					m_spApproxOperator->num_rows(), m_spApproxOperator->num_cols());

		// 	apply iterator: c = B*d
			if(!step(m_spApproxOperator, c, d))
			{
				UG_LOG("ERROR in '"<<name()<<"::apply': Step Routine failed.\n");
				return false;
			}

		//	apply scaling
			const number kappa = damping()->damping(c, d, m_spApproxOperator);
			if(kappa != 1.0){
				c *= kappa;
			}

		//	Correction is always consistent
			#ifdef 	UG_PARALLEL
			if(!c.change_storage_type(PST_CONSISTENT))
				UG_THROW(name() << "::apply': Cannot change "
						"parallel storage type of correction to consistent.");
			#endif

		//	we're done
			return true;
		}

	///	compute new correction c = B*d and update defect d:= d - L*c
	/**
	 * This method implements the virtual method of the ILinearIterator-interface.
	 * Basically, the request is forwarded to the 'apply'-method and then the
	 * update of the defect is computed afterwards.
	 *
	 * \param[out]		c		correction
	 * \param[in, out]	d		defect on entry, updated defect on exit
	 * \returns			bool	success flag
	 */
		virtual bool apply_update_defect(vector_type& c, vector_type& d)
		{
		//	compute new correction
			if(!apply(c, d)) return false;

		// 	update defect d := d - A*c
			m_spDefectOperator->apply_sub(d, c);

		//	we're done
			return true;
		}

		virtual void set_approximation(SmartPtr<MatrixOperator<matrix_type,vector_type> > approx)
		{
			UG_COND_THROW(!approx.valid(), "");
			m_spApproxOperator = approx;
			if(!preprocess(m_spApproxOperator))
			{
				UG_THROW("ERROR in '"<<name()<<"::init': Preprocess failed.\n");
				return;
			}
			m_bOtherApproxOperator = true;
		}

	/// virtual destructor
		virtual ~IPreconditioner() {};

		///	underlying matrix based operator for calculation of defect
		SmartPtr<MatrixOperator<matrix_type, vector_type> > defect_operator()
		{
			return m_spDefectOperator;
		}

		///	underlying matrix based operator used for the preconditioner
		SmartPtr<MatrixOperator<matrix_type, vector_type> > approx_operator()
		{
			return m_spApproxOperator;
		}

	protected:
	///	underlying matrix based operator for calculation of defect
		SmartPtr<ILinearOperator<vector_type> > m_spDefectOperator;
	///	underlying matrix based operator used for the preconditioner
		SmartPtr<MatrixOperator<matrix_type, vector_type> > m_spApproxOperator;


	/// init flag indicating if init has been called
		bool m_bInit;

		bool m_bOtherApproxOperator;

};


}  // end namespace ug
#endif /* __H__LIB_ALGEBRA__OPERATOR__INTERFACE__PRECONDITIONER__ */
