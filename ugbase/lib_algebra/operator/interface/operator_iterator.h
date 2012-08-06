/*
 * operator_iterator.h
 *
 *  Created on: 22.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_ALGEBRA__OPERATOR__INTERFACE__OPERATOR_ITERATOR__
#define __H__LIB_ALGEBRA__OPERATOR__INTERFACE__OPERATOR_ITERATOR__

#include "lib_algebra/operator/interface/operator.h"
#include "lib_algebra/operator/damping.h"
#include "lib_algebra/operator/debug_writer.h"

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
template <typename X, typename Y = X>
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
		virtual bool init(SmartPtr<ILinearOperator<Y,X> > J, const Y& u) = 0;

	///	initialize for linear operator L
	/**
	 * This method passes the operator L that used as underlying by this
	 * operator. In addition some preparation step can be made.
	 *
	 * \param[in]	L		linear operator to use as underlying
	 * \returns		bool	success flag
	 */
		virtual bool init(SmartPtr<ILinearOperator<Y,X> > L) = 0;

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

	///	sets a scaling for the correction
	/**
	 * Sets a scaling for the correction, i.e., once the correction has been
	 * computed, c = B*d, the correction is scaled by a factor, c := s*c, where
	 * s is provided by the passed scaling class. Note, that the scaling factor
	 * may depend on the defect and correction. The internal update of the defect
	 * as done in apply_update_defect must be performed with respect to the
	 * scaled correction.
	 */
		void set_damp(SmartPtr<IDamping<X,Y> > spScaling) {
			m_spDamping = spScaling;
		}

	///	sets the damping to a constant factor
		void set_damp(number factor) {
			m_spDamping = new ConstantDamping<X,Y>(factor);
		}

	///	returns the scaling
		SmartPtr<IDamping<X,Y> > damping() {return m_spDamping;}

	///	clone
		virtual SmartPtr<ILinearIterator<X,Y> > clone() = 0;

	/// virtual destructor
		virtual ~ILinearIterator() {};

	///	constructor
		ILinearIterator() {set_damp(1.0);};

	protected:
	///	the scaling
		SmartPtr<IDamping<X,Y> > m_spDamping;
};

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

	protected:
		using ILinearIterator<vector_type>::damping;
		using DebugWritingObject<TAlgebra>::set_debug;
		using DebugWritingObject<TAlgebra>::debug_writer;
		using DebugWritingObject<TAlgebra>::write_debug;

	public:
	///	default constructor
		IPreconditioner() :
			m_spOperator(NULL), m_bInit(false)
		{};

	///	constructor setting debug writer
		IPreconditioner(SmartPtr<IDebugWriter<algebra_type> > spDebugWriter) :
			DebugWritingObject<TAlgebra>(spDebugWriter),
			m_spOperator(NULL), m_bInit(false)
		{};

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
		virtual bool preprocess(MatrixOperator<matrix_type, vector_type>& mat) = 0;

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
		virtual bool step(MatrixOperator<matrix_type, vector_type>& mat, vector_type& c, const vector_type& d)  = 0;

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
		virtual bool init(SmartPtr<ILinearOperator<vector_type> > L)
		{
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
			m_spOperator = Op;

		//	Check that matrix exists
			if(m_spOperator.invalid())
				UG_THROW(name() << "::init': Passed Operator is invalid.");

		//	Preprocess
			if(!preprocess(*m_spOperator))
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
			if(d.size() != m_spOperator->num_rows())
				UG_THROW("Vector [size= "<<d.size()<<"] and Row [size= "
				               <<m_spOperator->num_rows()<<"] sizes have to match!");
			if(c.size() != m_spOperator->num_cols())
				UG_THROW("Vector [size= "<<c.size()<<"] and Column [size= "
				               <<m_spOperator->num_cols()<<"] sizes have to match!");
			if(d.size() != c.size())
				UG_THROW("Vector [d size= "<<d.size()<<", c size = "
				               << c.size() << "] sizes have to match!");

		// 	apply iterator: c = B*d
			if(!step(*m_spOperator, c, d))
			{
				UG_LOG("ERROR in '"<<name()<<"::apply': Step Routine failed.\n");
				return false;
			}

		//	apply scaling
			const number kappa = damping()->damping(c, d, m_spOperator);
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
			if(!m_spOperator->matmul_minus(d, c))
			{
				UG_LOG("ERROR in '"<<name()<<"::apply_update_defect': "
						"Cannot execute matmul_minus to compute d:=d-A*c.\n");
				return false;
			}

		//	we're done
			return true;
		}

	/// virtual destructor
		virtual ~IPreconditioner() {};

	protected:
	///	underlying matrix based operator
		SmartPtr<MatrixOperator<matrix_type, vector_type> > m_spOperator;

	/// init flag indicating if init has been called
		bool m_bInit;
};

} // end namespace ug

#endif /* __H__LIB_ALGEBRA__OPERATOR__INTERFACE__OPERATOR_ITERATOR__ */

