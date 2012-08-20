/*
 * operator_inverse.h
 *
 *  Created on: 22.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_ALGEBRA__OPERATOR__INTERFACE__OPERATOR_INVERSE__
#define __H__LIB_ALGEBRA__OPERATOR__INTERFACE__OPERATOR_INVERSE__

#include "common/util/smart_pointer.h"
#include "common/profiler/profiler.h"

#include "lib_algebra/operator/interface/operator.h"
#include "lib_algebra/operator/interface/operator_iterator.h"

#include "lib_algebra/operator/debug_writer.h"
#include "lib_algebra/operator/convergence_check.h"

#define PROFILE_LS
#ifdef PROFILE_LS
	#define LS_PROFILE_FUNC()		PROFILE_FUNC()
	#define LS_PROFILE_BEGIN(name)	PROFILE_BEGIN(name)
	#define LS_PROFILE_END()		PROFILE_END()
#else
	#define LS_PROFILE_FUNC()
	#define LS_PROFILE_BEGIN(name)
	#define LS_PROFILE_END()
#endif

namespace ug{

///////////////////////////////////////////////////////////////////////////////
// Inverse of a (Nonlinear-) Operator
///////////////////////////////////////////////////////////////////////////////

// \todo: The prepare methods seems obsolet, since it could be handled in the
//		  apply method. It depends on how expensive the prepare is, and how
//		  often we would apply the operator to an already apply one (?)
/// describes an inverse mapping X->Y
/**
 * This class is the base class for the inversion of the operator given in form
 * of the IOperator interface class. Given a operator N(u), the basic usage
 * of this class is to invert this operator, i.e. to compute
 *
 *     N(u) = 0, i.e. u := N^{-1}(0).
 *
 * This application has been split up into three steps:
 *
 * 1. init(N): This method initializes the inverse operator. The inverse operator
 * 			   is initialized the way that, its application will be the inverse
 * 			   application of the operator N passed in by this function. The
 * 			   prepare method can only be called, when this method has been
 * 			   called once.
 *
 * 2. prepare(u): This method prepares the function u, before it can be used in
 * 				  the apply method. Here, typically dirichlet values are set.
 *
 * 3. apply(u):	This method performs the inversion. Before this method is called
 * 				the init and prepare methods have to be called.
 *
 * This splitting has been made, since initialization and preparation may be
 * computationally expansive. Thus, the user of this class has the choice
 * when to call this initialization/preparation. E.g. when the operator is
 * applied several times on the same vectors, those have only to be prepared
 * once and the init of the operator is only needed once.
 */
template <typename X, typename Y = X>
class IOperatorInverse
{
	public:
	///	Domain space
		typedef X domain_function_type;

	/// Range space
		typedef Y codomain_function_type;

	public:
	/// initializes the operator for the inversion of the operator N: Y -> X
	/**
	 * This method sets the (nonlinear) operator that should be inverted by
	 * this inverse operator. In addition preparations can be made to
	 * facilitate the application of the inverse.
	 *
	 * \param[in]	N		operator that is to be inverted
	 * \returns 	bool	success flag
	 */
		virtual bool init(SmartPtr<IOperator<Y,X> > N) = 0;

	/// prepares the function u for application
	/**
	 * This method prepares the function u before it can be used to find the
	 * solution of N(u) = 0. Typically, dirichlet values are set here.
	 *
	 * \param[in]	u		domain function
	 * \returns		bool	success flag
	 */
		virtual bool prepare(X& u) = 0;

	/// apply the operator, i.e. u = N^{-1}(0)
	/**
	 * This method inverts the operator N and returns the solution u = N^{-1}(0).
	 *
	 * \param[in,out]	u		domain function with solution at output
	 * \returns			bool	success flag
	 */
		virtual bool apply(X& u) = 0;

	/// virtual destructor
		virtual ~IOperatorInverse() {};
};


///////////////////////////////////////////////////////////////////////////////
// Inverse of a Linear Operator
///////////////////////////////////////////////////////////////////////////////

/// describes an inverse linear mapping X->Y
/**
 * This class is the base class for the inversion of linear operator given in
 * form of the ILinearOperator interface class. Given a operator L, the basic
 * usage of this class is to invert this operator, i.e. to compute the solution
 * u of
 *
 * 		L*u = f     i.e. u := L^{-1} f
 *
 * This application has been split up into three steps:
 *
 * 1. init():  This method initializes the inverse operator. The inverse operator
 * 			   is initialized the way that, its application will be the inverse
 * 			   application of the operator L passed in by this function. The
 * 			   prepare method can only be called, when this method has been
 * 			   called once.
 *
 * 3. apply():	This method performs the inversion. Before this method is called
 * 				the init and prepare methods have to be called.
 *
 * This splitting has been made, since initialization and preparation may be
 * computationally expansive. Thus, the user of this class has the choice
 * when to call this initialization/preparation. E.g. when the operator is
 * applied several times on the same vectors, those have only to be prepared
 * once and the init of the operator is only needed once.
 *
 * \tparam	X		domain space
 * \tparam	Y		range space
 */
template <typename X, typename Y = X>
class ILinearOperatorInverse
{
	public:
	///	Domain space
		typedef X domain_function_type;

	///	Range space
		typedef Y codomain_function_type;

	public:
	///	constructor setting convergence check to (100, 1e-12, 1e-12, true)
		ILinearOperatorInverse()
			: m_spLinearOperator(NULL),
			  m_spConvCheck(new StdConvCheck<X>(100, 1e-12, 1e-12, true))
		{}

	///	Default constructor
		ILinearOperatorInverse(SmartPtr<IConvergenceCheck<X> > spConvCheck)
			: m_spLinearOperator(NULL),
			  m_spConvCheck(spConvCheck)
		{}

	/// virtual destructor
		virtual ~ILinearOperatorInverse() {};

	///	returns the name of the operator inverse
	/**
	 * This method returns the name of the inverse operator. This function is
	 * typically needed, when the inverse operator is used inside of another and
	 * some debug output should be printed
	 *
	 * \returns 	const char* 	name of inverse operator
	 */
		virtual const char* name() const = 0;

	/// initializes for the inverse for a linear operator
	/**
	 * This method passes the operator L that is inverted by this operator. In
	 * addition some preparation step can be made.
	 *
	 * \param[in]	L		linear operator to invert
	 * \returns		bool	success flag
	 */
		virtual bool init(SmartPtr<ILinearOperator<Y,X> > L)
		{
		//	remember operator
			m_spLinearOperator = L;
			return true;
		}

	/// initializes for the inverse for a linearized operator at linearization point u
	/**
	 * This method passes the linear operator J(u) that should be inverted by
	 * this operator. As second argument the linearization point is passed.
	 * This is needed e.g. for the geometric multigrid method, that inverts
	 * a linearized operator based on coarser grid operators, that have to be
	 * initialized based on the linearization point.
	 *
	 * \param[in]	J		linearized operator to invert
	 * \param[in]	u		linearization point
	 * \returns		bool	success flag
	 */
		virtual bool init(SmartPtr<ILinearOperator<Y,X> > J, const Y& u)
		{
		//	remember operator
			m_spLinearOperator = J;
			return true;
		}

	///	applies inverse operator, i.e. returns u = A^{-1} f
	/**
	 * This method applies the inverse operator, i.e. u = A^{-1} f. The
	 * domain function f remains unchanged.
	 * Note, that this method can always be implemented by creating a copy of
	 * f and calling apply_return_defect with this copy.
	 *
	 * \param[in]	f		right-hand side
	 * \param[out]	u		solution
	 * \returns		bool	success flag
	 */
		virtual bool apply(Y& u, const X& f) = 0;

	///	applies inverse operator, i.e. returns u = A^{-1} f and returns defect d := f - A*u
	/**
	 * This method applies the inverse operator, i.e. u = A^{-1} f. The
	 * domain function f is changed in the way, that the defect d := f - A*u
	 * is returned in the function. This is always useful, when the inverting
	 * algorithm can (or must) update the defect during computation (this is
	 * e.g. the case for the geometric multigrid method).
	 * Note, that this method can always be implemented by calling apply and
	 * then computing d := f - A*u.
	 *
	 * \param[in,out]	f		right-hand side
	 * \param[out]		u		solution
	 * \returns			bool	success flag
	 */
		virtual bool apply_return_defect(Y& u, X& f) = 0;

	///	returns the convergence check
		ConstSmartPtr<IConvergenceCheck<X> > convergence_check() const {return m_spConvCheck;}

	///	returns the convergence check
		SmartPtr<IConvergenceCheck<X> > convergence_check() {return m_spConvCheck;}

	/// returns the current defect
		number defect() const {return convergence_check()->defect();}

	/// returns the current number of steps
		int step() const {return convergence_check()->step();}

	/// returns the current relative reduction
		number reduction() const {return convergence_check()->reduction();}

	///	returns the standard offset for output
		virtual int standard_offset() const {return 3;}

	///	set the convergence check
		void set_convergence_check(SmartPtr<IConvergenceCheck<X> > spConvCheck)
		{
			m_spConvCheck = spConvCheck;
			m_spConvCheck->set_offset(standard_offset());
		};

	///	returns the current Operator this Inverse Operator is initialized for
		SmartPtr<ILinearOperator<Y,X> > linear_operator()
		{
			if(m_spLinearOperator.invalid())
				UG_THROW(name() << ": Linear Operator that should be "
				               	 "inverted has not been set.");

			return m_spLinearOperator;
		}

	protected:
	/// Operator that is inverted by this Inverse Operator
		SmartPtr<ILinearOperator<Y,X> > m_spLinearOperator;

	///	smart pointer holding the convergence check
		SmartPtr<IConvergenceCheck<X> > m_spConvCheck;
};


///////////////////////////////////////////////////////////////////////////////
// Inverse of a Linear Operator using a ILinearIterator as preconditioner
///////////////////////////////////////////////////////////////////////////////

/// describes an inverse linear mapping X->X
/**
 * This a useful derived class from ILinearOperatorInverse, that uses a
 * ILinearIterator in order to precondition the solution process. This is
 * used e.g. in LinearSolver, CG and BiCGStab.
 *
 * \tparam	X		domain and range space
 */
template <typename X>
class IPreconditionedLinearOperatorInverse
	: public ILinearOperatorInverse<X>,
	  public VectorDebugWritingObject<X>
{
	public:
	///	Domain space
		typedef X domain_function_type;

	///	Range space
		typedef X codomain_function_type;

	///	Base class
		typedef ILinearOperatorInverse<X,X> base_type;

	protected:
		using base_type::name;
		using base_type::linear_operator;
		using base_type::apply_return_defect;
		using VectorDebugWritingObject<X>::write_debug;

	public:
	///	Empty constructor
		IPreconditionedLinearOperatorInverse()
			: m_bRecompute(false), m_spPrecond(NULL)
		{}

	///	constructor setting the preconditioner
		IPreconditionedLinearOperatorInverse(SmartPtr<ILinearIterator<X,X> > spPrecond)
			: m_bRecompute(false), m_spPrecond(spPrecond)
		{}

	///	constructor setting the preconditioner
		IPreconditionedLinearOperatorInverse(SmartPtr<ILinearIterator<X,X> > spPrecond,
		                                     SmartPtr<IConvergenceCheck<X> > spConvCheck)
			: 	base_type(spConvCheck),
				m_bRecompute(false), m_spPrecond(spPrecond)
		{}

	///	sets the preconditioner
		void set_preconditioner(SmartPtr<ILinearIterator<X, X> > spPrecond)
		{
			m_spPrecond = spPrecond;
		}

	///	returns the preconditioner
		SmartPtr<ILinearIterator<X, X> > preconditioner(){return m_spPrecond;}

	///	initializes the solver for an operator
		virtual bool init(SmartPtr<ILinearOperator<X,X> > J, const X& u)
		{
			if(!base_type::init(J, u)) return false;

			LS_PROFILE_BEGIN(LS_InitPrecond);
			if(m_spPrecond.valid())
				if(!m_spPrecond->init(J, u))
					UG_THROW(name() << "::init: Cannot init Preconditioner "
													"Operator for Operator J.");
			LS_PROFILE_END();

			return true;
		}

	///	initializes the solver for an operator
		virtual bool init(SmartPtr<ILinearOperator<X,X> > L)
		{
			if(!base_type::init(L)) return false;

			LS_PROFILE_BEGIN(LS_InitPrecond);
			if(m_spPrecond.valid())
				if(!m_spPrecond->init(L))
					UG_THROW(name() <<"::prepare: Cannot init Preconditioner "
														"Operator for Operator L.");
			LS_PROFILE_END();

			return true;
		}

		virtual bool apply(X& x, const X& b)
		{
		//	copy defect
			X bTmp; bTmp.resize(b.size()); bTmp = b;

		//	solve on copy of defect
			bool bRes = apply_return_defect(x, bTmp);

		//	write updated defect
			write_debug(bTmp, "LS_UpdatedDefectEnd.vec");

		//	compute defect again, for debug purpose
			if(m_bRecompute)
			{
			//	recompute defect
				bTmp = b; linear_operator()->apply_sub(bTmp, x);
				number norm = bTmp.norm();

			//	print norm of recomputed defect
				UG_LOG("%%%% DEBUG "<<name()<<": (Re)computed defect has norm: "
				       <<norm<<"\n");

			//	write true end defect
				write_debug(bTmp, "LS_TrueDefectEnd.vec");
			}

		//	return
			return bRes;
		}

	///	for debug: computes norm again after whole calculation of apply
		void set_compute_fresh_defect_when_finished(bool bRecompute)
		{
			m_bRecompute = bRecompute;
		}

	protected:
	///	flag if fresh defect should be computed when finish for debug purpose
		bool m_bRecompute;

	///	Iterator used in the iterative scheme to compute the correction and update the defect
		SmartPtr<ILinearIterator<X,X> > m_spPrecond;

};

///////////////////////////////////////////////////////////////////////////////
// Inverse of a Matrix-based Linear Operator
///////////////////////////////////////////////////////////////////////////////

/// describes an inverse linear mapping X->Y based on a matrix
/**
 * This class is the base class for the inversion of linear matrix-based operator
 * given in form of the IMatrixOperator interface class. Given a operator L,
 * the basic usage of this class is to invert this operator, i.e. to compute
 * the solution u of
 *
 * 		L*u = f     i.e. u := L^{-1} f
 *
 *
 * \tparam	X		domain space (i.e. a vector corresponding to the matrix)
 * \tparam	Y		range space (i.e. a vector corresponding to the matrix)
 * \tparam	M		matrix type used to represent linear mapping
 */
template <typename M, typename X, typename Y = X>
class IMatrixOperatorInverse
	: public virtual ILinearOperatorInverse<X,Y>
{
	public:
	///	Domain space
		typedef X domain_function_type;

	///	Range space
		typedef Y codomain_function_type;

	///	Matrix type
		typedef M matrix_type;

	public:
	///	initializes this inverse operator for a matrix-based operator
	/**
	 * This method passes the operator A that is inverted by this operator. In
	 * addition some preparation step can be made.
	 *
	 * \param[in]	A		linear matrix-basewd operator to invert
	 * \returns		bool	success flag
	 */
		virtual bool init(SmartPtr<MatrixOperator<M,Y,X> > A) = 0;

	/// applies the inverse operator, i.e. returns u = A^{-1} * f
	/**
	 * This method applies the inverse operator.
	 *
	 * \param[out]	u		solution
	 * \param[in]	f		right-hand side
	 * \returns		bool	success flag
	 */
		virtual bool apply(Y& u, const X& f) = 0;

	/// applies the inverse operator and updates the defect
	/**
	 * This method applies the inverse operator and updates the defect, i.e.
	 * returns u = A^{-1} * f and in f the last defect d:= f - A*u is returned.
	 *
	 * \param[out]	u		solution
	 * \param[in]	f		right-hand side on entry, defect on exit
	 * \returns		bool	success flag
	 */
		virtual bool apply_return_defect(Y& u, X& f) = 0;

	/// virtual destructor
		virtual ~IMatrixOperatorInverse() {};

	public:
	///	initializes this inverse operator for a linear operator
	/**
	 * This method implements the ILinearOperatorInverse interface method.
	 * Basically, the request is forwarded to the matrix-based init method,
	 * if the the operator is matrix-based. If the operator is not matrix-based
	 * this inverse can not be used and false is returned
	 *
	 * \param[in]	A		linear matrix-based operator to invert
	 * \param[in]	u		linearization point
	 * \returns		bool	success flag
	 */
		virtual bool init(SmartPtr<ILinearOperator<Y,X> > A, const Y&u)
		{
		//	forget about u and forward request.
			return init(A);
		}

	///	initializes this inverse operator for a linear operator
	/**
	 * This method implements the ILinearOperatorInverse interface method.
	 * Basically, the request is forwarded to the matrix-based init method,
	 * if the the operator is matrix-based. If the operator is not matrix-based
	 * this inverse can not be used and false is returned
	 *
	 * \param[in]	A		linear matrix-based operator to invert
	 * \returns		bool	success flag
	 */
		virtual bool init(SmartPtr<ILinearOperator<Y,X> > A)
		{
		//	cast operator
			SmartPtr<MatrixOperator<M,Y,X> > op =
									A.template cast_dynamic<MatrixOperator<M,Y,X> >();

		//	check if correct types are present
			if(op.invalid())
				UG_THROW("IMatrixOperatorInverse::init:"
						" Passed operator is not matrix-based.");

		//	forward request
			return init(op);
		}
};


} // end namespace ug

#endif /* __H__LIB_ALGEBRA__OPERATOR__INTERFACE__OPERATOR_INVERSE__ */
