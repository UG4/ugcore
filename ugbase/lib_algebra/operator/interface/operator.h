/*
 * operator.h
 *
 *  Created on: 22.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_ALGEBRA__OPERATOR__INTERFACE__OPERATOR__
#define __H__LIB_ALGEBRA__OPERATOR__INTERFACE__OPERATOR__

namespace ug{

///////////////////////////////////////////////////////////////////////////////
// (Nonlinear-) Operator
///////////////////////////////////////////////////////////////////////////////

/// describes a mapping X->Y
/**
 * This class is the base class for all mappings between two spaces. The domain
 * space and the codomain space are passed as template parameters. In particular,
 * the mapping can be nonlinear. For linear (or linearized) mappings see
 * ILinearizedOperator. The basic usage of this class is to provide the
 * computation of:
 *
 * 		d := N(u),
 *
 * where d is from the codomain space, u a function of the domain space and
 * N() is a (nonlinear-) mapping.
 *
 * This application is splitted into three methods, that have to be called
 * in the correct order:
 *
 * 1. init(): This method initializes the operator for application. It has
 * 			  to be called once before one of the other two methods can be
 * 			  invoked. There is no need to call this method more than once, but
 * 			  sometimes - due to parameter change - this is desirable and
 * 			  can be done.
 *
 * 2. prepare(): This method is used to prepare the in- and output vector used
 * 				 later in apply. It can be called, after the init method has
 * 				 been called at least once. The prepare method is e.g. used
 * 				 to set dirichlet values.
 *
 * 3. apply(): This method can be called when the operator has been initialized
 * 			   by a call of init and with two functions, that have been prepare
 * 			   using the prepare method. It maps the function from the domain
 * 			   space to the range space.
 *
 * This splitting has been made, since initialization and preparation may be
 * computationally expansive. Thus, the user of this class has the choice
 * when to call this initialization/preparation. E.g. when the operator is
 * applied several times on the same vectors, those have only to be prepared
 * once and the init of the operator is only needed once.
 *
 * \tparam	X 	Domain space function
 * \tparam	Y	Range space function
 */
template <typename X, typename Y = X>
class IOperator
{
	public:
	///	Domain space
		typedef X domain_function_type;

	///	Range space
		typedef Y codomain_function_type;

	public:
	/// initializes the operator
	/**
	 * This method initializes the operator. It must be called before any of
	 * the other methods are called.
	 *
	 * \returns 	bool	success flag
	 */
		virtual void init() = 0;

	/// prepares domain and codomain functions for application
	/**
	 * This method prepares the in- and output functions for the application
	 * and has to be called once before the apply method can be invoked with
	 * the functions used here.
	 *
	 * \param[in]	u		domain function
	 * \param[out]	d		codomain function
	 * \returns 	bool	flag if preparation successful
	 */
		virtual void prepare(Y& d, X& u) = 0;

	///	computes the nonlinear mapping d := N(u)
	/**
	 * This method maps a function from the domain space to the range space.
	 * Note, that is must be called with functions, that have previously been
	 * prepared using the 'prepare'-method and that the operator must have been
	 * initialized using the 'init'-method
	 *
	 * \param[in]	u		domain function
	 * \param[out]	d		codomain function
	 * \returns 	bool	flag if application successful
	 */
		virtual void apply(Y& d, const X& u) = 0;

	///	virtual destructor
		virtual ~IOperator() {};
};

///////////////////////////////////////////////////////////////////////////////
// Linearized Operator
///////////////////////////////////////////////////////////////////////////////

/// describes a linear mapping X->Y
/**
 * This class is the base class for all linear mappings between two spaces.
 * The domain space and the codomain space are passed as template parameters.
 * The mapping must be linear. For nonlinear mappings see IOperator. The basic
 * usage of this class is to provide the computation of:
 *
 * 		f := L*u,    (resp.  d := J(u) * c in iterative schemes),
 *
 * where f (resp. d) is from the codomain space, u (resp. c) a function of
 * the domain space and L is a linear mapping (resp. J(u) a linearized mapping)
 *
 * This application is splitted into two steps, that have to be called in the
 * correct order:
 *
 * 1. init() or init(u):
 * 		Theses methods initialize the operator for application. One of these
 * 		methods has to be called once before one of the other two methods can be
 * 		invoked. There is no need to init the operator more than once, but
 * 		sometimes - due to parameter change - this is desirable and can be done.
 *
 * 2. apply() or apply_sub():
 * 		These methods can be called when the operator has been initialized
 * 		by a call of init. These function perform the linear mapping, where in
 * 		the case of apply_sub() the result is subtracted from the input function.
 *
 * This splitting has been made, since initialization may be computationally
 * expansive. Thus, the user of this class has the choice when to call this
 * initialization. E.g. when the operator is applied several times, the init
 * of the operator is only needed once.
 *
 * \tparam	X 	Domain space function
 * \tparam	Y	Range space function
 */
template <typename X, typename Y = X>
class ILinearOperator : public IOperator<X,Y>
{
	public:
	///	Domain space
		typedef X domain_function_type;

	///	Range space
		typedef Y codomain_function_type;

	public:
	///	init operator depending on a function u
	/**
	 * This method initializes the operator. Once initialized the 'apply'-method
	 * can be called. The function u is passed here, since the linear operator
	 * may be the linearization of some non-linear operator. Thus, the operator
	 * depends on the linearization point.
	 * If the operator is not a linearization, this method can be implemented
	 * by simply calling init() and forgetting about the linearization point.
	 *
	 * \param[in]	u		function (linearization point)
	 * \returns 	bool	success flag
	 */
		virtual void init(const X& u) = 0;

	///	init operator
	/**
	 * This method initializes the operator. Once initialized the 'apply'-method
	 * can be called.
	 * If the operator is a linearization this function returns false.
	 *
	 * \returns 	bool	success flag
	 */
		virtual void init() = 0;

	///	default implementation for IOperator interface
		virtual void prepare(Y& d, X& u) {}

	// 	applies the operator
	/**
	 * This method applies the operator, i.e. f = L*u (or d = J(u)*c in
	 * iterative schemes). Note, that the operator must have been initialized
	 * once before this method can be used.
	 *
	 * \param[in]	u		domain function
	 * \param[out]	f		codomain function
	 * \returns		bool	success flag
	 */
		virtual void apply(Y& f, const X& u) = 0;

	// 	applies the operator and subtracts the result from the input
	/**
	 * This method applies the operator and subracts the result from the input
	 * codomain function, i.e. f -= L*u (or d -= J(u)*c in iterative schemes).
	 * Note, that the operator must have been initialized once before this
	 * method can be used.
	 *
	 * \param[in]		u		domain function
	 * \param[in,out]	f		codomain function
	 * \returns			bool	success flag
	 */
		virtual void apply_sub(Y& f, const X& u) = 0;

	/// virtual	destructor
		virtual ~ILinearOperator() {};
};


///////////////////////////////////////////////////////////////////////////////
// Matrix based linear operator
///////////////////////////////////////////////////////////////////////////////

template <typename M, typename X, typename Y = X>
class MatrixOperator :	public virtual ILinearOperator<X,Y>,
						public M
{
	public:
	// 	Domain space
		typedef X domain_function_type;

	// 	Range space
		typedef Y codomain_function_type;

	// 	Matrix type
		typedef M matrix_type;

	public:
	// 	Init Operator J(u)
		virtual void init(const X& u) {}

	// 	Init Operator L
		virtual void init() {}

	// 	Apply Operator f = L*u (e.g. d = J(u)*c in iterative scheme)
		virtual void apply(Y& f, const X& u) {matrix_type::apply(f,u);}

	// 	Apply Operator, i.e. f = f - L*u;
		virtual void apply_sub(Y& f, const X& u) {matrix_type::matmul_minus(f,u);}

	// 	Access to matrix
		virtual M& get_matrix() {return *this;};
};

} // end namespace ug

#endif /* __H__LIB_ALGEBRA__OPERATOR__INTERFACE__OPERATOR__ */
