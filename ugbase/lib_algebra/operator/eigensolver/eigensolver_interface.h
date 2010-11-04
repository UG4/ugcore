/*
 * eigensolver_interface.h
 *
 *  Created on: 01.11.2010
 *      Author: Martin Rupp
 */

#ifndef __H__LIB_ALGEBRA__OPERATOR__EIGENSOLVER_INTERFACE__
#define __H__LIB_ALGEBRA__OPERATOR__EIGENSOLVER_INTERFACE__


namespace ug{
///////////////////////////////////////////////////////////
// Operator
///////////////////////////////////////////////////////////

// describes a mapping X->Y
template <typename X, typename Y>
class IEigensolver
{
	public:
	// 	Domain space
		typedef X domain_function_type;

	// 	Range space
		typedef Y codomain_function_type;

	public:
	// 	Init Operator
		virtual bool init() = 0;

	// 	Prepare out function
		virtual bool prepare(Y& d, X& u) = 0;

	// 	Apply Operator, i.e. d := N(u);
		virtual bool apply(Y& d, const X& u) = 0;

	// 	Destructor
		virtual ~IOperator() {};
};

} // end namespace ug

#endif /* __H__LIB_ALGEBRA__OPERATOR__EIGENSOLVER_INTERFACE__ */
