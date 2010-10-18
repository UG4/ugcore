/*
 * user_data.h
 *
 *  Created on: 13.10.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__USER_DATA__
#define __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__USER_DATA__

#include <boost/function.hpp>

namespace ug {


template <int dim>
class IUserNumberProvider
{
	public:
	//	Functor Type
		typedef boost::function<void (number& n, const MathVector<dim>& x, number& time)> functor_type;

	/// provides the functor
		virtual functor_type get_functor() const = 0;

	///	virtual destructor
		virtual ~IUserNumberProvider(){}
};



template <int dim>
class IUserVectorProvider
{
	public:
	//	Functor Type
		typedef boost::function<void (MathVector<dim>& v, const MathVector<dim>& x, number& time)> functor_type;

	/// provides the functor
		virtual functor_type get_functor() const = 0;

	///	virtual destructor
		virtual ~IUserVectorProvider(){}
};


template <int dim>
class IUserMatrixProvider
{
	public:
	//	Functor Type
		typedef boost::function<void (MathMatrix<dim,dim>& D, const MathVector<dim>& x, number& time)> functor_type;

	/// provides the functor
		virtual functor_type get_functor() const = 0;

	///	virtual destructor
		virtual ~IUserMatrixProvider(){}
};



template <int dim>
class IBoundaryNumberProvider
{
	public:
	//	Functor Type
		typedef boost::function<bool (number& n, const MathVector<dim>& x, number& time)> functor_type;

	/// provides the functor
		virtual functor_type get_functor() const = 0;

	///	virtual destructor
		virtual ~IBoundaryNumberProvider(){}
};



}

#endif /* __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__USER_DATA__ */
