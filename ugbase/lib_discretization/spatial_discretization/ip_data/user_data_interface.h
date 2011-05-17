/*
 * user_data_interface.h
 *
 *  Created on: 13.10.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__IP_DATA__USER_DATA_INTERFACE__
#define __H__UG__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__IP_DATA__USER_DATA_INTERFACE__

#include "common/common.h"
#include "common/math/ugmath.h"
#include <boost/function.hpp>

#include "ip_data.h"
#include "data_linker.h"

namespace ug {

/**
 * \brief User Data
 *
 * User Data that can be used in assembling routines.
 *
 * \defgroup lib_disc_user_data User Data
 * \ingroup lib_discretization
 */

/// \addtogroup lib_disc_user_data
/// @{

/// User Data Interface
/**
 * Specifies and provides the User Data Functor
 * \tparam		TData	provided data type
 * \tparam		dim		world dimension
 */
template<typename TData, int dim>
class IUserData
	: public IPData<TData, dim>
{
	public:
	///	Functor type for evaluation function
		typedef boost::function<void (TData& value,
		                              const MathVector<dim>& x,
		                              number time)>
		functor_type;

	/// provides the functor
		virtual functor_type get_functor() const = 0;
};

/// User Function Interface
/**
 * Specifies and provides the User Function Functor
 * \tparam		TData	provided data type
 * \tparam		TPos	Position type of world position
 */
template<typename TData, int dim, typename TDataIn>
class IUserFunction
	: public DataLinkerEqualData<TData, dim, TDataIn>
{
	public:
	///	virtual operator
		virtual void operator() (TData& out, int numArgs, ...) = 0;

	/// virtual destructor
		virtual ~IUserFunction() {}
};

/// scalar boundary user data
/**
 * Provides a functor to evaluate scalar user data at boundary. The difference
 * to IUserNumber is, that the functor returns true iff position is dirichlet
 * and false else.
 *
 * \tparam 	dim		World dimension
 */
template <typename TData, int dim>
class IBoundaryData
{
	public:
	///	Functor Type
		typedef boost::function<bool (TData& value,
		                              const MathVector<dim>& x,
		                              number time)>
		functor_type;

	/// provides the functor
		virtual functor_type get_functor() const = 0;

	///	virtual destructor
		virtual ~IBoundaryData(){}
};

/// @}

} /// end namespace ug

#endif /* __H__UG__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__IP_DATA__USER_DATA_INTERFACE__ */
