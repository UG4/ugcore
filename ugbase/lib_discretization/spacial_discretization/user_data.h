/*
 * user_data.h
 *
 *  Created on: 13.10.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__USER_DATA__
#define __H__UG__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__USER_DATA__

#include "common/common.h"
#include "common/math/ugmath.h"
#include <boost/function.hpp>

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

/// scalar user data
/**
 * Provides a functor to evaluate scalar user data.
 * \tparam 	dim		World dimension
 */
template <int dim>
class IUserNumber
{
	public:
	///	Functor Type for evaluation function
		typedef boost::function<void (number& value,
		                              const MathVector<dim>& x,
		                              number& time)>
		functor_type;

	/// provides the functor
		virtual functor_type get_functor() const = 0;

	///	virtual destructor
		virtual ~IUserNumber(){}
};

/// vector user data
/**
 * Provides a functor to evaluate vector user data.
 * \tparam 	dim		World dimension
 */
template <int dim>
class IUserVector
{
	public:
	///	Functor Type
		typedef boost::function<void (MathVector<dim>& v,
		                              const MathVector<dim>& x,
		                              number& time)>
		functor_type;

	/// provides the functor
		virtual functor_type get_functor() const = 0;

	///	virtual destructor
		virtual ~IUserVector(){}
};

/// matrix user data
/**
 * Provides a functor to evaluate matrix user data.
 * \tparam 	dim		World dimension
 */
template <int dim>
class IUserMatrix
{
	public:
	///	Functor Type
		typedef boost::function<void (MathMatrix<dim,dim>& D,
		                              const MathVector<dim>& x,
		                              number& time)>
		functor_type;

	/// provides the functor
		virtual functor_type get_functor() const = 0;

	///	virtual destructor
		virtual ~IUserMatrix(){}
};

/// scalar boundary user data
/**
 * Provides a functor to evaluate scalar user data at boundary. The difference
 * to IUserNumber is, that the functor returns true iff position is dirichlet
 * and false else.
 *
 * \tparam 	dim		World dimension
 */
template <int dim>
class IBoundaryNumberProvider
{
	public:
	///	Functor Type
		typedef boost::function<bool (number& value,
		                              const MathVector<dim>& x,
		                              number& time)>
		functor_type;

	/// provides the functor
		virtual functor_type get_functor() const = 0;

	///	virtual destructor
		virtual ~IBoundaryNumberProvider(){}
};


///////////////////////////////////////
///////////////////////////////////////
// Const Data
///////////////////////////////////////
///////////////////////////////////////

/// constant scalar user data
template <int dim>
class ConstUserNumber
	: public IUserNumber<dim>
{
	public:
	///	Functor Type
		typedef typename IUserNumber<dim>::functor_type functor_type;

	///	return functor
		virtual functor_type get_functor() const {return boost::ref(*this);}

	public:
	///	creates empty user number
		ConstUserNumber() {m_Number = 0.0;}

	///	set constant value
		void set(number val)
		{
			m_Number = val;
		}

	///	print current setting
		void print() const
		{
			UG_LOG("ConstUserNumber:" << m_Number << "\n");
		}

	///	evaluate
		void operator() (number& c, const MathVector<dim>& x,
		                 number time = 0.0) const
		{
			c = m_Number;
		}

	protected:
		number m_Number;
};

/// constant vector user data
template <int dim>
class ConstUserVector
	: public IUserVector<dim>
{
	public:
	///	Functor Type
		typedef typename IUserVector<dim>::functor_type functor_type;

	///	return functor
		virtual functor_type get_functor() const {return boost::ref(*this);}

	public:
	///	Constructor
		ConstUserVector() {set_all_entries(0.0);}

	///	set all vector entries
		void set_all_entries(number val) { m_Vector = val;}

	///	set i'th vector entry
		void set_entry(size_t i, number val)
		{
			m_Vector[i] = val;
		}

	///	print current setting
		void print() const
		{
			UG_LOG("ConstUserVector:" << m_Vector << "\n");
		}

	/// evaluate
		void operator() (MathVector<dim>& v, const MathVector<dim>& x,
		                 number time = 0.0) const
		{
			v = m_Vector;
		}

	protected:
		MathVector<dim> m_Vector;
};

/// constant matrix user data
template <int dim>
class ConstUserMatrix
	: public IUserMatrix<dim>
{
	public:
	///	Functor Type
		typedef typename IUserMatrix<dim>::functor_type functor_type;

	///	return functor
		virtual functor_type get_functor() const {return boost::ref(*this);}

	public:
	///	Constructor
		ConstUserMatrix() {set_diag_tensor(1.0);}

	///	set diagonal of matrix to a vector
		void set_diag_tensor(number val)
		{
			for(size_t i = 0; i < dim; ++i){
				for(size_t j = 0; j < dim; ++j){
					m_Tensor[i][j] = 0;
				}
				m_Tensor[i][i] = val;
			}
		}

	///	sets all entries of the matrix
		void set_all_entries(number val)
		{
			for(size_t i = 0; i < dim; ++i){
				for(size_t j = 0; j < dim; ++j){
					m_Tensor[i][j] = val;
				}
			}
		}

	///	sets a single entry
		void set_entry(size_t i, size_t j, number val)
		{
			m_Tensor[i][j] = val;
		}

	///	print current setting
		void print() const
		{
			UG_LOG("ConstUserMatrix:\n" << m_Tensor << "\n");
		}

	///	evaluate
		void operator() (MathMatrix<dim, dim>& D, const MathVector<dim>& x,
		                 number time = 0.0) const
		{
			D = m_Tensor;
		}

	protected:
		MathMatrix<dim, dim> m_Tensor;
};

/// constant dirichlet boundary scalar user data
template <int dim>
class ConstBoundaryNumber : public IBoundaryNumberProvider<dim>
{
	public:
	///	Functor Type
		typedef typename IBoundaryNumberProvider<dim>::functor_type functor_type;

	///	return functor
		virtual functor_type get_functor() const {return boost::ref(*this);}

	public:
	///	Constructor
		ConstBoundaryNumber() {m_Number = 0.0;}

	///	set value
		void set(number val)
		{
			m_Number = val;
		}

	///	print current setting
		void print() const
		{
			UG_LOG("ConstBoundaryNumber:" << m_Number << "\n");
		}

	///	evaluate and return true for dirichlet value
		bool operator() (number& c, const MathVector<dim>& x,
		                 number time = 0.0) const
		{
			c = m_Number;
			return true;
		}

	protected:
		number m_Number;
};

/// @}

} /// end namespace ug

#endif /* __H__UG__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__USER_DATA__ */
