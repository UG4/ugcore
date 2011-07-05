/*
 * const_user_data.h
 *
 *  Created on: 13.10.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__IP_DATA__CONST_USER_DATA__
#define __H__UG__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__IP_DATA__CONST_USER_DATA__

#include "common/common.h"
#include "common/math/ugmath.h"
#include <boost/function.hpp>

#include "user_data_interface.h"

namespace ug {

/// \addtogroup lib_disc_user_data
/// @{

/// constant scalar user data
template <int dim>
class ConstUserNumber
	: public IUserData<number, dim>
{
	///	Base class type
		typedef IUserData<number, dim> base_type;

		using base_type::num_series;
		using base_type::num_ip;
		using base_type::ip;
		using base_type::time;
		using base_type::value;

	public:
	///	Functor Type
		typedef typename base_type::functor_type functor_type;

	///	return functor
		virtual functor_type get_functor() const {return boost::ref(*this);}

	public:
	///	creates empty user number
		ConstUserNumber() { set(0.0);}

	///	set constant value
		void set(number val) {m_Number = val;}

	///	print current setting
		void print() const {UG_LOG("ConstUserNumber:" << m_Number << "\n");}

	///	evaluate
		void operator() (number& c, const MathVector<dim>& x,
		                 number time = 0.0) const
		{
			c = m_Number;
		}

	///	implement as a IPData
		virtual bool compute(bool bDeriv = false)
		{
			for(size_t s = 0; s < num_series(); ++s)
				for(size_t i = 0; i < num_ip(s); ++i)
					value(s,i) = m_Number;
			return true;
		}

	///	returns if data is constant
		virtual bool constant_data() const {return true;}

	protected:
		number m_Number;
};

/// constant vector user data
template <int dim>
class ConstUserVector
	: public IUserData<MathVector<dim>, dim>
{
	/// Base class type
		typedef IUserData<MathVector<dim>, dim> base_type;

		using base_type::num_series;
		using base_type::num_ip;
		using base_type::ip;
		using base_type::time;
		using base_type::value;

	public:
	///	Functor Type
		typedef typename base_type::functor_type functor_type;

	///	return functor
		virtual functor_type get_functor() const {return boost::ref(*this);}

	public:
	///	Constructor
		ConstUserVector() {set_all_entries(0.0);}

	///	set all vector entries
		void set_all_entries(number val) { m_Vector = val;}

	///	set i'th vector entry
		void set_entry(size_t i, number val){m_Vector[i] = val;}

	///	print current setting
		void print() const {UG_LOG("ConstUserVector:" << m_Vector << "\n");}

	/// evaluate
		void operator() (MathVector<dim>& v, const MathVector<dim>& x,
		                 number time = 0.0) const
		{
			v = m_Vector;
		}

	///	returns if data is constant
		virtual bool constant_data() const {return true;}

	///	implement as a IPData
		virtual bool compute(bool bDeriv = false)
		{
			for(size_t s = 0; s < num_series(); ++s)
				for(size_t i = 0; i < num_ip(s); ++i)
					value(s,i) = m_Vector;
			return true;
		}

	protected:
		MathVector<dim> m_Vector;
};

/// constant matrix user data
template <int dim>
class ConstUserMatrix
	: public IUserData<MathMatrix<dim, dim>, dim>
{
	/// Base class type
		typedef IUserData<MathMatrix<dim, dim>, dim> base_type;

		using base_type::num_series;
		using base_type::num_ip;
		using base_type::ip;
		using base_type::time;
		using base_type::value;

	public:
	///	Functor Type
		typedef typename base_type::functor_type functor_type;

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
		void set_entry(size_t i, size_t j, number val){m_Tensor[i][j] = val;}

	///	print current setting
		void print() const{UG_LOG("ConstUserMatrix:\n" << m_Tensor << "\n");}

	///	evaluate
		void operator() (MathMatrix<dim, dim>& D, const MathVector<dim>& x,
		                 number time = 0.0) const
		{
			D = m_Tensor;
		}

	///	returns if data is constant
		virtual bool constant_data() const {return true;}

	///	implement as a IPData
		virtual bool compute(bool bDeriv = false)
		{
			for(size_t s = 0; s < num_series(); ++s)
				for(size_t i = 0; i < num_ip(s); ++i)
					value(s,i) = m_Tensor;
			return true;
		}

	protected:
		MathMatrix<dim, dim> m_Tensor;
};

/// constant dirichlet boundary scalar user data
template <int dim>
class ConstBoundaryNumber : public IBoundaryData<number, dim>
{
	public:
	///	Functor Type
		typedef typename IBoundaryData<number, dim>::functor_type functor_type;

	///	return functor
		virtual functor_type get_functor() const {return boost::ref(*this);}

	public:
	///	Constructor
		ConstBoundaryNumber() {set(0.0);}

	///	set value
		void set(number val) {m_Number = val;}

	///	print current setting
		void print() const{UG_LOG("ConstBoundaryNumber:" << m_Number << "\n");}

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

#endif /* __H__UG__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__IP_DATA__CONST_USER_DATA__ */
