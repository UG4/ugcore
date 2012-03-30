/*
 * const_user_data.h
 *
 *  Created on: 13.10.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__IP_DATA__CONST_USER_DATA__
#define __H__UG__LIB_DISC__SPATIAL_DISC__IP_DATA__CONST_USER_DATA__

#include "common/common.h"
#include "common/math/ugmath.h"

#include <boost/function.hpp>
#include "ip_data.h"

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

/// constant scalar user data
template <int dim>
class ConstUserNumber
	: public IPData<number, dim>,
	  public boost::function<void (number& res, const MathVector<dim>& x,number time)>
{
	///	Base class type
		typedef IPData<number, dim> base_type;

	///	Functor type
		typedef boost::function<void (number& res, const MathVector<dim>& x,number time)> func_type;

		using base_type::num_series;
		using base_type::num_ip;
		using base_type::ip;
		using base_type::time;
		using base_type::value;

	public:
	///	creates empty user number
		ConstUserNumber() : func_type(boost::ref(*this)) {set(0.0);}

	///	creates user number with value
		ConstUserNumber(number val) : func_type(boost::ref(*this)) {set(val);}

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

	///	callback, invoked when data storage changed
		virtual void value_storage_changed(const size_t seriesID)
		{
			for(size_t i = 0; i < num_ip(seriesID); ++i)
				value(seriesID,i) = m_Number;
		}

	protected:
		number m_Number;
};

/// constant vector user data
template <int dim>
class ConstUserVector
	: public IPData<MathVector<dim>, dim>,
	  public boost::function<void (MathVector<dim>& res, const MathVector<dim>& x,number time)>
{
	/// Base class type
		typedef IPData<MathVector<dim>, dim> base_type;

	///	Functor type
		typedef boost::function<void (MathVector<dim>& res, const MathVector<dim>& x,number time)> func_type;

		using base_type::num_series;
		using base_type::num_ip;
		using base_type::ip;
		using base_type::time;
		using base_type::value;

	public:
	///	Constructor
		ConstUserVector() : func_type(boost::ref(*this)) {set_all_entries(0.0);}

	///	creates user number with value
		ConstUserVector(number val) : func_type(boost::ref(*this)) {set_all_entries(val);}

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

	///	callback, invoked when data storage changed
		virtual void value_storage_changed(const size_t seriesID)
		{
			for(size_t i = 0; i < num_ip(seriesID); ++i)
				value(seriesID,i) = m_Vector;
		}

	protected:
		MathVector<dim> m_Vector;
};

/// constant matrix user data
template <int dim>
class ConstUserMatrix
	: public IPData<MathMatrix<dim, dim>, dim>,
	  public boost::function<void (MathMatrix<dim, dim>& res, const MathVector<dim>& x,number time)>
{
	/// Base class type
		typedef IPData<MathMatrix<dim, dim>, dim> base_type;

	///	Functor type
		typedef boost::function<void (MathMatrix<dim, dim>& res, const MathVector<dim>& x,number time)> func_type;

		using base_type::num_series;
		using base_type::num_ip;
		using base_type::ip;
		using base_type::time;
		using base_type::value;

	public:
	///	Constructor
		ConstUserMatrix() : func_type(boost::ref(*this)) {set_diag_tensor(1.0);}

	///	Constructor setting the diagonal
		ConstUserMatrix(number val) : func_type(boost::ref(*this)) {set_diag_tensor(val);}

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

	///	callback, invoked when data storage changed
		virtual void value_storage_changed(const size_t seriesID)
		{
			for(size_t i = 0; i < num_ip(seriesID); ++i)
				value(seriesID,i) = m_Tensor;
		}

	protected:
		MathMatrix<dim, dim> m_Tensor;
};

/// constant dirichlet boundary scalar user data
template <int dim>
class ConstBoundaryNumber
	 : public boost::function<bool (number& res, const MathVector<dim>& x,number time)>
{
	/// functor type
		typedef boost::function<bool (number& res, const MathVector<dim>& x,number time)> func_type;

	public:
	///	Constructor
		ConstBoundaryNumber() : func_type(boost::ref(*this)) {set(0.0);}

	///	Constructor
		ConstBoundaryNumber(number val) : func_type(boost::ref(*this)) {set(val);}

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

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__IP_DATA__CONST_USER_DATA__ */
