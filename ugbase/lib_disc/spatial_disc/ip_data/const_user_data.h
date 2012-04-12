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
	: public IPData<number, dim>
{
	///	Base class type
		typedef IPData<number, dim> base_type;

		using base_type::num_series;
		using base_type::num_ip;
		using base_type::ip;
		using base_type::time;
		using base_type::value;

	public:
	///	creates empty user number
		ConstUserNumber() {set(0.0);}

	///	creates user number with value
		ConstUserNumber(number val) {set(val);}

	///	set constant value
		void set(number val) {m_Number = val;}

	///	print current setting
		void print() const {UG_LOG("ConstUserNumber:" << m_Number << "\n");}

	///	evaluate
		void operator() (number& c, const MathVector<dim>& x,
		                 number time, int si) const
		{
			c = m_Number;
		}

	///	implement as a IPData
		virtual void compute(bool bDeriv = false)
		{
			for(size_t s = 0; s < num_series(); ++s)
				for(size_t i = 0; i < num_ip(s); ++i)
					value(s,i) = m_Number;
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
	: public IPData<MathVector<dim>, dim>
{
	/// Base class type
		typedef IPData<MathVector<dim>, dim> base_type;

		using base_type::num_series;
		using base_type::num_ip;
		using base_type::ip;
		using base_type::time;
		using base_type::value;

	public:
	///	Constructor
		ConstUserVector() {set_all_entries(0.0);}

	///	creates user number with value
		ConstUserVector(number val) {set_all_entries(val);}

	///	set all vector entries
		void set_all_entries(number val) { m_Vector = val;}

	///	set i'th vector entry
		void set_entry(size_t i, number val){m_Vector[i] = val;}

	///	print current setting
		void print() const {UG_LOG("ConstUserVector:" << m_Vector << "\n");}

	/// evaluate
		void operator() (MathVector<dim>& v, const MathVector<dim>& x,
		                 number time, int si) const
		{
			v = m_Vector;
		}

	///	returns if data is constant
		virtual bool constant_data() const {return true;}

	///	implement as a IPData
		virtual void compute(bool bDeriv = false)
		{
			for(size_t s = 0; s < num_series(); ++s)
				for(size_t i = 0; i < num_ip(s); ++i)
					value(s,i) = m_Vector;
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
	: public IPData<MathMatrix<dim, dim>, dim>
{
	/// Base class type
		typedef IPData<MathMatrix<dim, dim>, dim> base_type;

		using base_type::num_series;
		using base_type::num_ip;
		using base_type::ip;
		using base_type::time;
		using base_type::value;

	public:
	///	Constructor
		ConstUserMatrix() {set_diag_tensor(1.0);}

	///	Constructor setting the diagonal
		ConstUserMatrix(number val) {set_diag_tensor(val);}

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
		                 number time, int si) const
		{
			D = m_Tensor;
		}

	///	returns if data is constant
		virtual bool constant_data() const {return true;}

	///	implement as a IPData
		virtual void compute(bool bDeriv = false)
		{
			for(size_t s = 0; s < num_series(); ++s)
				for(size_t i = 0; i < num_ip(s); ++i)
					value(s,i) = m_Tensor;
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

/// @}

} /// end namespace ug

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__IP_DATA__CONST_USER_DATA__ */
