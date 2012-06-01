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

///////////////////////////////////////////////////////////////////////////////
// Base class for const data
///////////////////////////////////////////////////////////////////////////////

template <typename TData, int dim, typename TImpl>
class StdConstUserData
	: 	public IPData<TData,dim>
{
	public:
		StdConstUserData() {}

		virtual void operator() (TData& value,
		                         const MathVector<dim>& globIP,
		                         number time, int si) const
		{
			getImpl().evaluate(value);
		}

		virtual void operator() (TData vValue[],
		                         const MathVector<dim> vGlobIP[],
		                         number time, int si, const size_t nip) const
		{
			for(size_t ip = 0; ip < nip; ++ip)
				getImpl().evaluate(vValue[ip]);
		}

		////////////////
		// one value
		////////////////

		virtual void operator() (TData& value,
		                         const MathVector<dim>& globIP,
		                         number time, int si,
		                         LocalVector& u,
		                         GeometricObject* elem,
		                         const MathVector<dim> vCornerCoords[],
		                         const MathVector<1>& locIP) const
		{
			getImpl().evaluate(value);
		}

		virtual void operator() (TData& value,
		                         const MathVector<dim>& globIP,
		                         number time, int si,
		                         LocalVector& u,
		                         GeometricObject* elem,
		                         const MathVector<dim> vCornerCoords[],
		                         const MathVector<2>& locIP) const
		{
			getImpl().evaluate(value);
		}

		virtual void operator() (TData& value,
		                         const MathVector<dim>& globIP,
		                         number time, int si,
		                         LocalVector& u,
		                         GeometricObject* elem,
		                         const MathVector<dim> vCornerCoords[],
		                         const MathVector<3>& locIP) const
		{
			getImpl().evaluate(value);
		}

		////////////////
		// vector of values
		////////////////

		virtual void operator()(TData vValue[],
		                        const MathVector<dim> vGlobIP[],
		                        number time, int si,
		                        LocalVector& u,
		                        GeometricObject* elem,
		                        const MathVector<dim> vCornerCoords[],
		                        const MathVector<1> vLocIP[],
		                        const size_t nip) const
		{
			for(size_t ip = 0; ip < nip; ++ip)
				getImpl().evaluate(vValue[ip]);
		}

		virtual void operator()(TData vValue[],
		                        const MathVector<dim> vGlobIP[],
		                        number time, int si,
		                        LocalVector& u,
		                        GeometricObject* elem,
		                        const MathVector<dim> vCornerCoords[],
		                        const MathVector<2> vLocIP[],
		                        const size_t nip) const
		{
			for(size_t ip = 0; ip < nip; ++ip)
				getImpl().evaluate(vValue[ip]);
		}

		virtual void operator()(TData vValue[],
		                        const MathVector<dim> vGlobIP[],
		                        number time, int si,
		                        LocalVector& u,
		                        GeometricObject* elem,
		                        const MathVector<dim> vCornerCoords[],
		                        const MathVector<3> vLocIP[],
		                        const size_t nip) const
		{
			for(size_t ip = 0; ip < nip; ++ip)
				getImpl().evaluate(vValue[ip]);
		}

	protected:
	///	access to implementation
		TImpl& getImpl() {return static_cast<TImpl&>(*this);}

	///	const access to implementation
		const TImpl& getImpl() const {return static_cast<const TImpl&>(*this);}
};


///////////////////////////////////////////////////////////////////////////////
// Implementations
///////////////////////////////////////////////////////////////////////////////

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
	: public StdConstUserData<number, dim, ConstUserNumber<dim> >
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
		inline void evaluate (number& value) const
		{
			value = m_Number;
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
	: public StdConstUserData<MathVector<dim>, dim, ConstUserVector<dim> >
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
		inline void evaluate (MathVector<dim>& value) const
		{
			value = m_Vector;
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
	: public StdConstUserData<MathMatrix<dim, dim>, dim, ConstUserMatrix<dim> >
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
		inline void evaluate (MathMatrix<dim, dim>& value) const
		{
			value = m_Tensor;
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

/// constant tensor user data
template <int TRank, int dim>
class ConstUserTensor
	: public StdConstUserData<MathTensor<TRank, dim>, dim, ConstUserTensor<TRank,dim> >
{
	/// Base class type
		typedef IPData<MathTensor<TRank, dim>, dim> base_type;

		using base_type::num_series;
		using base_type::num_ip;
		using base_type::ip;
		using base_type::time;
		using base_type::value;

	public:
	///	Constructor
		ConstUserTensor() {set(0.0);}

	///	Constructor setting the diagonal
		ConstUserTensor(number val) {set(val);}

	///	set diagonal of matrix to a vector
		void set(number val)
		{
			m_Tensor.set(val);
		}

	///	print current setting
		void print() const{UG_LOG("ConstUserTensor:\n" << m_Tensor << "\n");}

	///	evaluate
		inline void evaluate (MathTensor<TRank, dim>& value) const
		{
			value = m_Tensor;
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
		MathTensor<TRank, dim> m_Tensor;
};

/// creates user data of desired type
template <typename TData, int dim>
SmartPtr<IPData<TData,dim> > CreateConstUserData(number val, TData dummy);

template <int dim>
inline SmartPtr<IPData<number,dim> > CreateConstUserData(number val, number)
{
	return CreateSmartPtr(new ConstUserNumber<dim>(val));
};

template <int dim>
SmartPtr<IPData<MathVector<dim>,dim> > CreateConstUserData(number val, MathVector<dim>)
{
	return CreateSmartPtr(new ConstUserVector<dim>(val));
}

template <int dim>
SmartPtr<IPData<MathMatrix<dim,dim>,dim> > CreateConstUserData(number val, MathMatrix<dim,dim>)
{
	return CreateSmartPtr(new ConstUserMatrix<dim>(val));
}

template <int dim>
SmartPtr<IPData<MathTensor<4,dim>,dim> > CreateConstUserData(number val, MathTensor<4,dim>)
{
	return CreateSmartPtr(new ConstUserTensor<4,dim>(val));
}

/// @}

} /// end namespace ug

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__IP_DATA__CONST_USER_DATA__ */
