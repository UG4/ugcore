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

#include "lib_disc/spatial_disc/ip_data/std_ip_data.h"

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
	: public StdConstIPData<ConstUserNumber<dim>, number, dim>
{
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
		inline void evaluate (number& value) const {value = m_Number;}

	protected:
		number m_Number;
};

/// constant vector user data
template <int dim>
class ConstUserVector
	: public StdConstIPData<ConstUserVector<dim>, MathVector<dim>, dim>
{
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
		inline void evaluate (MathVector<dim>& value) const{value = m_Vector;}

	protected:
		MathVector<dim> m_Vector;
};

/// constant matrix user data
template <int dim>
class ConstUserMatrix
	: public StdConstIPData<ConstUserMatrix<dim>, MathMatrix<dim, dim>, dim>
{
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
		inline void evaluate (MathMatrix<dim, dim>& value) const{value = m_Tensor;}

	protected:
		MathMatrix<dim, dim> m_Tensor;
};

/// constant tensor user data
template <int TRank, int dim>
class ConstUserTensor
	: public StdConstIPData<ConstUserTensor<TRank,dim>, MathTensor<TRank, dim>, dim>
{
	public:
	///	Constructor
		ConstUserTensor() {set(0.0);}

	///	Constructor setting the diagonal
		ConstUserTensor(number val) {set(val);}

	///	set diagonal of matrix to a vector
		void set(number val) {m_Tensor.set(val);}

	///	print current setting
		void print() const{UG_LOG("ConstUserTensor:\n" << m_Tensor << "\n");}

	///	evaluate
		inline void evaluate (MathTensor<TRank, dim>& value) const{value = m_Tensor;}

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
