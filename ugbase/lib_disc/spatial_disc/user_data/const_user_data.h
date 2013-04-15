/*
 * const_user_data.h
 *
 *  Created on: 13.10.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__USER_DATA__CONST_USER_DATA__
#define __H__UG__LIB_DISC__SPATIAL_DISC__USER_DATA__CONST_USER_DATA__

#include "common/common.h"
#include "common/math/ugmath.h"

#include "std_user_data.h"

namespace ug {


///////////////////////////////////////////////////////////////////////////////
// Base class for Constant Data
///////////////////////////////////////////////////////////////////////////////

/**
 * This class is a base class for all Constant user data. The data thus does not
 * depend neither on space, time or subset nor on the a computed solution.
 * In order to use the interface, the deriving class must implement the method:
 *
 * inline void evaluate(TData& data) const
 *
 */
template <typename TImpl, typename TData, int dim>
class StdConstData
	: 	public StdUserData<StdConstData<TImpl,TData,dim>, TData, dim>
{
	public:
		virtual void operator() (TData& value,
								 const MathVector<dim>& globIP,
								 number time, int si) const
		{
			getImpl().evaluate(value);
		}

		virtual void operator()(TData vValue[],
								const MathVector<dim> vGlobIP[],
								number time, int si, const size_t nip) const
		{
			for(size_t ip = 0; ip < nip; ++ip)
				getImpl().evaluate(vValue[ip]);
		}

		template <int refDim>
		inline void evaluate(TData vValue[],
		                     const MathVector<dim> vGlobIP[],
		                     number time, int si,
		                     GeometricObject* elem,
		                     const MathVector<dim> vCornerCoords[],
		                     const MathVector<refDim> vLocIP[],
		                     const size_t nip,
		                     LocalVector* u,
		                     const MathMatrix<refDim, dim>* vJT = NULL) const
		{
			for(size_t ip = 0; ip < nip; ++ip)
				getImpl().evaluate(vValue[ip]);
		}

	///	implement as a UserData
		virtual void compute(LocalVector* u, GeometricObject* elem,
		                     const MathVector<dim> vCornerCoords[], bool bDeriv = false)
		{
			for(size_t s = 0; s < this->num_series(); ++s)
				for(size_t ip = 0; ip < this->num_ip(s); ++ip)
					getImpl().evaluate(this->value(s,ip));
		}

	///	callback, invoked when data storage changed
		virtual void value_storage_changed(const size_t seriesID)
		{
			for(size_t ip = 0; ip < this->num_ip(seriesID); ++ip)
				getImpl().evaluate(this->value(seriesID,ip));
		}

	///	returns if data is constant
		virtual bool constant() const {return true;}

	///	returns if grid function is needed for evaluation
		virtual bool requires_grid_fct() const {return false;}

	///	returns if provided data is continuous over geometric object boundaries
		virtual bool continuous() const {return true;}

	protected:
	///	access to implementation
		TImpl& getImpl() {return static_cast<TImpl&>(*this);}

	///	const access to implementation
		const TImpl& getImpl() const {return static_cast<const TImpl&>(*this);}
};

///////////////////////////////////////////////////////////////////////////////
// Constant UserData
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
	: public StdConstData<ConstUserNumber<dim>, number, dim>
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
/**
 * Constant vector user data that can be used in assembling routines.
 *
 * \param dim the dimensionality of the vector itself (for ex. 2 for vectors of two components)
 * \param worldDim the dimensionality of the space embedding the grid (for ex. 3 for 3d PDE problems)
 */
template <int dim, int worldDim = dim>
class ConstUserVector
	: public StdConstData<ConstUserVector<dim, worldDim>, MathVector<dim>, worldDim>
{
	public:
	///	Constructor: no arguments, zero entries
		ConstUserVector() {set_all_entries(0.0);}

	///	Constructor: set all the entries to the given value
		ConstUserVector(number val) {set_all_entries(val);}

	///	Constructor: initialize with a given std::vector
		ConstUserVector(const std::vector<number>& val) {set_vector(val);}

	///	set all vector entries
		void set_all_entries(number val) { m_Vector = val;}

	///	set i'th vector entry
		void set_entry(size_t i, number val){m_Vector[i] = val;}
	
	/// set from a given vector:
		void set_vector(const std::vector<number>& val)
		{
			if(val.size() != dim) UG_THROW("Size mismatch in ConstUserVector");
			for(size_t i = 0; i < dim; i++) m_Vector[i] = val[i];
		}

	///	print current setting
		void print() const {UG_LOG("ConstUserVector:" << m_Vector << "\n");}

	/// evaluate
		inline void evaluate (MathVector<dim>& value) const{value = m_Vector;}

	protected:
		MathVector<dim> m_Vector;
};

/// constant matrix user data
/**
 * Constant matrix user data that can be used in assembling routines.
 *
 * \param N the row size of the matrix
 * \param M the column size of the matrix
 * \param worldDim the dimensionality of the space embedding the grid (for ex. 3 for 3d PDE problems)
 */
template <int N, int M = N, int worldDim = N>
class ConstUserMatrix
	: public StdConstData<ConstUserMatrix<N, M, worldDim>, MathMatrix<N, M>, worldDim>
{
	public:
	///	Constructor
		ConstUserMatrix() {set_diag_tensor(1.0);}

	///	Constructor setting the diagonal
		ConstUserMatrix(number val) {set_diag_tensor(val);}

	///	set diagonal of matrix to a vector
		void set_diag_tensor(number val)
		{
			for(size_t i = 0; i < N; ++i){
				for(size_t j = 0; j < M; ++j){
					m_Tensor[i][j] = 0;
				}
				m_Tensor[i][i] = val;
			}
		}

	///	sets all entries of the matrix
		void set_all_entries(number val)
		{
			for(size_t i = 0; i < N; ++i){
				for(size_t j = 0; j < M; ++j){
					m_Tensor[i][j] = val;
				}
			}
		}

	///	sets a single entry
		void set_entry(size_t i, size_t j, number val){m_Tensor[i][j] = val;}

	///	print current setting
		void print() const{UG_LOG("ConstUserMatrix:\n" << m_Tensor << "\n");}

	///	evaluate
		inline void evaluate (MathMatrix<N, M>& value) const{value = m_Tensor;}

	protected:
		MathMatrix<N, M> m_Tensor;
};

/// constant tensor user data
template <int TRank, int dim>
class ConstUserTensor
	: public StdConstData<ConstUserTensor<TRank,dim>, MathTensor<TRank, dim>, dim>
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
SmartPtr<CplUserData<TData,dim> > CreateConstUserData(number val, TData dummy);

template <int dim>
inline SmartPtr<CplUserData<number,dim> > CreateConstUserData(number val, number)
{
	return CreateSmartPtr(new ConstUserNumber<dim>(val));
};

template <int dim>
SmartPtr<CplUserData<MathVector<dim>,dim> > CreateConstUserData(number val, MathVector<dim>)
{
	return CreateSmartPtr(new ConstUserVector<dim>(val));
}

template <int dim>
SmartPtr<CplUserData<MathMatrix<dim,dim>,dim> > CreateConstUserData(number val, MathMatrix<dim,dim>)
{
	return CreateSmartPtr(new ConstUserMatrix<dim>(val));
}

template <int dim>
SmartPtr<CplUserData<MathTensor<4,dim>,dim> > CreateConstUserData(number val, MathTensor<4,dim>)
{
	return CreateSmartPtr(new ConstUserTensor<4,dim>(val));
}

/// @}

} /// end namespace ug

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__USER_DATA__CONST_USER_DATA__ */
