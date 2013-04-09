/*
 * std_const_data.h
 *
 *  Created on: 03.07.2012
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__STD_CONST_DATA__
#define __H__UG__LIB_DISC__SPATIAL_DISC__STD_CONST_DATA__

#include "std_user_data.h"

namespace ug{

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
	: 	public StdUserData<		StdConstData<TImpl,TData,dim>,
								UserData<TData,dim>,
								TData,dim>
{
	public:
		inline void evaluate (TData& value,
		                      const MathVector<dim>& globIP,
		                      number time, int si) const
		{
			getImpl().evaluate(value);
		}

		inline void evaluate (TData vValue[],
		                      const MathVector<dim> vGlobIP[],
		                      number time, int si, const size_t nip) const
		{
			for(size_t ip = 0; ip < nip; ++ip)
				getImpl().evaluate(vValue[ip]);
		}

		template <int refDim>
		inline void evaluate (TData& value,
		                      const MathVector<dim>& globIP,
		                      number time, int si,
		                      LocalVector& u,
		                      GeometricObject* elem,
		                      const MathVector<dim> vCornerCoords[],
		                      const MathVector<refDim>& locIP) const
		{
			getImpl().evaluate(value);
		}

		template <int refDim>
		inline void evaluate(TData vValue[],
		                     const MathVector<dim> vGlobIP[],
		                     number time, int si,
		                     LocalVector& u,
		                     GeometricObject* elem,
		                     const MathVector<dim> vCornerCoords[],
		                     const MathVector<refDim> vLocIP[],
		                     const size_t nip,
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

} // namespace ug

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__STD_CONST_DATA__ */
