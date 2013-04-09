/*
 * std_pos_data.h
 *
 *  Created on: 03.07.2012
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__STD_POS_DATA__
#define __H__UG__LIB_DISC__SPATIAL_DISC__STD_POS_DATA__

#include "std_user_data.h"

namespace ug{

///////////////////////////////////////////////////////////////////////////////
// Base class for Position-Time-UserData
///////////////////////////////////////////////////////////////////////////////

/**
 * This class is a base class for all position and time dependent user data.
 * The data thus does not on the a computed solution.
 * In order to use the interface, the deriving class must implement the method:
 *
 * inline TRet evaluate(TData& D, const MathVector<dim>& x, number time, int si) const
 *
 */
template <typename TImpl, typename TData, int dim, typename TRet = void>
class StdPositionData
	: 	public StdUserData<		StdPositionData<TImpl,TData,dim,TRet>,
								CplUserData<TData,dim,TRet>,
								TData,dim,TRet>
{
	public:
		inline TRet evaluate (TData& value,
		                      const MathVector<dim>& globIP,
		                      number time, int si) const
		{
			return this->getImpl().evaluate(value, globIP, time, si);
		}

		inline void evaluate (TData vValue[],
		                      const MathVector<dim> vGlobIP[],
		                      number time, int si, const size_t nip) const
		{
			for(size_t ip = 0; ip < nip; ++ip)
				this->getImpl().evaluate(vValue[ip], vGlobIP[ip], time, si);
		}

		template <int refDim>
		inline TRet evaluate (TData& value,
		                      const MathVector<dim>& globIP,
		                      number time, int si,
		                      LocalVector& u,
		                      GeometricObject* elem,
		                      const MathVector<dim> vCornerCoords[],
		                      const MathVector<refDim>& locIP) const
		{
			return this->getImpl().evaluate(value,globIP,time,si);
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
				this->getImpl().evaluate(vValue[ip],vGlobIP[ip],time,si);
		}

	///	implement as a UserData
		virtual void compute(LocalVector* u, GeometricObject* elem,
		                     const MathVector<dim> vCornerCoords[], bool bDeriv = false)
		{
			const number t = this->time();
			const int si = this->subset();

			for(size_t s = 0; s < this->num_series(); ++s)
				for(size_t ip = 0; ip < this->num_ip(s); ++ip)
					this->getImpl().evaluate(this->value(s,ip), this->ip(s, ip), t, si);
		}

	///	returns if data is constant
		virtual bool constant() const {return false;}

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

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__STD_POS_DATA__ */
