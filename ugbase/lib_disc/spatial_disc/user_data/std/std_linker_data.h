/*
 * std_linker_data.h
 *
 *  Created on: 03.07.2012
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__STD_LINKER_DATA__
#define __H__UG__LIB_DISC__SPATIAL_DISC__STD_LINKER_DATA__

#include "std_user_data.h"
#include "../data_linker.h"

namespace ug{

///////////////////////////////////////////////////////////////////////////////
// Base class for Linker UserData
///////////////////////////////////////////////////////////////////////////////

/**
 * This class is a base class for linker user data.
 * In order to use the interface, the deriving class must implement the method:
 *
 * template \<int refDim\>
 * inline TRet evaluate(TData vValue[],
						const MathVector\<dim\> vGlobIP[],
						number time, int si,
						LocalVector& u,
						GeometricObject* elem,
						const MathVector\<dim\> vCornerCoords[],
						const MathVector\<refDim\> vLocIP[],
						const size_t nip,
						const MathMatrix\<refDim, dim\>* vJT = NULL) const
 *
 */
template <typename TImpl, typename TData, int dim>
class StdDataLinker
	: 	public StdUserData<		StdDataLinker<TImpl,TData,dim>,
								DataLinker<TData,dim>,
								TData,dim>
{
	public:
		inline void evaluate (TData& value,
		                      const MathVector<dim>& globIP,
		                      number time, int si) const
		{
			getImpl().evaluate(value,globIP,time,si);
		}

		inline void evaluate (TData vValue[],
		                      const MathVector<dim> vGlobIP[],
		                      number time, int si, const size_t nip) const
		{
			for(size_t ip = 0; ip < nip; ++ip)
				getImpl().evaluate(vValue[ip],vGlobIP[ip],time,si);
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
			getImpl().template evaluate<refDim>(value,globIP,time,si,u,elem,
			                                    vCornerCoords,locIP);
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
			getImpl().template evaluate<refDim>(vValue,vGlobIP,time,si,u,elem,
			                                    vCornerCoords,vLocIP,nip, vJT);
		}


	///	returns that a grid function is needed for evaluation
		virtual bool requires_grid_fct() const
		{
			for(size_t i = 0; i < this->m_vspICplUserData.size(); ++i)
				if(this->m_vspUserDataInfo[i]->requires_grid_fct())
					return true;
			return false;
		}

	///	returns if provided data is continuous over geometric object boundaries
		virtual bool continuous() const
		{
			bool bRet = true;
			for(size_t i = 0; i < this->m_vspICplUserData.size(); ++i)
				bRet &= this->m_vspUserDataInfo[i]->continuous();
			return bRet;
		}

	///	sets the associated function pattern
		virtual void set_function_pattern(const FunctionPattern& fctPatt)
		{
			for(size_t i = 0; i < this->m_vspICplUserData.size(); ++i)
				this->m_vspUserDataInfo[i]->set_function_pattern(fctPatt);
		}

	protected:
	///	access to implementation
		TImpl& getImpl() {return static_cast<TImpl&>(*this);}

	///	const access to implementation
		const TImpl& getImpl() const {return static_cast<const TImpl&>(*this);}
};

} // namespace ug

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__STD_LINKER_DATA__ */
