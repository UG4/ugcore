/*
 * std_user_data.h
 *
 *  Created on: 03.07.2012
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__STD_USER_DATA__
#define __H__UG__LIB_DISC__SPATIAL_DISC__STD_USER_DATA__

#include "../user_data.h"

namespace ug{

///////////////////////////////////////////////////////////////////////////////
// Wrapper Class for UserData
///////////////////////////////////////////////////////////////////////////////

/**
 * This class is used as a wrapper class for user data in order to ease the
 * implementation of the virtual evaluation operators through templated
 * methods evaluate.
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
template <typename TImpl, typename TBase, typename TData, int dim, typename TRet = void>
class StdUserData : public TBase
{
	public:
	///	returns value for a global position
		virtual TRet operator() (TData& value,
								 const MathVector<dim>& globIP,
								 number time, int si) const
		{
			return getImpl().evaluate(value, globIP, time, si);
		}

	///	returns value for local and global position
	///	\{
		virtual TRet operator() (TData& value,
		                         const MathVector<dim>& globIP,
		                         number time, int si,
		                         LocalVector& u,
		                         GeometricObject* elem,
		                         const MathVector<dim> vCornerCoords[],
		                         const MathVector<1>& locIP) const
		{
			return getImpl().template evaluate<1>(value,globIP,time,si,u,
			                                      elem,vCornerCoords,locIP);
		}

		virtual TRet operator() (TData& value,
		                         const MathVector<dim>& globIP,
		                         number time, int si,
		                         LocalVector& u,
		                         GeometricObject* elem,
		                         const MathVector<dim> vCornerCoords[],
		                         const MathVector<2>& locIP) const
		{
			return getImpl().template evaluate<2>(value,globIP,time,si,u,
			                                      elem,vCornerCoords,locIP);
		}

		virtual TRet operator() (TData& value,
		                         const MathVector<dim>& globIP,
		                         number time, int si,
		                         LocalVector& u,
		                         GeometricObject* elem,
		                         const MathVector<dim> vCornerCoords[],
		                         const MathVector<3>& locIP) const
		{
			return getImpl().template evaluate<3>(value,globIP,time,si,u,
			                                      elem,vCornerCoords,locIP);
		}
	///	\}

	///	returns value for global positions
		virtual void operator()(TData vValue[],
								const MathVector<dim> vGlobIP[],
								number time, int si, const size_t nip) const
		{
				getImpl().evaluate(vValue,vGlobIP,time,si,nip);
		}

	///	returns values for local and global positions
	///	\{
		virtual void operator()(TData vValue[],
		                        const MathVector<dim> vGlobIP[],
		                        number time, int si,
		                        LocalVector& u,
		                        GeometricObject* elem,
		                        const MathVector<dim> vCornerCoords[],
		                        const MathVector<1> vLocIP[],
		                        const size_t nip,
		                        const MathMatrix<1, dim>* vJT = NULL) const
		{
			getImpl().template evaluate<1>(vValue,vGlobIP,time,si,u,elem,
										   vCornerCoords,vLocIP,nip, vJT);
		}

		virtual void operator()(TData vValue[],
		                        const MathVector<dim> vGlobIP[],
		                        number time, int si,
		                        LocalVector& u,
		                        GeometricObject* elem,
		                        const MathVector<dim> vCornerCoords[],
		                        const MathVector<2> vLocIP[],
		                        const size_t nip,
		                        const MathMatrix<2, dim>* vJT = NULL) const
		{
			getImpl().template evaluate<2>(vValue,vGlobIP,time,si,u,elem,
										   vCornerCoords,vLocIP,nip, vJT);
		}

		virtual void operator()(TData vValue[],
		                        const MathVector<dim> vGlobIP[],
		                        number time, int si,
		                        LocalVector& u,
		                        GeometricObject* elem,
		                        const MathVector<dim> vCornerCoords[],
		                        const MathVector<3> vLocIP[],
		                        const size_t nip,
		                        const MathMatrix<3, dim>* vJT = NULL) const
		{
			getImpl().template evaluate<3>(vValue,vGlobIP,time,si,u,elem,
										   vCornerCoords,vLocIP,nip, vJT);
		}

	///	\}
	protected:
	///	access to implementation
		TImpl& getImpl() {return static_cast<TImpl&>(*this);}

	///	const access to implementation
		const TImpl& getImpl() const {return static_cast<const TImpl&>(*this);}
};

} // namespace ug

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__STD_USER_DATA__ */
