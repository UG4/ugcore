/*
 * proj_gauss_seidel_impl.h
 *
 *  Created on: 10.10.2013
 *      Author: raphaelprohl
 */

#ifndef __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__PROJ_GAUSS_SEIDEL_IMPL__
#define __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__PROJ_GAUSS_SEIDEL_IMPL__

#include "proj_gauss_seidel.h"

#define PROFILE_PROJ_GS
#ifdef PROFILE_PROJ_GS
	#define PROJ_GS_PROFILE_FUNC()		PROFILE_FUNC()
	#define PROJ_GS_PROFILE_BEGIN(name)	PROFILE_BEGIN_GROUP(name, "projGS")
	#define PROJ_GS_PROFILE_END()		PROFILE_END()
#else
	#define PROJ_GS_PROFILE_FUNC()
	#define PROJ_GS_PROFILE_BEGIN(name)
	#define PROJ_GS_PROFILE_END()
#endif

namespace ug{

template <typename TDomain, typename TAlgebra>
bool
ProjGaussSeidel<TDomain, TAlgebra>::
init(SmartPtr<ILinearOperator<vector_type> > J, const vector_type& u)
{
	PROFILE_FUNC_GROUP("projGS");
	return false;
}

template <typename TDomain, typename TAlgebra>
bool
ProjGaussSeidel<TDomain, TAlgebra>::
init(SmartPtr<ILinearOperator<vector_type> > L)
{
	PROFILE_FUNC_GROUP("projGS");
	return false;
}

template <typename TDomain, typename TAlgebra>
bool
ProjGaussSeidel<TDomain, TAlgebra>::
apply(vector_type &c, const vector_type& d)
{
	PROFILE_FUNC_GROUP("projGS");
	return false;
}

template <typename TDomain, typename TAlgebra>
bool
ProjGaussSeidel<TDomain, TAlgebra>::
apply_update_defect(vector_type &c, vector_type& d)
{
	PROFILE_FUNC_GROUP("projGS");
	return false;
}

template <typename TDomain, typename TAlgebra>
SmartPtr<ILinearIterator<typename TAlgebra::vector_type> >
ProjGaussSeidel<TDomain, TAlgebra>::
clone()
{
	SmartPtr<ProjGaussSeidel<TDomain, TAlgebra> > clone(
		new ProjGaussSeidel<TDomain, TAlgebra>());

	return clone;
}

} // end namespace ug

#endif /* __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__PROJ_GAUSS_SEIDEL_IMPL__ */
