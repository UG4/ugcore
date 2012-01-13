/*
 * sqp_elem_util.h
 *
 *  Created on: 12.01.2012
 *      Author: Raphael Prohl
 */

#ifndef __H__UG__LIB_DISC__OPERATOR__NON_LINEAR_OPERATOR__SQP_METHOD__SQP_ELEM_UTIL__
#define __H__UG__LIB_DISC__OPERATOR__NON_LINEAR_OPERATOR__SQP_METHOD__SQP_ELEM_UTIL__

// extern includes
#include <iostream>
#include <vector>

// other ug4 modules
#include "common/common.h"
#include "lib_grid/lg_base.h"

// intern headers
#include "lib_disc/reference_element/reference_element.h"
#include "lib_disc/common/function_group.h"
#include "lib_disc/common/local_algebra.h"
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h"
#include "lib_disc/spatial_disc/ip_data/data_evaluator.h"


#define PROFILE_SQP_ELEM
#ifdef PROFILE_SQP_ELEM
	#define SQP_ELEM_PROFILE_FUNC()		PROFILE_FUNC()
	#define SQP_ELEM_PROFILE_BEGIN(name)	PROFILE_BEGIN(name)
	#define SQP_ELEM_PROFILE_END()		PROFILE_END()
#else
	#define SQP_ELEM_PROFILE_FUNC()
	#define SQP_ELEM_PROFILE_BEGIN(name)
	#define SQP_ELEM_PROFILE_END()
#endif

namespace ug{

template <	typename TElem,
			typename TDoFDistribution,
			typename TAlgebra>
bool
CheckTolerance(const std::vector<IElemDisc*>& vElemDisc,
       	const IDoFDistribution<TDoFDistribution>& dofDistr,
       	int si, bool bNonRegularGrid,
       	const typename TAlgebra::vector_type& u)
{
	// 	check if at least on element exist, else return
		if(dofDistr.template num<TElem>(si) == 0) return true;

	//	get current time and vector
		//const number time = vSol.time(0);
		//const typename TAlgebra::vector_type& u = vSol.solution(0);

	//	create data evaluator
		DataEvaluator Eval;

	//	prepare for given elem discs
		if(!Eval.set_elem_discs(vElemDisc, dofDistr.get_function_pattern(), bNonRegularGrid, true))
		{
			UG_LOG("ERROR in 'CheckTolerance': "
					"Cannot init DataEvaluator with IElemDiscs.\n");
			return false;
		}

	//	set time-independent
	//	LocalVectorTimeSeries locTimeSeries; //time series of local vectors
	//	bool bNeedLocTimeSeries = Eval.set_time_dependent(true, time, &locTimeSeries);

	// 	local indices and local algebra
		LocalIndices ind; LocalVector locU;

	//	prepare element discs
		if(!Eval.template prepare_elem_loop<TElem>(ind, 0.0, true)) //time, true))
		{
			UG_LOG("ERROR in 'CheckTolerance': "
					"Cannot prepare element loop.\n");
			return false;
		}

	//	get element iterator
		typename geometry_traits<TElem>::const_iterator iter, iterBegin, iterEnd;
		iterBegin = dofDistr.template begin<TElem>(si);
		iterEnd = dofDistr.template end<TElem>(si);

		// 	Loop over all elements
		for(iter = iterBegin; iter != iterEnd; ++iter)
		{
		// 	get Element
			TElem* elem = *iter;

		// 	get global indices
			dofDistr.indices(elem, ind, Eval.use_hanging());

		// 	adapt local algebra
			locU.resize(ind);

		// 	read local values of u
			GetLocalVector(locU, u);

		//	read local values of time series
		//	if(bNeedLocTimeSeries)
		//	{
		//		locTimeSeries.read_values(vSol, ind);
		//		locTimeSeries.read_times(vSol);
		//	}

		// 	check tolerance
			if(!Eval.sqp_check_tolerance_elem(elem, locU))
			{
				UG_LOG("ERROR in 'CheckTolerance': "
						"Cannot check sqp tolerance.\n");
				return false;
			}

		}

	//	we're done
	return true;
}

template <	typename TElem,
			typename TDoFDistribution,
			typename TAlgebra>
bool
UpdateSQPVariables(const std::vector<IElemDisc*>& vElemDisc,
       	const IDoFDistribution<TDoFDistribution>& dofDistr,
       	int si, bool bNonRegularGrid,
       	const typename TAlgebra::vector_type& u)
{
	// 	check if at least on element exist, else return
		if(dofDistr.template num<TElem>(si) == 0) return true;

	//	get current time and vector
		//const number time = vSol.time(0);
		//const typename TAlgebra::vector_type& u = vSol.solution(0);

	//	create data evaluator
		DataEvaluator Eval;

	//	prepare for given elem discs
		if(!Eval.set_elem_discs(vElemDisc, dofDistr.get_function_pattern(), bNonRegularGrid, true))
		{
			UG_LOG("ERROR in 'UpdateSQPVariables': "
					"Cannot init DataEvaluator with IElemDiscs.\n");
			return false;
		}

	//	set time-independent
	//	LocalVectorTimeSeries locTimeSeries; //time series of local vectors
	//	bool bNeedLocTimeSeries = Eval.set_time_dependent(true, time, &locTimeSeries);

	// 	local indices and local algebra
		LocalIndices ind; LocalVector locU;

	//	prepare element discs
		if(!Eval.template prepare_elem_loop<TElem>(ind, 0.0, true)) //time, true))
		{
			UG_LOG("ERROR in 'UpdateSQPVariables': "
					"Cannot prepare element loop.\n");
			return false;
		}

	//	get element iterator
		typename geometry_traits<TElem>::const_iterator iter, iterBegin, iterEnd;
		iterBegin = dofDistr.template begin<TElem>(si);
		iterEnd = dofDistr.template end<TElem>(si);

		// 	Loop over all elements
		for(iter = iterBegin; iter != iterEnd; ++iter)
		{
		// 	get Element
			TElem* elem = *iter;

		// 	get global indices
			dofDistr.indices(elem, ind, Eval.use_hanging());

		// 	adapt local algebra
			locU.resize(ind);

		// 	read local values of u
			GetLocalVector(locU, u);

		//	read local values of time series
		//	if(bNeedLocTimeSeries)
		//	{
		//		locTimeSeries.read_values(vSol, ind);
		//		locTimeSeries.read_times(vSol);
		//	}

		// 	update variables
			if(!Eval.sqp_variables_update_elem(elem, locU))
			{
				UG_LOG("ERROR in 'UpdateSQPVariables': "
						"Cannot update sqp variables.\n");
				return false;
			}

		}

	//	we're done
	return true;
}


}

#endif /* __H__UG__LIB_DISC__OPERATOR__NON_LINEAR_OPERATOR__SQP_METHOD__SQP_ELEM_UTIL__ */

