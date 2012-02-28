/*
 * sqp_impl.h
 *
 *  Created on: 11.01.2012
 *      Author: Raphael Prohl
 */

#ifndef __H__UG__LIB_DISC__OPERATOR__NON_LINEAR_OPERATOR__SQP_METHOD__SQP_IMPL__
#define __H__UG__LIB_DISC__OPERATOR__NON_LINEAR_OPERATOR__SQP_METHOD__SQP_IMPL__

#include <iostream>
#include <sstream>

#include "sqp.h"
#include "lib_disc/function_spaces/grid_function_util.h"

#include "lib_disc/common/groups_util.h"
#include "lib_disc/spatial_disc/domain_disc.h"
#include "sqp_elem_util.h"

#define PROFILE_SQP
#ifdef PROFILE_SQP
	#define SQP_PROFILE_FUNC()		PROFILE_FUNC()
	#define SQP_PROFILE_BEGIN(name)	PROFILE_BEGIN(name)
	#define SQP_PROFILE_END()		PROFILE_END()
#else
	#define SQP_PROFILE_FUNC()
	#define SQP_PROFILE_BEGIN(name)
	#define SQP_PROFILE_END()
#endif

namespace ug{

template <typename TDomain, typename TAlgebra>
bool
SQPMethod<TDomain, TAlgebra>::
init()
{
	return true;
}


template <typename TDomain, typename TAlgebra>
bool
SQPMethod<TDomain, TAlgebra>::
prepare()
{
//	Check if Tolerance Check has been set
	if(m_toleranceCheck <= 0.0)
	{
		UG_LOG("ERROR in 'SQPMethod::prepare': Tolerance Check"
				"must be set to a valid value.\n");
		return false;
	}

	return true;
}

template <typename TDomain, typename TAlgebra>
bool
SQPMethod<TDomain, TAlgebra>::
check_tolerance(const vector_type& u, const SurfaceDoFDistribution& dd, base_type& domDisc)
{

	//soll auf ElemDisc (FE1NonlinearElasticity) - Methode zugreifen,
	//die elementweise über die attached Data läuft

	//typedef typename DomainDiscretization<TDomain,TDoFDistribution,TAlgebra> domDisc;
	//DomainDiscretization<TDomain,TDoFDistribution,TAlgebra> domDisc;

//	update the elem discs
	//if(!DomainDiscretization<TDomain,TDoFDistributxion,TAlgebra>::update_disc_items())
	if(!domDisc.update_disc_items())
	//	if(!update_disc_items())
		UG_THROW_FATAL("Cannot update disc items.");

//	Union of Subsets
	SubsetGroup unionSubsets;
	std::vector<SubsetGroup> vSSGrp;

//	create list of all subsets
	const ISubsetHandler& sh = dd.function_pattern().subset_handler();
	if(!CreateSubsetGroups(vSSGrp, unionSubsets, domDisc.m_vElemDisc, sh))
			UG_THROW_FATAL("ERROR in 'SQPMethod':"
				" Can not Subset Groups and union.\n");

//	loop subsets
	for(size_t i = 0; i < unionSubsets.num_subsets(); ++i)
	{
	//	get subset
		const int si = unionSubsets[i];

	//	get dimension of the subset
		const int dim = unionSubsets.dim(i);

	//	request if subset is regular grid
		bool bNonRegularGrid = !unionSubsets.regular_grid(i);

	//	overrule by regular grid if required
		if(domDisc.m_bForceRegGrid) bNonRegularGrid = false;

	//	Elem Disc on the subset
		std::vector<IElemDisc*> vSubsetElemDisc;

	//	get all element discretizations that work on the subset
		GetElemDiscOnSubset(vSubsetElemDisc, domDisc.m_vElemDisc, vSSGrp, si);

	//	success flag
		bool bSuc = true;

	//	assemble on suitable elements
		switch(dim)
		{
		case 1:
			bSuc &= CheckTolerance<Edge,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, u);
			break;
		case 2:
			bSuc &= CheckTolerance<Triangle,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, u);
			bSuc &= CheckTolerance<Quadrilateral,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, u);
			break;
		case 3:
			bSuc &= CheckTolerance<Tetrahedron,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, u);
			bSuc &= CheckTolerance<Pyramid,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, u);
			bSuc &= CheckTolerance<Prism,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, u);
			bSuc &= CheckTolerance<Hexahedron,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, u);
			break;
		default:
			UG_THROW_FATAL("ERROR in 'SQPMethod::check_tolerance':"
						"Dimension "<<dim<<" (subset="<<si<<") not supported.\n");
		}

	//	check success
		if(!bSuc)
			UG_THROW_FATAL("ERROR in 'SQPMethod::check_tolerance':"
							" Assembling of elements of Dimension " << dim << " in "
							" subset "<<si<< " failed.\n");
	}

	return true;
}

template <typename TDomain, typename TAlgebra>
bool
SQPMethod<TDomain, TAlgebra>::
update_variables(const vector_type& u, const SurfaceDoFDistribution& dd, base_type& domDisc)
{

//soll auf ElemDisc (FE1NonlinearElasticity) - Methode zugreifen,
//die elementweise über die attached Data läuft

//	update the elem discs
	if(!domDisc.update_disc_items())
		UG_THROW_FATAL("Cannot update disc items.");

	//typedef typename DomainDiscretization<TDomain,TDoFDistribution,TAlgebra> domDisc;
	//DomainDiscretization<TDomain,TDoFDistribution,TAlgebra> domDisc;

//	Union of Subsets
	SubsetGroup unionSubsets;
	std::vector<SubsetGroup> vSSGrp;

//	create list of all subsets
	const ISubsetHandler& sh = dd.get_function_pattern().subset_handler();
	if(!CreateSubsetGroups(vSSGrp, unionSubsets, domDisc.m_vElemDisc, sh))
		UG_THROW_FATAL("ERROR in 'SQPMethod':"
				" Can not Subset Groups and union.\n");

//	loop subsets
	for(size_t i = 0; i < unionSubsets.num_subsets(); ++i)
	{
	//	get subset
		const int si = unionSubsets[i];

	//	get dimension of the subset
		const int dim = unionSubsets.dim(i);

	//	request if subset is regular grid
		bool bNonRegularGrid = !unionSubsets.regular_grid(i);

	//	overrule by regular grid if required
		if(domDisc.m_bForceRegGrid) bNonRegularGrid = false;

	//	Elem Disc on the subset
		std::vector<IElemDisc*> vSubsetElemDisc;

	//	get all element discretizations that work on the subset
		GetElemDiscOnSubset(vSubsetElemDisc, domDisc.m_vElemDisc, vSSGrp, si);

	//	success flag
		bool bSuc = true;

	//	assemble on suitable elements
		switch(dim)
		{
		case 1:
			bSuc &= UpdateSQPVariables<Edge,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, u);
			break;
		case 2:
			bSuc &= UpdateSQPVariables<Triangle,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, u);
			bSuc &= UpdateSQPVariables<Quadrilateral,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, u);
			break;
		case 3:
			bSuc &= UpdateSQPVariables<Tetrahedron,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, u);
			bSuc &= UpdateSQPVariables<Pyramid,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, u);
			bSuc &= UpdateSQPVariables<Prism,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, u);
			bSuc &= UpdateSQPVariables<Hexahedron,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, u);
			break;
		default:
			UG_THROW_FATAL("ERROR in 'SQPMethod::update_variables':"
						"Dimension "<<dim<<" (subset="<<si<<") not supported.\n");
		}

	//	check success
		if(!bSuc)
			UG_THROW_FATAL("ERROR in 'SQPMethod::update_variables':"
							" Assembling of elements of Dimension " << dim << " in "
							" subset "<<si<< " failed.\n");
	}

	return true;
}

}

#endif /* __H__UG__LIB_DISC__OPERATOR__NON_LINEAR_OPERATOR__SQP_METHOD__SQP_IMPL__ */

