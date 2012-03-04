/*
 * error_indicator.h
 *
 *  Created on: 30.03.2011
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__FUNCTION_SPACE__ERROR_INDICATOR__
#define __H__UG__LIB_DISC__FUNCTION_SPACE__ERROR_INDICATOR__

#include <vector>

#include "common/common.h"
#include "common/util/provider.h"
#include "lib_grid/algorithms/refinement/refiner_interface.h"
#include "lib_disc/common/geometry_util.h"
#include "lib_disc/reference_element/reference_element_util.h"

namespace ug{

template <typename TFunction>
void ComputeGradient(TFunction& u, size_t fct,
                     MultiGrid::AttachmentAccessor<
                     typename TFunction::element_type,
                     ug::Attachment<number> >& aaError)
{
	static const int dim = TFunction::dim;
	typedef typename TFunction::const_element_iterator const_iterator;
	typedef typename TFunction::element_type element_type;

//	get position accessor
	typename TFunction::domain_type::position_accessor_type& aaPos
			= u.domain()->position_accessor();

//	some storage
	MathMatrix<dim, dim> JTInv;
	std::vector<MathVector<dim> > vLocalGrad;
	std::vector<MathVector<dim> > vGlobalGrad;
	std::vector<MathVector<dim> > vCorner;

//	get iterator over elements
	const_iterator iter = u.template begin<element_type>();
	const_iterator iterEnd = u.template end<element_type>();

//	loop elements
	for(; iter != iterEnd; ++iter)
	{
	//	get the element
		element_type* elem = *iter;

	//	reference object type
		ReferenceObjectID roid = elem->reference_object_id();

	//	get trial space
		const DimLocalShapeFunctionSet<dim>& lsfs =
				LocalShapeFunctionSetProvider::get<dim>(roid, LFEID(LFEID::LAGRANGE, 1));

	//	create a reference mapping
		DimReferenceMapping<dim, dim>& map
			= ReferenceMappingProvider::get<dim, dim>(roid);

	//	get local Mid Point
		MathVector<dim> localIP = ReferenceElementCenter<dim>(roid);

	//	number of shape functions
		const size_t numSH = lsfs.num_sh();
		vLocalGrad.resize(numSH);
		vGlobalGrad.resize(numSH);

	//	evaluate reference gradient at local midpoint
		lsfs.grads(&vLocalGrad[0], localIP);

	//	get corners of element
		CollectCornerCoordinates(vCorner, *elem, aaPos);

	//	update mapping
		map.update(&vCorner[0]);

	//	compute jacobian
		map.jacobian_transposed_inverse(JTInv, localIP);

	//	compute size (volume) of element
		const number elemSize = ElementSize<dim>(roid, &vCorner[0]);

	//	compute gradient at mid point by summing contributions of all shape fct
		MathVector<dim> MidGrad; VecSet(MidGrad, 0.0);
		for(size_t sh = 0 ; sh < numSH; ++sh)
		{
		//	get global Gradient
			MatVecMult(vGlobalGrad[sh], JTInv, vLocalGrad[sh]);

		//	get vertex
			VertexBase* vert = elem->vertex(sh);

		//	get of of vertex
			std::vector<MultiIndex<2> > ind;
			u.inner_multi_indices(vert, fct, ind);

		//	scale global gradient
			vGlobalGrad[sh] *= BlockRef(u[ind[0][0]], ind[0][1]);

		//	sum up
			MidGrad += vGlobalGrad[sh];
		}

	//	write result in array storage
		aaError[elem] = VecTwoNorm(MidGrad) * pow(elemSize, 2./dim);
	}
}

/// marks elements according to an attached error value field
/**
 * This function marks elements for refinement. The passed error attachment
 * is used as a weight for the amount of the error an each element. All elements
 * that have an indicated error with s* max <= err <= max are marked for refinement.
 * Here, max is the maximum error measured, s is a scaling quantity choosen by
 * the user. In addition, all elements with an error smaller than TOL
 * (user defined) are not refined.
 *
 * \param[in, out]	refiner		Refiner, elements marked on exit
 * \param[in]		u			Grid Function
 * \param[in]		TOL			Minimum error, such that an element is marked
 * \param[in]		scale		scaling factor indicating lower bound for marking
 * \param[in]		aaError		Error value attachment to elements
 */
template <typename TFunction>
void MarkElements(IRefiner& refiner,
                  TFunction& u,
                  number TOL, number scale,
                  MultiGrid::AttachmentAccessor<typename TFunction::element_type,
                  	  	  	  	  	  	  	    ug::Attachment<number> >& aaError)
{
	typedef typename TFunction::element_type element_type;
	typedef typename TFunction::const_element_iterator const_iterator;

//	reset maximum of error
	number max = 0.0;

//	get element iterator
	const_iterator iter = u.template begin<element_type>();
	const_iterator iterEnd = u.template end<element_type>();

//	loop all elements to find the maximum of the error
	for( ;iter != iterEnd; ++iter)
	{
	//	get element
		element_type* elem = *iter;

	//	search for maximum
		if(aaError[elem] > max)
			max = aaError[elem];
	}

	UG_LOG("  +++  Gradient Error Indicator  +++\n");
#ifdef UG_PARALLEL
	if(pcl::GetNumProcesses() > 1){
		pcl::ProcessCommunicator com;
		number maxLocal = max;
		com.allreduce(&maxLocal, &max, 1, PCL_DT_DOUBLE, PCL_RO_MAX);
		UG_LOG("  +++ Max Error on Proc " << pcl::GetProcRank() <<
			   " is " << maxLocal << ".\n");
	}
#endif

	UG_LOG("  +++ Max Error is " << max << ".\n");

//	check if something to do
	if(max <= TOL) return;

//	Compute minimum
	number min = max*scale;
	if(min < TOL) min = TOL;

	UG_LOG("  +++ Refining all elements with error >= " << min << ".\n");

//	reset counter
	int numMarked = 0;

	iter = u.template begin<element_type>();
	iterEnd = u.template end<element_type>();

//	loop elements for marking
	for(; iter != iterEnd; ++iter)
	{
	//	get element
		element_type* elem = *iter;

	//	check if element error is in range
		if(aaError[elem] >= min)
		{
		//	mark element and increase counter
			refiner.mark(elem);
			numMarked++;
		}
	}

#ifdef UG_PARALLEL
	if(pcl::GetNumProcesses() > 1){
		pcl::ProcessCommunicator com;
		int numMarkedLocal = numMarked;
		com.allreduce(&numMarkedLocal, &numMarked, 1, PCL_DT_INT, PCL_RO_SUM);
		UG_LOG("  +++ " << numMarkedLocal << " Elements marked on Proc "
			   << pcl::GetProcRank() << ".\n");
	}
#endif

	UG_LOG("  +++ " << numMarked << " Elements marked for refinement.\n");
}

template <typename TDomain, typename TDD, typename TAlgebra>
void MarkForRefinement_GradientIndicator(IRefiner& refiner,
#ifdef UG_PARALLEL
                                         ParallelGridFunction<GridFunction<TDomain, TDD, TAlgebra> >& u,
#else
                                         GridFunction<TDomain, TDD, TAlgebra>& u,
#endif
                                         const char* fctName,
                                         number TOL, number scale)
{
//	types
	typedef GridFunction<TDomain, TDD, TAlgebra> TFunction;
	typedef typename TFunction::domain_type::grid_type grid_type;
	typedef typename TFunction::element_type element_type;

//	function id
	const size_t fct = u.fct_id_by_name(fctName);

//	get multigrid
	SmartPtr<grid_type> pMG = u.domain()->grid();

// 	attach error field
	typedef Attachment<number> ANumber;
	ANumber aError;
	pMG->template attach_to<element_type>(aError);
	MultiGrid::AttachmentAccessor<element_type, ANumber> aaError(*pMG, aError);

// 	Compute error on elements
	ComputeGradient(u, fct, aaError);

// 	Mark elements for refinement
	MarkElements(refiner, u, TOL, scale, aaError);

// 	detach error field
	pMG->template detach_from<element_type>(aError);
};

} // end namespace ug

#endif /* __H__UG__LIB_DISC__FUNCTION_SPACE__ERROR_INDICATOR__ */
