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
#include "lib_disc/reference_element/reference_mapping_provider.h"
#include "lib_disc/spatial_disc/disc_util/fvcr_geom.h"

namespace ug{

/// marks elements according to an attached error value field
/**
 * This function marks elements for refinement. The passed error attachment
 * is used as a weight for the amount of the error an each element. All elements
 * that have an indicated error with s* max <= err <= max are marked for refinement.
 * Here, max is the maximum error measured, s is a scaling quantity chosen by
 * the user. In addition, all elements with an error smaller than TOL
 * (user defined) are not refined.
 *
 * \param[in, out]	refiner		Refiner, elements marked on exit
 * \param[in]		dd			dof distribution
 * \param[in]		TOL			Minimum error, such that an element is marked
 * \param[in]		scale		scaling factor indicating lower bound for marking
 * \param[in]		aaError		Error value attachment to elements
 */
template <typename TElem>
void MarkElements(MultiGrid::AttachmentAccessor<TElem, ug::Attachment<number> >& aaError,
                  IRefiner& refiner,
                  ConstSmartPtr<DoFDistribution> dd,
                  number TOL,
                  number refineFrac, number coarseFrac,
				  int maxLevel)
{
	typedef typename DoFDistribution::traits<TElem>::const_iterator const_iterator;

//	reset maximum of error
	number max = 0.0, min = std::numeric_limits<number>::max();
	number totalErr = 0.0;
	int numElem = 0;

//	get element iterator
	const_iterator iter = dd->template begin<TElem>();
	const_iterator iterEnd = dd->template end<TElem>();

//	loop all elements to find the maximum of the error
	for( ;iter != iterEnd; ++iter)
	{
	//	get element
		TElem* elem = *iter;

	//	search for maximum and minimum
		if(aaError[elem] > max)	max = aaError[elem];
		if(aaError[elem] < min)	min = aaError[elem];

	//	sum up total error
		totalErr += aaError[elem];
		numElem += 1;
	}

#ifdef UG_PARALLEL
	number maxLocal = max, minLocal = min, totalErrLocal = totalErr;
	int numElemLocal = numElem;
	if(pcl::NumProcs() > 1){
		pcl::ProcessCommunicator com;
		max = com.allreduce(maxLocal, PCL_RO_MAX);
		min = com.allreduce(minLocal, PCL_RO_MIN);
		totalErr = com.allreduce(totalErrLocal, PCL_RO_SUM);
		numElem = com.allreduce(numElemLocal, PCL_RO_SUM);
	}
#endif
	UG_LOG("  +++++  Gradient Error Indicator on "<<numElem<<" Elements +++++\n");
#ifdef UG_PARALLEL
	if(pcl::NumProcs() > 1){
		UG_LOG("  +++ Element Errors on Proc " << pcl::ProcRank() <<
			   ": maximum=" << max << ", minimum="<<min<<", sum="<<totalErr<<".\n");
	}
#endif

	UG_LOG("  +++ Element Errors: maximum=" << max << ", minimum="<<min<<
	       ", sum="<<totalErr<<".\n");

//	check if total error is smaller than tolerance. If that is the case we're done
	if(totalErr < TOL)
	{
		UG_LOG("  +++ Total error "<<totalErr<<" smaller than TOL ("<<TOL<<"). done.");
		return;
	}

//	Compute minimum
	number minErrToRefine = max * refineFrac;
	UG_LOG("  +++ Refining elements if error greater " << refineFrac << "*" <<max<<
			" = "<< minErrToRefine << ".\n");
	number maxErrToCoarse = min * (1+coarseFrac);
	if(maxErrToCoarse < TOL/numElem) maxErrToCoarse = TOL/numElem;
	UG_LOG("  +++ Coarsening elements if error smaller "<< maxErrToCoarse << ".\n");

//	reset counter
	int numMarkedRefine = 0, numMarkedCoarse = 0;

	iter = dd->template begin<TElem>();
	iterEnd = dd->template end<TElem>();

//	loop elements for marking
	for(; iter != iterEnd; ++iter)
	{
	//	get element
		TElem* elem = *iter;

	//	marks for refinement
		if(aaError[elem] >= minErrToRefine)
			if(dd->multi_grid()->get_level(elem) <= maxLevel)
			{
				refiner.mark(elem, RM_REFINE);
				numMarkedRefine++;
			}

	//	marks for coarsening
		if(aaError[elem] <= maxErrToCoarse)
		{
			refiner.mark(elem, RM_COARSEN);
			numMarkedCoarse++;
		}
	}

#ifdef UG_PARALLEL
	if(pcl::NumProcs() > 1){
		UG_LOG("  +++ Marked for refinement on Proc "<<pcl::ProcRank()<<": " << numMarkedRefine << " Elements.\n");
		UG_LOG("  +++ Marked for coarsening on Proc "<<pcl::ProcRank()<<": " << numMarkedCoarse << " Elements.\n");
		pcl::ProcessCommunicator com;
		int numMarkedRefineLocal = numMarkedRefine, numMarkedCoarseLocal = numMarkedCoarse;
		numMarkedRefine = com.allreduce(numMarkedRefineLocal, PCL_RO_SUM);
		numMarkedCoarse = com.allreduce(numMarkedCoarseLocal, PCL_RO_SUM);
	}
#endif

	UG_LOG("  +++ Marked for refinement: " << numMarkedRefine << " Elements.\n");
	UG_LOG("  +++ Marked for coarsening: " << numMarkedCoarse << " Elements.\n");
}

/// marks elements according to an attached error value field
/**
 * This function marks elements for refinement. The passed error attachment
 * is used as a weight for the amount of the error an each element. All elements
 * that have an indicated error > refineTol are marked for refinement and
 * elements with an error < coarsenTol are marked for coarsening
 *
 * \param[in, out]	refiner		Refiner, elements marked on exit
 * \param[in]		dd			dof distribution
 * \param[in]		refTol		all elements with error > refTol are marked for refinement.
 * 								If refTol is negative, no element will be marked for refinement.
 * \param[in]		coarsenTol	all elements with error < coarsenTol are marked for coarsening.
 * 								If coarsenTol is negative, no element will be marked for coarsening.
 * \param[in]		aaError		Error value attachment to elements
 */
template <typename TElem>
void MarkElementsAbsolute(MultiGrid::AttachmentAccessor<TElem, ug::Attachment<number> >& aaError,
						  IRefiner& refiner,
						  ConstSmartPtr<DoFDistribution> dd,
						  number refTol,
						  number coarsenTol,
						  int minLevel,
						  int maxLevel)
{
	typedef typename DoFDistribution::traits<TElem>::const_iterator const_iterator;

	int numMarkedRefine = 0, numMarkedCoarse = 0;
	const_iterator iter = dd->template begin<TElem>();
	const_iterator iterEnd = dd->template end<TElem>();

//	loop elements for marking
	for(; iter != iterEnd; ++iter){
		TElem* elem = *iter;

	//	marks for refinement
		if((refTol >= 0)
			&& (aaError[elem] > refTol)
			&& (dd->multi_grid()->get_level(elem) < maxLevel))
		{
			refiner.mark(elem, RM_REFINE);
			numMarkedRefine++;
		}

	//	marks for coarsening
		if((coarsenTol >= 0)
			&& (aaError[elem] < coarsenTol)
			&& (dd->multi_grid()->get_level(elem) > minLevel))
		{
			refiner.mark(elem, RM_COARSEN);
			numMarkedCoarse++;
		}
	}

#ifdef UG_PARALLEL
	if(pcl::NumProcs() > 1){
		UG_LOG("  +++ Marked for refinement on Proc "<<pcl::ProcRank()<<": " << numMarkedRefine << " Elements.\n");
		UG_LOG("  +++ Marked for coarsening on Proc "<<pcl::ProcRank()<<": " << numMarkedCoarse << " Elements.\n");
		pcl::ProcessCommunicator com;
		int numMarkedRefineLocal = numMarkedRefine, numMarkedCoarseLocal = numMarkedCoarse;
		numMarkedRefine = com.allreduce(numMarkedRefineLocal, PCL_RO_SUM);
		numMarkedCoarse = com.allreduce(numMarkedCoarseLocal, PCL_RO_SUM);
	}
#endif

	UG_LOG("  +++ Marked for refinement: " << numMarkedRefine << " Elements.\n");
	UG_LOG("  +++ Marked for coarsening: " << numMarkedCoarse << " Elements.\n");
}


template <typename TFunction>
void ComputeGradientLagrange1(TFunction& u, size_t fct,
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
		const LocalShapeFunctionSet<dim>& lsfs =
				LocalFiniteElementProvider::get<dim>(roid, LFEID(LFEID::LAGRANGE, dim, 1));

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
			std::vector<DoFIndex> ind;
			u.inner_dof_indices(vert, fct, ind);

		//	scale global gradient
			vGlobalGrad[sh] *= DoFRef(u, ind[0]);

		//	sum up
			MidGrad += vGlobalGrad[sh];
		}

	//	write result in array storage
		aaError[elem] = VecTwoNorm(MidGrad) * pow(elemSize, 2./dim);
	}
}

template <typename TFunction>
void ComputeGradientCrouzeixRaviart(TFunction& u, size_t fct,
                     MultiGrid::AttachmentAccessor<
                     typename TFunction::element_type,
                     ug::Attachment<number> >& aaError)
{
	static const int dim = TFunction::dim;
	typedef typename TFunction::const_element_iterator const_iterator;
	typedef typename TFunction::domain_type domain_type;
	typedef typename domain_type::grid_type grid_type;
	typedef typename TFunction::element_type element_type;
	typedef typename element_type::side side_type;
	
	typename grid_type::template traits<side_type>::secure_container sides;

//	get position accessor
	typename TFunction::domain_type::position_accessor_type& aaPos
			= u.domain()->position_accessor();
			
	grid_type& grid = *u.domain()->grid();

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
	
	//  get sides of element		
		grid.associated_elements_sorted(sides, elem );

	//	reference object type
		ReferenceObjectID roid = elem->reference_object_id();

	//	get trial space
		const LocalShapeFunctionSet<dim>& lsfs =
				LocalFiniteElementProvider::get<dim>(roid, LFEID(LFEID::CROUZEIX_RAVIART, dim, 1));

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

		//	get of of vertex
			std::vector<DoFIndex> ind;
			u.inner_dof_indices(sides[sh], fct, ind);

		//	scale global gradient
			vGlobalGrad[sh] *= DoFRef(u, ind[0]);

		//	sum up
			MidGrad += vGlobalGrad[sh];
		}

	//	write result in array storage
		aaError[elem] = VecTwoNorm(MidGrad) * pow(elemSize, 2./dim);
	}
}

template <int dim> struct face_type_traits
{
    typedef void face_type0;
	typedef void face_type1;
};

template <> struct face_type_traits<1>
{
    typedef ReferenceVertex face_type0;
	typedef ReferenceVertex face_type1;
};

template <> struct face_type_traits<2>
{
    typedef ReferenceEdge face_type0;
	typedef ReferenceEdge face_type1;
};

template <> struct face_type_traits<3>
{
    typedef ReferenceTriangle face_type0;
	typedef ReferenceQuadrilateral face_type1;
};

template <typename TFunction>
void ComputeGradientPiecewiseConstant(TFunction& u, size_t fct,
                     MultiGrid::AttachmentAccessor<
                     typename TFunction::element_type,
                     ug::Attachment<number> >& aaError)
{
	static const int dim = TFunction::dim;
	typedef typename TFunction::const_element_iterator const_iterator;
	typedef typename TFunction::domain_type domain_type;
	typedef typename domain_type::grid_type grid_type;
	typedef typename TFunction::element_type element_type;
	typedef typename element_type::side side_type;
	
	typename grid_type::template traits<side_type>::secure_container sides;
	
	typedef typename face_type_traits<dim>::face_type0 face_type0;
	typedef typename face_type_traits<dim>::face_type1 face_type1;

//	get position accessor
	typename TFunction::domain_type::position_accessor_type& aaPos
			= u.domain()->position_accessor();
			
	grid_type& grid = *u.domain()->grid();

//	some storage
	MathMatrix<dim, dim> JTInv;
	
	std::vector<MathVector<dim> > vCorner;
	std::vector<MathVector<dim> > sideCorners;
	//MathVector<dim> sideCoPos[dim+1];
	std::vector<MathVector<dim> > vSideCoPos;

//	get iterator over elements
	const_iterator iter = u.template begin<element_type>();
	const_iterator iterEnd = u.template end<element_type>();

//	loop elements
	for(; iter != iterEnd; ++iter)
	{
	//	get the element
		element_type* elem = *iter;
		
		MathVector<dim> vGlobalGrad=0;
	
	//  get sides of element		
		grid.associated_elements_sorted(sides, elem );

	//	reference object type
		ReferenceObjectID roid = elem->reference_object_id();

		const DimReferenceElement<dim>& rRefElem
				= ReferenceElementProvider::get<dim>(roid);

	//	get corners of element
		CollectCornerCoordinates(vCorner, *elem, aaPos);

	//	compute size (volume) of element
		const number elemSize = ElementSize<dim>(roid, &vCorner[0]);
		
		typename grid_type::template traits<element_type>::secure_container assoElements;
		
	// assemble element-wise finite volume gradient
		for (size_t s=0;s<sides.size();s++){
			grid.associated_elements(assoElements,sides[s]);
			// face value is average of associated elements
			number faceValue = 0;
			size_t numOfAsso = assoElements.size();
			for (size_t i=0;i<numOfAsso;i++){
				std::vector<DoFIndex> ind;
				u.inner_dof_indices(assoElements[i], fct, ind);
				faceValue+=DoFRef(u, ind[0]);
			}
			faceValue/=(number)numOfAsso;
			MathVector<dim> normal;
			size_t numSideCo = rRefElem.num(dim-1,s,0);

			//alt:
			//for (size_t j=0;j<numSideCo;j++)
			//	sideCoPos[j] = vCorner[rRefElem.id(dim-1,s,0,j)];

			//neu:
			for (size_t i = 0; i < numSideCo; ++i)
				vSideCoPos.push_back(vCorner[rRefElem.id(dim-1, s, 0, i)]);

			// faces have dim corners in 1d, 2d
			// in 3d they have dim corners (triangle) or dim+1 corners (quadrilateral)
			if ((int)numSideCo==dim)
				//alt:
				//ElementNormal<face_type0,dim>(normal,sideCoPos);
				ElementNormal<face_type0,dim>(normal,&vSideCoPos[0]);
			else
				//ElementNormal<face_type1,dim>(normal,sideCoPos);
				ElementNormal<face_type0,dim>(normal,&vSideCoPos[0]);

			for (int d=0;d<dim;d++){
				vGlobalGrad[d] += faceValue * normal[d];
			}
		}
		vGlobalGrad/=(number)elemSize;

	//	write result in array storage
		aaError[elem] = VecTwoNorm(vGlobalGrad) * pow(elemSize, 2./dim);
	}
}


template <typename TDomain, typename TAlgebra>
void MarkForAdaption_GradientIndicator(IRefiner& refiner,
                                       GridFunction<TDomain, TAlgebra>& u,
                                       const char* fctName,
                                       number TOL,
                                       number refineFrac, number coarseFrac,
                                       int maxLevel)
{
//	types
	typedef GridFunction<TDomain, TAlgebra> TFunction;
	typedef typename TFunction::domain_type::grid_type grid_type;
	typedef typename TFunction::element_type element_type;
	const int dim = TFunction::dim;

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
	if (u.local_finite_element_id(fct) == LFEID(LFEID::LAGRANGE, dim, 1))
		ComputeGradientLagrange1(u, fct, aaError);
	else if (u.local_finite_element_id(fct) == LFEID(LFEID::CROUZEIX_RAVIART, dim, 1))
		ComputeGradientCrouzeixRaviart(u, fct, aaError);
	else if (u.local_finite_element_id(fct) == LFEID(LFEID::PIECEWISE_CONSTANT, dim, 0))
		ComputeGradientPiecewiseConstant(u,fct,aaError);
	else{
		UG_THROW("Non-supported finite element type " << u.local_finite_element_id(fct) << "\n");
	}
	
// 	Mark elements for refinement
	MarkElements<element_type> (aaError, refiner, u.dof_distribution(), TOL, refineFrac, coarseFrac, maxLevel);

// 	detach error field
	pMG->template detach_from<element_type>(aError);
};


template <typename TDomain, typename TAlgebra>
void MarkForAdaption_AbsoluteGradientIndicator(IRefiner& refiner,
                                       GridFunction<TDomain, TAlgebra>& u,
                                       const char* fctName,
                                       number refTol,
                                       number coarsenTol,
                                       int minLvl,
                                       int maxLevel)
{
//	types
	typedef GridFunction<TDomain, TAlgebra> TFunction;
	typedef typename TFunction::domain_type::grid_type grid_type;
	typedef typename TFunction::element_type element_type;
	const int dim = TFunction::dim;

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
	if (u.local_finite_element_id(fct) == LFEID(LFEID::LAGRANGE, dim, 1))
		ComputeGradientLagrange1(u, fct, aaError);
	else if (u.local_finite_element_id(fct) == LFEID(LFEID::CROUZEIX_RAVIART, dim, 1))
		ComputeGradientCrouzeixRaviart(u, fct, aaError);
	else if (u.local_finite_element_id(fct) == LFEID(LFEID::PIECEWISE_CONSTANT, dim, 0))
		ComputeGradientPiecewiseConstant(u,fct,aaError);
	else{
		UG_THROW("Non-supported finite element type " << u.local_finite_element_id(fct) << "\n");
	}

// 	Mark elements for refinement
	MarkElementsAbsolute<element_type> (aaError, refiner, u.dof_distribution (), refTol, coarsenTol, minLvl, maxLevel);

// 	detach error field
	pMG->template detach_from<element_type>(aError);
};


template <typename TFunction>
void computeGradientJump(TFunction& u,
                     MultiGrid::AttachmentAccessor<
                     typename TFunction::element_type,
                     ug::Attachment<number> >& aaGrad,
					 MultiGrid::AttachmentAccessor<
                     typename TFunction::element_type,
                     ug::Attachment<number> >& aaError)
{
	typedef typename TFunction::domain_type domain_type;
	typedef typename domain_type::grid_type grid_type;
	typedef typename TFunction::element_type element_type;
	typedef typename element_type::side side_type;
	typedef typename TFunction::template traits<side_type>::const_iterator side_iterator;
	
	grid_type& grid = *u.domain()->grid();

	//	get iterator over elements
	side_iterator iter = u.template begin<side_type>();
	side_iterator iterEnd = u.template end<side_type>();
	//	loop elements
	for(; iter != iterEnd; ++iter)
	{
		//	get the element
		side_type* side = *iter;
		typename grid_type::template traits<element_type>::secure_container neighElements;
		grid.associated_elements(neighElements,side);
		if (neighElements.size()!=2) continue;
		number localJump = std::abs(aaGrad[neighElements[0]]-aaGrad[neighElements[1]]);
		for (size_t i=0;i<2;i++)
			if (aaError[neighElements[i]]<localJump) aaError[neighElements[i]]=localJump;
	}
}

// indicator is based on elementwise computation of jump between elementwise gradients
// the value in an element is the maximum jump between the element gradient and gradient
// in elements with common face
template <typename TDomain, typename TAlgebra>
void MarkForAdaption_GradientJumpIndicator(IRefiner& refiner,
                                       GridFunction<TDomain, TAlgebra>& u,
                                       const char* fctName,
                                       number TOL,
                                       number refineFrac, number coarseFrac,
                                       int maxLevel)
{
//	types
	typedef GridFunction<TDomain, TAlgebra> TFunction;
	typedef typename TFunction::domain_type::grid_type grid_type;
	typedef typename TFunction::element_type element_type;
	const int dim = TFunction::dim;

//	function id
	const size_t fct = u.fct_id_by_name(fctName);

//	get multigrid
	SmartPtr<grid_type> pMG = u.domain()->grid();

// 	attach error field
	typedef Attachment<number> ANumber;
	ANumber aGrad;
	pMG->template attach_to<element_type>(aGrad);
	MultiGrid::AttachmentAccessor<element_type, ANumber> aaGrad(*pMG, aGrad);
	
	ANumber aError;
	pMG->template attach_to<element_type>(aError);
	MultiGrid::AttachmentAccessor<element_type, ANumber> aaError(*pMG, aError);

// 	Compute error on elements
	if (u.local_finite_element_id(fct) == LFEID(LFEID::LAGRANGE, dim, 1))
		ComputeGradientLagrange1(u, fct, aaGrad);
	else if (u.local_finite_element_id(fct) == LFEID(LFEID::CROUZEIX_RAVIART, dim, 1))
		ComputeGradientCrouzeixRaviart(u, fct, aaGrad);
	else if (u.local_finite_element_id(fct) == LFEID(LFEID::PIECEWISE_CONSTANT, dim, 0))
		ComputeGradientPiecewiseConstant(u,fct,aaGrad);
	else{
		UG_THROW("Non-supported finite element type " << u.local_finite_element_id(fct) << "\n");
	}
	computeGradientJump(u,aaGrad,aaError);
	
// 	Mark elements for refinement
	MarkElements<element_type> (aaError, refiner, u.dof_distribution(), TOL, refineFrac, coarseFrac, maxLevel);

// 	detach error field
	pMG->template detach_from<element_type>(aError);
};

template <typename TDomain, typename TAlgebra>
void MarkForAdaption_AbsoluteGradientJumpIndicator(IRefiner& refiner,
                                       GridFunction<TDomain, TAlgebra>& u,
                                       const char* fctName,
                                       number refTol,
                                       number coarsenTol,
                                       int minLvl,
                                       int maxLevel)
{
//	types
	typedef GridFunction<TDomain, TAlgebra> TFunction;
	typedef typename TFunction::domain_type::grid_type grid_type;
	typedef typename TFunction::element_type element_type;
	const int dim = TFunction::dim;

//	function id
	const size_t fct = u.fct_id_by_name(fctName);

//	get multigrid
	SmartPtr<grid_type> pMG = u.domain()->grid();

// 	attach error field
	typedef Attachment<number> ANumber;
	ANumber aGrad;
	pMG->template attach_to<element_type>(aGrad);
	MultiGrid::AttachmentAccessor<element_type, ANumber> aaGrad(*pMG, aGrad);
	
	ANumber aError;
	pMG->template attach_to<element_type>(aError);
	MultiGrid::AttachmentAccessor<element_type, ANumber> aaError(*pMG, aError);

// 	Compute error on elements
	if (u.local_finite_element_id(fct) == LFEID(LFEID::LAGRANGE, dim, 1))
		ComputeGradientLagrange1(u, fct, aaGrad);
	else if (u.local_finite_element_id(fct) == LFEID(LFEID::CROUZEIX_RAVIART, dim, 1))
		ComputeGradientCrouzeixRaviart(u, fct, aaGrad);
	else if (u.local_finite_element_id(fct) == LFEID(LFEID::PIECEWISE_CONSTANT, dim, 0))
		ComputeGradientPiecewiseConstant(u,fct,aaGrad);
	else{
		UG_THROW("Non-supported finite element type " << u.local_finite_element_id(fct) << "\n");
	}
	computeGradientJump(u,aaGrad,aaError);
	
// 	Mark elements for refinement
	MarkElementsAbsolute<element_type> (aaError, refiner, u.dof_distribution(), refTol, coarsenTol, minLvl, maxLevel);

// 	detach error field
	pMG->template detach_from<element_type>(aError);
};

} // end namespace ug

#endif /* __H__UG__LIB_DISC__FUNCTION_SPACE__ERROR_INDICATOR__ */
