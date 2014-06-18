/*
 * error_indicator.h
 *
 *  Created on: 30.03.2011
 *      Author: andreasvogel
 *		Additions by Sebastian Reiter
 */

#ifndef __H__UG__LIB_DISC__FUNCTION_SPACE__ERROR_INDICATOR__
#define __H__UG__LIB_DISC__FUNCTION_SPACE__ERROR_INDICATOR__

#include <vector>

#include "error_indicator_util.h"
#include "gradient_evaluators.h"
#include "common/common.h"
#include "common/util/provider.h"
#include "lib_grid/algorithms/refinement/refiner_interface.h"
#include "lib_disc/common/geometry_util.h"
#include "lib_disc/reference_element/reference_element_util.h"
#include "lib_disc/reference_element/reference_mapping_provider.h"
#include "lib_disc/spatial_disc/disc_util/fvcr_geom.h"
#include "integrate.h"
#include "common/profiler/profiler.h"

#ifdef UG_PARALLEL
 	#include "lib_grid/parallelization/util/compol_attachment_reduce.h"
#endif

namespace ug{


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
			Vertex* vert = elem->vertex(sh);

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

			for (size_t i = 0; i < numSideCo; ++i)
				vSideCoPos.push_back(vCorner[rRefElem.id(dim-1, s, 0, i)]);

			// faces have dim corners in 1d, 2d
			// in 3d they have dim corners (triangle) or dim+1 corners (quadrilateral)
			if ((int)numSideCo==dim)
				ElementNormal<face_type0,dim>(normal,&vSideCoPos[0]);
			else
				ElementNormal<face_type1,dim>(normal,&vSideCoPos[0]);

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
	PROFILE_FUNC();
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
	PROFILE_FUNC();
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
	PROFILE_FUNC();
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
	PROFILE_FUNC();
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

/**
 * \param[in]	maxL2Error	errors below maxL2Error are considered fine.
 */
template <typename TDomain, typename TAlgebra>
void MarkForAdaption_L2ErrorExact(IRefiner& refiner,
                                   SmartPtr<GridFunction<TDomain, TAlgebra> > u,
                                   SmartPtr<UserData<number, TDomain::dim> > spExactSol,
                                   const char* cmp,
                                   number minL2Error,
                                   number maxL2Error,
                                   number refFrac,
                                   int minLvl, int maxLvl,
                                   number time, int quadOrder)
{
	PROFILE_FUNC();
	using namespace std;
//	types
	typedef GridFunction<TDomain, TAlgebra> TFunction;
	typedef typename TFunction::domain_type::grid_type grid_t;
	typedef typename TFunction::element_type elem_t;
	const int dim = TFunction::dim;

//	function id
	const size_t fct = u->fct_id_by_name(cmp);
	UG_COND_THROW(fct >= u->num_fct(),
				  "Function space does not contain a function with name " << cmp);

//	get multigrid and position accessor
	grid_t& mg = *u->domain()->grid();

	typename TFunction::domain_type::position_accessor_type&
		aaPos = u->domain()->position_accessor();

// 	attach error field
	typedef Attachment<number> ANumber;
	ANumber aError;
	mg.template attach_to<elem_t>(aError);
	MultiGrid::AttachmentAccessor<elem_t, ANumber> aaError(mg, aError);

//	create integrand and perform integration
	SmartPtr<IIntegrand<number, dim> > spIntegrand
		= make_sp(new L2ErrorIntegrand<TFunction>(spExactSol, u, fct, time));

	number l2Error =
		Integrate<dim, dim>(u->template begin<elem_t>(), u->template end<elem_t>(),
				  			aaPos, *spIntegrand, quadOrder, "best", &aaError);

	#ifdef UG_PARALLEL
		pcl::ProcessCommunicator com;
		l2Error = com.allreduce(l2Error, PCL_RO_SUM);
	#endif

	l2Error = sqrt(l2Error);

	UG_LOG("maxError " << maxL2Error << ", l2Error " << l2Error << std::endl);

	if(l2Error > maxL2Error){
		typedef typename TFunction::template traits<elem_t>::const_iterator	ElemIter;
		size_t numElemsActive = 0;	// number of elements which may be refined
		size_t numElemsTotal = 0;	// total number of elements
		number maxElemError = 0;	// error in elements which may be refined
		number minElemError = numeric_limits<number>::max();
		number fixedError = 0;		// error in elements which can't be refined any further
		for(ElemIter iter = u->template begin<elem_t>(); iter != u->template end<elem_t>(); ++iter){
			++numElemsTotal;
			if(mg.get_level(*iter) < maxLvl){
				++numElemsActive;
				maxElemError = max(maxElemError, aaError[*iter]);
				minElemError = min(minElemError, aaError[*iter]);
			}
			else{
				fixedError += aaError[*iter];
			}
		}

		#ifdef UG_PARALLEL
			pcl::ProcessCommunicator com;
			minElemError = com.allreduce(minElemError, PCL_RO_MIN);
			maxElemError = com.allreduce(maxElemError, PCL_RO_MAX);
			fixedError = com.allreduce(fixedError, PCL_RO_MAX);
		#endif

	//	note that aaError contains error-squares
		maxElemError = minElemError + (maxElemError - minElemError) * sq(refFrac);
		number refThreshold = maxElemError;
		if(maxL2Error > 0){
		//	here we try to reduce the number of elements marked in the last steps.
		//	note that each element stores error-squares...
			//refThreshold = max(maxElemError, (sq(maxL2Error) - fixedError) / (number)numElemsActive);
			//refThreshold = max(maxElemError, sq(maxL2Error)/ (number)numElemsTotal);
		}

		MarkElementsAbsolute(aaError, refiner, u->dof_distribution(), refThreshold, -1,
						 minLvl, maxLvl);
	}

//	coarsening
	// number coarsenThreshold = -1;
	// if(minL2Error > 0)
	// 	coarsenThreshold = minL2Error * minL2Error / (number)numElemsActive;

// 	detach error field
	mg.template detach_from<elem_t>(aError);
};


/**	Exchages error-values and nbr-element-numbers between distributed sides and
 * adjusts those values in rim-shadows and rim-shadowing elements on the fly
 * (e.g. constrained and constraining elements).
 *
 * \param[in] u				The grid-function on which error-indication is based
 * \param[in] aSideError	ANumber attachment on side_t:
 *							The total error accumulated in element-sides.
 * \param[in] aNumElems		ANumber attachment on side_t:
 *							The number of elements which have contributed to aSideError.
 */
template <class side_t, class TFunction>
void ExchangeAndAdjustSideErrors(TFunction& u, ANumber aSideError, ANumber aNumElems)
{
	typedef typename TFunction::template traits<side_t>::const_iterator side_iter_t;
	typedef typename TFunction::domain_type::grid_type	grid_t;

	grid_t& g = *u.domain()->grid();
	Grid::AttachmentAccessor<side_t, ANumber> aaSideError(g, aSideError);
	Grid::AttachmentAccessor<side_t, ANumber> aaNumElems(g, aNumElems);

//	in a parallel environment we now have to sum the gradients over parallel
//	interfaces
//	since there may be constrained sides which locally do not have a constraining
//	side we do this before adding the constrained values to their constraining objects.
	#ifdef UG_PARALLEL
		typedef typename GridLayoutMap::Types<side_t>::Layout layout_t;
		DistributedGridManager& dgm = *g.distributed_grid_manager();
		GridLayoutMap& glm = dgm.grid_layout_map();
		pcl::InterfaceCommunicator<layout_t > icom;

	//	sum all copies at the h-master attachment
		ComPol_AttachmentReduce<layout_t, ANumber> compolSumErr(g, aSideError, PCL_RO_SUM);
		ComPol_AttachmentReduce<layout_t, ANumber> compolSumNum(g, aNumElems, PCL_RO_SUM);
		icom.exchange_data(glm, INT_H_SLAVE, INT_H_MASTER, compolSumErr);
		icom.exchange_data(glm, INT_H_SLAVE, INT_H_MASTER, compolSumNum);
		icom.communicate();

	//	copy the sum from the master to all of its slave-copies
		ComPol_CopyAttachment<layout_t, ANumber> compolCopyErr(g, aSideError);
		ComPol_CopyAttachment<layout_t, ANumber> compolCopyNum(g, aNumElems);
		icom.exchange_data(glm, INT_H_MASTER, INT_H_SLAVE, compolCopyErr);
		icom.exchange_data(glm, INT_H_MASTER, INT_H_SLAVE, compolCopyNum);
		icom.communicate();

	//	since we're copying from vmasters to vslaves later on, we have to make
	//	sure, that also all v-masters contain the correct values.
	//todo: communication to vmasters may not be necessary here...
	//		it is performed to make sure that all surface-rim-children
	//		contain their true value.
		icom.exchange_data(glm, INT_V_SLAVE, INT_V_MASTER, compolCopyErr);
		icom.exchange_data(glm, INT_V_SLAVE, INT_V_MASTER, compolCopyNum);

		icom.communicate();
	#endif

//	if we're currently operating on a surface function, we have to adjust
//	errors between shadowed and shadowing sides (constraining and constrained sides).
	if(u.grid_level().type() == GridLevel::SURFACE){
	//todo: avoid iteration over the whole grid!
	//todo: reduce communication
		const SurfaceView* surfView = u.approx_space()->surface_view().get();

		// for(side_iter_t iter = u.template begin<side_t>(SurfaceView::SHADOW_RIM);
		// 	iter != u.template end<side_t>(SurfaceView::SHADOW_RIM); ++iter)
		typedef typename Grid::traits<side_t>::iterator grid_side_iter_t;
		for(grid_side_iter_t iter = g.template begin<side_t>();
			iter != g.template end<side_t>(); ++iter)
		{
			side_t* s = *iter;
			if(!surfView->surface_state(s).partially_contains(SurfaceView::SHADOW_RIM))
				continue;

			const size_t numChildren = g.template num_children<side_t>(s);
			if(numChildren > 0){
				number w = 1. / number(numChildren);
				for(size_t i_child = 0; i_child < numChildren; ++i_child){
					side_t* c = g.template get_child<side_t>(s, i_child);
					aaSideError[s] += w * aaSideError[c];
					aaNumElems[s] += w * aaNumElems[c];
				}
			}
		}
	//	those elems which have local children now contain their true value (vslaves...)
		#ifdef UG_PARALLEL
			icom.exchange_data(glm, INT_V_SLAVE, INT_V_MASTER, compolCopyErr);
			icom.exchange_data(glm, INT_V_SLAVE, INT_V_MASTER, compolCopyNum);
			icom.communicate();
		#endif

		for(grid_side_iter_t iter = g.template begin<side_t>();
			iter != g.template end<side_t>(); ++iter)
		{
			side_t* s = *iter;
			if(!surfView->surface_state(s).partially_contains(SurfaceView::SHADOW_RIM))
				continue;

			const size_t numChildren = g.template num_children<side_t>(s);
			if(numChildren > 0){
				number w = 1. / number(numChildren);
				for(size_t i_child = 0; i_child < numChildren; ++i_child){
					side_t* c = g.template get_child<side_t>(s, i_child);
					aaSideError[c] = w * aaSideError[s];
					aaNumElems[c] = w * aaNumElems[s];
				}
			}
		}

	//	those elems which have a local parent now contain their true value (vmasters...)
		#ifdef UG_PARALLEL
		//	copy from v-masters to v-slaves, since there may be constrained
		//	sides which locally do not have a constraining element. Note that
		//	constrained V-Masters always have a local constraining element and thus
		//	contain the correct value.
			icom.exchange_data(glm, INT_V_MASTER, INT_V_SLAVE, compolCopyErr);
			icom.exchange_data(glm, INT_V_MASTER, INT_V_SLAVE, compolCopyNum);
			icom.communicate();
		#endif
	}
}

/** Calculates the jump between normal components of element-gradients on element sides,
 * then calculates the integral over those jumps on each side and finally sums those
 * integrals into aaError.*/
template <typename TGradientEvaluator, typename TFunction>
void EvaluateGradientJump_SideIntegral(TFunction& u, size_t fct,
					                   MultiGrid::AttachmentAccessor<
					                   		typename TFunction::element_type,
					                     	ug::Attachment<number> >& aaError,
					                   bool addErrSquareToAAError = false)
{
	using std::max;

	static const int dim = TFunction::dim;
	typedef typename TFunction::domain_type::grid_type	grid_t;
	typedef typename TFunction::const_element_iterator const_iterator;
	typedef typename TFunction::element_type elem_t;
	typedef typename elem_t::side side_t;
	typedef MathVector<dim>	vector_t;

//	get position accessor
	typename TFunction::domain_type::position_accessor_type& aaPos
			= u.domain()->position_accessor();

//	we need an attachment on the sides to gather gradient-jumps in parallel.
//	Note that a serial version could be created without this additional attachment.
	grid_t& g = *u.domain()->grid();
	ANumber aSideError;
	ANumber aNumElems;
	g.template attach_to_dv<side_t>(aSideError, 0);
	g.template attach_to_dv<side_t>(aNumElems, 0);
	Grid::AttachmentAccessor<side_t, ANumber> aaSideError(g, aSideError);
	Grid::AttachmentAccessor<side_t, ANumber> aaNumElems(g, aNumElems);

//	loop elements
	TGradientEvaluator gradEvaluator(&u, fct);
	const_iterator iterEnd = u.template end<elem_t>();
	for(const_iterator iter = u.template begin<elem_t>();
		iter != iterEnd; ++iter)
	{
	//	get the element
		elem_t* elem = *iter;
		vector_t elemGrad = gradEvaluator.evaluate(elem);

	//	iterate over all sides and add the normal component of the vector
	//	to the sides normGrad attachment
		for(size_t i = 0; i < elem->num_sides(); ++i){
			vector_t outerNorm = CalculateOuterNormal(elem, i, aaPos);
			number ng = VecDot(elemGrad, outerNorm);
			side_t* s = g.get_side(elem, i);
			aaSideError[s] += ng;
			++aaNumElems[s];
		}
	}

	ExchangeAndAdjustSideErrors<side_t>(u, aSideError, aNumElems);

//	now let's iterate over all elements again, this time integrating the jump
//	over the side and summing everything in aaError. Note that each element receives
//	50% of a sides error
	typename Grid::traits<side_t>::secure_container	sides;
	Grid::edge_traits::secure_container edges;
	for(const_iterator iter = u.template begin<elem_t>();
		iter != iterEnd; ++iter)
	{
		elem_t* elem = *iter;
		number err = 0;
		g.associated_elements(sides, elem);
		for(size_t i = 0; i < sides.size(); ++i){
			side_t* s = sides[i];
			if(aaNumElems[s] > 1){
				g.associated_elements(edges, s);
				number a = CalculateVolume(s, aaPos);
				number hs = ElementDiameter(g, aaPos, s);
				err += hs * sq(a * aaSideError[s]);
			}
		}

		if(addErrSquareToAAError)
			aaError[elem] += 0.5 * err;
		else
			aaError[elem] = sqrt(0.5 * err);
	}

	g.template detach_from<side_t>(aSideError);
	g.template detach_from<side_t>(aNumElems);
}


/**	calculates the length of the gradient in each element and then
 * fills aaError depending on the difference in length between neighbored elements.*/
template <typename TGradientEvaluator, typename TFunction>
void EvaluateGradientJump_Norm(TFunction& u, size_t fct,
			                     MultiGrid::AttachmentAccessor<
				                     typename TFunction::element_type,
				                     ug::Attachment<number> >& aaError)
{
	using std::max;

	static const int dim = TFunction::dim;
	typedef typename TFunction::domain_type::grid_type	grid_t;
	typedef typename TFunction::const_element_iterator const_iterator;
	typedef typename TFunction::element_type elem_t;
	typedef typename elem_t::side side_t;
	typedef MathVector<dim>	vector_t;

//	get position accessor
	typename TFunction::domain_type::position_accessor_type& aaPos
			= u.domain()->position_accessor();

//	some storage
	typename Grid::traits<side_t>::secure_container sides;

//	we need attachments on the sides to gather gradient-jumps in parallel.
//	Note that a serial version could be created without this additional attachment.
	grid_t& g = *u.domain()->grid();
	ANumber aSideError;
	ANumber aNumElems;
	g.template attach_to_dv<side_t>(aSideError, 0);
	g.template attach_to_dv<side_t>(aNumElems, 0);
	Grid::AttachmentAccessor<side_t, ANumber> aaSideError(g, aSideError);
	Grid::AttachmentAccessor<side_t, ANumber> aaNumElems(g, aNumElems);

//	loop elements and evaluate gradient
	TGradientEvaluator gradEvaluator(&u, fct);
	const_iterator iterEnd = u.template end<elem_t>();
	for(const_iterator iter = u.template begin<elem_t>();
		iter != iterEnd; ++iter)
	{
	//	get the element
		elem_t* elem = *iter;
		vector_t elemGrad = gradEvaluator.evaluate(elem);

		number elemContrib = VecTwoNorm(elemGrad);
		aaError[elem] = elemContrib;
		g.associated_elements(sides, elem);
		for(size_t i = 0; i < sides.size(); ++i){
			side_t* s = sides[i];
			aaSideError[s] += elemContrib;
			++aaNumElems[s];
		}
	}

	ExchangeAndAdjustSideErrors<side_t>(u, aSideError, aNumElems);

//	finally store the highest difference between an element-error and
//	averaged errors in associated sides in each element-error.
	for(const_iterator iter = u.template begin<elem_t>();
		iter != iterEnd; ++iter)
	{
		elem_t* elem = *iter;
		const number elemErr = aaError[elem];
		number err = 0;
		g.associated_elements(sides, elem);
		for(size_t i = 0; i < sides.size(); ++i){
			side_t* s = sides[i];
			if(aaNumElems[s] > 0)
				err = max(err, fabs(elemErr - aaSideError[s] / aaNumElems[s]));
		}
		//aaError[elem] = 2. * err * pow(CalculateVolume(elem, aaPos), number(dim-1)/number(dim));
		aaError[elem] = 2. * err * pow(CalculateVolume(elem, aaPos), 2./number(dim));
		//aaError[elem] = 2. * err * pow(CalculateVolume(elem, aaPos), 2./ (1. + 0.5*number(dim)));
		//aaError[elem] = 2. * err * CalculateVolume(elem, aaPos);
	}

	g.template detach_from<side_t>(aSideError);
	g.template detach_from<side_t>(aNumElems);
}

/** This gradient jump error indicator runs in parallel environments and
 * works with h-nodes.
 *
 * \param[in]	refFrac		given minElemError and maxElemError, all elements with
 *							an error > minElemError+refFrac*(maxElemError-minElemError)
 *							will be refined.
 * \param[in]	jumpType	defines the type of jump that shall be evaluated:
 *								- norm: the difference between gradient norms on
 *										neighboring elements is regarded.
 *								- sideInt:	The integral over the normal component
 *											of the gradient on each side is regarded.
 *
 * \todo: Add coarsenFrac and apply coarsen-marks, too.
 */
template <typename TDomain, typename TAlgebra>
void MarkForAdaption_GradientJump(IRefiner& refiner,
                                   SmartPtr<GridFunction<TDomain, TAlgebra> > u,
                                   const char* cmp,
                                   number refFrac,
                                   int minLvl, int maxLvl,
                                   std::string jumpType)
{
	PROFILE_FUNC();
	using namespace std;
//	types
	typedef GridFunction<TDomain, TAlgebra> TFunction;
	typedef typename TFunction::domain_type::grid_type grid_t;
	typedef typename TFunction::element_type elem_t;
	typedef typename TFunction::template traits<elem_t>::const_iterator	ElemIter;

//	function id
	const size_t fct = u->fct_id_by_name(cmp);
	UG_COND_THROW(fct >= u->num_fct(),
				  "Function space does not contain a function with name " << cmp);

//	get multigrid and position accessor
	grid_t& mg = *u->domain()->grid();

// 	attach error field
	typedef Attachment<number> ANumber;
	ANumber aError;
	mg.template attach_to<elem_t>(aError);
	MultiGrid::AttachmentAccessor<elem_t, ANumber> aaError(mg, aError);
	
//todo:	the type of the used gradient evaluator should depend on the used grid function.
	typedef GradientEvaluator_LagrangeP1<TFunction>	LagrangeP1Evaluator;
	if(jumpType == std::string("norm"))
		EvaluateGradientJump_Norm<LagrangeP1Evaluator>(*u, fct, aaError);
	else if(jumpType == std::string("sideInt"))
		EvaluateGradientJump_SideIntegral<LagrangeP1Evaluator>(*u, fct, aaError);
	else{
		UG_THROW("Unsupported jumpType in MarkForAdaption_GradientJump: "
				 "Valid values are: norm, sideInt");
	}


	number maxElemError = 0;	// error in elements which may be refined
	number minElemError = numeric_limits<number>::max();
	for(ElemIter iter = u->template begin<elem_t>(); iter != u->template end<elem_t>(); ++iter){
		if(mg.get_level(*iter) < maxLvl){
			maxElemError = max(maxElemError, aaError[*iter]);
			minElemError = min(minElemError, aaError[*iter]);
		}
	}

	#ifdef UG_PARALLEL
		pcl::ProcessCommunicator com;
		minElemError = com.allreduce(minElemError, PCL_RO_MIN);
		maxElemError = com.allreduce(maxElemError, PCL_RO_MAX);
	#endif

//	note that aaError contains error-squares
	number refTol = minElemError + (maxElemError - minElemError) * refFrac;
	MarkElementsAbsolute(aaError, refiner, u->dof_distribution(), refTol, -1,
					 	 minLvl, maxLvl);

// 	detach error field
	mg.template detach_from<elem_t>(aError);
};


template <typename TFunction>
void EvaluateResidualErrorP1(SmartPtr<TFunction> u,
	                           SmartPtr<UserData<number, TFunction::dim> > f,
	                           const char* cmp,
	                           number time,
	                           int quadOrder, std::string quadType,
	                           MultiGrid::AttachmentAccessor<
				                     typename TFunction::element_type,
				                     ug::Attachment<number> >& aaError)
{
	using namespace std;
//	types
	typedef typename TFunction::domain_type::grid_type grid_t;
	typedef typename TFunction::element_type elem_t;
	typedef typename TFunction::template traits<elem_t>::const_iterator	ElemIter;
	const int dim = TFunction::dim;

//	function id
	const size_t fct = u->fct_id_by_name(cmp);
	UG_COND_THROW(fct >= u->num_fct(),
				  "Function space does not contain a function with name " << cmp);

//	get multigrid and position accessor
	grid_t& mg = *u->domain()->grid();
	typename TFunction::domain_type::position_accessor_type&
		aaPos = u->domain()->position_accessor();
	
//	evaluate L2-Norm of f and store element contributions in aaError
	SmartPtr<IIntegrand<number, dim> > spIntegrand
		= make_sp(new UserDataIntegrand<number, TFunction>(f, u, time));

	Integrate<dim, dim>(u->template begin<elem_t>(), u->template end<elem_t>(),
			  			aaPos, *spIntegrand, quadOrder, quadType, &aaError);

//	we have to square contributions in aaError and multiply them with the
//	square of the element-diameter
	for(ElemIter iter = u->template begin<elem_t>();
		iter != u->template end<elem_t>(); ++iter)
	{
		elem_t* elem = *iter;
		aaError[elem] *= aaError[elem] * ElementDiameterSq(mg, aaPos, elem);
	}

//	now evaluate contributions of the integral over gradient jumps on the element sides
//	their squares will be added to aaError.
	typedef GradientEvaluator_LagrangeP1<TFunction>	LagrangeP1Evaluator;
	EvaluateGradientJump_SideIntegral<LagrangeP1Evaluator>(*u, fct, aaError, true);

//	finally take the square root of aaError
	for(ElemIter iter = u->template begin<elem_t>();
		iter != u->template end<elem_t>(); ++iter)
	{
		aaError[*iter] = sqrt(aaError[*iter]);
	}
};

/** A classical residual error for the poisson problem on linear shape functions
 * works with h-nodes.
 *
 * Evaluates the element residual error sqrt{hT^2 * ||RT||^2 + 0.5*sum(hS*||RS||^2)}
 * where the sum contains all sides S of an element T. RT denotes the residuum
 * and RS the gradient-jump over sides of T.
 *
 * \param[in]	refTol		Threshold value. Only elements with a higher residual
 *							error than refTol are refined.
 *
 * \todo: Add coarsenFrac and apply coarsen-marks, too.
 */
template <typename TDomain, typename TAlgebra>
void MarkForAdaption_ResidualErrorP1Absolute(IRefiner& refiner,
                                   SmartPtr<GridFunction<TDomain, TAlgebra> > u,
                                   SmartPtr<UserData<number, TDomain::dim> > f,
                                   const char* cmp,
                                   number time,
                                   number refTol,
                                   int minLvl, int maxLvl,
                                   int quadOrder, std::string quadType)
{
	PROFILE_FUNC();
	using namespace std;
//	types
	typedef GridFunction<TDomain, TAlgebra>	TFunction;
	typedef typename TFunction::domain_type::grid_type grid_t;
	typedef typename TFunction::element_type elem_t;
	typedef typename TFunction::template traits<elem_t>::const_iterator	ElemIter;

// 	attach error field
	grid_t& mg = *u->domain()->grid();
	typedef Attachment<number> ANumber;
	ANumber aError;
	mg.template attach_to<elem_t>(aError);
	MultiGrid::AttachmentAccessor<elem_t, ANumber> aaError(mg, aError);
	
//	Evaluate the residual error
	EvaluateResidualErrorP1(u, f, cmp, time, quadOrder, quadType, aaError);
	MarkElementsAbsolute(aaError, refiner, u->dof_distribution(), refTol, -1,
					 	 minLvl, maxLvl);

// 	detach error field
	mg.template detach_from<elem_t>(aError);
};

/** A classical residual error for the poisson problem on linear shape functions
 * works with h-nodes.
 *
 * Evaluates the element residual error sqrt{hT^2 * ||RT||^2 + 0.5*sum(hS*||RS||^2)}
 * where the sum contains all sides S of an element T. RT denotes the residuum
 * and RS the gradient-jump over sides of T.
 *
* \param[in]	refFrac		given minElemError and maxElemError, all elements with
 *							an error > minElemError+refFrac*(maxElemError-minElemError)
 *							will be refined.
 *
 * \todo: Add coarsenFrac and apply coarsen-marks, too.
 */
template <typename TDomain, typename TAlgebra>
void MarkForAdaption_ResidualErrorP1Relative(IRefiner& refiner,
                                   SmartPtr<GridFunction<TDomain, TAlgebra> > u,
                                   SmartPtr<UserData<number, TDomain::dim> > f,
                                   const char* cmp,
                                   number time,
                                   number refFrac,
                                   int minLvl, int maxLvl,
                                   int quadOrder, std::string quadType)
{
	PROFILE_FUNC();
	using namespace std;
//	types
	typedef GridFunction<TDomain, TAlgebra> TFunction;
	typedef typename TFunction::domain_type::grid_type grid_t;
	typedef typename TFunction::element_type elem_t;
	typedef typename TFunction::template traits<elem_t>::const_iterator	ElemIter;

// 	attach error field
	grid_t& mg = *u->domain()->grid();
	typedef Attachment<number> ANumber;
	ANumber aError;
	mg.template attach_to<elem_t>(aError);
	MultiGrid::AttachmentAccessor<elem_t, ANumber> aaError(mg, aError);
	
//	Evaluate the residual error
	EvaluateResidualErrorP1(u, f, cmp, time, quadOrder, quadType, aaError);

//	find min- and max-values
	number maxElemError = 0;	// error in elements which may be refined
	number minElemError = numeric_limits<number>::max();
	for(ElemIter iter = u->template begin<elem_t>();
		iter != u->template end<elem_t>(); ++iter)
	{
		if(mg.get_level(*iter) < maxLvl){
			 number err = aaError[*iter] = sqrt(aaError[*iter]);
			 maxElemError = max(maxElemError, err);
			 minElemError = min(minElemError, err);
		}
	}

	#ifdef UG_PARALLEL
	 	pcl::ProcessCommunicator com;
	 	minElemError = com.allreduce(minElemError, PCL_RO_MIN);
	 	maxElemError = com.allreduce(maxElemError, PCL_RO_MAX);
	#endif

	number refTol = minElemError + (maxElemError - minElemError) * refFrac;
	MarkElementsAbsolute(aaError, refiner, u->dof_distribution(), refTol, -1,
					 	 minLvl, maxLvl);

// 	detach error field
	mg.template detach_from<elem_t>(aError);
};


template <typename TDomain, typename TAlgebra>
void MarkForAdaption_GradientAverage(IRefiner& refiner,
                                   SmartPtr<GridFunction<TDomain, TAlgebra> > u,
                                   const char* cmp,
                                   number refFrac,
                                   int minLvl, int maxLvl)
{
	PROFILE_FUNC();
	using namespace std;
//	types
	typedef GridFunction<TDomain, TAlgebra> TFunction;
	static const int dim = TFunction::dim;
	typedef typename TFunction::domain_type::grid_type grid_t;
	typedef typename TFunction::element_type elem_t;
	typedef typename TFunction::template traits<elem_t>::const_iterator	ElemIter;
	typedef MathVector<dim>	vector_t;

//	get position accessor
	typename TFunction::domain_type::position_accessor_type& aaPos
			= u->domain()->position_accessor();

//	function id
	const size_t fct = u->fct_id_by_name(cmp);
	UG_COND_THROW(fct >= u->num_fct(),
				  "Function space does not contain a function with name " << cmp);

//	get multigrid and position accessor
	grid_t& mg = *u->domain()->grid();

// 	attach error field
	typedef Attachment<number> ANumber;
	ANumber aError;
	mg.template attach_to<elem_t>(aError);
	MultiGrid::AttachmentAccessor<elem_t, ANumber> aaErrorElem(mg, aError);
	
//	attach gradient to vertices and initialize it with zero
	typedef Attachment<vector_t> AGrad;
	AGrad aGrad;
	mg.attach_to_vertices(aGrad);
	MultiGrid::AttachmentAccessor<Vertex, AGrad> aaGradVrt(mg, aGrad);
	vector_t zeroVec;
	VecSet(zeroVec, 0);
	SetAttachmentValues(aaGradVrt, mg.template begin<Vertex>(),
						mg.template end<Vertex>(), zeroVec);

//	attach contributor-counter to vertices. Required for parallel execution!
	ANumber aNumContribs;
	mg.attach_to_vertices(aNumContribs);
	MultiGrid::AttachmentAccessor<Vertex, ANumber> aaNumContribsVrt(mg, aNumContribs);
	SetAttachmentValues(aaNumContribsVrt, mg.template begin<Vertex>(),
						mg.template end<Vertex>(), 0);

//	average the gradients in the grids vertices
//todo:	the type of the used gradient evaluator should depend on the used grid function.
	typedef GradientEvaluator_LagrangeP1<TFunction>	LagrangeP1Evaluator;

//	loop elements and evaluate gradient
	LagrangeP1Evaluator gradEvaluator(u.get(), fct);
	Grid::vertex_traits::secure_container	vrts;
	
	ElemIter iterElemEnd = u->template end<elem_t>();
	for(ElemIter iter = u->template begin<elem_t>();
		iter != iterElemEnd; ++iter)
	{
	//	get the element
		elem_t* elem = *iter;
		vector_t elemGrad = gradEvaluator.evaluate(elem);

		mg.associated_elements(vrts, elem);
		for(size_t i = 0; i < vrts.size(); ++i){
			Vertex* v = vrts[i];
			aaGradVrt[v] += elemGrad;
			++aaNumContribsVrt[v];
		}
	}

//	in a parallel environment we now have to sum the gradients over parallel
//	interfaces
	#ifdef UG_PARALLEL
		typedef typename GridLayoutMap::Types<Vertex>::Layout layout_t;
		DistributedGridManager& dgm = *mg.distributed_grid_manager();
		GridLayoutMap& glm = dgm.grid_layout_map();
		pcl::InterfaceCommunicator<layout_t> icom;

	//	sum all copies at the h-master attachment
		ComPol_AttachmentReduce<layout_t, AGrad> compolSumGrad(mg, aGrad, PCL_RO_SUM);
		ComPol_AttachmentReduce<layout_t, ANumber> compolSumNum(mg, aNumContribs, PCL_RO_SUM);
		icom.exchange_data(glm, INT_H_SLAVE, INT_H_MASTER, compolSumGrad);
		icom.exchange_data(glm, INT_H_SLAVE, INT_H_MASTER, compolSumNum);
		icom.communicate();

	//	copy the sum from the master to all of its slave-copies
		ComPol_CopyAttachment<layout_t, AGrad> compolCopyGrad(mg, aGrad);
		ComPol_CopyAttachment<layout_t, ANumber> compolCopyNum(mg, aNumContribs);
		icom.exchange_data(glm, INT_H_MASTER, INT_H_SLAVE, compolCopyGrad);
		icom.exchange_data(glm, INT_H_MASTER, INT_H_SLAVE, compolCopyNum);
		icom.communicate();

	//todo: communication to vmasters may not be necessary here...
	//		it is performed to make sure that all surface-rim-children
	//		contain their true value.
		icom.exchange_data(glm, INT_V_SLAVE, INT_V_MASTER, compolCopyGrad);
		icom.exchange_data(glm, INT_V_SLAVE, INT_V_MASTER, compolCopyNum);

		icom.communicate();
	#endif

//	if we're currently operating on a surface function, we have to adjust
//	errors between shadowed and shadowing vertices
	if(u->grid_level().type() == GridLevel::SURFACE){
	//todo: avoid iteration over the whole grid!
	//todo: reduce communication
		const SurfaceView* surfView = u->approx_space()->surface_view().get();

		// for(side_iter_t iter = u.template begin<side_t>(SurfaceView::SHADOW_RIM);
		// 	iter != u.template end<side_t>(SurfaceView::SHADOW_RIM); ++iter)
		typedef Grid::traits<Vertex>::iterator grid_side_iter_t;
		for(grid_side_iter_t iter = mg.template begin<Vertex>();
			iter != mg.template end<Vertex>(); ++iter)
		{
			Vertex* s = *iter;
			if(!surfView->surface_state(s).partially_contains(SurfaceView::SHADOW_RIM))
				continue;

			Vertex* c = mg.get_child_vertex(s);
			if(c){
				aaGradVrt[s] += aaGradVrt[c];
				aaNumContribsVrt[s] += aaNumContribsVrt[c];
			}
		}
	//	those elems which have local children now contain their true value (vslaves...)
		#ifdef UG_PARALLEL
			icom.exchange_data(glm, INT_V_SLAVE, INT_V_MASTER, compolCopyGrad);
			icom.exchange_data(glm, INT_V_SLAVE, INT_V_MASTER, compolCopyNum);
			icom.communicate();
		#endif

	//	copy values up to children
		for(grid_side_iter_t iter = mg.template begin<Vertex>();
			iter != mg.template end<Vertex>(); ++iter)
		{
			Vertex* s = *iter;
			if(!surfView->surface_state(s).partially_contains(SurfaceView::SHADOW_RIM))
				continue;

			Vertex* c = mg.get_child_vertex(s);
			if(c){
				aaGradVrt[c] = aaGradVrt[s];
				aaNumContribsVrt[c] = aaNumContribsVrt[s];
			}
		}

	//	we'll also average values in h-nodes now. If a parent is locally available,
	//	it's associated values are correct at this point. Communication from vmaster
	//	to vslaves is performed afterwards.
		typedef MultiGrid::traits<ConstrainedVertex>::iterator constr_vrt_iter;
		for(constr_vrt_iter iter = mg.template begin<ConstrainedVertex>();
			iter != mg.template end<ConstrainedVertex>(); ++iter)
		{
			ConstrainedVertex* v = *iter;
			GridObject* p = v->get_constraining_object();
			if(p){
				VecSet(aaGradVrt[v], 0);
				aaNumContribsVrt[v] = 0;
				mg.associated_elements(vrts, p);
				int numConstr = 0;
				for(size_t i = 0; i < vrts.size(); ++i){
					if(!surfView->surface_state(vrts[i]).partially_contains(SurfaceView::SURFACE_RIM))
						continue;
					aaGradVrt[v] += aaGradVrt[vrts[i]];
					aaNumContribsVrt[v] += aaNumContribsVrt[vrts[i]];
					++numConstr;
				}
				aaGradVrt[v] /= (number)numConstr;
				aaNumContribsVrt[v] /= (number)numConstr;
			}
		}

	//	those elems which have a local parent now contain their true value (vmasters...)
		#ifdef UG_PARALLEL
		//	copy from v-masters to v-slaves, since there may be constrained
		//	sides which locally do not have a constraining element. Note that
		//	constrained V-Masters always have a local constraining element and thus
		//	contain the correct value.
			icom.exchange_data(glm, INT_V_MASTER, INT_V_SLAVE, compolCopyGrad);
			icom.exchange_data(glm, INT_V_MASTER, INT_V_SLAVE, compolCopyNum);
			icom.communicate();
		#endif
	}

//	evaluate integral over difference of element gradients and averaged
//	vertex gradients
	for(ElemIter iter = u->template begin<elem_t>();
		iter != iterElemEnd; ++iter)
	{
	//	get the element
		elem_t* elem = *iter;
		vector_t elemGrad = gradEvaluator.evaluate(elem);

		vector_t vrtAvrgGrad;
		VecSet(vrtAvrgGrad, 0);
		mg.associated_elements(vrts, elem);
		for(size_t i = 0; i < vrts.size(); ++i){
			Vertex* v = vrts[i];
			vector_t vg = aaGradVrt[v];
			vg /= aaNumContribsVrt[v];
			vrtAvrgGrad += vg;
		}
		vrtAvrgGrad /= (number)vrts.size();
		aaErrorElem[elem] = CalculateVolume(elem, aaPos)
					  * VecDistance(elemGrad, vrtAvrgGrad)
					  / (number)(dim+1);
	}

//	mark for adaption
	number maxElemError = 0;	// error in elements which may be refined
	number minElemError = numeric_limits<number>::max();
	for(ElemIter iter = u->template begin<elem_t>(); iter != u->template end<elem_t>(); ++iter){
		if(mg.get_level(*iter) < maxLvl){
			maxElemError = max(maxElemError, aaErrorElem[*iter]);
			minElemError = min(minElemError, aaErrorElem[*iter]);
		}
	}

	#ifdef UG_PARALLEL
		pcl::ProcessCommunicator com;
		minElemError = com.allreduce(minElemError, PCL_RO_MIN);
		maxElemError = com.allreduce(maxElemError, PCL_RO_MAX);
	#endif

//	note that aaErrorElem contains error-squares
	number refTol = minElemError + (maxElemError - minElemError) * refFrac;
	MarkElementsAbsolute(aaErrorElem, refiner, u->dof_distribution(), refTol, -1,
					 	 minLvl, maxLvl);

// 	detach error field
	mg.detach_from_vertices(aNumContribs);
	mg.detach_from_vertices(aGrad);
	mg.template detach_from<elem_t>(aError);
};

} // end namespace ug

#endif /* __H__UG__LIB_DISC__FUNCTION_SPACE__ERROR_INDICATOR__ */
