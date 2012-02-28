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

namespace ug{

template <typename TElem, typename TFunction>
bool ComputeGradient(TFunction& u,
                     MultiGrid::AttachmentAccessor<
                     	 typename geometry_traits<TElem>::geometric_base_object,
                     	 ug::Attachment<number> >& aaError)
{
//	get reference element
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;
	const ref_elem_type& refElem = Provider<ref_elem_type>::get();

//	get reference dimension
	static const int dim = ref_elem_type::dim;

//	get world dimension
	static const int worldDim = TFunction::domain_type::dim;

//	get position accessor
	typename TFunction::domain_type::position_accessor_type& aaPos
			= u.domain()->position_accessor();

//	get trial space
	const LocalShapeFunctionSet<ref_elem_type>& TrialSpace =
			LocalShapeFunctionSetProvider::
				get<ref_elem_type>
				(LFEID(LFEID::LAGRANGE, 1));

//	create a reference mapping
	ReferenceMapping<ref_elem_type, worldDim> mapping;

//	number of shape functions
	size_t num_sh = (size_t)ref_elem_type::num_corners;

//	local IP
	MathVector<dim> localIP;

//	some storage
	std::vector<MathVector<dim> > localGrad(num_sh);
	std::vector<MathVector<worldDim> > globalGrad(num_sh);
	std::vector<MathVector<worldDim> > vCorner(num_sh);
	MathMatrix<dim, dim> JTInv;

//	compute local midpoint
	VecSet(localIP, 0.0);
	for(size_t i = 0; i < num_sh; ++i)
		localIP += refElem.corner(i);
	localIP *= 1./(num_sh);

//	evaluate reference gradient at local midpoint
	TrialSpace.grads(&localGrad[0], localIP);

	typedef typename TFunction::template traits<TElem>::const_iterator const_iterator;

//	loop subsets
	for(int si = 0; si < u.num_subsets(); ++si)
	{
	//	get iterator over elements
		const_iterator iter = u.template begin<TElem>(si);
		const_iterator iterEnd = u.template end<TElem>(si);

	//	loop elements
		for(; iter != iterEnd; ++iter)
		{
		//	get the element
			TElem* elem = *iter;

		//	get corners of element
			CollectCornerCoordinates(vCorner, *elem, aaPos);

		//	update mapping
			mapping.update(&vCorner[0]);

		//	compute jacobian
			mapping.jacobian_transposed_inverse(JTInv, localIP);

		//	compute size (volume) of element
			const number elemSize = ElementSize<ref_elem_type, dim>(&vCorner[0]);

		//	compute gradient at mid point by summing contributions of all shape fct
			MathVector<dim> MidGrad; VecSet(MidGrad, 0.0);
			for(size_t sh = 0 ; sh < num_sh; ++sh)
			{
			//	get global Gradient
				MatVecMult(globalGrad[sh], JTInv,localGrad[sh]);

			//	get vertex
				VertexBase* vert = elem->vertex(sh);

			//	get of of vertex
				//\todo: this is for fct=0 only
				std::vector<MultiIndex<2> > ind;
				u.inner_multi_indices(vert, 0, ind);

			//	scale global gradient
				globalGrad[sh] *= (u.get_dof_value(ind[0][0], ind[0][1]));

			//	sum up
				MidGrad += globalGrad[sh];
			}

		//	write result in array storage
			aaError[elem] = VecTwoNorm(MidGrad) * pow(elemSize, 2./dim);
		}
	}

//	we're done
	return true;
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
template <typename TElemBase, typename TFunction>
bool MarkElements(IRefiner& refiner,
                  TFunction& u,
                  number TOL, number scale,
                  MultiGrid::AttachmentAccessor<TElemBase, ug::Attachment<number> >& aaError)
{
//	reset maximum of error
	number max = 0.0;

	typedef typename TFunction::template traits<TElemBase>::const_iterator const_iterator;

//	loop subsets
	for(int si = 0; si < u.num_subsets(); ++si)
	{
	//	get element iterator
		const_iterator iter = u.template begin<TElemBase>(si);
		const_iterator iterEnd = u.template end<TElemBase>(si);

	//	loop all elements to find the maximum of the error
		for( ;iter != iterEnd; ++iter)
		{
		//	get element
			TElemBase* elem = *iter;

		//	search for maximum
			if(aaError[elem] > max)
				max = aaError[elem];
		}
	}

#ifdef UG_PARALLEL
	pcl::ProcessCommunicator com;
	number maxLocal = max;
	com.allreduce(&maxLocal, &max, 1, PCL_DT_DOUBLE, PCL_RO_MAX);
	UG_LOG("Max Error on Proc " << pcl::GetProcRank() << " is " << maxLocal << ".");
#endif

	UG_LOG("Max Error is " << max << ". ");

//	check if something to do
	if(max <= TOL) return false;

//	Compute minimum
	number min = max*scale;
	if(min < TOL) min = TOL;

	UG_LOG("Refining all elements with error >= " << min << ". ");

//	reset counter
	int numMarked = 0;

//	loop subsets
	for(int si = 0; si < u.num_subsets(); ++si)
	{
		const_iterator iter = u.template begin<TElemBase>(si);
		const_iterator iterEnd = u.template end<TElemBase>(si);

	//	loop elements for marking
		for(; iter != iterEnd; ++iter)
		{
		//	get element
			TElemBase* elem = *iter;

		//	check if element error is in range
			if(aaError[elem] >= min)
			{
			//	mark element and increase counter
				refiner.mark(elem);
				numMarked++;
			}
		}
	}

#ifdef UG_PARALLEL
	int numMarkedLocal = numMarked;
	com.allreduce(&numMarkedLocal, &numMarked, 1, PCL_DT_INT, PCL_RO_SUM);
	UG_LOG(numMarkedLocal << " Elements marked on Proc " << pcl::GetProcRank() << ".");
#endif

	UG_LOG(numMarked << " Elements marked for refinement.\n");

//	we're done
	return true;
}


template <typename TFunction>
void MarkForRefinement_GradientIndicator_DIM(IRefiner& refiner,
                                             TFunction& u,
                                             number TOL, number scale,
                                             Int2Type<1>)
{
//	get domain
	typename TFunction::domain_type& domain = u.approximation_space().domain();

//	get multigrid
	typename TFunction::domain_type::grid_type& mg = domain.grid();

// 	attach error field
	typedef Attachment<number> ANumber;
	ANumber aError;
	mg.attach_to_edges(aError);
	MultiGrid::EdgeAttachmentAccessor<ANumber> aaError(mg, aError);

// 	Compute error on elements
	ComputeGradient<Edge, TFunction>(u, aaError);

// 	Mark elements for refinement
	if(!MarkElements<EdgeBase, TFunction>(refiner, u, TOL, scale, aaError))
		UG_LOG("No element marked. Not refining the grid.\n");

// 	detach error field
	mg.detach_from_edges(aError);
};

template <typename TFunction>
void MarkForRefinement_GradientIndicator_DIM(IRefiner& refiner,
                                             TFunction& u,
                                             number TOL, number scale,
                                             Int2Type<2>)
{
//	get domain
	SmartPtr<typename TFunction::domain_type> domain = u.domain();

//	get multigrid
	SmartPtr<typename TFunction::domain_type::grid_type> pMG = domain->grid();

// 	attach error field
	typedef Attachment<number> ANumber;
	ANumber aError;
	pMG->attach_to_faces(aError);
	MultiGrid::FaceAttachmentAccessor<ANumber> aaError(*pMG, aError);

// 	Compute error on elements
	ComputeGradient<Triangle, TFunction>(u, aaError);
	ComputeGradient<Quadrilateral, TFunction>(u, aaError);

// 	Mark elements for refinement
	if(!MarkElements<Face, TFunction>(refiner, u, TOL, scale, aaError))
		UG_LOG("No element marked. Not refining the grid.\n");

// 	detach error field
	pMG->detach_from_faces(aError);
};

template <typename TFunction>
void MarkForRefinement_GradientIndicator_DIM(IRefiner& refiner,
                                             TFunction& u,
                                             number TOL, number scale,
                                             Int2Type<3>)
{
//	get domain
	typename TFunction::domain_type& domain = u.approximation_space().domain();

//	get multigrid
	typename TFunction::domain_type::grid_type& mg = domain.grid();

// 	attach error field
	typedef Attachment<number> ANumber;
	ANumber aError;
	mg.attach_to_volumes(aError);
	MultiGrid::VolumeAttachmentAccessor<ANumber> aaError(mg, aError);

// 	Compute error on elements
	ComputeGradient<Tetrahedron, TFunction>(u, aaError);
	ComputeGradient<Hexahedron, TFunction>(u, aaError);

// 	Mark elements for refinement
	if(!MarkElements<Volume, TFunction>(refiner, u, TOL, scale, aaError))
		UG_LOG("No element marked. Not refining the grid.\n");

// 	detach error field
	mg.detach_from_faces(aError);
};

template <typename TFunction>
void MarkForRefinement_GradientIndicator(IRefiner& refiner,
                                         TFunction& u,
                                         number TOL, number scale)
{
	MarkForRefinement_GradientIndicator_DIM(refiner, u, TOL, scale,
	                                    Int2Type<TFunction::domain_type::dim>());
}

} // end namespace ug

#endif /* __H__UG__LIB_DISC__FUNCTION_SPACE__ERROR_INDICATOR__ */
