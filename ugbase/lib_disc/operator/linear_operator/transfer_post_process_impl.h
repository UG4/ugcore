/*
 * prolongation_operator_impl.h
 *
 *  Created on: 04.12.2009
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__TRANSFER_POST_PROCESS_IMPL__
#define __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__TRANSFER_POST_PROCESS_IMPL__

#include "transfer_post_process.h"
#include "lib_disc/common/groups_util.h"
#include "lib_disc/common/geometry_util.h"

namespace ug{

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// 	AverageComponent
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template <typename TDomain, typename TAlgebra>
void AverageComponent<TDomain, TAlgebra>::set_levels(GridLevel level)
{
	m_level = level;

	if(m_level.type() != GridLevel::LEVEL)
		UG_THROW("AverageComponent<TDomain, TAlgebra>::set_levels:"
				" Can only project between level dof distributions, but "<<m_level);
}

template <typename TDomain, typename TAlgebra>
void AverageComponent<TDomain, TAlgebra>::init()
{
	if(!m_spApproxSpace.valid())
		UG_THROW("AverageComponent<TDomain, TAlgebra>::init: "
				"Approximation Space not set. Cannot init Projection.");

//	read functions
	try{
		m_fctGrp.clear();
		m_fctGrp.set_function_pattern(m_spApproxSpace->dof_distribution_info());
		m_fctGrp.add(TokenizeString(m_symbFct));
	}
	UG_CATCH_THROW("AverageComponent: Cannot parse functions for p.");

	m_bInit = true;
}

template <typename TDomain, typename TAlgebra>
void AverageComponent<TDomain, TAlgebra>::
post_process(vector_type& u)
{
	PROFILE_FUNC_GROUP("gmg");
	const int dim = TDomain::dim;
	const DoFDistribution& dd = *m_spApproxSpace->level_dof_distribution(m_level.level());

	std::vector<MultiIndex<2> > vMultInd;

//  iterators
	typedef typename DoFDistribution::template dim_traits<dim>::const_iterator const_iterator;
	typedef typename DoFDistribution::template dim_traits<dim>::geometric_base_object Element;
	const_iterator iter, iterBegin, iterEnd;

//	check piecewise-constant
	for(size_t f = 0; f < m_fctGrp.size(); f++)
	{
		const size_t fct = m_fctGrp[f];
		if(dd.local_finite_element_id(fct) != LFEID(LFEID::PIECEWISE_CONSTANT, 0))
			UG_THROW("Only implemented for piecewise constant.");
	}

//	compute integral of components
	std::vector<number> vIntegral(m_fctGrp.size(), 0.0);
	std::vector<number> vArea(m_fctGrp.size(), 0.0);
	std::vector<MathVector<dim> > vCorner;

//  loop subsets on fine level
	for(int si = 0; si < dd.num_subsets(); ++si)
	{
		iterBegin = dd.template begin<Element>(si);
		iterEnd = dd.template end<Element>(si);

	//  loop vertices for fine level subset
		for(iter = iterBegin; iter != iterEnd; ++iter)
		{
		//	get element
			Element* elem = *iter;

		//	get corners of element
			CollectCornerCoordinates(vCorner, *elem, *m_spApproxSpace->domain());

		//	size of element
			const number elemSize = ElementSize<dim>(elem->reference_object_id(), &vCorner[0]);

		//	loop all components
			for(size_t f = 0; f < m_fctGrp.size(); f++)
			{
			//	get fct index
				const size_t fct = m_fctGrp[f];

			//	check that fct defined on subset
				if(!dd.is_def_in_subset(fct, si)) continue;

			//  get global indices
				dd.inner_multi_indices(elem, fct, vMultInd);

			//	sum up value
				for(size_t i = 0; i < vMultInd.size(); ++i)
				{
					vArea[f] += elemSize;
					vIntegral[f] += DoFRef(u, vMultInd[i]) * elemSize;
				}
			}
		}
	}

//  loop subsets on fine level
	for(int si = 0; si < dd.num_subsets(); ++si)
	{
		iterBegin = dd.template begin<Element>(si);
		iterEnd = dd.template end<Element>(si);

	//  loop vertices for fine level subset
		for(iter = iterBegin; iter != iterEnd; ++iter)
		{
		//	get element
			Element* elem = *iter;

		//	loop all components
			for(size_t f = 0; f < m_fctGrp.size(); f++)
			{
			//	get fct index
				const size_t fct = m_fctGrp[f];

			//	check that fct defined on subset
				if(!dd.is_def_in_subset(fct, si)) continue;

			//  get global indices
				dd.inner_multi_indices(elem, fct, vMultInd);

			//	sum up value
				for(size_t i = 0; i < vMultInd.size(); ++i)
				{
					DoFRef(u, vMultInd[i]) -= vIntegral[f] / vArea[f];
				}
			}
		}
	}

}

template <typename TDomain, typename TAlgebra>
SmartPtr<ITransferPostProcess<TAlgebra> >
AverageComponent<TDomain, TAlgebra>::clone()
{
	SmartPtr<AverageComponent> op(new AverageComponent(m_spApproxSpace, m_symbFct.c_str()));
	return op;
}

} // end namespace ug

#endif /* __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__TRANSFER_POST_PROCESS_IMPL__ */
