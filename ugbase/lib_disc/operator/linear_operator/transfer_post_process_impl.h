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
#include "lib_disc/function_spaces/integrate.h"

namespace ug{

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// 	AverageComponent
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template <typename TDomain, typename TAlgebra>
void AverageComponent<TDomain, TAlgebra>::
post_process(SmartPtr<vector_type> spU)
{
	PROFILE_FUNC_GROUP("gmg");

	if(m_vCmp.empty())
		return;

	SmartPtr<GridFunction<TDomain, TAlgebra> > spGF
					= spU.template cast_dynamic<GridFunction<TDomain, TAlgebra> >();

	if(spGF.invalid())
		UG_THROW("AverageComponent: expects correction to be a GridFunction.");

	ConstSmartPtr<DoFDistributionInfo> ddinfo =
								spGF->approx_space()->dof_distribution_info();

//	compute integral of components
	const number area = Integral(1.0, spGF);
	std::vector<number> vIntegral(m_vCmp.size(), 0.0);
	for(size_t f = 0; f < m_vCmp.size(); f++)
		vIntegral[f] = Integral(spGF, m_vCmp[f].c_str());

//	subtract value
	for(size_t f = 0; f < m_vCmp.size(); f++)
	{
		const number sub = vIntegral[f] / area;
		const size_t fct = spGF->fct_id_by_name(m_vCmp[f].c_str());

		if(ddinfo->max_fct_dofs(fct, VERTEX)) subtract_value<VertexBase>(spGF, fct, sub);
		if(ddinfo->max_fct_dofs(fct, EDGE)) subtract_value<EdgeBase>(spGF, fct, sub);
		if(ddinfo->max_fct_dofs(fct, FACE)) subtract_value<Face>(spGF, fct, sub);
		if(ddinfo->max_fct_dofs(fct, VOLUME)) subtract_value<Volume>(spGF, fct, sub);
	}
}

template <typename TDomain, typename TAlgebra>
template <typename TBaseElem>
void AverageComponent<TDomain, TAlgebra>::
subtract_value(SmartPtr<GridFunction<TDomain, TAlgebra> > spGF, size_t fct, number sub)
{
	typedef typename GridFunction<TDomain, TAlgebra>::template traits<TBaseElem>::const_iterator iter_type;

	iter_type iter = spGF->template begin<TBaseElem>();
	iter_type iterEnd = spGF->template end<TBaseElem>();

//  loop elems
	std::vector<DoFIndex> vMultInd;
	for(; iter != iterEnd; ++iter)
	{
	//	get element
		TBaseElem* elem = *iter;

	//  get global indices
		spGF->inner_dof_indices(elem, fct, vMultInd);

	//	sum up value
		for(size_t i = 0; i < vMultInd.size(); ++i)
		{
			DoFRef(*spGF, vMultInd[i]) -= sub;
		}
	}
}

template <typename TDomain, typename TAlgebra>
SmartPtr<ITransferPostProcess<TAlgebra> >
AverageComponent<TDomain, TAlgebra>::clone()
{
	return SmartPtr<AverageComponent>(new AverageComponent(m_vCmp));
}

} // end namespace ug

#endif /* __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__TRANSFER_POST_PROCESS_IMPL__ */
