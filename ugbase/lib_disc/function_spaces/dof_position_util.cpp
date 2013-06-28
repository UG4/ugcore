/*
 * dof_position_util.cpp
 *
 *  Created on: 11.01.2012
 *      Author: andreasvogel
 */

#include "dof_position_util.h"
#include "grid_function_impl.h"
#include "approximation_space.h"
#include "lib_disc/domain.h"

namespace ug{

////////////////////////////////////////////////////////////////////////////////
//	Extract Positions
////////////////////////////////////////////////////////////////////////////////

template <typename TDomain>
void ExtractPositionsVertex(ConstSmartPtr<TDomain> domain,
                            ConstSmartPtr<DoFDistribution> dd,
                            std::vector<MathVector<TDomain::dim> >& vPos,
                            const std::vector<int>* pvMapGlobalToPatch)
{
//	get position accessor
	const typename TDomain::position_accessor_type& aaPos = domain->position_accessor();

//	iterator
	typename DoFDistribution::traits<VertexBase>::const_iterator iter, iterEnd;

//	algebra indices vector
	std::vector<size_t> ind;

//	get iterators
	iter = dd->begin<VertexBase>();
	iterEnd = dd->end<VertexBase>();

//	loop all vertices
	for(;iter != iterEnd; ++iter)
	{
	//	get vertex
		VertexBase* v = *iter;

	//	load indices associated with vertex
		dd->inner_algebra_indices(v, ind);

	//	write position
		for(size_t i = 0; i < ind.size(); ++i)
		{
			size_t index = ind[i];
			if(pvMapGlobalToPatch){
				if((*pvMapGlobalToPatch)[index] < 0) continue;
				index = (*pvMapGlobalToPatch)[index];
			}
			vPos[index] = aaPos[v];
		}
	}
}

template <typename TDomain, typename TBaseElem>
void ExtractPositionsElem(ConstSmartPtr<TDomain> domain,
                          ConstSmartPtr<DoFDistribution> dd,
                          std::vector<MathVector<TDomain::dim> >& vPos,
                          const std::vector<int>* pvMapGlobalToPatch)
{
	typename DoFDistribution::traits<TBaseElem>::const_iterator iter, iterEnd;

//	vector for positions
	std::vector<MathVector<TDomain::dim> > vElemPos;

//	algebra indices vector
	std::vector<MultiIndex<2> > ind;

//	loop all subsets
	for(int si = 0; si < dd->num_subsets(); ++si)
	{
	//	get iterators
		iter = dd->begin<TBaseElem>(si);
		iterEnd = dd->end<TBaseElem>(si);

	//	loop all elements
		for(;iter != iterEnd; ++iter)
		{
		//	get vertex
			TBaseElem* elem = *iter;

		//	loop all functions
			for(size_t fct = 0; fct < dd->num_fct(); ++fct)
			{
			//	skip non-used function
				if(!dd->is_def_in_subset(fct,si)) continue;

			//	load indices associated with element function
				dd->inner_multi_indices(elem, fct, ind);

			//	load positions associated with element and function
				InnerDoFPosition(vElemPos, elem, *(const_cast<TDomain*>(domain.get())),
				                 dd->local_finite_element_id(fct));

			//	check correct size
				UG_ASSERT(ind.size() == vElemPos.size(), "Num MultiIndex ("<<ind.size()
						  <<") and Num Position ("<<vElemPos.size()<<") must match."
						 "GeomObject dim="<<geometry_traits<TBaseElem>::BASE_OBJECT_ID);

			//	write position
				for(size_t sh = 0; sh < ind.size(); ++sh)
				{
					size_t index = ind[sh][0];
					if(pvMapGlobalToPatch){
						if((*pvMapGlobalToPatch)[index] < 0) continue;
						index = (*pvMapGlobalToPatch)[index];
					}
					vPos[index] = vElemPos[sh];
				}
			}
		}
	}
}

template <typename TDomain>
void ExtractPositions(ConstSmartPtr<TDomain> domain,
                      ConstSmartPtr<DoFDistribution> dd,
                      std::vector<MathVector<TDomain::dim> >& vPos,
                      const std::vector<int>* pvMapGlobalToPatch)
{
//	number of total dofs
	int nr = dd->num_indices();

//	resize positions
	vPos.resize(nr);

//	extract for all element types
	if(dd->max_dofs(VERTEX)) ExtractPositionsVertex<TDomain>(domain, dd, vPos, pvMapGlobalToPatch);
	if(dd->max_dofs(EDGE)) ExtractPositionsElem<TDomain, EdgeBase>(domain, dd, vPos, pvMapGlobalToPatch);
	if(dd->max_dofs(FACE)) ExtractPositionsElem<TDomain, Face>(domain, dd, vPos, pvMapGlobalToPatch);
	if(dd->max_dofs(VOLUME)) ExtractPositionsElem<TDomain, Volume>(domain, dd, vPos, pvMapGlobalToPatch);
}

template void ExtractPositions(ConstSmartPtr<Domain1d> domain, ConstSmartPtr<DoFDistribution> dd, std::vector<MathVector<Domain1d::dim> >& vPos, const std::vector<int>* pvMapGlobalToPatch);
template void ExtractPositions(ConstSmartPtr<Domain2d> domain, ConstSmartPtr<DoFDistribution> dd, std::vector<MathVector<Domain2d::dim> >& vPos, const std::vector<int>* pvMapGlobalToPatch);
template void ExtractPositions(ConstSmartPtr<Domain3d> domain, ConstSmartPtr<DoFDistribution> dd, std::vector<MathVector<Domain3d::dim> >& vPos, const std::vector<int>* pvMapGlobalToPatch);

////////////////////////////////////////////////////////////////////////////////
//	Extract (Positions, Index) Pairs
////////////////////////////////////////////////////////////////////////////////

template <typename TDomain>
void ExtractPositionsVertex(ConstSmartPtr<TDomain> domain,
                            ConstSmartPtr<DoFDistribution> dd,
                            std::vector<std::pair<MathVector<TDomain::dim>, size_t> >& vPosPair)
{
//	get position accessor
	const typename TDomain::position_accessor_type& aaPos = domain->position_accessor();

//	resize positions
	vPosPair.resize(dd->num_indices());

	typedef DoFDistribution::traits<VertexBase>::const_iterator const_iterator;

//	loop all vertices
	const_iterator iter = dd->begin<VertexBase>();
	const_iterator iterEnd = dd->end<VertexBase>();

//	algebra indices vector
	std::vector<size_t> ind;

	for(;iter != iterEnd; ++iter)
	{
	//	get vertex
		VertexBase* v = *iter;

	//	load indices associated with vertex
		dd->inner_algebra_indices(v, ind);

	//	write position and index
		for(size_t i = 0; i < ind.size(); ++i)
		{
			const size_t index = ind[i];
			vPosPair[index].first = aaPos[v];
			vPosPair[index].second = index;
		}
	}
}


template <typename TDomain, typename TBaseElem>
void ExtractPositionsElem(ConstSmartPtr<TDomain> domain,
                          ConstSmartPtr<DoFDistribution> dd,
                          std::vector<std::pair<MathVector<TDomain::dim>, size_t> >& vPosPair)
{
	typename DoFDistribution::traits<TBaseElem>::const_iterator iter, iterEnd;

//	vector for positions
	std::vector<MathVector<TDomain::dim> > vElemPos;

//	algebra indices vector
	std::vector<MultiIndex<2> > ind;

//	loop all subsets
	for(int si = 0; si < dd->num_subsets(); ++si)
	{
	//	get iterators
		iter = dd->begin<TBaseElem>(si);
		iterEnd = dd->end<TBaseElem>(si);

	//	loop all elements
		for(;iter != iterEnd; ++iter)
		{
		//	get vertex
			TBaseElem* elem = *iter;

		//	loop all functions
			for(size_t fct = 0; fct < dd->num_fct(); ++fct)
			{
			//	skip non-used function
				if(!dd->is_def_in_subset(fct,si)) continue;

			//	load indices associated with element function
				dd->inner_multi_indices(elem, fct, ind);

			//	load positions associated with element and function
				InnerDoFPosition(vElemPos, elem, *(const_cast<TDomain*>(domain.get())),
				                 dd->local_finite_element_id(fct));

			//	check correct size
				UG_ASSERT(ind.size() == vElemPos.size(), "Num MultiIndex ("<<ind.size()
						  <<") and Num Position ("<<vElemPos.size()<<") must match."
						 "GeomObject dim="<<geometry_traits<TBaseElem>::BASE_OBJECT_ID);

			//	write position
				for(size_t sh = 0; sh < ind.size(); ++sh)
				{
					const size_t index = ind[sh][0];
					vPosPair[index].first = vElemPos[sh];
					vPosPair[index].second = index;

				}
			}
		}
	}
}

template <typename TDomain>
void ExtractPositions(ConstSmartPtr<TDomain> domain,
                      ConstSmartPtr<DoFDistribution> dd,
                      std::vector<std::pair<MathVector<TDomain::dim>, size_t> >& vPosPair)
{
//	number of total dofs
	int nr = dd->num_indices();

//	resize positions
	vPosPair.resize(nr);

//	extract for all element types
	if(dd->max_dofs(VERTEX)) ExtractPositionsVertex<TDomain>(domain, dd, vPosPair);
	if(dd->max_dofs(EDGE)) ExtractPositionsElem<TDomain, EdgeBase>(domain, dd, vPosPair);
	if(dd->max_dofs(FACE)) ExtractPositionsElem<TDomain, Face>(domain, dd, vPosPair);
	if(dd->max_dofs(VOLUME)) ExtractPositionsElem<TDomain, Volume>(domain, dd, vPosPair);
}

template void ExtractPositions(ConstSmartPtr<Domain1d> domain, ConstSmartPtr<DoFDistribution> dd, std::vector<std::pair<MathVector<Domain1d::dim>, size_t> >& vPos);
template void ExtractPositions(ConstSmartPtr<Domain2d> domain, ConstSmartPtr<DoFDistribution> dd, std::vector<std::pair<MathVector<Domain2d::dim>, size_t> >& vPos);
template void ExtractPositions(ConstSmartPtr<Domain3d> domain, ConstSmartPtr<DoFDistribution> dd, std::vector<std::pair<MathVector<Domain3d::dim>, size_t> >& vPos);

////////////////////////////////////////////////////////////////////////////////
//	Extract (Positions, Index) Pairs for a single component
////////////////////////////////////////////////////////////////////////////////


template <typename TDomain, typename TBaseElem>
void ExtractPositionsElem(ConstSmartPtr<TDomain> domain,
                          ConstSmartPtr<DoFDistribution> dd,
                          const size_t fct,
                          std::vector<std::pair<MathVector<TDomain::dim>, size_t> >& vPosPair)
{
	typename DoFDistribution::traits<TBaseElem>::const_iterator iter, iterEnd;

//	vector for positions
	std::vector<MathVector<TDomain::dim> > vElemPos;

//	algebra indices vector
	std::vector<MultiIndex<2> > ind;

//	a pair
	std::pair<MathVector<TDomain::dim>, size_t> pair;

//	loop all subsets
	for(int si = 0; si < dd->num_subsets(); ++si)
	{
	//	get iterators
		iter = dd->begin<TBaseElem>(si);
		iterEnd = dd->end<TBaseElem>(si);

	//	skip non-used function
		if(!dd->is_def_in_subset(fct,si)) continue;

	//	loop all elements
		for(;iter != iterEnd; ++iter)
		{
		//	get vertex
			TBaseElem* elem = *iter;

		//	load indices associated with element function
			dd->inner_multi_indices(elem, fct, ind);

		//	load positions associated with element and function
			InnerDoFPosition(vElemPos, elem, *(const_cast<TDomain*>(domain.get())),
							 dd->local_finite_element_id(fct));

		//	check correct size
			UG_ASSERT(ind.size() == vElemPos.size(), "Num MultiIndex ("<<ind.size()
					  <<") and Num Position ("<<vElemPos.size()<<") must match."
					 "GeomObject dim="<<geometry_traits<TBaseElem>::BASE_OBJECT_ID);

		//	write position
			for(size_t sh = 0; sh < ind.size(); ++sh)
			{
				const size_t index = ind[sh][0];

				pair.first = vElemPos[sh];
				pair.second = index;

				vPosPair.push_back(pair);
			}
		}
	}
}

template <typename TDomain>
void ExtractPositions(ConstSmartPtr<TDomain> domain,
                      ConstSmartPtr<DoFDistribution> dd,
                      const size_t fct,
                      std::vector<std::pair<MathVector<TDomain::dim>, size_t> >& vPosPair)
{
//	resize positions
	vPosPair.clear();

//	extract for all element types
	if(dd->max_dofs(VERTEX)) ExtractPositionsElem<TDomain, VertexBase>(domain, dd, fct, vPosPair);
	if(dd->max_dofs(EDGE)) ExtractPositionsElem<TDomain, EdgeBase>(domain, dd, fct, vPosPair);
	if(dd->max_dofs(FACE)) ExtractPositionsElem<TDomain, Face>(domain, dd, fct, vPosPair);
	if(dd->max_dofs(VOLUME)) ExtractPositionsElem<TDomain, Volume>(domain, dd, fct, vPosPair);
}

template void ExtractPositions(ConstSmartPtr<Domain1d> domain, ConstSmartPtr<DoFDistribution> dd, const size_t fct, std::vector<std::pair<MathVector<Domain1d::dim>, size_t> >& vPos);
template void ExtractPositions(ConstSmartPtr<Domain2d> domain, ConstSmartPtr<DoFDistribution> dd, const size_t fct, std::vector<std::pair<MathVector<Domain2d::dim>, size_t> >& vPos);
template void ExtractPositions(ConstSmartPtr<Domain3d> domain, ConstSmartPtr<DoFDistribution> dd, const size_t fct, std::vector<std::pair<MathVector<Domain3d::dim>, size_t> >& vPos);

} // end namespace ug
