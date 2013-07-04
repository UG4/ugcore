/*
 * dof_position_util.cpp
 *
 *  Created on: 11.01.2012
 *      Author: andreasvogel
 */

#include "dof_position_util.h"
#include "approximation_space.h"
#include "lib_disc/domain.h"
#include "lib_disc/domain_util.h"
#include "lib_disc/domain_traits.h"

#include "lib_disc/local_finite_element/local_finite_element_provider.h"
#include "lib_disc/local_finite_element/local_dof_set.h"
#include "lib_disc/reference_element/reference_mapping_provider.h"
#include "lib_disc/reference_element/reference_element_util.h"

namespace ug{

////////////////////////////////////////////////////////////////////////////////
//	DoFPositions
////////////////////////////////////////////////////////////////////////////////

template <int refDim, int dim>
bool InnerDoFPosition(std::vector<MathVector<dim> >& vPos, const ReferenceObjectID roid,
                      const std::vector<MathVector<dim> >& vVertPos, const LFEID& lfeID)
{
//	\TODO: This is a lousy quick-hack. Remove, when dof pos on lower dim elemens
	//		can be handeled correctly for non-lagrangian spaces.
/*	if(lfeID == LFEID(LFEID::CROUZEIX_RAVIART, dim, 1))
	{
		vPos.clear();
		if(lfeID.dim() != refDim+1) return true;

		MathVector<dim> center;
		VecSet(center, 0.0);
		for(size_t co = 0; co < vVertPos.size(); ++co)
			VecAppend(center, vVertPos[co]);
		VecScale(center, center, 1./(vVertPos.size()));

		vPos.push_back(center);
		return true;
	}
	if(lfeID == LFEID(LFEID::PIECEWISE_CONSTANT, dim, 0))
	{
		vPos.clear();
		if(lfeID.dim() != refDim) return true;

		MathVector<dim> center;
		VecSet(center, 0.0);
		for(size_t co = 0; co < vVertPos.size(); ++co)
			VecAppend(center, vVertPos[co]);
		VecScale(center, center, 1./(vVertPos.size()));

		vPos.push_back(center);
		return true;
	}
	if(lfeID == LFEID(LFEID::NEDELEC, dim, 1))
	{
		vPos.clear();
		if(refDim != 1) return true;

		MathVector<dim> center;
		VecSet(center, 0.0);
		for(size_t co = 0; co < vVertPos.size(); ++co)
			VecAppend(center, vVertPos[co]);
		VecScale(center, center, 1./(vVertPos.size()));

		vPos.push_back(center);
		return true;
	}
*/

//	get local dof set
	const DimLocalDoFSet<refDim>& lds = LocalFiniteElementProvider::get_dofs<refDim>(roid, lfeID);

//	create a reference mapping
	 DimReferenceMapping<refDim,dim>& map =
					ReferenceMappingProvider::get<refDim,dim>(roid, vVertPos);

//	clear pos
	vPos.clear();

//	bool flag if position is exact, or no exact position available for shapes
	bool bExact = true;

//	loop all shape functions
	for(size_t sh = 0; sh < lds.num_sh(); ++sh)
	{
	//	check if dof in interior
		if(lds.local_dof(sh).dim() != refDim) continue;

	//	get local position
		MathVector<refDim> locPos;
		bExact &= lds.position(sh, locPos);

	//	map to global position
		MathVector<dim> globPos;
		map.local_to_global(globPos, locPos);

	//	add
		vPos.push_back(globPos);
	}

//	return if positions are given exactly
	return bExact;
};

template <int dim>
bool InnerDoFPositionVertex(std::vector<MathVector<dim> >& vPos, const ReferenceObjectID roid,
                            const std::vector<MathVector<dim> >& vVertPos, const LFEID& lfeID)
{
	UG_ASSERT(vVertPos.size() == 1, "Vertex should have only on inner Vertex");

//	get local dof set
	const CommonLocalDoFSet& lds = LocalFiniteElementProvider::get_dofs(lfeID);

//	clear pos
	vPos.clear();

//	loop all shape functions
	for(int sh = 0; sh < lds.num_dof(ROID_VERTEX); ++sh)
		vPos.push_back(vVertPos[0]);

//	return if positions are given exactly
	return true;
};


template <int dim>
bool InnerDoFPosition(std::vector<MathVector<dim> >& vPos, const ReferenceObjectID roid,
                      const std::vector<MathVector<dim> >& vCornerCoord, const LFEID& lfeID)
{
	switch(ReferenceElementDimension(roid))
	{
		case VERTEX: return InnerDoFPositionVertex<dim>(vPos, roid, vCornerCoord, lfeID);
		case EDGE:   return InnerDoFPosition<1,dim>(vPos, roid, vCornerCoord, lfeID);
		case FACE:   return InnerDoFPosition<2,dim>(vPos, roid, vCornerCoord, lfeID);
		case VOLUME: return InnerDoFPosition<3,dim>(vPos, roid, vCornerCoord, lfeID);
		default: UG_THROW("Base Object type not found.");
	}
}

template <typename TDomain>
bool InnerDoFPosition(std::vector<MathVector<TDomain::dim> >& vPos, GeometricObject* elem,
                      const TDomain& domain, const LFEID& lfeID)
{
//	reference object id
	const ReferenceObjectID roid = elem->reference_object_id();

//	get the vertices
	std::vector<MathVector<TDomain::dim> > vVertPos;
	switch(elem->base_object_id())
	{
		case VERTEX: CollectCornerCoordinates(vVertPos, *static_cast<VertexBase*>(elem), domain, true); break;
		case EDGE: CollectCornerCoordinates(vVertPos, *static_cast<EdgeBase*>(elem), domain, true); break;
		case FACE: CollectCornerCoordinates(vVertPos, *static_cast<Face*>(elem), domain, true); break;
		case VOLUME: CollectCornerCoordinates(vVertPos, *static_cast<Volume*>(elem), domain, true); break;
		default: UG_THROW("Base Object type not found.");
	}

//	forward
	return InnerDoFPosition<TDomain::dim>(vPos, roid, vVertPos, lfeID);
}



template <int refDim, int dim>
bool DoFPosition(std::vector<MathVector<dim> >& vPos, const ReferenceObjectID roid,
                 const std::vector<MathVector<dim> >& vVertPos, const LFEID& lfeID)
{
//	\TODO: This is a lousy quick-hack. Remove, when dof pos on lower dim elemens
	//		can be handeled correctly for non-lagrangian spaces.
/*	if(lfeID.type() == LFEID::CROUZEIX_RAVIART)
	{
		vPos.clear();

		// case, that we are on a side
		if(lfeID.dim() == refDim+1)
		{
			MathVector<dim> center;
			VecSet(center, 0.0);
			for(size_t co = 0; co < vVertPos.size(); ++co)
				VecAppend(center, vVertPos[co]);
			VecScale(center, center, 1./(vVertPos.size()));
			vPos.push_back(center);
			return true;
		}

		// case, that we are on elements with lower dim, than a side
		if(lfeID.dim() > refDim+1) return true;

		// case, that we are on element with higher dim than side (i.e. the volume element itself)
		if(lfeID.dim() == refDim && lfeID.dim() == dim)
		{
			const ReferenceElement& refElem = ReferenceElementProvider::get(roid);

			for(size_t side = 0; side < refElem.num(dim-1); ++side)
			{
				MathVector<dim> center;
				VecSet(center, 0.0);
				for(size_t co = 0; co < refElem.num(dim-1, side, VERTEX); ++co)
					VecAppend(center, vVertPos[refElem.id(dim-1, side, VERTEX, co)]);
				VecScale(center, center, 1./(refElem.num(dim-1)));
				vPos.push_back(center);
			}
			return true;
		}

		// other cases should never happen
		UG_THROW("Special case for Crouzeix-Raviart: case should not happen.");
	}
	if(lfeID.type() == LFEID::PIECEWISE_CONSTANT)
	{
		vPos.clear();
		if(lfeID.dim() != refDim) return true;

		MathVector<dim> center;
		VecSet(center, 0.0);
		for(size_t co = 0; co < vVertPos.size(); ++co)
			VecAppend(center, vVertPos[co]);
		VecScale(center, center, 1./(vVertPos.size()));

		vPos.push_back(center);
		return true;
	}
*/
//	get local shape function set
	const DimLocalDoFSet<refDim>& lds
					= LocalFiniteElementProvider::get_dofs<refDim>(roid, lfeID);

//	create a reference mapping
	const DimReferenceMapping<refDim,dim>& map =
					ReferenceMappingProvider::get<refDim,dim>(roid, vVertPos);

//	clear pos
	vPos.resize(lds.num_sh());

//	bool flag if position is exact, or no exact position available for shapes
	bool bExact = true;

//	loop all shape functions
	for(size_t sh = 0; sh < lds.num_sh(); ++sh)
	{
	//	get local position
		MathVector<refDim> locPos;
		bExact &= lds.position(sh, locPos);

	//	map to global position
		map.local_to_global(vPos[sh], locPos);
	}

//	return if positions are given exactly
	return bExact;
};

template <int dim>
bool DoFPositionVertex(std::vector<MathVector<dim> >& vPos, const ReferenceObjectID roid,
                       const std::vector<MathVector<dim> >& vVertPos, const LFEID& lfeID)
{
//	get local dof set
	const CommonLocalDoFSet& lds = LocalFiniteElementProvider::get_dofs(lfeID);

//	clear pos
	vPos.clear();

//	loop all shape functions
	for(size_t co = 0; co < vVertPos.size(); ++co)
		for(int sh = 0; sh < lds.num_dof(ROID_VERTEX); ++sh)
		vPos.push_back(vVertPos[co]);

//	return if positions are given exactly
	return true;
};

template <int dim>
bool DoFPosition(std::vector<MathVector<dim> >& vPos, const ReferenceObjectID roid,
                 const std::vector<MathVector<dim> >& vCornerCoord, const LFEID& lfeID)
{
	switch(ReferenceElementDimension(roid))
	{
		case VERTEX: return DoFPositionVertex<dim>(vPos, roid, vCornerCoord, lfeID);
		case EDGE:   return DoFPosition<1,dim>(vPos, roid, vCornerCoord, lfeID);
		case FACE:   return DoFPosition<2,dim>(vPos, roid, vCornerCoord, lfeID);
		case VOLUME: return DoFPosition<3,dim>(vPos, roid, vCornerCoord, lfeID);
		default: UG_THROW("Base Object type not found.");
	}
}

template <typename TDomain>
bool DoFPosition(std::vector<MathVector<TDomain::dim> >& vPos, GeometricObject* elem,
                 const TDomain& domain, const LFEID& lfeID)
{
//	reference object id
	const ReferenceObjectID roid = elem->reference_object_id();

//	get the vertices
	std::vector<MathVector<TDomain::dim> > vVertPos;
	switch(elem->base_object_id())
	{
		case VERTEX: CollectCornerCoordinates(vVertPos, *static_cast<VertexBase*>(elem), domain, true); break;
		case EDGE: CollectCornerCoordinates(vVertPos, *static_cast<EdgeBase*>(elem), domain, true); break;
		case FACE: CollectCornerCoordinates(vVertPos, *static_cast<Face*>(elem), domain, true); break;
		case VOLUME: CollectCornerCoordinates(vVertPos, *static_cast<Volume*>(elem), domain, true); break;
		default: UG_THROW( "Base Object type not found.");
	}

//	forward
	return DoFPosition<TDomain::dim>(vPos, roid, vVertPos, lfeID);
}

#ifdef UG_DIM_1
template bool InnerDoFPosition<Domain1d>(std::vector<MathVector<1> >& vPos, GeometricObject* elem, const Domain1d& domain, const LFEID& lfeID);
template bool InnerDoFPosition<1>(std::vector<MathVector<1> >& vPos, const ReferenceObjectID roid,const std::vector<MathVector<1> >& vCornerCoord, const LFEID& lfeID);
template bool DoFPosition<Domain1d>(std::vector<MathVector<1> >& vPos, GeometricObject* elem, const Domain1d& domain, const LFEID& lfeID);
template bool DoFPosition<1>(std::vector<MathVector<1> >& vPos, const ReferenceObjectID roid,const std::vector<MathVector<1> >& vCornerCoord, const LFEID& lfeID);
#endif
#ifdef UG_DIM_2
template bool InnerDoFPosition<Domain2d>(std::vector<MathVector<2> >& vPos, GeometricObject* elem, const Domain2d& domain, const LFEID& lfeID);
template bool InnerDoFPosition<2>(std::vector<MathVector<2> >& vPos, const ReferenceObjectID roid,const std::vector<MathVector<2> >& vCornerCoord, const LFEID& lfeID);
template bool DoFPosition<Domain2d>(std::vector<MathVector<2> >& vPos, GeometricObject* elem, const Domain2d& domain, const LFEID& lfeID);
template bool DoFPosition<2>(std::vector<MathVector<2> >& vPos, const ReferenceObjectID roid,const std::vector<MathVector<2> >& vCornerCoord, const LFEID& lfeID);
#endif
#ifdef UG_DIM_3
template bool InnerDoFPosition<Domain3d>(std::vector<MathVector<3> >& vPos, GeometricObject* elem, const Domain3d& domain, const LFEID& lfeID);
template bool InnerDoFPosition<3>(std::vector<MathVector<3> >& vPos, const ReferenceObjectID roid,const std::vector<MathVector<3> >& vCornerCoord, const LFEID& lfeID);
template bool DoFPosition<Domain3d>(std::vector<MathVector<3> >& vPos, GeometricObject* elem, const Domain3d& domain, const LFEID& lfeID);
template bool DoFPosition<3>(std::vector<MathVector<3> >& vPos, const ReferenceObjectID roid,const std::vector<MathVector<3> >& vCornerCoord, const LFEID& lfeID);
#endif


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
