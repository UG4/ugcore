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
bool InnerDoFPositionElem(std::vector<MathVector<dim> >& vPos, const ReferenceObjectID roid,
                      const std::vector<MathVector<dim> >& vVertPos, const LFEID& lfeID)
{
//	get local dof set
	const DimLocalDoFSet<refDim>& lds =
					LocalFiniteElementProvider::get_dofs<refDim>(roid, lfeID);

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
		case EDGE:   return InnerDoFPositionElem<1,dim>(vPos, roid, vCornerCoord, lfeID);
		case FACE:   return InnerDoFPositionElem<2,dim>(vPos, roid, vCornerCoord, lfeID);
		case VOLUME: return InnerDoFPositionElem<3,dim>(vPos, roid, vCornerCoord, lfeID);
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
bool DoFPositionElem(std::vector<MathVector<dim> >& vPos, const ReferenceObjectID roid,
                 const std::vector<MathVector<dim> >& vVertPos, const LFEID& lfeID)
{
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
		case EDGE:   return DoFPositionElem<1,dim>(vPos, roid, vCornerCoord, lfeID);
		case FACE:   return DoFPositionElem<2,dim>(vPos, roid, vCornerCoord, lfeID);
		case VOLUME: return DoFPositionElem<3,dim>(vPos, roid, vCornerCoord, lfeID);
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


////////////////////////////////////////////////////////////////////////////////
//	ShapesAtGlobalPosition
////////////////////////////////////////////////////////////////////////////////

template <int dim>
void ShapesAtGlobalPositionVertex(std::vector<std::vector<number> >& vvShape,
                                  const std::vector<MathVector<dim> >& vGlobPos,
                                  const LFEID& lfeID)
{
//	get local position of DoF
	std::vector<MathVector<0> > vLocPos(vGlobPos.size(), 0.0);

//	evaluate coarse shape fct at fine local point
	try{
		const LocalShapeFunctionSet<0>& lsfs =
				LocalFiniteElementProvider::get<0>(ROID_VERTEX, lfeID);
		lsfs.shapes(vvShape, vLocPos);
	}
	UG_CATCH_THROW("ShapesAtGlobalPosition: Cannot evalute shapes.")
}

template <int refDim, int dim>
void ShapesAtGlobalPositionElem(std::vector<std::vector<number> >& vvShape,
                                const std::vector<MathVector<dim> >& vGlobPos,
                                const ReferenceObjectID roid,
                                const std::vector<MathVector<dim> >& vCornerCoord,
                                const LFEID& lfeID)
{
//	get local position of DoF
	std::vector<MathVector<refDim> > vLocPos(vGlobPos.size(), 0.0);
	try{
		DimReferenceMapping<refDim, dim>& map =
				ReferenceMappingProvider::get<refDim, dim>(roid, vCornerCoord);
		map.global_to_local(vLocPos, vGlobPos);
	}
	UG_CATCH_THROW("ShapesAtGlobalPosition: Cannot find elem-local Positions.");

//	evaluate coarse shape fct at fine local point
	try{
		const LocalShapeFunctionSet<refDim>& lsfs =
				LocalFiniteElementProvider::get<refDim>(roid, lfeID);
		lsfs.shapes(vvShape, vLocPos);
	}
	UG_CATCH_THROW("ShapesAtGlobalPosition: Cannot evalute shapes.")
}

template <int dim>
void ShapesAtGlobalPosition(std::vector<std::vector<number> >& vvShape,
                           const std::vector<MathVector<dim> >& vGlobPos,
                           const ReferenceObjectID roid,
                           const std::vector<MathVector<dim> >& vCornerCoord,
                           const LFEID& lfeID)
{
	switch(ReferenceElementDimension(roid))
	{
		case VERTEX: return ShapesAtGlobalPositionVertex<dim>(vvShape, vGlobPos, lfeID);
		case EDGE:   return ShapesAtGlobalPositionElem<1,dim>(vvShape, vGlobPos, roid, vCornerCoord, lfeID);
		case FACE:   return ShapesAtGlobalPositionElem<2,dim>(vvShape, vGlobPos, roid, vCornerCoord, lfeID);
		case VOLUME: return ShapesAtGlobalPositionElem<3,dim>(vvShape, vGlobPos, roid, vCornerCoord, lfeID);
		default: UG_THROW("Base Object type not found.");
	}
}

template <typename TDomain>
void ShapesAtGlobalPosition(std::vector<std::vector<number> >& vvShape,
                           const std::vector<MathVector<TDomain::dim> >& vGlobPos,
                           GeometricObject* elem, const TDomain& domain, const LFEID& lfeID)
{
	const int baseDim = elem->base_object_id();
	if(baseDim == VERTEX)
		return ShapesAtGlobalPositionVertex<TDomain::dim>(vvShape, vGlobPos, lfeID);

//	get the vertices
	std::vector<MathVector<TDomain::dim> > vCornerCoord;
	switch(baseDim)
	{
		case EDGE: CollectCornerCoordinates(vCornerCoord, *static_cast<EdgeBase*>(elem), domain, true); break;
		case FACE: CollectCornerCoordinates(vCornerCoord, *static_cast<Face*>(elem), domain, true); break;
		case VOLUME: CollectCornerCoordinates(vCornerCoord, *static_cast<Volume*>(elem), domain, true); break;
		default: UG_THROW( "Base Object type not found.");
	}

//	reference object id
	const ReferenceObjectID roid = elem->reference_object_id();

//	forward
	return ShapesAtGlobalPosition<TDomain::dim>(vvShape, vGlobPos, roid, vCornerCoord, lfeID);
}


////////////////////////////////////////////////////////////////////////////////
//	Extract Positions
////////////////////////////////////////////////////////////////////////////////

template <typename TDomain>
void ExtractPositionsVertex(ConstSmartPtr<TDomain> domain,
                            ConstSmartPtr<DoFDistribution> dd,
                            std::vector<MathVector<TDomain::dim> >& vPos)
{
//	get position accessor
	const typename TDomain::position_accessor_type& aaPos = domain->position_accessor();

//	iterator
	typename DoFDistribution::traits<VertexBase>::const_iterator iter, iterEnd;

//	algebra indices vector
	std::vector<size_t> ind;

//	get iterators
	iter = dd->begin<VertexBase>(SurfaceView::ALL);
	iterEnd = dd->end<VertexBase>(SurfaceView::ALL);

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
			const size_t index = ind[i];
			vPos[index] = aaPos[v];
		}
	}
}

template <typename TDomain, typename TBaseElem>
void ExtractPositionsElem(ConstSmartPtr<TDomain> domain,
                          ConstSmartPtr<DoFDistribution> dd,
                          std::vector<MathVector<TDomain::dim> >& vPos)
{
	typename DoFDistribution::traits<TBaseElem>::const_iterator iter, iterEnd;

//	vector for positions
	std::vector<MathVector<TDomain::dim> > vElemPos;

//	algebra indices vector
	std::vector<DoFIndex> ind;

//	loop all subsets
	for(int si = 0; si < dd->num_subsets(); ++si)
	{
	//	get iterators
		iter = dd->begin<TBaseElem>(si, SurfaceView::ALL);
		iterEnd = dd->end<TBaseElem>(si, SurfaceView::ALL);

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
				dd->inner_dof_indices(elem, fct, ind);

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
					vPos[index] = vElemPos[sh];
				}
			}
		}
	}
}

template <typename TDomain>
void ExtractPositions(ConstSmartPtr<TDomain> domain,
                      ConstSmartPtr<DoFDistribution> dd,
                      std::vector<MathVector<TDomain::dim> >& vPos)
{
//	number of total dofs
	int nr = dd->num_indices();

//	resize positions
	vPos.resize(nr);

//	extract for all element types
	if(dd->max_dofs(VERTEX)) ExtractPositionsVertex<TDomain>(domain, dd, vPos);
	if(dd->max_dofs(EDGE)) ExtractPositionsElem<TDomain, EdgeBase>(domain, dd, vPos);
	if(dd->max_dofs(FACE)) ExtractPositionsElem<TDomain, Face>(domain, dd, vPos);
	if(dd->max_dofs(VOLUME)) ExtractPositionsElem<TDomain, Volume>(domain, dd, vPos);
}


template <typename TDomain, typename TBaseElem>
void ExtractAlgebraIndices2(ConstSmartPtr<TDomain> domain,
                  ConstSmartPtr<DoFDistribution> dd,
                  std::vector<size_t>& fctIndex)
{
	typename DoFDistribution::traits<TBaseElem>::const_iterator iter, iterEnd;

//	algebra indices vector
	std::vector<size_t> ind;

//	loop all subsets
	for(int si = 0; si < dd->num_subsets(); ++si)
	{
	//	get iterators
		iter = dd->begin<TBaseElem>(si, SurfaceView::ALL);
		iterEnd = dd->end<TBaseElem>(si, SurfaceView::ALL);

	//	loop all elements
		for(;iter != iterEnd; ++iter)
		{
		//	get element
			TBaseElem* elem = *iter;

		//	loop all functions
			for(size_t fct = 0; fct < dd->num_fct(); ++fct)
			{
			//	skip non-used function
				if(!dd->is_def_in_subset(fct,si)) continue;

			//	load indices associated with element function
				dd->inner_algebra_indices_for_fct(elem, ind, true, fct);

				for(size_t sh = 0; sh < ind.size(); ++sh)
				{
					const size_t index = ind[sh];
					if(index < fctIndex.size())
						fctIndex[index] = fct;
				}
			}
		}
	}
}


template <typename TDomain>
void ExtractAlgebraIndices(ConstSmartPtr<TDomain> domain,
                      ConstSmartPtr<DoFDistribution> dd,
                      std::vector<size_t> &fctIndex)
{
//	number of total dofs
	int nr = dd->num_indices();

//	resize positions
	fctIndex.resize(nr);

//	extract for all element types
	if(dd->max_dofs(VERTEX)) ExtractAlgebraIndices2<TDomain, VertexBase>(domain, dd, fctIndex);
	if(dd->max_dofs(EDGE)) ExtractAlgebraIndices2<TDomain, EdgeBase>(domain, dd, fctIndex);
	if(dd->max_dofs(FACE)) ExtractAlgebraIndices2<TDomain, Face>(domain, dd, fctIndex);
	if(dd->max_dofs(VOLUME)) ExtractAlgebraIndices2<TDomain, Volume>(domain, dd, fctIndex);
}


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
	const_iterator iter = dd->begin<VertexBase>(SurfaceView::ALL);
	const_iterator iterEnd = dd->end<VertexBase>(SurfaceView::ALL);

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
	std::vector<DoFIndex> ind;

//	loop all subsets
	for(int si = 0; si < dd->num_subsets(); ++si)
	{
	//	get iterators
		iter = dd->begin<TBaseElem>(si, SurfaceView::ALL);
		iterEnd = dd->end<TBaseElem>(si, SurfaceView::ALL);

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
				dd->inner_dof_indices(elem, fct, ind);

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
	std::vector<DoFIndex> ind;

//	a pair
	std::pair<MathVector<TDomain::dim>, size_t> pair;

//	loop all subsets
	for(int si = 0; si < dd->num_subsets(); ++si)
	{
	//	get iterators
		iter = dd->begin<TBaseElem>(si, SurfaceView::ALL);
		iterEnd = dd->end<TBaseElem>(si, SurfaceView::ALL);

	//	skip non-used function
		if(!dd->is_def_in_subset(fct,si)) continue;

	//	loop all elements
		for(;iter != iterEnd; ++iter)
		{
		//	get vertex
			TBaseElem* elem = *iter;

		//	load indices associated with element function
			dd->inner_dof_indices(elem, fct, ind);

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

////////////////////////////////////////////////////////////////////////////////
//	Checks correct DoF Positions
////////////////////////////////////////////////////////////////////////////////

template <typename TDomain, typename TBaseElem>
bool CheckDoFElem(ConstSmartPtr<TDomain> domain,
                  ConstSmartPtr<DoFDistribution> dd,
                  std::vector<MathVector<TDomain::dim> >& vPos)
{
	typename DoFDistribution::traits<TBaseElem>::const_iterator iter, iterEnd;

	bool bRes = true;

//	vector for positions
	std::vector<MathVector<TDomain::dim> > vElemPos;

//	algebra indices vector
	std::vector<DoFIndex> ind;

//	loop all subsets
	for(int si = 0; si < dd->num_subsets(); ++si)
	{
	//	get iterators
		iter = dd->begin<TBaseElem>(si, SurfaceView::ALL);
		iterEnd = dd->end<TBaseElem>(si, SurfaceView::ALL);

	//	loop all elements
		for(;iter != iterEnd; ++iter)
		{
		//	get element
			TBaseElem* elem = *iter;

		//	loop all functions
			for(size_t fct = 0; fct < dd->num_fct(); ++fct)
			{
			//	skip non-used function
				if(!dd->is_def_in_subset(fct,si)) continue;

			//	load indices associated with element function
				dd->inner_dof_indices(elem, fct, ind);

			//	load positions associated with element and function
				InnerDoFPosition(vElemPos, elem, *(const_cast<TDomain*>(domain.get())),
				                 dd->local_finite_element_id(fct));

				bool bWrite = false;
				for(size_t sh = 0; sh < ind.size(); ++sh)
				{
					size_t index = ind[sh][0];

					if(vPos[index] != MathVector<TDomain::dim>(-1)){
							if(VecDistance(vPos[index], vElemPos[sh]) < 1e-10) continue;

						if(!bWrite)
							UG_LOG(" **** inner_multi_index (start) ******\n")
						bWrite = true;
						bRes = false;
						UG_LOG("CheckDoFPositions "<<sh<<": inner_dof_indices: index: "
						       <<index<<" at "<<vElemPos[sh]<<", but previously: "
						       <<vPos[index]<<"\n");
					}

					vPos[index] = vElemPos[sh];

				}

				if(bWrite){
					std::vector<MathVector<TDomain::dim> > vVertPos;
					CollectCornerCoordinates(vVertPos, elem, *domain, true);
					UG_LOG("From Elem "<<elem<<" ("<<elem->reference_object_id()<<"):\n")
					for(size_t i = 0; i < vVertPos.size(); ++i)
						UG_LOG("with corner "<<i<<": "<<vVertPos[i]<<"\n");
					UG_LOG(" **** inner_multi_index (end) ******\n")
				}

				 /////////////////////////
				//	load indices associated with element function
					dd->dof_indices(elem, fct, ind);

				//	load positions associated with element and function
					DoFPosition(vElemPos, elem, *(const_cast<TDomain*>(domain.get())),
					                 dd->local_finite_element_id(fct));

				//	write position
					bWrite = false;
					for(size_t sh = 0; sh < ind.size(); ++sh)
					{
						size_t index = ind[sh][0];

						if(vPos[index] != MathVector<TDomain::dim>(-1)){
							if(VecDistance(vPos[index], vElemPos[sh]) < 1e-10) continue;

							if(!bWrite)
								UG_LOG(" **** multi_index (start) ******\n")
							bWrite = true;
							bRes = false;
							UG_LOG("CheckDoFPositions "<<sh<<": dof_indices: index: "
								   <<index<<" at "<<vElemPos[sh]<<", but previously: "
								   <<vPos[index]<<"\n");
						}

						vPos[index] = vElemPos[sh];
					}

					if(bWrite){
						std::vector<MathVector<TDomain::dim> > vVertPos;
						CollectCornerCoordinates(vVertPos, elem, *domain, true);
						UG_LOG("From Elem "<<elem<<" ("<<elem->reference_object_id()<<"):\n")
						for(size_t i = 0; i < vVertPos.size(); ++i)
							UG_LOG("with corner "<<i<<": "<<vVertPos[i]<<"\n");
						UG_LOG(" **** multi_index (end) ******\n")
					}
			}
		}
	}

	return bRes;
}

template <typename TDomain>
bool CheckDoFPositions(ConstSmartPtr<TDomain> domain,
                       ConstSmartPtr<DoFDistribution> dd)
{
//	number of total dofs
	int nr = dd->num_indices();
	std::vector<MathVector<TDomain::dim> > vPos(nr, -1);

//	extract for all element types
	bool bRes = true;
	if(dd->max_dofs(VERTEX)) bRes &= CheckDoFElem<TDomain, VertexBase>(domain, dd, vPos);
	if(dd->max_dofs(EDGE)) bRes &= CheckDoFElem<TDomain, EdgeBase>(domain, dd, vPos);
	if(dd->max_dofs(FACE)) bRes &= CheckDoFElem<TDomain, Face>(domain, dd, vPos);
	if(dd->max_dofs(VOLUME)) bRes &= CheckDoFElem<TDomain, Volume>(domain, dd, vPos);
	return bRes;
}

#ifdef UG_DIM_1
template bool InnerDoFPosition<Domain1d>(std::vector<MathVector<1> >& vPos, GeometricObject* elem, const Domain1d& domain, const LFEID& lfeID);
template bool InnerDoFPosition<1>(std::vector<MathVector<1> >& vPos, const ReferenceObjectID roid,const std::vector<MathVector<1> >& vCornerCoord, const LFEID& lfeID);
template bool DoFPosition<Domain1d>(std::vector<MathVector<1> >& vPos, GeometricObject* elem, const Domain1d& domain, const LFEID& lfeID);
template bool DoFPosition<1>(std::vector<MathVector<1> >& vPos, const ReferenceObjectID roid,const std::vector<MathVector<1> >& vCornerCoord, const LFEID& lfeID);
template void ExtractPositions(ConstSmartPtr<Domain1d> domain, ConstSmartPtr<DoFDistribution> dd, std::vector<MathVector<Domain1d::dim> >& vPos);
template void ExtractPositions(ConstSmartPtr<Domain1d> domain, ConstSmartPtr<DoFDistribution> dd, std::vector<std::pair<MathVector<Domain1d::dim>, size_t> >& vPos);
template void ExtractPositions(ConstSmartPtr<Domain1d> domain, ConstSmartPtr<DoFDistribution> dd, const size_t fct, std::vector<std::pair<MathVector<Domain1d::dim>, size_t> >& vPos);
template bool CheckDoFPositions(ConstSmartPtr<Domain1d> domain, ConstSmartPtr<DoFDistribution> dd);
template void ExtractAlgebraIndices<Domain1d>(ConstSmartPtr<Domain1d> domain, ConstSmartPtr<DoFDistribution> dd, std::vector<size_t> &fctIndex);
template void ShapesAtGlobalPosition<1>(std::vector<std::vector<number> >& vvShape, const std::vector<MathVector<1> >& vGlobPos, const ReferenceObjectID roid,const std::vector<MathVector<1> >& vCornerCoord, const LFEID& lfeID);
template void ShapesAtGlobalPosition<Domain1d>(std::vector<std::vector<number> >& vvShape, const std::vector<MathVector<1> >& vGlobPos,GeometricObject* elem, const Domain1d& domain, const LFEID& lfeID);
#endif
#ifdef UG_DIM_2
template bool InnerDoFPosition<Domain2d>(std::vector<MathVector<2> >& vPos, GeometricObject* elem, const Domain2d& domain, const LFEID& lfeID);
template bool InnerDoFPosition<2>(std::vector<MathVector<2> >& vPos, const ReferenceObjectID roid,const std::vector<MathVector<2> >& vCornerCoord, const LFEID& lfeID);
template bool DoFPosition<Domain2d>(std::vector<MathVector<2> >& vPos, GeometricObject* elem, const Domain2d& domain, const LFEID& lfeID);
template bool DoFPosition<2>(std::vector<MathVector<2> >& vPos, const ReferenceObjectID roid,const std::vector<MathVector<2> >& vCornerCoord, const LFEID& lfeID);
template void ExtractPositions(ConstSmartPtr<Domain2d> domain, ConstSmartPtr<DoFDistribution> dd, std::vector<MathVector<Domain2d::dim> >& vPos);
template void ExtractPositions(ConstSmartPtr<Domain2d> domain, ConstSmartPtr<DoFDistribution> dd, std::vector<std::pair<MathVector<Domain2d::dim>, size_t> >& vPos);
template void ExtractPositions(ConstSmartPtr<Domain2d> domain, ConstSmartPtr<DoFDistribution> dd, const size_t fct, std::vector<std::pair<MathVector<Domain2d::dim>, size_t> >& vPos);
template bool CheckDoFPositions(ConstSmartPtr<Domain2d> domain, ConstSmartPtr<DoFDistribution> dd);
template void ExtractAlgebraIndices<Domain2d>(ConstSmartPtr<Domain2d> domain, ConstSmartPtr<DoFDistribution> dd, std::vector<size_t> &fctIndex);
template void ShapesAtGlobalPosition<2>(std::vector<std::vector<number> >& vvShape, const std::vector<MathVector<2> >& vGlobPos, const ReferenceObjectID roid,const std::vector<MathVector<2> >& vCornerCoord, const LFEID& lfeID);
template void ShapesAtGlobalPosition<Domain2d>(std::vector<std::vector<number> >& vvShape, const std::vector<MathVector<2> >& vGlobPos,GeometricObject* elem, const Domain2d& domain, const LFEID& lfeID);
#endif
#ifdef UG_DIM_3
template bool InnerDoFPosition<Domain3d>(std::vector<MathVector<3> >& vPos, GeometricObject* elem, const Domain3d& domain, const LFEID& lfeID);
template bool InnerDoFPosition<3>(std::vector<MathVector<3> >& vPos, const ReferenceObjectID roid,const std::vector<MathVector<3> >& vCornerCoord, const LFEID& lfeID);
template bool DoFPosition<Domain3d>(std::vector<MathVector<3> >& vPos, GeometricObject* elem, const Domain3d& domain, const LFEID& lfeID);
template bool DoFPosition<3>(std::vector<MathVector<3> >& vPos, const ReferenceObjectID roid,const std::vector<MathVector<3> >& vCornerCoord, const LFEID& lfeID);
template void ExtractPositions(ConstSmartPtr<Domain3d> domain, ConstSmartPtr<DoFDistribution> dd, std::vector<MathVector<Domain3d::dim> >& vPos);
template void ExtractPositions(ConstSmartPtr<Domain3d> domain, ConstSmartPtr<DoFDistribution> dd, std::vector<std::pair<MathVector<Domain3d::dim>, size_t> >& vPos);
template void ExtractPositions(ConstSmartPtr<Domain3d> domain, ConstSmartPtr<DoFDistribution> dd, const size_t fct, std::vector<std::pair<MathVector<Domain3d::dim>, size_t> >& vPos);
template bool CheckDoFPositions(ConstSmartPtr<Domain3d> domain, ConstSmartPtr<DoFDistribution> dd);
template void ExtractAlgebraIndices<Domain3d>(ConstSmartPtr<Domain3d> domain, ConstSmartPtr<DoFDistribution> dd, std::vector<size_t> &fctIndex);
template void ShapesAtGlobalPosition<3>(std::vector<std::vector<number> >& vvShape, const std::vector<MathVector<3> >& vGlobPos, const ReferenceObjectID roid,const std::vector<MathVector<3> >& vCornerCoord, const LFEID& lfeID);
template void ShapesAtGlobalPosition<Domain3d>(std::vector<std::vector<number> >& vvShape, const std::vector<MathVector<3> >& vGlobPos,GeometricObject* elem, const Domain3d& domain, const LFEID& lfeID);
#endif

} // end namespace ug
