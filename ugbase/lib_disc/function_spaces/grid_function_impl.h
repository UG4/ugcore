/*
 * grid_function_impl.h
 *
 *  Created on: 13.06.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__FUNCTION_SPACE__GRID_FUNCTION_IMPL__
#define __H__UG__LIB_DISC__FUNCTION_SPACE__GRID_FUNCTION_IMPL__

#include "grid_function.h"

#include "lib_algebra/algebra_type.h"
#include "lib_disc/local_finite_element/local_shape_function_set.h"
#include "lib_disc/local_finite_element/local_dof_set.h"
#include "lib_disc/domain_util.h"

#ifdef UG_PARALLEL
	#include "pcl/pcl.h"
	#include "lib_algebra/parallelization/parallelization.h"
#endif

namespace ug{

////////////////////////////////////////////////////////////////////////////////
//	DoFPositions
////////////////////////////////////////////////////////////////////////////////

template <typename TElem, typename TDomain>
bool
InnerDoFPosition(std::vector<MathVector<TDomain::dim> >& vPos,
                 TElem* elem, TDomain& domain, LFEID lfeID)
{
//	reference element
	typedef typename reference_element_traits<TElem>::reference_element_type
			reference_element_type;

//	reference element dimension
	static const int refDim = reference_element_type::dim;

//	physical world dimension
	static const int dim = TDomain::dim;

//	reference object id
	static const ReferenceObjectID roid = reference_element_type::REFERENCE_OBJECT_ID;

//	vector for the vertex positions
	std::vector<MathVector<dim> > vVertPos(reference_element_type::num_corners);

//	get the vertices
	CollectCornerCoordinates(vVertPos, *elem, domain, true);

//	create a reference mapping
	ReferenceMapping<reference_element_type, dim> map(&(vVertPos[0]));

//	get local shape function set
	const DimLocalShapeFunctionSet<refDim>& lsfs
		= LocalShapeFunctionSetProvider::get<refDim>(roid, lfeID);

//	get local dof set
	const ILocalDoFSet& lds = LocalDoFSetProvider::get(roid, lfeID);

//	typedef local position type
	typedef typename DimLocalShapeFunctionSet<refDim>::position_type
		local_pos_type;

//	clear pos
	vPos.clear();

//	bool flag if position is exact, or no exact position available for shapes
	bool bExact = true;

//	loop all shape functions
	for(size_t sh = 0; sh < lsfs.num_sh(); ++sh)
	{
	//	check if dof in interior
		if(lds.local_dof(sh).dim() != refDim) continue;

	//	get local position
		local_pos_type locPos;
		bExact &= lsfs.position(sh, locPos);

	//	map to global position
		MathVector<dim> globPos;
		map.local_to_global(globPos, locPos);

	//	add
		vPos.push_back(globPos);
	}

//	return if positions are given exactly
	return bExact;
};

template <typename TDomain>
bool
DoFPositionOnVertex(std::vector<MathVector<TDomain::dim> >& vPos, VertexBase* elem, TDomain& domain, LFEID lfeID)
{
//\todo: handle finite element spaces with more than on DoF in Vertex

//	get the vertices
	CollectCornerCoordinates(vPos, *elem, domain, true);

//	return if positions are given exactly
	return true;
};

template <typename TDomain>
bool InnerDoFPosition(std::vector<MathVector<TDomain::dim> >& vPos, GeometricObject* elem, TDomain& domain, LFEID lfeID)
{
	switch(elem->base_object_type_id())
	{
		case VERTEX: return InnerDoFPosition(vPos, static_cast<VertexBase*>(elem), domain, lfeID);
		case EDGE:   return InnerDoFPosition(vPos, static_cast<EdgeBase*>(elem), domain, lfeID);
		case FACE:   return InnerDoFPosition(vPos, static_cast<Face*>(elem), domain, lfeID);
		case VOLUME: return InnerDoFPosition(vPos, static_cast<Volume*>(elem), domain, lfeID);
	}
	throw(UGFatalError("Base Object type not found."));
}

template <typename TDomain>
bool InnerDoFPosition(std::vector<MathVector<TDomain::dim> >& vPos, VertexBase* elem, TDomain& domain, LFEID lfeID)
{
	return  DoFPositionOnVertex(vPos, elem, domain, lfeID);
}

template <typename TDomain>
bool InnerDoFPosition(std::vector<MathVector<TDomain::dim> >& vPos, EdgeBase* elem, TDomain& domain, LFEID lfeID)
{
	switch(elem->shared_pipe_section())
	{
		case SPSEDGE_EDGE: 			   return InnerDoFPosition(vPos, static_cast<Edge*>(elem), domain, lfeID);
		case SPSEDGE_CONSTRAINED_EDGE: return InnerDoFPosition(vPos, static_cast<ConstrainedEdge*>(elem), domain, lfeID);
		case SPSEDGE_CONSTRAINING_EDGE:return InnerDoFPosition(vPos, static_cast<ConstrainingEdge*>(elem), domain, lfeID);
	}
	throw(UGFatalError("Edge type not found."));
}

template <typename TDomain>
bool InnerDoFPosition(std::vector<MathVector<TDomain::dim> >& vPos, Face* elem, TDomain& domain, LFEID lfeID)
{
	switch(elem->shared_pipe_section())
	{
		case SPSFACE_TRIANGLE: return InnerDoFPosition(vPos, static_cast<Triangle*>(elem), domain, lfeID);
		case SPSFACE_CONSTRAINED_TRIANGLE: return InnerDoFPosition(vPos, static_cast<ConstrainedTriangle*>(elem), domain, lfeID);
		case SPSFACE_CONSTRAINING_TRIANGLE: return InnerDoFPosition(vPos, static_cast<ConstrainingTriangle*>(elem), domain, lfeID);
		case SPSFACE_QUADRILATERAL: return InnerDoFPosition(vPos, static_cast<Quadrilateral*>(elem), domain, lfeID);
		case SPSFACE_CONSTRAINED_QUADRILATERAL: return InnerDoFPosition(vPos, static_cast<ConstrainedQuadrilateral*>(elem), domain, lfeID);
		case SPSFACE_CONSTRAINING_QUADRILATERAL: return InnerDoFPosition(vPos, static_cast<ConstrainingQuadrilateral*>(elem), domain, lfeID);
	}
	throw(UGFatalError("Face type not found."));
}

template <typename TDomain>
bool InnerDoFPosition(std::vector<MathVector<TDomain::dim> >& vPos, Volume* elem, TDomain& domain, LFEID lfeID)
{
	switch(elem->shared_pipe_section())
	{
		case SPSVOL_TETRAHEDRON: return InnerDoFPosition(vPos, static_cast<Tetrahedron*>(elem), domain, lfeID);
		case SPSVOL_PYRAMID: return InnerDoFPosition(vPos, static_cast<Pyramid*>(elem), domain, lfeID);
		case SPSVOL_PRISM: return InnerDoFPosition(vPos, static_cast<Prism*>(elem), domain, lfeID);
		case SPSVOL_HEXAHEDRON: return InnerDoFPosition(vPos, static_cast<Hexahedron*>(elem), domain, lfeID);
	}
	throw(UGFatalError("Volume type not found."));
}




template <typename TElem, typename TDomain>
bool DoFPosition(std::vector<MathVector<TDomain::dim> >& vPos, TElem* elem, TDomain& domain, LFEID lfeID)
{
//	reference element
	typedef typename reference_element_traits<TElem>::reference_element_type
			reference_element_type;

//	physical world dimension
	static const int dim = TDomain::dim;

//	reference object id
	static const ReferenceObjectID roid = reference_element_type::REFERENCE_OBJECT_ID;

//	vector for the vertex positions
	std::vector<MathVector<TDomain::dim> > vVertPos(reference_element_type::num_corners);

//	get the vertices
	CollectCornerCoordinates(vVertPos, *elem, domain, true);

//	create a reference mapping
	ReferenceMapping<reference_element_type, dim> map(&(vVertPos[0]));

//	get local shape function set
	const DimLocalShapeFunctionSet<dim>& lsfs
		= LocalShapeFunctionSetProvider::get<reference_element_type>(roid, lfeID);

//	typedef local position type
	typedef typename DimLocalShapeFunctionSet<dim>::position_type local_pos_type;

//	clear pos
	vPos.resize(lsfs.num_sh());

//	bool flag if position is exact, or no exact position available for shapes
	bool bExact = true;

//	loop all shape functions
	for(size_t sh = 0; sh < lsfs.num_sh(); ++sh)
	{
	//	get local position
		local_pos_type locPos;
		bExact &= lsfs.position(sh, locPos);

	//	map to global position
		map.local_to_global(vPos[sh], locPos);
	}

//	return if positions are given exactly
	return bExact;
};

template <typename TDomain>
bool DoFPosition(std::vector<MathVector<TDomain::dim> >& vPos, GeometricObject* elem, TDomain& domain, LFEID lfeID)
{
	switch(elem->base_object_type_id())
	{
		case VERTEX: return DoFPosition(vPos, static_cast<VertexBase*>(elem), domain, lfeID);
		case EDGE:   return DoFPosition(vPos, static_cast<EdgeBase*>(elem), domain, lfeID);
		case FACE:   return DoFPosition(vPos, static_cast<Face*>(elem), domain, lfeID);
		case VOLUME: return DoFPosition(vPos, static_cast<Volume*>(elem), domain, lfeID);
	}
	throw(UGFatalError("Base Object type not found."));
}

template <typename TDomain>
bool DoFPosition(std::vector<MathVector<TDomain::dim> >& vPos, VertexBase* elem, TDomain& domain, LFEID lfeID)
{
	return  DoFPositionOnVertex(vPos, elem, domain, lfeID);
}

template <typename TDomain>
bool DoFPosition(std::vector<MathVector<TDomain::dim> >& vPos, EdgeBase* elem, TDomain& domain, LFEID lfeID)
{
	switch(elem->shared_pipe_section())
	{
		case SPSEDGE_EDGE: 			   return DoFPosition(vPos, static_cast<Edge*>(elem), domain, lfeID);
		case SPSEDGE_CONSTRAINED_EDGE: return DoFPosition(vPos, static_cast<ConstrainedEdge*>(elem), domain, lfeID);
		case SPSEDGE_CONSTRAINING_EDGE:return DoFPosition(vPos, static_cast<ConstrainingEdge*>(elem), domain, lfeID);
	}
	throw(UGFatalError("Edge type not found."));
}

template <typename TDomain>
bool DoFPosition(std::vector<MathVector<TDomain::dim> >& vPos, Face* elem, TDomain& domain, LFEID lfeID)
{
	switch(elem->shared_pipe_section())
	{
		case SPSFACE_TRIANGLE: return DoFPosition(vPos, static_cast<Triangle*>(elem), domain, lfeID);
		case SPSFACE_CONSTRAINED_TRIANGLE: return DoFPosition(vPos, static_cast<ConstrainedTriangle*>(elem), domain, lfeID);
		case SPSFACE_CONSTRAINING_TRIANGLE: return DoFPosition(vPos, static_cast<ConstrainingTriangle*>(elem), domain, lfeID);
		case SPSFACE_QUADRILATERAL: return DoFPosition(vPos, static_cast<Quadrilateral*>(elem), domain, lfeID);
		case SPSFACE_CONSTRAINED_QUADRILATERAL: return DoFPosition(vPos, static_cast<ConstrainedQuadrilateral*>(elem), domain, lfeID);
		case SPSFACE_CONSTRAINING_QUADRILATERAL: return DoFPosition(vPos, static_cast<ConstrainingQuadrilateral*>(elem), domain, lfeID);
	}
	throw(UGFatalError("Face type not found."));
}

template <typename TDomain>
bool DoFPosition(std::vector<MathVector<TDomain::dim> >& vPos, Volume* elem, TDomain& domain, LFEID lfeID)
{
	switch(elem->shared_pipe_section())
	{
		case SPSVOL_TETRAHEDRON: return DoFPosition(vPos, static_cast<Tetrahedron*>(elem), domain, lfeID);
		case SPSVOL_PYRAMID: return DoFPosition(vPos, static_cast<Pyramid*>(elem), domain, lfeID);
		case SPSVOL_PRISM: return DoFPosition(vPos, static_cast<Prism*>(elem), domain, lfeID);
		case SPSVOL_HEXAHEDRON: return DoFPosition(vPos, static_cast<Hexahedron*>(elem), domain, lfeID);
	}
	throw(UGFatalError("Volume type not found."));
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template <typename TDomain, typename TDD, typename TAlgebra>
GridFunction<TDomain, TDD, TAlgebra>::
GridFunction(SmartPtr<approximation_space_type> approxSpace,
             SmartPtr<TDD> spDoFDistr)
 : IDDGridFunction<TDD>(spDoFDistr), m_spDomain(approxSpace->domain())
{
	check_algebra();
	resize_values(num_indices());
#ifdef UG_PARALLEL
//	set layouts
	copy_layouts_into_vector();

//	set storage type
	this->set_storage_type(PST_UNDEFINED);
#endif
};

template <typename TDomain, typename TDD, typename TAlgebra>
GridFunction<TDomain, TDD, TAlgebra>::
GridFunction(SmartPtr<approximation_space_type> approxSpace)
	: IDDGridFunction<TDD>(approxSpace->surface_dof_distribution()),
	  m_spDomain(approxSpace->domain())
{
	check_algebra();
	resize_values(num_indices());
#ifdef UG_PARALLEL
//	set layouts
	copy_layouts_into_vector();

//	set storage type
	this->set_storage_type(PST_UNDEFINED);
#endif
};

template <typename TDomain, typename TDD, typename TAlgebra>
void
GridFunction<TDomain, TDD, TAlgebra>::check_algebra()
{
//	get blocksize of algebra
	const int blockSize = algebra_type::blockSize;

//	a)	If blocksize fixed and > 1, we need grouping.
	if(blockSize > 1 && !this->m_spDD->grouped())
	{
		UG_THROW_FATAL("Fixed block algebra needs grouped dofs.");
	}
//	b) 	If blocksize flexible, we group
	else if (blockSize == AlgebraType::VariableBlockSize
			&& !this->m_spDD->grouped())
	{
		UG_THROW_FATAL("Variable block algebra needs grouped dofs.");
	}
//	c)	If blocksize == 1, we do not group. This will allow us to handle
//		this case for any problem.
	else if (blockSize == 1 && this->m_spDD->grouped())
	{
		UG_THROW_FATAL("block 1x1 algebra needs non-grouped dofs.");
	}
}

template <typename TDomain, typename TDD, typename TAlgebra>
template <typename TElem>
bool
GridFunction<TDomain, TDD, TAlgebra>::
dof_positions(TElem* elem, size_t fct, std::vector<MathVector<dim> >& vPos) const
{
	return DoFPosition(vPos, elem, *domain(), this->local_finite_element_id(fct));
};

template <typename TDomain, typename TDD, typename TAlgebra>
template <typename TElem>
bool
GridFunction<TDomain, TDD, TAlgebra>::
inner_dof_positions(TElem* elem, size_t fct, std::vector<MathVector<dim> >& vPos) const
{
	return InnerDoFPosition(vPos, elem, *domain(), this->local_finite_element_id(fct));
};

template <typename TDomain, typename TDD, typename TAlgebra>
void
GridFunction<TDomain, TDD, TAlgebra>::
clone_pattern(const this_type& v)
{
// 	copy approximation space
	m_spDomain = v.m_spDomain;

//	assign dof distribution (resizes vector)
	this->m_spDD = v.m_spDD;

//	resize the vector
	resize_values(num_indices());

#ifdef UG_PARALLEL
//	set layouts
	copy_layouts_into_vector();

//	copy storage type
	vector_type::copy_storage_type(v);
#endif
};

template <typename TDomain, typename TDD, typename TAlgebra>
void
GridFunction<TDomain, TDD, TAlgebra>::
resize_values(size_t s, number defaultValue)
{
//	remember old values
	const size_t oldSize = vector_type::size();

//	resize vector
	vector_type::resize(s);

//	set vector to zero-values
	for(size_t i = oldSize; i < s; ++i)
		this->operator[](i) = defaultValue;
}

template <typename TDomain, typename TDD, typename TAlgebra>
void
GridFunction<TDomain, TDD, TAlgebra>::
permute_values(const std::vector<size_t>& vIndNew)
{
//	check sizes
	if(vIndNew.size() != this->size())
		UG_THROW_FATAL("GridFunction::permute_values: For a permutation the"
				" index set must have same cardinality as vector.");

// \todo: avoid tmp vector, only copy values into new vector and use that one
//	create tmp vector
	vector_type vecTmp; vecTmp.resize(this->size());

//	loop indices and copy values
	for(size_t i = 0; i < vIndNew.size(); ++i)
		vecTmp[vIndNew[i]] = this->operator[](i);

//	copy tmp vector into this vector
	this->assign(vecTmp);
}

template <typename TDomain, typename TDD, typename TAlgebra>
void
GridFunction<TDomain, TDD, TAlgebra>::
copy_values(const std::vector<std::pair<size_t, size_t> >& vIndexMap,bool bDisjunct)
{
//	disjunct case
	if(bDisjunct)
		for(size_t i = 0; i < vIndexMap.size(); ++i)
			this->operator[](vIndexMap[i].second)
				= this->operator[](vIndexMap[i].first);

//	other case not implemented
	else UG_THROW_FATAL("Not implemented.");
}

template <typename TDomain, typename TDD, typename TAlgebra>
void GridFunction<TDomain, TDD, TAlgebra>::assign(const vector_type& v)
{
//	check size
	if(v.size() != vector_type::size())
		UG_THROW_FATAL("GridFunction: Assigned vector has incorrect size.");

//	assign vector
	*(dynamic_cast<vector_type*>(this)) = v;

#ifdef UG_PARALLEL
//	set layouts
	copy_layouts_into_vector();

//	copy storage type
	vector_type::copy_storage_type(v);
#endif
}

template <typename TDomain, typename TDD, typename TAlgebra>
void GridFunction<TDomain, TDD, TAlgebra>::assign(const this_type& v)
{
// 	copy approximation space
	m_spDomain = v.m_spDomain;

//	assign dof distribution (resizes vector)
	this->m_spDD = v.m_spDD;

//	resize the vector
	resize_values(num_indices());

//  copy values
	*(dynamic_cast<vector_type*>(this)) = *dynamic_cast<const vector_type*>(&v);

#ifdef UG_PARALLEL
//	set layouts
	copy_layouts_into_vector();

//	copy storage type
	vector_type::copy_storage_type(v);
#endif
}

} // end namespace ug

#endif /* __H__UG__LIB_DISC__FUNCTION_SPACE__GRID_FUNCTION_IMPL__ */
