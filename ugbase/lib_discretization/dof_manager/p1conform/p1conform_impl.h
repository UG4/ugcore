/*
 * p1conform_impl.h
 *
 *  Created on: 13.07.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__DOF_MANAGER__P1CONFORM_IMPL__
#define __H__LIB_DISCRETIZATION__DOF_MANAGER__P1CONFORM_IMPL__

#include <vector>

#include "./p1conform.h"

namespace ug{

///////////////////////////////////////
// P1ConformDoFDistribution
///////////////////////////////////////

///////////// LocalIndex access /////////////////

template<typename TElem>
bool
P1ConformDoFDistribution::
has_dofs_on() const
{
//	get base obj type
	uint type = geometry_traits<TElem>::BASE_OBJECT_TYPE_ID;

//	only in case of a Vertex, we have a DoF
	if(type == VERTEX) return true;
	else return false;
}

template<typename TElem>
size_t
P1ConformDoFDistribution::
num_indices(int si, const FunctionGroup& fctGrp) const
{
//	get reference element type
	typedef typename reference_element_traits<TElem>::reference_element_type
			ref_elem_type;

//	compile-time number of DoFs
	static const size_t numCo = ref_elem_type::num_corners;

//	count number of function on this subset
	size_t numFct = 0;
	for(size_t fct = 0; fct < fctGrp.num_fct(); ++fct)
	{
		if(is_def_in_subset(fctGrp[fct], si))
			numFct++;
	}

//	return number of algebraic indices
	return numFct * numCo;
}

template<typename TElem>
void
P1ConformDoFDistribution::
indices(TElem* elem, LocalIndices& ind, bool bHang) const
{
//	get reference element type
	typedef typename reference_element_traits<TElem>::reference_element_type
			ref_elem_type;

//	compile-time number of DoFs
	static const size_t numCo = ref_elem_type::num_corners;

//	resize the number of functions
	ind.resize_fct(num_fct());
	for(size_t fct = 0; fct < num_fct(); ++fct)
		ind.resize_dof(fct, 0);

//	add normal dofs
	for(size_t i = 0; i < numCo; ++i)
	{
	//	get vertex
		VertexBase* vrt = GetVertex(elem, i);

	//	get subset index
		int si = m_pISubsetHandler->get_subset_index(vrt);
		UG_ASSERT(si >= 0, "Invalid subset index " << si);

	//	read algebra index
		const size_t firstindex = first_index(vrt, si);

	//	loop all functions
		for(size_t fct = 0; fct < num_fct(); ++fct)
		{
		//	check if function is defined on the subset
			if(!is_def_in_subset(fct, si)) continue;

		//	compute index
			const size_t index = firstindex + m_vvOffsets[si][fct];

		//	add dof to local indices
			ind.push_back_index(fct, index);
		}
	}

//	If no hanging dofs are required, we're done
	if(!bHang) return;

// 	Handle Hanging DoFs on Natural edges
//	collect all edges
	std::vector<EdgeBase*> vEdges;
	CollectEdgesSorted(vEdges, *(m_pISubsetHandler->get_assigned_grid()), elem);

//	loop all edges
	for(size_t ed = 0; ed < vEdges.size(); ++ed)
	{
	//	only constraining edges are of interest
		ConstrainingEdge* edge = dynamic_cast<ConstrainingEdge*>(vEdges[ed]);
		if(edge == NULL) continue;

	//	loop constraining vertices
		for(size_t i = 0; i != edge->num_constrained_vertices(); ++i)
		{
		//	get vertex
			VertexBase* vrt = edge->constrained_vertex(i);

		//	get subset index
			int si = m_pISubsetHandler->get_subset_index(vrt);
			UG_ASSERT(si >= 0, "Invalid subset index " << si);

		//	read algebra index
			const size_t firstindex = first_index(vrt, si);

		//	loop functions
			for(size_t fct = 0; fct < num_fct(); ++fct)
			{
			//	check that function is defined on subset
				if(!is_def_in_subset(fct, si)) continue;

			//	compute algebra index
				const size_t index = firstindex + m_vvOffsets[si][fct];

			//	increase number of indices
				ind.push_back_index(fct, index);
			}
		}
	}

// 	Handle Hanging DoFs on Natural faces

//	Collect all faces
	std::vector<Face*> vFaces;
	CollectFacesSorted(vFaces, *(m_pISubsetHandler->get_assigned_grid()), elem);

//	loop faces
	for(size_t fa = 0; fa < vFaces.size(); ++fa)
	{
	//	only constraining quads are of interest
		ConstrainingQuadrilateral* quad =
				dynamic_cast<ConstrainingQuadrilateral*>(vFaces[fa]);
		if(quad == NULL) continue;

	//	loop hanging vertices
		for(size_t i = 0; i < quad->num_constrained_vertices(); ++i)
		{
		//	get vertex
			VertexBase* vrt = quad->constrained_vertex(i);

		//	get subset index
			int si = m_pISubsetHandler->get_subset_index(vrt);
			UG_ASSERT(si >= 0, "Invalid subset index " << si);

		//	read algebraic index
			const size_t firstindex = first_index(vrt, si);

		//	loop functions
			for(size_t fct = 0; fct < num_fct(); ++fct)
			{
			//	check that function is defined on subset
				if(!is_def_in_subset(fct, si)) continue;

			//	compute algebra index
				const size_t index = firstindex + m_vvOffsets[si][fct];

			//	increase number of algebraic indices
				ind.push_back_index(fct, index);
			}
		}
	}

//	we're done
	return;
}

///////////// Multi index access /////////////////

template<typename TElem>
size_t
P1ConformDoFDistribution::
num_multi_indices(TElem* elem, size_t fct) const
{
//	get reference element type
	typedef typename reference_element_traits<TElem>::reference_element_type
				reference_element_type;

//	number of multi indices for each function is the number of corners
	return reference_element_type::num_corners;
}

template<typename TElem>
size_t
P1ConformDoFDistribution::
num_inner_multi_indices(TElem* elem, size_t fct) const
{
//	get reference element type
	typedef typename reference_element_traits<TElem>::reference_element_type
				reference_element_type;

//	if elem is a vertex, we have a DoF else no DoF
	if(reference_element_type::REFERENCE_OBJECT_ID == ROID_VERTEX)
		return 1;
	else
		return 0;
}

template<typename TElem>
size_t
P1ConformDoFDistribution::
get_multi_indices(TElem* elem, size_t fct, multi_index_vector_type& ind) const
{
//	get reference element type
	typedef typename reference_element_traits<TElem>::reference_element_type
			ref_elem_type;

//	compile-time number of DoFs
	static const size_t numDoFs = ref_elem_type::num_corners;

//	resize indices
	ind.resize(numDoFs);

//	add indices
	for(size_t i = 0; i < numDoFs; ++i)
	{
	//	get vertex
		VertexBase* vrt = GetVertex(elem, i);

	//	get subset index
		int si = m_pISubsetHandler->get_subset_index(vrt);
		UG_ASSERT(si >= 0, "Invalid subset index " << si);

	//	fill algebra index
		ind[i][0] = first_index(vrt, si) + m_vvOffsets[si][fct];
		ind[i][1] = 0;
	}

//	return number of DoFs
	return numDoFs;
}

template<typename TElem>
size_t
P1ConformDoFDistribution::
get_inner_multi_indices(TElem* elem, size_t fct, multi_index_vector_type& ind) const
{
//	get reference element type
	typedef typename reference_element_traits<TElem>::reference_element_type
				reference_element_type;

//	if elem is a vertex, we have a DoF else no DoF
	if(reference_element_type::REFERENCE_OBJECT_ID == ROID_VERTEX)
	{
	//	get Vertex itself
		VertexBase* vrt = GetVertex(elem, 0);

	//	get subset index
		int si = m_pISubsetHandler->get_subset_index(vrt);
		UG_ASSERT(si >= 0, "Invalid subset index " << si);

	//	resize indices
		ind.resize(1);

	// 	fill algebra index
		ind[0][0] = first_index(vrt, si) + m_vvOffsets[si][fct];
		ind[0][1] = 0;

	//	return number of indices
		return 1;
	}
	else
	{
	//	clear, since no indices given
		ind.clear();

	//	return number of indices
		return 0;
	}
}

///////////// Algebra index access /////////////////

template<typename TElem>
size_t
P1ConformDoFDistribution::
num_algebra_indices(TElem* elem, size_t fct) const
{
//	get reference element type
	typedef typename reference_element_traits<TElem>::reference_element_type
				reference_element_type;

//	number of multi indices for each function is the number of corners
	return reference_element_type::num_corners;
}

template<typename TElem>
size_t
P1ConformDoFDistribution::
num_inner_algebra_indices(TElem* elem, size_t fct) const
{
//	get reference element type
	typedef typename reference_element_traits<TElem>::reference_element_type
				reference_element_type;

//	if elem is a vertex, we have a DoF else no DoF
	if(reference_element_type::REFERENCE_OBJECT_ID == ROID_VERTEX)
		return 1;
	else
		return 0;
}

template<typename TElem>
void
P1ConformDoFDistribution::
get_algebra_indices(TElem* elem, algebra_index_vector_type& ind) const
{
//	get reference element
	typedef typename reference_element_traits<TElem>::reference_element_type
			ref_elem_type;

//	clear indices
	ind.clear();

//	fill vector of algebraic indices
	for(size_t fct = 0; fct < num_fct(); ++fct)
	{
		for(size_t i = 0; i < (size_t)ref_elem_type::num_corners; ++i)
		{
		//	get vertex
			VertexBase* vrt = GetVertex(elem, i);

		//	get subset index
			int si = m_pISubsetHandler->get_subset_index(vrt);
			UG_ASSERT(si >= 0, "Invalid subset index " << si);

		//	\todo: can this happen ???
			if(!is_def_in_subset(fct, si)) continue;

		//	get algebra index
			const size_t index = first_index(vrt, si) + m_vvOffsets[si][fct];

		//	write algebra index
			ind.push_back(index);
		}
	}
}

template<typename TElem>
void
P1ConformDoFDistribution::
get_inner_algebra_indices(TElem* elem, algebra_index_vector_type& ind) const
{
//	clear indices
	ind.clear();

//	get base obj type
	static const uint type = geometry_traits<TElem>::BASE_OBJECT_TYPE_ID;

//	only in case of vertex, we have DoFs
	if(type == VERTEX)
	{
	//	get vertex
		VertexBase* vrt = GetVertex(elem, 0);

	//	get subset index
		const int si = m_pISubsetHandler->get_subset_index(vrt);
		UG_ASSERT(si >= 0, "Invalid subset index " << si);

	//	get first algebra index
		const size_t firstIndex = first_index(vrt, si);

	//	loop functions
		for(size_t fct = 0; fct < num_fct(); ++fct)
		{
		//	\todo: Can this happen ???
			if(!is_def_in_subset(fct, si)) continue;

		//	get algebra index
			const size_t index = firstIndex + m_vvOffsets[si][fct];

		//	write algebra index
			ind.push_back(index);
		}
	}
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// GroupedP1ConformDoFDistribution
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

template<typename TElem>
bool
GroupedP1ConformDoFDistribution::
has_dofs_on() const
{
//	get base obj type
	uint type = geometry_traits<TElem>::BASE_OBJECT_TYPE_ID;

//	only in case of a Vertex, we have a DoF
	if(type == VERTEX) return true;
	else return false;
}

template<typename TElem>
size_t
GroupedP1ConformDoFDistribution::
num_indices(int si, const FunctionGroup& fctGrp) const
{
//	get reference element type
	typedef typename reference_element_traits<TElem>::reference_element_type
			ref_elem_type;

//	compile-time number of DoFs
	static const size_t numCo = ref_elem_type::num_corners;

	for(size_t fct = 0; fct < fctGrp.num_fct(); ++fct)
	{
		if(is_def_in_subset(fctGrp[fct], si))
			return numCo;
	}

	return 0;
}

template<typename TElem>
void
GroupedP1ConformDoFDistribution::
indices(TElem* elem, LocalIndices& ind, bool bHang) const
{
//	get reference element type
	typedef typename reference_element_traits<TElem>::reference_element_type
			ref_elem_type;

//	compile-time number of DoFs
	static const size_t numCo = ref_elem_type::num_corners;

//	resize the number of functions
	ind.resize_fct(num_fct());
	for(size_t fct = 0; fct < num_fct(); ++fct)
		ind.resize_dof(fct, 0);

//	add normal dofs
	for(size_t i = 0; i < numCo; ++i)
	{
	//	get vertex
		VertexBase* vrt = GetVertex(elem, i);

	//	get subset index
		int si = m_pISubsetHandler->get_subset_index(vrt);
		UG_ASSERT(si >= 0, "Invalid subset index " << si);

	//	read algebra index
		const size_t firstindex = alg_index(vrt, si);

	//	loop all functions
		for(size_t fct = 0; fct < num_fct(); ++fct)
		{
		//	check if function is defined on the subset
			if(!is_def_in_subset(fct, si)) continue;

		//	compute index
			const size_t index = firstindex;
			const size_t comp = m_vvOffsets[si][fct];

		//	add dof to local indices
			ind.push_back_multi_index(fct, index, comp);
		}
	}

//	If no hanging dofs are required, we're done
	if(!bHang) return;

// 	Handle Hanging DoFs on Natural edges
//	collect all edges
	std::vector<EdgeBase*> vEdges;
	CollectEdgesSorted(vEdges, *(m_pISubsetHandler->get_assigned_grid()), elem);

//	loop all edges
	for(size_t ed = 0; ed < vEdges.size(); ++ed)
	{
	//	only constraining edges are of interest
		ConstrainingEdge* edge = dynamic_cast<ConstrainingEdge*>(vEdges[ed]);
		if(edge == NULL) continue;

	//	loop constraining vertices
		for(size_t i = 0; i != edge->num_constrained_vertices(); ++i)
		{
		//	get vertex
			VertexBase* vrt = edge->constrained_vertex(i);

		//	get subset index
			int si = m_pISubsetHandler->get_subset_index(vrt);
			UG_ASSERT(si >= 0, "Invalid subset index " << si);

		//	read algebra index
			const size_t firstindex = alg_index(vrt, si);

		//	loop functions
			for(size_t fct = 0; fct < num_fct(); ++fct)
			{
			//	check that function is defined on subset
				if(!is_def_in_subset(fct, si)) continue;

			//	compute index
				const size_t index = firstindex;
				const size_t comp = m_vvOffsets[si][fct];

			//	add dof to local indices
				ind.push_back_multi_index(fct, index, comp);
			}
		}
	}

// 	Handle Hanging DoFs on Natural faces

//	Collect all faces
	std::vector<Face*> vFaces;
	CollectFacesSorted(vFaces, *(m_pISubsetHandler->get_assigned_grid()), elem);

//	loop faces
	for(size_t fa = 0; fa < vFaces.size(); ++fa)
	{
	//	only constraining quads are of interest
		ConstrainingQuadrilateral* quad =
				dynamic_cast<ConstrainingQuadrilateral*>(vFaces[fa]);
		if(quad == NULL) continue;

	//	loop hanging vertices
		for(size_t i = 0; i < quad->num_constrained_vertices(); ++i)
		{
		//	get vertex
			VertexBase* vrt = quad->constrained_vertex(i);

		//	get subset index
			int si = m_pISubsetHandler->get_subset_index(vrt);
			UG_ASSERT(si >= 0, "Invalid subset index " << si);

		//	read algebraic index
			const size_t firstindex = alg_index(vrt, si);

		//	loop functions
			for(size_t fct = 0; fct < num_fct(); ++fct)
			{
			//	check that function is defined on subset
				if(!is_def_in_subset(fct, si)) continue;

			//	compute index
				const size_t index = firstindex;
				const size_t comp = m_vvOffsets[si][fct];

			//	add dof to local indices
				ind.push_back_multi_index(fct, index, comp);
			}
		}
	}

//	we're done
	return;
}

///////////// Multi Index access /////////////////

template<typename TElem>
size_t
GroupedP1ConformDoFDistribution::
num_multi_indices(TElem* elem, size_t fct) const
{
//	get reference element type
	typedef typename reference_element_traits<TElem>::reference_element_type
				reference_element_type;

//	number of multi indices for each function is the number of corners
	return reference_element_type::num_corners;
}

template<typename TElem>
size_t
GroupedP1ConformDoFDistribution::
num_inner_multi_indices(TElem* elem, size_t fct) const
{
//	get reference element type
	typedef typename reference_element_traits<TElem>::reference_element_type
				reference_element_type;

//	if elem is a vertex, we have a DoF else no DoF
	if(reference_element_type::REFERENCE_OBJECT_ID == ROID_VERTEX)
		return 1;
	else
		return 0;
}

template<typename TElem>
size_t
GroupedP1ConformDoFDistribution::
get_multi_indices(TElem* elem, size_t fct, multi_index_vector_type& ind) const
{
//	get reference element type
	typedef typename reference_element_traits<TElem>::reference_element_type
			ref_elem_type;

//	compile-time number of DoFs
	static const size_t numDoFs = ref_elem_type::num_corners;

//	resize indices
	ind.resize(numDoFs);

//	add indices
	for(size_t i = 0; i < numDoFs; ++i)
	{
	//	get vertex
		VertexBase* vrt = GetVertex(elem, i);

	//	get subset index
		int si = m_pISubsetHandler->get_subset_index(vrt);
		UG_ASSERT(si >= 0, "Invalid subset index " << si);

	//	fill algebra index
		ind[i][0] = alg_index(vrt, si);
		ind[i][1] = fct;
	}

//	return number of DoFs
	return numDoFs;
}

template<typename TElem>
size_t
GroupedP1ConformDoFDistribution::
get_inner_multi_indices(TElem* elem, size_t fct, multi_index_vector_type& ind) const
{
//	get reference element type
	typedef typename reference_element_traits<TElem>::reference_element_type
				reference_element_type;

//	if elem is a vertex, we have a DoF else no DoF
	if(reference_element_type::REFERENCE_OBJECT_ID == ROID_VERTEX)
	{
	//	get Vertex itself
		VertexBase* vrt = GetVertex(elem, 0);

	//	get subset index
		int si = m_pISubsetHandler->get_subset_index(vrt);
		UG_ASSERT(si >= 0, "Invalid subset index " << si);

	//	resize indices
		ind.resize(1);

	// 	fill algebra index
		ind[0][0] = alg_index(vrt, si);
		ind[0][1] = fct;

	//	return number of indices
		return 1;
	}
	else
	{
	//	clear, since no indices given
		ind.clear();

	//	return number of indices
		return 0;
	}
}


///////////// Algebra Index access /////////////////

/// number of algebra indices on element for a function (Element + Closure of Element)
template<typename TElem>
size_t
GroupedP1ConformDoFDistribution::
num_algebra_indices(TElem* elem, size_t fct) const
{
//	get reference element type
	typedef typename reference_element_traits<TElem>::reference_element_type
				reference_element_type;

//	number of multi indices for each function is the number of corners
	return reference_element_type::num_corners;
}

/// number of algebras indices on element for a function (only inner part of Element)
template<typename TElem>
size_t
GroupedP1ConformDoFDistribution::
num_inner_algebra_indices(TElem* elem, size_t fct) const
{
//	get reference element type
	typedef typename reference_element_traits<TElem>::reference_element_type
				reference_element_type;

//	if elem is a vertex, we have a DoF else no DoF
	if(reference_element_type::REFERENCE_OBJECT_ID == ROID_VERTEX)
		return 1;
	else
		return 0;
}

template<typename TElem>
void
GroupedP1ConformDoFDistribution::
get_algebra_indices(TElem* elem, algebra_index_vector_type& ind) const
{
//	get reference element
	typedef typename reference_element_traits<TElem>::reference_element_type
			ref_elem_type;

//	clear indices
	ind.clear();

// 	if no functions, return
	const int elem_si = m_pISubsetHandler->get_subset_index(elem);
	if(num_fct(elem_si) == 0) return;

//	fill vector of algebraic indices
	for(size_t i = 0; i < (size_t)ref_elem_type::num_corners; ++i)
	{
	//	get vertex
		VertexBase* vrt = GetVertex(elem, i);;

	//	get subset index
		int si = m_pISubsetHandler->get_subset_index(vrt);
		UG_ASSERT(si >= 0, "Invalid subset index " << si);

	//	get algebra indices
		const size_t index = alg_index(vrt, si);

	//	write algebra index
		ind.push_back(index);
	}
}

template<typename TElem>
void
GroupedP1ConformDoFDistribution::
get_inner_algebra_indices(TElem* elem, algebra_index_vector_type& ind) const
{
//	clear indices
	ind.clear();

//	get base obj type
	static const uint type = geometry_traits<TElem>::BASE_OBJECT_TYPE_ID;

	//	only in case of vertex, we have DoFs
	if(type == VERTEX)
	{
	//	get vertex
		VertexBase* vrt = GetVertex(elem, 0);

	//	get subset index
		int si = m_pISubsetHandler->get_subset_index(vrt);
		UG_ASSERT(si >= 0, "Invalid subset index " << si);

	//	if no functions given in this subset, nothing to do
		if(num_fct(si) == 0) return;

	//	get algebra index
		const size_t index = alg_index(vrt, si);

	//	write algebra index
		ind.push_back(index);
	}
}


} // end namespace ug

#endif /* __H__LIB_DISCRETIZATION__DOF_MANAGER__P1CONFORM_IMPL__ */
