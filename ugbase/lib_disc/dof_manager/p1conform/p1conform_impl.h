/*
 * p1conform_impl.h
 *
 *  Created on: 13.07.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__DOF_MANAGER__P1CONFORM_IMPL__
#define __H__UG__LIB_DISC__DOF_MANAGER__P1CONFORM_IMPL__

#include <vector>

#include "./p1conform.h"

namespace ug{

///////////////////////////////////////
// P1ConformDoFDistribution
///////////////////////////////////////

template<typename TElem>
void
P1DoFDistribution::
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
		const size_t firstindex = first_index(vrt);

	//	loop all functions
		for(size_t fct = 0; fct < num_fct(); ++fct)
		{
		//	check if function is defined on the subset
			if(!is_def_in_subset(fct, si)) continue;

			if(!m_bGrouped)
			{
			//	compute index
				const size_t index = firstindex + m_vvOffsets[si][fct];

			//	add dof to local indices
				ind.push_back_index(fct, index);
			}
			else
			{
			//	compute index
				const size_t index = firstindex;
				const size_t comp = m_vvOffsets[si][fct];

			//	add dof to local indices
				ind.push_back_multi_index(fct, index, comp);
			}
		}
	}

//	If no hanging dofs are required, we're done
	if(!bHang) return;

// 	Handle Hanging DoFs on Natural edges
//	collect all edges
	std::vector<EdgeBase*> vEdges;
	CollectEdgesSorted(vEdges, *(m_pISubsetHandler->grid()), elem);

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
			const size_t firstindex = first_index(vrt);

		//	loop functions
			for(size_t fct = 0; fct < num_fct(); ++fct)
			{
			//	check that function is defined on subset
				if(!is_def_in_subset(fct, si)) continue;

				if(!m_bGrouped)
				{
				//	compute index
					const size_t index = firstindex + m_vvOffsets[si][fct];

				//	add dof to local indices
					ind.push_back_index(fct, index);
				}
				else
				{
				//	compute index
					const size_t index = firstindex;
					const size_t comp = m_vvOffsets[si][fct];

				//	add dof to local indices
					ind.push_back_multi_index(fct, index, comp);
				}
			}
		}
	}

// 	Handle Hanging DoFs on Natural faces

//	Collect all faces
	std::vector<Face*> vFaces;
	CollectFacesSorted(vFaces, *(m_pISubsetHandler->grid()), elem);

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
			const size_t firstindex = first_index(vrt);

		//	loop functions
			for(size_t fct = 0; fct < num_fct(); ++fct)
			{
			//	check that function is defined on subset
				if(!is_def_in_subset(fct, si)) continue;

				if(!m_bGrouped)
				{
				//	compute index
					const size_t index = firstindex + m_vvOffsets[si][fct];

				//	add dof to local indices
					ind.push_back_index(fct, index);
				}
				else
				{
				//	compute index
					const size_t index = firstindex;
					const size_t comp = m_vvOffsets[si][fct];

				//	add dof to local indices
					ind.push_back_multi_index(fct, index, comp);
				}
			}
		}
	}

//	we're done
	return;
}

template<typename TElem>
size_t
P1DoFDistribution::
multi_indices(TElem* elem, size_t fct, multi_index_vector_type& ind) const
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

	//	get first index
		const size_t firstindex = first_index(vrt);

	//	fill algebra index
		if(!m_bGrouped)
		{
			ind[i][0] = firstindex + m_vvOffsets[si][fct];
			ind[i][1] = 0;
		}
		else
		{
			ind[i][0] = firstindex;
			ind[i][1] = m_vvOffsets[si][fct];
		}
	}

//	return number of DoFs
	return numDoFs;
}

template<typename TElem>
size_t
P1DoFDistribution::
inner_multi_indices(TElem* elem, size_t fct, multi_index_vector_type& ind) const
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

	//	get first index
		const size_t firstindex = first_index(vrt);

	//	fill algebra index
		if(!m_bGrouped)
		{
			ind[0][0] = firstindex + m_vvOffsets[si][fct];
			ind[0][1] = 0;
		}
		else
		{
			ind[0][0] = firstindex;
			ind[0][1] = m_vvOffsets[si][fct];
		}

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

template<typename TElem>
size_t
P1DoFDistribution::
algebra_indices(TElem* elem, algebra_index_vector_type& ind) const
{
//	get reference element
	typedef typename reference_element_traits<TElem>::reference_element_type
			ref_elem_type;

//	clear indices
	ind.clear();

//	fill vector of algebraic indices
	for(size_t i = 0; i < (size_t)ref_elem_type::num_corners; ++i)
	{
	//	get vertex
		VertexBase* vrt = GetVertex(elem, i);

	//	get subset index
		int si = m_pISubsetHandler->get_subset_index(vrt);
		UG_ASSERT(si >= 0, "Invalid subset index " << si);

	//	get first index
		const size_t firstindex = first_index(vrt);

		if(!m_bGrouped)
		{
			for(size_t fct = 0; fct < num_fct(); ++fct)
			{
			//	\todo: can this happen ???
				if(!is_def_in_subset(fct, si)) continue;

				const size_t index = firstindex + m_vvOffsets[si][fct];
				ind.push_back(index);
			}
		}
		else
		{
			if(num_fct(si) > 0)
				ind.push_back(firstindex);
		}
	}

//	return number of indices
	return ind.size();
}

template<typename TElem>
size_t
P1DoFDistribution::
inner_algebra_indices(TElem* elem, algebra_index_vector_type& ind) const
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
		const size_t firstIndex = first_index(vrt);

		if(!m_bGrouped)
		{
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
		else
		{
			if(num_fct(si) > 0)
				ind.push_back(firstIndex);
		}
	}

//	return number of indices
	return ind.size();
}

} // end namespace ug

#endif /* __H__UG__LIB_DISC__DOF_MANAGER__P1CONFORM_IMPL__ */
