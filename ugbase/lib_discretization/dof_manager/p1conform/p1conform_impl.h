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
void
P1ConformDoFDistribution::
update_indices(TElem* elem, LocalIndices& ind, bool withHanging) const
{
//	get reference element type
	typedef typename reference_element_traits<TElem>::reference_element_type
			ref_elem_type;

//	compile-time number of DoFs
	static const size_t numCo = ref_elem_type::num_corners;

//	CASE: no Hanging DoFs
	if(!withHanging)
	{
	//	update algebraic indices
		for(size_t i = 0; i < numCo; ++i)
		{
		//	get vertex
			VertexBase* vrt = GetVertex(elem, i);

		//	get subset index
			int si = m_pISubsetHandler->get_subset_index(vrt);

		//	read algebra index
			const size_t index = first_index(vrt, si);

		//	loop all functions
			size_t numFct = 0;
			for(size_t fct = 0; fct < ind.num_fct(); ++fct)
			{
			//	\todo: can this happen ???
				if(!is_def_in_subset(ind.unique_id(fct), si)) continue;

			//	set algebraic index
				ind.set_index(i + numFct * numCo,
				              index + m_vvOffsets[si][ind.unique_id(fct)]);

			//	increase number of functions
				numFct++;
			}
		}
	}

//	CASE: Hanging DoFs
	else
	{
	//	clear indices
		ind.clear();

	//	reset alg dof counter
		size_t algDof = 0;

	// 	handle Natural DoFs
		for(size_t i = 0; i < numCo; ++i)
		{
		//	get natural vertex
			VertexBase* vrt = GetVertex(elem, i);

		//	get subset index
			int si = m_pISubsetHandler->get_subset_index(vrt);

		//	read algebra index
			const size_t index = first_index(vrt, si);

		//	loop functions
			for(size_t fct = 0; fct < ind.num_fct(); ++fct)
			{
			//	\todo: can this happen ???
				if(!is_def_in_subset(ind.unique_id(fct), si)) continue;

			//	compute algebra index
				const size_t theIndex = index + m_vvOffsets[si][ind.unique_id(fct)];
				UG_ASSERT(theIndex < m_numDoFs, "Adding index " << theIndex <<
				          ", but only " << m_numDoFs << " present.");

			//	increase number of indices
				ind.set_num_indices(algDof+1);

			//	set algebra index
				ind.set_index(algDof, theIndex);

			//	set mapping (fct, dof) -> (alg index, alg comp)
				LocalIndices::multi_index_type dof_ind;
				dof_ind[0] = algDof;
				dof_ind[1] = 0;
				ind.add_dof(fct, dof_ind);

			//	increase number of algebraic indices
				algDof++;
			}
		}

	// 	Handle Hanging DoFs on Natural edges
		{
		//	collect all edges
			std::vector<EdgeBase*> vEdges;
			EdgeBase* ed = dynamic_cast<EdgeBase*>(elem);
			if(ed != NULL)
				CollectEdgesSorted(vEdges,
				                   *(m_pISubsetHandler->get_assigned_grid()), ed);

			Face* face = dynamic_cast<Face*>(elem);
			if(face != NULL)
				CollectEdgesSorted(vEdges,
				                   *(m_pISubsetHandler->get_assigned_grid()), face);

			Volume* vol = dynamic_cast<Volume*>(elem);
			if(vol != NULL)
				CollectEdgesSorted(vEdges,
				                   *(m_pISubsetHandler->get_assigned_grid()), vol);

		//	loop all edges
			for(size_t i = 0; i < vEdges.size(); ++i)
			{
			//	only constraining edges are of interest
				ConstrainingEdge* edge = dynamic_cast<ConstrainingEdge*>(vEdges[i]);
				if(edge == NULL) continue;

			//	loop constraining vertices
				for(VertexBaseIterator iter = edge->constrained_vertices_begin();
						iter != edge->constrained_vertices_end(); ++iter)
				{
				//	get vertex
					VertexBase* vrt = *iter;

				//	get subset index
					int si = m_pISubsetHandler->get_subset_index(vrt);

				//	read algebra index
					const size_t index = first_index(vrt, si);

				//	loop functions
					for(size_t fct = 0; fct < ind.num_fct(); ++fct)
					{
					//	\todo: Can this happen ???
						if(!is_def_in_subset(ind.unique_id(fct), si)) continue;

					//	compute algebra index
						const size_t theIndex = index + m_vvOffsets[si][ind.unique_id(fct)];
						UG_ASSERT(theIndex < m_numDoFs, "Adding index " << theIndex <<
						          ", but only " << m_numDoFs << " present.");

					//	increase number of indices
						ind.set_num_indices(algDof+1);

					//	write algebra index
						ind.set_index(algDof, theIndex);


					//	set mapping (fct, dof) -> (alg index, alg comp)
						LocalIndices::multi_index_type dof_ind;
						dof_ind[0] = algDof;
						dof_ind[1] = 0;
						ind.add_dof(fct, dof_ind);

					//	increase number of algebraic indices
						algDof++;
					}
				}
			}
		}

	// 	Handle Hanging DoFs on Natural edges
		{
		//	Collect all faces
			std::vector<Face*> vFaces; vFaces.clear();
			Face* face = dynamic_cast<Face*>(elem);
			if(face != NULL)
				CollectFacesSorted(vFaces,
				                   *(m_pISubsetHandler->get_assigned_grid()), face);
			Volume* vol = dynamic_cast<Volume*>(elem);
			if(vol != NULL)
				CollectFacesSorted(vFaces,
				                   *(m_pISubsetHandler->get_assigned_grid()), vol);

		//	loop faces
			for(size_t i = 0; i < vFaces.size(); ++i)
			{
			//	only constraining quads are of interest
				ConstrainingQuadrilateral* quad =
						dynamic_cast<ConstrainingQuadrilateral*>(vFaces[i]);
				if(quad == NULL) continue;

			//	loop hanging vertices
				for(VertexBaseIterator iter = quad->constrained_vertices_begin();
						iter != quad->constrained_vertices_end(); ++iter)
				{
				//	get vertex
					VertexBase* vrt = *iter;

				//	get subset index
					int si = m_pISubsetHandler->get_subset_index(vrt);

				//	read algebraic index
					const size_t index = first_index(vrt, si);

				//	loop functions
					for(size_t fct = 0; fct < ind.num_fct(); ++fct)
					{
					// \todo: can this happen ???
						if(!is_def_in_subset(ind.unique_id(fct), si)) continue;

					//	compute algebra index
						const size_t theIndex = index + m_vvOffsets[si][ind.unique_id(fct)];
						UG_ASSERT(theIndex < m_numDoFs, "Adding index " << theIndex <<
						          ", but only " << m_numDoFs << " present.");

					//	increase number of algebraic indices
						ind.set_num_indices(algDof+1);

					//	set algebraic index
						ind.set_index(algDof, theIndex);


					//	set mapping (fct, dof) -> (alg index, alg comp)
						LocalIndices::multi_index_type dof_ind;
						dof_ind[0] = algDof;
						dof_ind[1] = 0;
						ind.add_dof(fct, dof_ind);

					//	increase number of algebraic DoFs
						algDof++;
					}
				}
			}
		}
	}

//	we're done
	return;
}

template<typename TElem>
void
P1ConformDoFDistribution::
update_inner_indices(TElem* elem, LocalIndices& ind) const
{
//	get reference element type
	typedef typename reference_element_traits<TElem>::reference_element_type
				reference_element_type;

//	only in case of a Vertex, we have a DoF
	if(reference_element_type::REFERENCE_OBJECT_ID == ROID_VERTEX)
		return update_indices(elem, ind);
	else
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

//	get reference element type
	typedef typename reference_element_traits<TElem>::reference_element_type
				reference_element_type;

//	only in case of vertex, we have DoFs
	if(reference_element_type::REFERENCE_OBJECT_ID == ROID_VERTEX)
	{
	//	get vertex
		VertexBase* vrt = GetVertex(elem, 0);

	//	get subset index
		int si = m_pISubsetHandler->get_subset_index(vrt);

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


///////////// LocalIndex access /////////////////

template<typename TElem>
void
GroupedP1ConformDoFDistribution::
update_indices(TElem* elem, LocalIndices& ind, bool withHanging) const
{
//	get reference element type
	typedef typename reference_element_traits<TElem>::reference_element_type
			ref_elem_type;

//	compile-time number of DoFs
	static const size_t numDoFs = ref_elem_type::num_corners;

	if(withHanging) throw(UGError("Not implemented"));

//	update all algebraic indices
	for(size_t i = 0; i < numDoFs; ++i)
	{
	//	get vertex
		VertexBase* vrt = GetVertex(elem, i);

	//	get subset index
		int si = m_pISubsetHandler->get_subset_index(vrt);

	//	read algebra index
		const size_t index = alg_index(vrt, si);

	//	write algebra index
		ind.set_index(i, index);
	}
}

template<typename TElem>
void
GroupedP1ConformDoFDistribution::
update_inner_indices(TElem* elem, LocalIndices& ind) const
{
//	get reference element type
	typedef typename reference_element_traits<TElem>::reference_element_type
				reference_element_type;

//	only in case of a Vertex, we have a DoF
	if(reference_element_type::REFERENCE_OBJECT_ID == ROID_VERTEX)
		return update_indices(elem, ind);
	else
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

//	get reference element type
	typedef typename reference_element_traits<TElem>::reference_element_type
				reference_element_type;

//	only in case of vertex, we have DoFs
	if(reference_element_type::REFERENCE_OBJECT_ID == ROID_VERTEX)
	{
	//	get vertex
		VertexBase* vrt = GetVertex(elem, 0);

	//	get subset index
		int si = m_pISubsetHandler->get_subset_index(vrt);

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
