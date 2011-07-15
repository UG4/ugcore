/*
 * dof_distribution_impl.h
 *
 *  Created on: 13.07.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__DOF_MANAGER__DOF_DISTRIBUTION_IMPL__
#define __H__LIB_DISCRETIZATION__DOF_MANAGER__DOF_DISTRIBUTION_IMPL__

#include "dof_distribution.h"

namespace ug{

template <typename TImpl>
void IDoFDistribution<TImpl>::manage_grid_function(IGridFunction<TImpl>& gf)
{
//	search for grid function in list
	typename std::vector<IGridFunction<TImpl>*>::iterator it;
	it = find(m_vManagedGridFunc.begin(), m_vManagedGridFunc.end(), &gf);

//	add if not found
	if(it == m_vManagedGridFunc.end())
		m_vManagedGridFunc.push_back(&gf);
}

template <typename TImpl>
void IDoFDistribution<TImpl>::unmanage_grid_function(IGridFunction<TImpl>& gf)
{
//	search for grid function in list
	typename std::vector<IGridFunction<TImpl>*>::iterator it;
	it = find(m_vManagedGridFunc.begin(), m_vManagedGridFunc.end(), &gf);

//	remove if found
	if(it != m_vManagedGridFunc.end())
		m_vManagedGridFunc.erase(it);
}

template <typename TImpl>
bool IDoFDistribution<TImpl>::permute_indices(std::vector<size_t>& vIndNew)
{
//	check size of permuting array. Must have same size as index set
	if(vIndNew.size() != num_dofs())
	{
		UG_LOG("ERROR in 'IDoFDistribution<TImpl>::permute_indices': Passed "
				" permutation does not have the size of the index set "
				<<num_dofs()<<", but has size "<<vIndNew.size()<<"\n");
		return false;
	}

//	swap indices as implemented
	if(!getImpl().permute_indices(vIndNew)) return false;

//	in parallel adjust also the layouts
#ifdef UG_PARALLEL
	PermuteIndicesInIndexLayout(m_slaveLayout, vIndNew);
	PermuteIndicesInIndexLayout(m_masterLayout, vIndNew);
	PermuteIndicesInIndexLayout(m_verticalSlaveLayout, vIndNew);
	PermuteIndicesInIndexLayout(m_verticalMasterLayout, vIndNew);
#endif

//	swap values of handled grid functions
	for(size_t i = 0; i < m_vManagedGridFunc.size(); ++i)
		m_vManagedGridFunc[i]->permute_values(vIndNew);

//	we're done
	return true;
}

template <typename TImpl>
void IDoFDistribution<TImpl>::grid_obj_added(GeometricObject* obj)
{
	switch(obj->base_object_type_id())
	{
		case VERTEX:grid_obj_added(static_cast<VertexBase*>(obj)); return;
		case EDGE: 	grid_obj_added(static_cast<EdgeBase*>(obj)); return;
		case FACE:	grid_obj_added(static_cast<Face*>(obj)); return;
		case VOLUME:grid_obj_added(static_cast<Volume*>(obj)); return;
	}
	throw(UGFatalError("GeomObject type not known."));
}

template <typename TImpl>
void IDoFDistribution<TImpl>::grid_obj_to_be_removed(GeometricObject* obj)
{
	switch(obj->base_object_type_id())
	{
		case VERTEX:grid_obj_to_be_removed(static_cast<VertexBase*>(obj)); return;
		case EDGE: 	grid_obj_to_be_removed(static_cast<EdgeBase*>(obj)); return;
		case FACE:	grid_obj_to_be_removed(static_cast<Face*>(obj)); return;
		case VOLUME:grid_obj_to_be_removed(static_cast<Volume*>(obj)); return;
	}
	throw(UGFatalError("GeomObject type not known."));
}
template <typename TImpl>
bool IDoFDistribution<TImpl>::indices_swaped(const std::vector<std::pair<size_t, size_t> >& vIndexMap,
							bool bDisjunct)
{
//	swap values of handled grid functions
	for(size_t i = 0; i < m_vManagedGridFunc.size(); ++i)
		m_vManagedGridFunc[i]->copy_values(vIndexMap, bDisjunct);

//	we're done
	return true;
}

template <typename TImpl>
void IDoFDistribution<TImpl>::num_indices_changed(size_t newSize)
{
//	swap values of handled grid functions
	for(size_t i = 0; i < m_vManagedGridFunc.size(); ++i)
		m_vManagedGridFunc[i]->resize_values(newSize);
}

////////////////////////////////////////////////////////////////////////////////
//	indices
////////////////////////////////////////////////////////////////////////////////

template <typename TImpl>
void IDoFDistribution<TImpl>::indices(GeometricObject* elem, LocalIndices& ind, bool bHang) const
{
	switch(elem->base_object_type_id())
	{
		case VERTEX: return indices(static_cast<VertexBase*>(elem), ind, bHang);
		case EDGE:   return indices(static_cast<EdgeBase*>(elem), ind, bHang);
		case FACE:   return indices(static_cast<Face*>(elem), ind, bHang);
		case VOLUME: return indices(static_cast<Volume*>(elem), ind, bHang);
	}
	throw(UGFatalError("Base Object type not found."));
}

template <typename TImpl>
void IDoFDistribution<TImpl>::indices(VertexBase* elem, LocalIndices& ind, bool bHang) const
{
	switch(elem->shared_pipe_section())
	{
		case SPSVRT_VERTEX: 		return indices(static_cast<Vertex*>(elem), ind, bHang);
		case SPSVRT_HANGING_VERTEX: return indices(static_cast<HangingVertex*>(elem), ind, bHang);
	}
	throw(UGFatalError("Vertex type not found."));
}

template <typename TImpl>
void IDoFDistribution<TImpl>::indices(EdgeBase* elem, LocalIndices& ind, bool bHang) const
{
	switch(elem->shared_pipe_section())
	{
		case SPSEDGE_EDGE: 			   return indices(static_cast<Edge*>(elem), ind, bHang);
		case SPSEDGE_CONSTRAINED_EDGE: return indices(static_cast<ConstrainedEdge*>(elem), ind, bHang);
		case SPSEDGE_CONSTRAINING_EDGE:return indices(static_cast<ConstrainingEdge*>(elem), ind, bHang);
	}
	throw(UGFatalError("Edge type not found."));
}

template <typename TImpl>
void IDoFDistribution<TImpl>::indices(Face* elem, LocalIndices& ind, bool bHang) const
{
	switch(elem->shared_pipe_section())
	{
		case SPSFACE_TRIANGLE: return indices(static_cast<Triangle*>(elem), ind, bHang);
		case SPSFACE_CONSTRAINED_TRIANGLE: return indices(static_cast<ConstrainedTriangle*>(elem), ind, bHang);
		case SPSFACE_CONSTRAINING_TRIANGLE: return indices(static_cast<ConstrainingTriangle*>(elem), ind, bHang);
		case SPSFACE_QUADRILATERAL: return indices(static_cast<Quadrilateral*>(elem), ind, bHang);
		case SPSFACE_CONSTRAINED_QUADRILATERAL: return indices(static_cast<ConstrainedQuadrilateral*>(elem), ind, bHang);
		case SPSFACE_CONSTRAINING_QUADRILATERAL: return indices(static_cast<ConstrainingQuadrilateral*>(elem), ind, bHang);
	}
	throw(UGFatalError("Face type not found."));
}

template <typename TImpl>
void IDoFDistribution<TImpl>::indices(Volume* elem, LocalIndices& ind, bool bHang) const
{
	switch(elem->shared_pipe_section())
	{
		case SPSVOL_TETRAHEDRON: return indices(static_cast<Tetrahedron*>(elem), ind, bHang);
		case SPSVOL_PYRAMID: return indices(static_cast<Pyramid*>(elem), ind, bHang);
		case SPSVOL_PRISM: return indices(static_cast<Prism*>(elem), ind, bHang);
		case SPSVOL_HEXAHEDRON: return indices(static_cast<Hexahedron*>(elem), ind, bHang);
	}
	throw(UGFatalError("Volume type not found."));
}

////////////////////////////////////////////////////////////////////////////////
//	multi_indices
////////////////////////////////////////////////////////////////////////////////

template <typename TImpl>
size_t IDoFDistribution<TImpl>::multi_indices(GeometricObject* elem, size_t fct, multi_index_vector_type& ind) const
{
	switch(elem->base_object_type_id())
	{
		case VERTEX: return multi_indices(static_cast<VertexBase*>(elem), fct, ind);
		case EDGE:   return multi_indices(static_cast<EdgeBase*>(elem), fct, ind);
		case FACE:   return multi_indices(static_cast<Face*>(elem), fct, ind);
		case VOLUME: return multi_indices(static_cast<Volume*>(elem), fct, ind);
	}
	throw(UGFatalError("Base Object type not found."));
}

template <typename TImpl>
size_t IDoFDistribution<TImpl>::multi_indices(VertexBase* elem, size_t fct, multi_index_vector_type& ind) const
{
	switch(elem->shared_pipe_section())
	{
		case SPSVRT_VERTEX: 		return multi_indices(static_cast<Vertex*>(elem), fct, ind);
		case SPSVRT_HANGING_VERTEX: return multi_indices(static_cast<HangingVertex*>(elem), fct, ind);
	}
	throw(UGFatalError("Vertex type not found."));
}

template <typename TImpl>
size_t IDoFDistribution<TImpl>::multi_indices(EdgeBase* elem, size_t fct, multi_index_vector_type& ind) const
{
	switch(elem->shared_pipe_section())
	{
		case SPSEDGE_EDGE: 			   return multi_indices(static_cast<Edge*>(elem), fct, ind);
		case SPSEDGE_CONSTRAINED_EDGE: return multi_indices(static_cast<ConstrainedEdge*>(elem), fct, ind);
		case SPSEDGE_CONSTRAINING_EDGE:return multi_indices(static_cast<ConstrainingEdge*>(elem), fct, ind);
	}
	throw(UGFatalError("Edge type not found."));
}

template <typename TImpl>
size_t IDoFDistribution<TImpl>::multi_indices(Face* elem, size_t fct, multi_index_vector_type& ind) const
{
	switch(elem->shared_pipe_section())
	{
		case SPSFACE_TRIANGLE: return multi_indices(static_cast<Triangle*>(elem), fct, ind);
		case SPSFACE_CONSTRAINED_TRIANGLE: return multi_indices(static_cast<ConstrainedTriangle*>(elem), fct, ind);
		case SPSFACE_CONSTRAINING_TRIANGLE: return multi_indices(static_cast<ConstrainingTriangle*>(elem), fct, ind);
		case SPSFACE_QUADRILATERAL: return multi_indices(static_cast<Quadrilateral*>(elem), fct, ind);
		case SPSFACE_CONSTRAINED_QUADRILATERAL: return multi_indices(static_cast<ConstrainedQuadrilateral*>(elem), fct, ind);
		case SPSFACE_CONSTRAINING_QUADRILATERAL: return multi_indices(static_cast<ConstrainingQuadrilateral*>(elem), fct, ind);
	}
	throw(UGFatalError("Face type not found."));
}

template <typename TImpl>
size_t IDoFDistribution<TImpl>::multi_indices(Volume* elem, size_t fct, multi_index_vector_type& ind) const
{
	switch(elem->shared_pipe_section())
	{
		case SPSVOL_TETRAHEDRON: return multi_indices(static_cast<Tetrahedron*>(elem), fct, ind);
		case SPSVOL_PYRAMID: return multi_indices(static_cast<Pyramid*>(elem), fct, ind);
		case SPSVOL_PRISM: return multi_indices(static_cast<Prism*>(elem), fct, ind);
		case SPSVOL_HEXAHEDRON: return multi_indices(static_cast<Hexahedron*>(elem), fct, ind);
	}
	throw(UGFatalError("Volume type not found."));
}


////////////////////////////////////////////////////////////////////////////////
//	inner_multi_indices
////////////////////////////////////////////////////////////////////////////////

template <typename TImpl>
size_t IDoFDistribution<TImpl>::inner_multi_indices(GeometricObject* elem, size_t fct, multi_index_vector_type& ind) const
{
	switch(elem->base_object_type_id())
	{
		case VERTEX: return inner_multi_indices(static_cast<VertexBase*>(elem), fct, ind);
		case EDGE:   return inner_multi_indices(static_cast<EdgeBase*>(elem), fct, ind);
		case FACE:   return inner_multi_indices(static_cast<Face*>(elem), fct, ind);
		case VOLUME: return inner_multi_indices(static_cast<Volume*>(elem), fct, ind);
	}
	throw(UGFatalError("Base Object type not found."));
}

template <typename TImpl>
size_t IDoFDistribution<TImpl>::inner_multi_indices(VertexBase* elem, size_t fct, multi_index_vector_type& ind) const
{
	switch(elem->shared_pipe_section())
	{
		case SPSVRT_VERTEX: 		return inner_multi_indices(static_cast<Vertex*>(elem), fct, ind);
		case SPSVRT_HANGING_VERTEX: return inner_multi_indices(static_cast<HangingVertex*>(elem), fct, ind);
	}
	throw(UGFatalError("Vertex type not found."));
}

template <typename TImpl>
size_t IDoFDistribution<TImpl>::inner_multi_indices(EdgeBase* elem, size_t fct, multi_index_vector_type& ind) const
{
	switch(elem->shared_pipe_section())
	{
		case SPSEDGE_EDGE: 			   return inner_multi_indices(static_cast<Edge*>(elem), fct, ind);
		case SPSEDGE_CONSTRAINED_EDGE: return inner_multi_indices(static_cast<ConstrainedEdge*>(elem), fct, ind);
		case SPSEDGE_CONSTRAINING_EDGE:return inner_multi_indices(static_cast<ConstrainingEdge*>(elem), fct, ind);
	}
	throw(UGFatalError("Edge type not found."));
}

template <typename TImpl>
size_t IDoFDistribution<TImpl>::inner_multi_indices(Face* elem, size_t fct, multi_index_vector_type& ind) const
{
	switch(elem->shared_pipe_section())
	{
		case SPSFACE_TRIANGLE: return inner_multi_indices(static_cast<Triangle*>(elem), fct, ind);
		case SPSFACE_CONSTRAINED_TRIANGLE: return inner_multi_indices(static_cast<ConstrainedTriangle*>(elem), fct, ind);
		case SPSFACE_CONSTRAINING_TRIANGLE: return inner_multi_indices(static_cast<ConstrainingTriangle*>(elem), fct, ind);
		case SPSFACE_QUADRILATERAL: return inner_multi_indices(static_cast<Quadrilateral*>(elem), fct, ind);
		case SPSFACE_CONSTRAINED_QUADRILATERAL: return inner_multi_indices(static_cast<ConstrainedQuadrilateral*>(elem), fct, ind);
		case SPSFACE_CONSTRAINING_QUADRILATERAL: return inner_multi_indices(static_cast<ConstrainingQuadrilateral*>(elem), fct, ind);
	}
	throw(UGFatalError("Face type not found."));
}

template <typename TImpl>
size_t IDoFDistribution<TImpl>::inner_multi_indices(Volume* elem, size_t fct, multi_index_vector_type& ind) const
{
	switch(elem->shared_pipe_section())
	{
		case SPSVOL_TETRAHEDRON: return inner_multi_indices(static_cast<Tetrahedron*>(elem), fct, ind);
		case SPSVOL_PYRAMID: return inner_multi_indices(static_cast<Pyramid*>(elem), fct, ind);
		case SPSVOL_PRISM: return inner_multi_indices(static_cast<Prism*>(elem), fct, ind);
		case SPSVOL_HEXAHEDRON: return inner_multi_indices(static_cast<Hexahedron*>(elem), fct, ind);
	}
	throw(UGFatalError("Volume type not found."));
}



////////////////////////////////////////////////////////////////////////////////
//	algebra_indices
////////////////////////////////////////////////////////////////////////////////

template <typename TImpl>
size_t IDoFDistribution<TImpl>::algebra_indices(GeometricObject* elem, algebra_index_vector_type& ind) const
{
	switch(elem->base_object_type_id())
	{
		case VERTEX: return algebra_indices(static_cast<VertexBase*>(elem), ind);
		case EDGE:   return algebra_indices(static_cast<EdgeBase*>(elem), ind);
		case FACE:   return algebra_indices(static_cast<Face*>(elem), ind);
		case VOLUME: return algebra_indices(static_cast<Volume*>(elem), ind);
	}
	throw(UGFatalError("Base Object type not found."));
}

template <typename TImpl>
size_t IDoFDistribution<TImpl>::algebra_indices(VertexBase* elem, algebra_index_vector_type& ind) const
{
	switch(elem->shared_pipe_section())
	{
		case SPSVRT_VERTEX: 		return algebra_indices(static_cast<Vertex*>(elem), ind);
		case SPSVRT_HANGING_VERTEX: return algebra_indices(static_cast<HangingVertex*>(elem), ind);
	}
	throw(UGFatalError("Vertex type not found."));
}

template <typename TImpl>
size_t IDoFDistribution<TImpl>::algebra_indices(EdgeBase* elem, algebra_index_vector_type& ind) const
{
	switch(elem->shared_pipe_section())
	{
		case SPSEDGE_EDGE: 			   return algebra_indices(static_cast<Edge*>(elem), ind);
		case SPSEDGE_CONSTRAINED_EDGE: return algebra_indices(static_cast<ConstrainedEdge*>(elem), ind);
		case SPSEDGE_CONSTRAINING_EDGE:return algebra_indices(static_cast<ConstrainingEdge*>(elem), ind);
	}
	throw(UGFatalError("Edge type not found."));
}

template <typename TImpl>
size_t IDoFDistribution<TImpl>::algebra_indices(Face* elem, algebra_index_vector_type& ind) const
{
	switch(elem->shared_pipe_section())
	{
		case SPSFACE_TRIANGLE: return algebra_indices(static_cast<Triangle*>(elem), ind);
		case SPSFACE_CONSTRAINED_TRIANGLE: return algebra_indices(static_cast<ConstrainedTriangle*>(elem), ind);
		case SPSFACE_CONSTRAINING_TRIANGLE: return algebra_indices(static_cast<ConstrainingTriangle*>(elem), ind);
		case SPSFACE_QUADRILATERAL: return algebra_indices(static_cast<Quadrilateral*>(elem), ind);
		case SPSFACE_CONSTRAINED_QUADRILATERAL: return algebra_indices(static_cast<ConstrainedQuadrilateral*>(elem), ind);
		case SPSFACE_CONSTRAINING_QUADRILATERAL: return algebra_indices(static_cast<ConstrainingQuadrilateral*>(elem), ind);
	}
	throw(UGFatalError("Face type not found."));
}

template <typename TImpl>
size_t IDoFDistribution<TImpl>::algebra_indices(Volume* elem, algebra_index_vector_type& ind) const
{
	switch(elem->shared_pipe_section())
	{
		case SPSVOL_TETRAHEDRON: return algebra_indices(static_cast<Tetrahedron*>(elem), ind);
		case SPSVOL_PYRAMID: return algebra_indices(static_cast<Pyramid*>(elem), ind);
		case SPSVOL_PRISM: return algebra_indices(static_cast<Prism*>(elem), ind);
		case SPSVOL_HEXAHEDRON: return algebra_indices(static_cast<Hexahedron*>(elem), ind);
	}
	throw(UGFatalError("Volume type not found."));
}

////////////////////////////////////////////////////////////////////////////////
//	inner_algebra_indices
////////////////////////////////////////////////////////////////////////////////

template <typename TImpl>
size_t IDoFDistribution<TImpl>::inner_algebra_indices(GeometricObject* elem, algebra_index_vector_type& ind) const
{
	switch(elem->base_object_type_id())
	{
		case VERTEX: return inner_algebra_indices(static_cast<VertexBase*>(elem), ind);
		case EDGE:   return inner_algebra_indices(static_cast<EdgeBase*>(elem), ind);
		case FACE:   return inner_algebra_indices(static_cast<Face*>(elem), ind);
		case VOLUME: return inner_algebra_indices(static_cast<Volume*>(elem), ind);
	}
	throw(UGFatalError("Base Object type not found."));
}

template <typename TImpl>
size_t IDoFDistribution<TImpl>::inner_algebra_indices(VertexBase* elem, algebra_index_vector_type& ind) const
{
	switch(elem->shared_pipe_section())
	{
		case SPSVRT_VERTEX: 		return inner_algebra_indices(static_cast<Vertex*>(elem), ind);
		case SPSVRT_HANGING_VERTEX: return inner_algebra_indices(static_cast<HangingVertex*>(elem), ind);
	}
	throw(UGFatalError("Vertex type not found."));
}

template <typename TImpl>
size_t IDoFDistribution<TImpl>::inner_algebra_indices(EdgeBase* elem, algebra_index_vector_type& ind) const
{
	switch(elem->shared_pipe_section())
	{
		case SPSEDGE_EDGE: 			   return inner_algebra_indices(static_cast<Edge*>(elem), ind);
		case SPSEDGE_CONSTRAINED_EDGE: return inner_algebra_indices(static_cast<ConstrainedEdge*>(elem), ind);
		case SPSEDGE_CONSTRAINING_EDGE:return inner_algebra_indices(static_cast<ConstrainingEdge*>(elem), ind);
	}
	throw(UGFatalError("Edge type not found."));
}

template <typename TImpl>
size_t IDoFDistribution<TImpl>::inner_algebra_indices(Face* elem, algebra_index_vector_type& ind) const
{
	switch(elem->shared_pipe_section())
	{
		case SPSFACE_TRIANGLE: return inner_algebra_indices(static_cast<Triangle*>(elem), ind);
		case SPSFACE_CONSTRAINED_TRIANGLE: return inner_algebra_indices(static_cast<ConstrainedTriangle*>(elem), ind);
		case SPSFACE_CONSTRAINING_TRIANGLE: return inner_algebra_indices(static_cast<ConstrainingTriangle*>(elem), ind);
		case SPSFACE_QUADRILATERAL: return inner_algebra_indices(static_cast<Quadrilateral*>(elem), ind);
		case SPSFACE_CONSTRAINED_QUADRILATERAL: return inner_algebra_indices(static_cast<ConstrainedQuadrilateral*>(elem), ind);
		case SPSFACE_CONSTRAINING_QUADRILATERAL: return inner_algebra_indices(static_cast<ConstrainingQuadrilateral*>(elem), ind);
	}
	throw(UGFatalError("Face type not found."));
}

template <typename TImpl>
size_t IDoFDistribution<TImpl>::inner_algebra_indices(Volume* elem, algebra_index_vector_type& ind) const
{
	switch(elem->shared_pipe_section())
	{
		case SPSVOL_TETRAHEDRON: return inner_algebra_indices(static_cast<Tetrahedron*>(elem), ind);
		case SPSVOL_PYRAMID: return inner_algebra_indices(static_cast<Pyramid*>(elem), ind);
		case SPSVOL_PRISM: return inner_algebra_indices(static_cast<Prism*>(elem), ind);
		case SPSVOL_HEXAHEDRON: return inner_algebra_indices(static_cast<Hexahedron*>(elem), ind);
	}
	throw(UGFatalError("Volume type not found."));
}


} // end namespace ug

#endif /* __H__LIB_DISCRETIZATION__DOF_MANAGER__DOF_DISTRIBUTION_IMPL__ */
