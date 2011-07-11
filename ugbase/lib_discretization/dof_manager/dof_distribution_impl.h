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
	uint type = obj->base_object_type_id();
	switch(type)
	{
		case VERTEX:grid_obj_added(reinterpret_cast<VertexBase*>(obj)); return;
		case EDGE: 	grid_obj_added(reinterpret_cast<EdgeBase*>(obj)); return;
		case FACE:	grid_obj_added(reinterpret_cast<Face*>(obj)); return;
		case VOLUME:grid_obj_added(reinterpret_cast<Volume*>(obj)); return;
	}
	throw(UGFatalError("GeomObject type not known."));
}

template <typename TImpl>
void IDoFDistribution<TImpl>::grid_obj_to_be_removed(GeometricObject* obj)
{
	uint type = obj->base_object_type_id();
	switch(type)
	{
		case VERTEX:grid_obj_to_be_removed(reinterpret_cast<VertexBase*>(obj)); return;
		case EDGE: 	grid_obj_to_be_removed(reinterpret_cast<EdgeBase*>(obj)); return;
		case FACE:	grid_obj_to_be_removed(reinterpret_cast<Face*>(obj)); return;
		case VOLUME:grid_obj_to_be_removed(reinterpret_cast<Volume*>(obj)); return;
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

} // end namespace ug

#endif /* __H__LIB_DISCRETIZATION__DOF_MANAGER__DOF_DISTRIBUTION_IMPL__ */
