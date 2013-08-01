/*
 * mg_dof_distribution_impl.h
 *
 *  Created on: 08.03.2012
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__DOF_MANAGER__MG_DOF_DISTRIBUTION_IMPL__
#define __H__UG__LIB_DISC__DOF_MANAGER__MG_DOF_DISTRIBUTION_IMPL__

#include <vector>
#include <set>
#include "mg_dof_distribution.h"
#include "lib_grid/tools/periodic_boundary_manager.h"

namespace ug{


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// MGDoFDistribution
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

template <typename TBaseElem>
GeometricObject* MGDoFDistribution::get_parent(TBaseElem* elem) const
{
	return multi_grid()->get_parent(elem);
}

template <typename TBaseElem>
TBaseElem* MGDoFDistribution::parent_if_copy(TBaseElem* elem) const
{
	GeometricObject* pParent = multi_grid()->get_parent(elem);
	TBaseElem* parent = dynamic_cast<TBaseElem*>(pParent);
	if(parent != NULL &&
		multi_grid()->template num_children<TBaseElem>(parent) == 1) return parent;
	else return NULL;
}

template <typename TBaseElem>
TBaseElem* MGDoFDistribution::parent_if_same_type(TBaseElem* elem) const
{
	GeometricObject* pParent = multi_grid()->get_parent(elem);
	return dynamic_cast<TBaseElem*>(pParent);
}

template <typename TBaseElem>
TBaseElem* MGDoFDistribution::child_if_copy(TBaseElem* elem) const
{
	if(multi_grid()->template num_children<TBaseElem>(elem) != 1) return NULL;
	return multi_grid()->template get_child<TBaseElem>(elem, 0);
}


inline size_t& MGDoFDistribution::obj_index(GeometricObject* obj)
{
	switch(obj->base_object_id())
	{
		case VERTEX: return obj_index(static_cast<VertexBase*>(obj));
		case EDGE:   return obj_index(static_cast<EdgeBase*>(obj));
		case FACE:   return obj_index(static_cast<Face*>(obj));
		case VOLUME: return obj_index(static_cast<Volume*>(obj));
		default: UG_THROW("Base Object type not found.");
	}
}

inline const size_t& MGDoFDistribution::obj_index(GeometricObject* obj) const
{
	switch(obj->base_object_id())
	{
		case VERTEX: return obj_index(static_cast<VertexBase*>(obj));
		case EDGE:   return obj_index(static_cast<EdgeBase*>(obj));
		case FACE:   return obj_index(static_cast<Face*>(obj));
		case VOLUME: return obj_index(static_cast<Volume*>(obj));
		default: UG_THROW("Base Object type not found.");
	}
};


template <typename TBaseObject>
bool MGDoFDistribution::
add(TBaseObject* obj, const ReferenceObjectID roid, const int si, LevInfo& li)
{
	UG_ASSERT(si >= 0, "Invalid subset index passed");

//	if no dofs on this subset for the roid, do nothing
	if(num_dofs(roid,si) == 0) return false;

	bool master = false;

	if(m_spMG->has_periodic_boundaries())
	{
		PeriodicBoundaryManager& pbm = *m_spMG->periodic_boundary_manager();
		// ignore slaves
		if(pbm.is_slave(obj))
			return false;

		if(pbm.is_master(obj))
		{
			master = true;
		}
	}

//	compute the number of indices needed on the Geometric object
	size_t numNewIndex = 1;
	if(!m_bGrouped) numNewIndex = num_dofs(roid,si);

// 	set first available index to the object. The first available index is the
//	first managed index plus the size of the index set. (If holes are in the
//	index set, this is not treated here, holes remain)
	obj_index(obj) = li.numIndex;

//	number of managed indices and the number of managed indices on the subset has
//	changed. Thus, increase the counters.
	li.numIndex += numNewIndex;
	li.vNumIndexOnSubset[si] += numNewIndex;

	// if obj is a master, assign all its slaves
	if(master) {
		typedef typename PeriodicBoundaryManager::Group<TBaseObject>::SlaveContainer SlaveContainer;
		typedef typename PeriodicBoundaryManager::Group<TBaseObject>::SlaveIterator SlaveIterator;
		SlaveContainer& slaves = *m_spMG->periodic_boundary_manager()->slaves(obj);
		size_t master_index = obj_index(obj);
		for(SlaveIterator iter = slaves.begin(); iter != slaves.end(); ++iter)
		{
			obj_index(*iter) = master_index;
		}
	}

	return true;
}

} // end namespace ug

#endif /* __H__UG__LIB_DISC__DOF_MANAGER__MG_DOF_DISTRIBUTION_IMPL__ */
