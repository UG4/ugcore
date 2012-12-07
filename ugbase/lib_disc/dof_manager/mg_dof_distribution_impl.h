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



template <typename TBaseObject, typename T>
void MGDoFDistribution::
add(TBaseObject* obj, const ReferenceObjectID roid, const int si,
    LevInfo<T>& li)
{
//	if no dofs on this subset for the roid, do nothing
	if(m_vvNumDoFsOnROID[roid][si] == 0) return;

	bool master = false;

	if(m_spMG->has_periodic_boundaries())
	{
		PeriodicBoundaryManager& pbm = *m_spMG->periodic_boundary_manager();
		// ignore slaves
		if(pbm.is_slave(obj))
			return;

		if(pbm.is_master(obj))
		{
			master = true;
		}
	}
//	compute the number of indices needed on the Geometric object
	size_t numNewIndex = 1;
	if(!m_bGrouped) numNewIndex = m_vvNumDoFsOnROID[roid][si];

// 	set first available index to the object. The first available index is the
//	first managed index plus the size of the index set. (If holes are in the
//	index set, this is not treated here, holes remain)
	obj_index(obj) = li.sizeIndexSet;
	
//	the size of the index set has changed. adjust counter
	li.sizeIndexSet += numNewIndex;

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
}

template <typename TBaseObject, typename T>
void MGDoFDistribution::
add_from_free(TBaseObject* obj, const ReferenceObjectID roid, const int si,
              LevInfo<T>& li)
{
//	if no dofs on this subset for the roid, do nothing
	if(m_vvNumDoFsOnROID[roid][si] == 0) return;

	bool master = false;

	if(m_spMG->has_periodic_boundaries())
	{
		PeriodicBoundaryManager& pbm = *m_spMG->periodic_boundary_manager();
		// ignore slaves
		if(pbm.is_slave(obj))
			return;

		if(pbm.is_master(obj))
		{
			master = true;
		}
	}

//	compute the number of indices needed on the Geometric object
	size_t numNewIndex = 1;
	if(!m_bGrouped) numNewIndex = m_vvNumDoFsOnROID[roid][si];

//	a) 	if no holes are in the index set,
	if(li.free_index_available())
	{
	//	get a free index (a hole) and use it
		obj_index(obj) = li.pop_free_index();
	}
	else
	{
	// 	set first available index to the object. The first available index is the
	//	first managed index plus the size of the index set. (If holes are in the
	//	index set, this is not treated here, wholes remain)
		obj_index(obj) = li.sizeIndexSet;

	//	the size of the index set has changed. adjust counter
		li.sizeIndexSet += numNewIndex;
	}

//	number of managed indices and the number of managed indices on the subset has
//	changed. Thus, increase the counters.
	li.numIndex += numNewIndex;
	li.vNumIndexOnSubset[si] += numNewIndex;

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
}

template <typename TBaseObject, typename T>
void MGDoFDistribution::
defragment(TBaseObject* obj, const ReferenceObjectID roid, const int si,
           LevInfo<T>& li, std::vector<std::pair<size_t, size_t> >& vReplaced)
{
	// todo handle periodic case
//	get old (current) index
	const size_t oldIndex = obj_index(obj);

// 	check if index must be replaced by lower one
	if(oldIndex < li.numIndex) return;

//	must have holes in the index set
	UG_ASSERT(li.free_index_available(), "Hole in index set, but no free index:"
	          <<" oldIndex: "<<oldIndex<<", li.numIndex: "<<li.numIndex<<", "<<
	          " roid: "<<roid<<", si: "<<si);

//	get new index from stack
	while(1)
	{
	//	get a free index (a hole) and use it
		const size_t newIndex = li.pop_free_index();

	//	check that index is admissible
		if(newIndex < li.numIndex)
		{
		//	set new index
			obj_index(obj) = newIndex;

		//	remember replacement
			vReplaced.push_back(std::pair<size_t,size_t>(oldIndex,newIndex));

		//	done
			break;
		}
	//	one free index less
		li.sizeIndexSet -= 1;

	//	else try next
		if(!li.free_index_available())
			UG_THROW("No more free index, but still need to defragment.");
	}

//	compute the number of indices needed on the Geometric object
	size_t numNewIndex = 1;
	if(!m_bGrouped) numNewIndex = m_vvNumDoFsOnROID[roid][si];

//	number of Indices stays the same, but size of index set is changed.
	li.sizeIndexSet -= numNewIndex;
}


template <typename TBaseObject, typename T>
void MGDoFDistribution::
erase(TBaseObject* obj, const ReferenceObjectID roid, const int si,
      LevInfo<T>& li)
{
	// todo handle periodic case
//	if no indices needed, we do nothing
	if(m_vvNumDoFsOnROID[roid][si] == 0) return;

//	store the index of the object, that will be erased as a available hole of the
//	index set
	bool bNonContained = li.push_free_index(obj_index(obj));
	if(!bNonContained) return;

//	compute number of indices on the geometric object
	size_t numNewIndex = 1;
	if(!m_bGrouped) numNewIndex = m_vvNumDoFsOnROID[roid][si];

//	number of managed indices has changed, thus decrease counter. Note, that the
//	size of the index set remains unchanged.
	li.numIndex -= numNewIndex;
	li.vNumIndexOnSubset[si] -= numNewIndex;
}

inline
void MGDoFDistribution::
copy(GeometricObject* objNew, GeometricObject* objOld)
{
//	check subsets
	UG_ASSERT(m_spMGSH->get_subset_index(objNew) ==
			  m_spMGSH->get_subset_index(objOld),
			  "Subset index "<<m_spMGSH->get_subset_index(objNew)<<
			  "of replacing obj must match the one of replaced obj "
			  <<m_spMGSH->get_subset_index(objOld));

//	simply copy the index
	obj_index(objNew) = obj_index(objOld);
}

template <typename T>
void MGDoFDistribution::add(GeometricObject* elem, const ReferenceObjectID roid,
                            const int si, LevInfo<T>& li)
{
	switch(elem->base_object_id())
	{
		case VERTEX: return add(static_cast<VertexBase*>(elem), roid, si, li);
		case EDGE: return add(static_cast<EdgeBase*>(elem), roid, si, li);
		case FACE: return add(static_cast<Face*>(elem), roid, si, li);
		case VOLUME: return add(static_cast<Volume*>(elem), roid, si, li);
		default: UG_THROW("Geometric Base element not found.");
	}
}

template <typename T>
void MGDoFDistribution::add_from_free(GeometricObject* elem, const ReferenceObjectID roid,
                                      const int si, LevInfo<T>& li)
{
	switch(elem->base_object_id())
	{
		case VERTEX: return add_from_free(static_cast<VertexBase*>(elem), roid, si, li);
		case EDGE: return add_from_free(static_cast<EdgeBase*>(elem), roid, si, li);
		case FACE: return add_from_free(static_cast<Face*>(elem), roid, si, li);
		case VOLUME: return add_from_free(static_cast<Volume*>(elem), roid, si, li);
		default: UG_THROW("Geometric Base element not found.");
	}
}

template <typename T>
void MGDoFDistribution::erase(GeometricObject* elem, const ReferenceObjectID roid,
                              const int si, LevInfo<T>& li)
{
	switch(elem->base_object_id())
	{
		case VERTEX: return erase(static_cast<VertexBase*>(elem), roid, si, li);
		case EDGE: return erase(static_cast<EdgeBase*>(elem), roid, si, li);
		case FACE: return erase(static_cast<Face*>(elem), roid, si, li);
		case VOLUME: return erase(static_cast<Volume*>(elem), roid, si, li);
		default: UG_THROW("Geometric Base element not found.");
	}
}

template <typename T>
void MGDoFDistribution::defragment(GeometricObject* elem, const ReferenceObjectID roid, const int si,
                                   LevInfo<T>& li, std::vector<std::pair<size_t, size_t> >& vReplaced)
{
	switch(elem->base_object_id())
	{
		case VERTEX: return defragment(static_cast<VertexBase*>(elem), roid, si, li, vReplaced);
		case EDGE: return defragment(static_cast<EdgeBase*>(elem), roid, si, li, vReplaced);
		case FACE: return defragment(static_cast<Face*>(elem), roid, si, li, vReplaced);
		case VOLUME: return defragment(static_cast<Volume*>(elem), roid, si, li, vReplaced);
		default: UG_THROW("Geometric Base element not found.");
	}
}


} // end namespace ug

#endif /* __H__UG__LIB_DISC__DOF_MANAGER__MG_DOF_DISTRIBUTION_IMPL__ */
