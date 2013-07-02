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
bool MGDoFDistribution::
add(TBaseObject* obj, const ReferenceObjectID roid, const int si,
    LevInfo<T>& li)
{
	if(si == -1){
		if(m_strictSubsetChecks){
			UG_THROW("Only elements which are assigned to a subset may be added"
					" to the dof manager.");
		}
		else{
			obj_index(obj) = NOT_YET_ASSIGNED;
			return false;
		}
	}

//	if no dofs on this subset for the roid, do nothing
	if(num_dofs(roid,si) == 0) return false;
	if(num_dofs(roid,si) == DoFDistributionInfo::NOT_SPECIFIED) {
		UG_THROW("DoFDistribution: Trying to add dof to element of type "<<roid<<
		         " on Subset "<<si<<", but not all function spaces seem to "
		         "implement the element type.");
	}

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

	return true;
}

template <typename TBaseObject, typename T>
bool MGDoFDistribution::
add(TBaseObject* obj, const ReferenceObjectID roid, const int si,
    LevInfo<T>& li,std::vector<std::pair<size_t,size_t> >& vReplaced)
{
	if(si == -1){
		if(m_strictSubsetChecks){
			UG_THROW("Only elements which are assigned to a subset may be added"
					" to the dof manager.");
		}
		else{
			obj_index(obj) = NOT_YET_ASSIGNED;
			return false;
		}
	}

//	if no dofs on this subset for the roid, do nothing
	if(num_dofs(roid,si) == 0) return false;
	if(num_dofs(roid,si) == DoFDistributionInfo::NOT_SPECIFIED) {
		UG_THROW("DoFDistribution: Trying to add dof to element of type "<<roid<<
		         " on Subset "<<si<<", but not all function spaces seem to "
		         "implement the element type.");
	}

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

//	get old (current) index
	const size_t oldIndex = obj_index(obj);

//	compute the number of indices needed on the Geometric object
	size_t numNewIndex = 1;
	if(!m_bGrouped) numNewIndex = num_dofs(roid,si);

	const size_t newIndex = li.sizeIndexSet;

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

	//	remember replacement
	vReplaced.push_back(std::pair<size_t,size_t>(oldIndex, newIndex));

	return true;
}

template <typename TBaseObject, typename T>
bool MGDoFDistribution::
add_from_free(TBaseObject* obj, const ReferenceObjectID roid, const int si,
              LevInfo<T>& li)
{
	if(si == -1){
		if(m_strictSubsetChecks){
			UG_THROW("Only elements which are assigned to a subset may be added"
					" to the dof manager.");
		}
		else{
			obj_index(obj) = NOT_YET_ASSIGNED;
			return false;
		}
	}

//	if no dofs on this subset for the roid, do nothing
	if(num_dofs(roid,si) == 0) return false;
	if(num_dofs(roid,si) == DoFDistributionInfo::NOT_SPECIFIED) {
		UG_THROW("DoFDistribution: Trying to add dof to element of type "<<roid<<
		         " on Subset "<<si<<", but not all function spaces seem to "
		         "implement the element type.");
	}

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

//	a) 	if no holes are in the index set,
	if(li.free_index_available(numNewIndex))
	{
	//	get a free index (a hole) and use it
		obj_index(obj) = li.pop_free_index(numNewIndex);
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
	return true;
}

template <typename TBaseObject, typename T>
bool MGDoFDistribution::
defragment(TBaseObject* obj, const ReferenceObjectID roid, const int si,
           LevInfo<T>& li, std::vector<std::pair<size_t, size_t> >& vReplaced)
{
	bool master = false;

	//	compute the number of indices needed on the Geometric object
	size_t numNewIndex=1;
	if(!m_bGrouped) numNewIndex= num_dofs(roid,si);

	if(m_spMG->has_periodic_boundaries())
	{
		PeriodicBoundaryManager& pbm = *m_spMG->periodic_boundary_manager();
		// ignore slaves
		if(pbm.is_slave(obj))
		{
			return true;
		}

		if(pbm.is_master(obj))
		{
			master = true;
		}
	}
//	get old (current) index
	const size_t oldIndex = obj_index(obj);

// 	check if index must be replaced by lower one
	if(oldIndex < li.numIndex) return true;

//	if no holes in index set return false
	if (li.free_index_available(numNewIndex)==false) return false;

//	get new index from stack
	while(1)
	{
	//	get a free index (a hole) and use it
		const size_t newIndex = li.pop_free_index(numNewIndex);

	//	check that index is admissible
		if(newIndex < li.numIndex)
		{
		//	set new index
			obj_index(obj) = newIndex;

		//  if obj is a master, also replace all its slave
			if(master)
			{
				typedef typename PeriodicBoundaryManager::Group<TBaseObject>::SlaveContainer SlaveContainer;
				typedef typename PeriodicBoundaryManager::Group<TBaseObject>::SlaveIterator SlaveIterator;
				SlaveContainer& slaves = *m_spMG->periodic_boundary_manager()->slaves(obj);
				for(SlaveIterator iter = slaves.begin(); iter != slaves.end(); ++iter)
				{
					obj_index(*iter) = newIndex;
				}
			}

		//	remember replacement
			vReplaced.push_back(std::pair<size_t,size_t>(oldIndex, newIndex));

		//	done
			break;
		}
	//	one free index less
		li.sizeIndexSet -= 1;

	//	else try next
		if(!li.free_index_available())
			UG_THROW("No more free index, but still need to defragment.");
	}

//	number of Indices stays the same, but size of index set is changed.
	li.sizeIndexSet -= numNewIndex;
	return true;
}


template <typename TBaseObject, typename T>
void MGDoFDistribution::
erase(TBaseObject* obj, const ReferenceObjectID roid, const int si,
      LevInfo<T>& li)
{
//	if no indices needed, we do nothing
	if(num_dofs(roid,si) == 0) return;

	// in periodic case simply ignore slaves
	if(m_spMG->has_periodic_boundaries())
		if(m_spMG->periodic_boundary_manager()->is_slave(obj))
			return;

//	compute number of indices on the geometric object
	size_t numNewIndex = 1;
	if(!m_bGrouped) numNewIndex = num_dofs(roid,si);

//	store the index of the object, that will be erased as a available hole of the
//	index set
	bool bNonContained = li.push_free_index(obj_index(obj), numNewIndex);

	obj_index(obj) = NOT_YET_ASSIGNED;

	if(!bNonContained) return;

//	number of managed indices has changed, thus decrease counter. Note, that the
//	size of the index set remains unchanged.
	li.numIndex -= numNewIndex;
	li.vNumIndexOnSubset[si] -= numNewIndex;
}

template <class TElemNew, class TElemOld>
inline
void MGDoFDistribution::
copy(TElemNew* objNew, TElemOld* objOld)
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
bool MGDoFDistribution::add(GeometricObject* elem, const ReferenceObjectID roid,
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
	return false;
}

template <typename T>
bool MGDoFDistribution::add(GeometricObject* elem, const ReferenceObjectID roid,
                            const int si, LevInfo<T>& li,std::vector<std::pair<size_t,size_t> >& vReplaced)
{
	switch(elem->base_object_id())
	{
		case VERTEX: return add(static_cast<VertexBase*>(elem), roid, si, li,vReplaced);
		case EDGE: return add(static_cast<EdgeBase*>(elem), roid, si, li,vReplaced);
		case FACE: return add(static_cast<Face*>(elem), roid, si, li,vReplaced);
		case VOLUME: return add(static_cast<Volume*>(elem), roid, si, li,vReplaced);
		default: UG_THROW("Geometric Base element not found.");
	}
	return false;
}

template <typename T>
bool MGDoFDistribution::add_from_free(GeometricObject* elem, const ReferenceObjectID roid,
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
	return false;
}

template <typename T>
void MGDoFDistribution::erase(GeometricObject* elem, const ReferenceObjectID roid,
                              const int si, LevInfo<T>& li)
{
	switch(elem->base_object_id())
	{
		case VERTEX: return erase(static_cast<VertexBase*>(elem), roid, si, li);
		case EDGE: 	 return erase(static_cast<EdgeBase*>(elem), roid, si, li);
		case FACE:	 return erase(static_cast<Face*>(elem), roid, si, li);
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
