/*
 * local_finite_element_provider_impl.h
 *
 *  Created on: 17.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__LOCAL_FINITE_ELEMENT__LOCAL_FINITE_ELEMENT_PROVIDER_IMPL__
#define __H__UG__LIB_DISC__LOCAL_FINITE_ELEMENT__LOCAL_FINITE_ELEMENT_PROVIDER_IMPL__

#include "common/common.h"
#include "local_finite_element_provider.h"

namespace ug{

///////////////////////////////////////////
// LocalFiniteElementProvider
///////////////////////////////////////////

template <int dim, typename TShape, typename TGrad>
std::map<LFEID, LocalFiniteElementProvider::LocalShapeFunctionSets<dim, TShape, TGrad> >&
LocalFiniteElementProvider::lsfs_map()
{
	typedef std::map<LFEID, LocalShapeFunctionSets<dim, TShape, TGrad> > Map;
	static Map map;
	return map;
};

template <int dim>
std::map<LFEID, LocalFiniteElementProvider::DimLocalDoFSets<dim> >&
LocalFiniteElementProvider::lds_map()
{
	typedef std::map<LFEID, DimLocalDoFSets<dim> > Map;
	static Map map;
	return map;
}

template <int dim, typename TShape, typename TGrad>
void LocalFiniteElementProvider::
register_set(const LFEID& id,
             ConstSmartPtr<LocalShapeFunctionSet<dim, TShape, TGrad> > set)
{
//	get type of map
	typedef std::map<LFEID, LocalShapeFunctionSets<dim, TShape, TGrad> > Map;
	Map& map = inst().lsfs_map<dim, TShape, TGrad>();
	LocalShapeFunctionSets<dim, TShape, TGrad>& vLSFS = map[id];
	const ReferenceObjectID roid = set->roid();

	if(vLSFS[roid].valid()){
		UG_THROW("LocalFiniteElementProvider::register_set(): "
				"Reference type already registered for trial space: "<<id<<" and "
				" Reference element type "<<roid<<".");
	} else {
		vLSFS[roid] = set;
	}

	if(m_mContSpace.find(id) == m_mContSpace.end()){
		m_mContSpace.insert(std::map<LFEID, bool>::value_type(id, set->continuous()));
	}else{
		if(m_mContSpace[id] != set->continuous())
			UG_THROW("LocalFiniteElementProvider::register_set(): "
					"Reference type says continuous:"<<set->continuous()<<", but "
					" other Reference element says "<<m_mContSpace[id]<<".");
	}

//	register also as DimLocalDoFSet
	register_set(id, set.template cast_dynamic<DimLocalDoFSet<dim> >());
}

template <int dim>
void LocalFiniteElementProvider::
register_set(const LFEID& id,
             ConstSmartPtr<DimLocalDoFSet<dim> > set)
{
//	get type of map
	typedef std::map<LFEID, DimLocalDoFSets<dim> > Map;
	Map& map = inst().lds_map<dim>();
	DimLocalDoFSets<dim>& vLDS = map[id];

	const ReferenceObjectID roid = set->roid();

	if(vLDS[roid].valid()){
		UG_THROW("LocalFiniteElementProvider::register_set(): "
				 "Reference type already registered for trial space: "<<id<<" and "
				 " Reference element type "<<roid<<".");
	} else {
		vLDS[roid] = set;
	}

//	register also as LocalDoFSet
	register_set(id, set.template cast_dynamic<LocalDoFSet>());
}


template <int dim, typename TShape, typename TGrad>
ConstSmartPtr<LocalShapeFunctionSet<dim, TShape, TGrad> >
LocalFiniteElementProvider::
getptr(ReferenceObjectID roid, const LFEID& id, bool bCreate)
{
//	init provider and get map
	typedef std::map<LFEID, LocalShapeFunctionSets<dim, TShape, TGrad> > Map;
	Map& map = inst().lsfs_map<dim, TShape, TGrad>();

//	search for identifier
	typename Map::const_iterator iter = map.find(id);
	if(iter == map.end() || (iter->second)[roid].invalid())
	{
	//	try to create the set
		if(bCreate){
			create_set(roid, id);
			return getptr<dim, TShape, TGrad>(roid, id, false);
		}
		return NULL;
	}

//	return shape function set
	return (iter->second)[roid];
}

template <int dim, typename TShape, typename TGrad>
const LocalShapeFunctionSet<dim, TShape, TGrad>&
LocalFiniteElementProvider::
get(ReferenceObjectID roid, const LFEID& id, bool bCreate)
{
	ConstSmartPtr<LocalShapeFunctionSet<dim, TShape, TGrad> > ptr =
			getptr<dim,TShape,TGrad>(roid, id, bCreate);

	if(ptr.valid()) return *ptr;
	else
		UG_THROW("LocalFiniteElementProvider: Local Shape Function Set not "
				 "found for "<<roid<<" (world dim: "<<dim<<") and type = "<<id<<
				 ". (This is usually due to: a) The function set is not implemented at "
				 " all, or b) The finite element space is discontinuous but the "
				 "evaluation is requested on a subelement, i.e. a grid object "
				 "with dimension less than the dimension where the finite element"
				 " space is defined.)");
}

template <int dim>
ConstSmartPtr<DimLocalDoFSet<dim> >
LocalFiniteElementProvider::
get_dof_ptr(ReferenceObjectID roid, const LFEID& id, bool bCreate)
{
//	init provider and get map
	typedef std::map<LFEID, DimLocalDoFSets<dim> > Map;
	Map& map = inst().lds_map<dim>();

//	search for identifier
	typename Map::const_iterator iter = map.find(id);
	if(iter == map.end() || (iter->second)[roid].invalid())
	{
	//	try to create the set
		if(bCreate){
			create_dof_set(roid, id);
			return get_dof_ptr<dim>(roid, id, false);
		}
		return NULL;
	}

//	return shape function set
	return (iter->second)[roid];
}

template <int dim>
const DimLocalDoFSet<dim>&
LocalFiniteElementProvider::
get_dofs(ReferenceObjectID roid, const LFEID& id, bool bCreate)
{
	ConstSmartPtr<DimLocalDoFSet<dim> > ptr =
			get_dof_ptr<dim>(roid, id, bCreate);

	if(ptr.valid()) return *ptr;
	else
		UG_THROW("LocalFiniteElementProvider: Local DoF Set not "
				 "found for "<<roid<<" (world dim: "<<dim<<") and type = "<<id);
}

}
#endif /* __H__UG__LIB_DISC__LOCAL_FINITE_ELEMENT__LOCAL_FINITE_ELEMENT_PROVIDER_IMPL__ */
