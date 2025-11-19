/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
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
	using Map = std::map<LFEID, LocalShapeFunctionSets<dim, TShape, TGrad> >;
	static Map map;
	return map;
};

template <int dim>
std::map<LFEID, LocalFiniteElementProvider::DimLocalDoFSets<dim> >&
LocalFiniteElementProvider::lds_map()
{
	using Map = std::map<LFEID, DimLocalDoFSets<dim> >;
	static Map map;
	return map;
}

template <int dim, typename TShape, typename TGrad>
void LocalFiniteElementProvider::
register_set(const LFEID& id,
             ConstSmartPtr<LocalShapeFunctionSet<dim, TShape, TGrad> > set)
{
//	get type of map
	using Map = std::map<LFEID, LocalShapeFunctionSets<dim, TShape, TGrad> >;
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
	using Map = std::map<LFEID, DimLocalDoFSets<dim> >;
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
	using Map = std::map<LFEID, LocalShapeFunctionSets<dim, TShape, TGrad> >;
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
		return nullptr;
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
	using Map = std::map<LFEID, DimLocalDoFSets<dim> >;
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
		return nullptr;
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
#endif