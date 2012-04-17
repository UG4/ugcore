/*
 * local_shape_function_set_impl.h
 *
 *  Created on: 17.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__LOCAL_FINITE_ELEMENT__LOCAL_SHAPE_FUNCTION_SET_IMPL__
#define __H__UG__LIB_DISC__LOCAL_FINITE_ELEMENT__LOCAL_SHAPE_FUNCTION_SET_IMPL__

#include "common/common.h"
#include "local_shape_function_set.h"

namespace ug{

///////////////////////////////////////////
// LocalShapeFunctionSetProvider
///////////////////////////////////////////

template <typename TRefElem>
std::map<LFEID, const LocalShapeFunctionSet<TRefElem>* >&
LocalShapeFunctionSetProvider::get_map()
{
//	get type of map
	typedef std::map<LFEID, const LocalShapeFunctionSet<TRefElem>* > Map;

//	create static map
	static Map sShapeFunctionSetMap;

//	return map
	return sShapeFunctionSetMap;
};

template <int dim>
std::vector<std::map<LFEID, const DimLocalShapeFunctionSet<dim>* > >&
LocalShapeFunctionSetProvider::get_dim_map()
{
//	get type of map
	typedef std::vector<std::map<LFEID, const DimLocalShapeFunctionSet<dim>*> > VecMap;

//	create static map
	static VecMap sShapeFunctionSetMap(NUM_REFERENCE_OBJECTS);

//	return map
	return sShapeFunctionSetMap;
};

template <int dim>
std::vector<DimLocalShapeFunctionSet<dim>*>&
LocalShapeFunctionSetProvider::get_dynamic_allocated_vector()
{
//	get type of map
	typedef std::vector<DimLocalShapeFunctionSet<dim>*> Vec;

//	create static map
	static Vec sShapeFunctionSetMap;

//	return map
	return sShapeFunctionSetMap;
};

template <typename TRefElem>
void
LocalShapeFunctionSetProvider::
clear_maps()
{
//	Reference Object type and dim
	static const ReferenceObjectID roid = TRefElem::REFERENCE_OBJECT_ID;
	static const int dim = TRefElem::dim;

//	get type and map
	typedef std::map<LFEID, const LocalShapeFunctionSet<TRefElem>* > Map;
	Map& map = get_map<TRefElem>();

//	clear
	map.clear();

//	get type and map
	typedef std::vector<std::map<LFEID, const DimLocalShapeFunctionSet<dim>* > > VecMap;
	typedef std::map<LFEID, const DimLocalShapeFunctionSet<dim>* > DimMap;
	VecMap& vDimMap = get_dim_map<dim>();
	DimMap& dimMap = vDimMap[roid];

//	clear
	dimMap.clear();
}

template <typename TRefElem>
void
LocalShapeFunctionSetProvider::
register_set(LFEID type, const LocalShapeFunctionSet<TRefElem>& set)
{
//	Reference Object type and dim
	static const ReferenceObjectID roid = TRefElem::REFERENCE_OBJECT_ID;
	static const int dim = TRefElem::dim;

//	get type of map
	typedef std::map<LFEID, const LocalShapeFunctionSet<TRefElem>* > Map;
	Map& map = get_map<TRefElem>();
	typedef typename Map::value_type MapPair;

//	insert into map
	MapPair pair = MapPair(type, &set);
	if(map.insert(pair).second == false)
		UG_THROW_FATAL("LocalShapeFunctionSetProvider::register_set(): "
				"Reference type already registered for trial space: "<<type<<" and "
				" Reference element type "<<roid<<".");

//	get type of map
	typedef std::vector<std::map<LFEID, const DimLocalShapeFunctionSet<dim>* > > VecMap;
	typedef std::map<LFEID, const DimLocalShapeFunctionSet<dim>* > DimMap;
	VecMap& vDimMap = get_dim_map<dim>();
	DimMap& dimMap = vDimMap[roid];

//	insert into map
	typedef typename DimMap::value_type DimMapPair;
	if(dimMap.insert(DimMapPair(type, &set)).second == false)
		UG_THROW_FATAL("LocalShapeFunctionSetProvider::register_set(): "
				"Reference type already registered for trial space: "<<type<<" and "
				" Reference element type "<<roid<<".");
}


template <typename TRefElem>
bool
LocalShapeFunctionSetProvider::
unregister_set(LFEID id)
{
//	Reference Object type and dim
	static const ReferenceObjectID roid = TRefElem::REFERENCE_OBJECT_ID;
	static const int dim = TRefElem::dim;

//	get type of map
	typedef std::map<LFEID, const LocalShapeFunctionSet<TRefElem>* > Map;

//	init provider and get map
	static Map& map = inst().get_map<TRefElem>();

//	erase element
	if(map.erase(id) != 1) return false;

//	get map
	typedef std::vector<std::map<LFEID, const DimLocalShapeFunctionSet<dim>* > > VecMap;
	typedef std::map<LFEID, const DimLocalShapeFunctionSet<dim>* > DimMap;
	VecMap& vDimMap = get_dim_map<dim>();
	DimMap& dimMap = vDimMap[roid];

//	erase element
	return (dimMap.erase(id) == 1);
}

template <typename TRefElem>
const LocalShapeFunctionSet<TRefElem>&
LocalShapeFunctionSetProvider::
get(LFEID id, bool bCreate)
{
//	get type of map
	typedef std::map<LFEID, const LocalShapeFunctionSet<TRefElem>* > Map;
	const static ReferenceObjectID roid = TRefElem::REFERENCE_OBJECT_ID;

//	init provider and get map
	static Map& map = inst().get_map<TRefElem>();

//	search for identifier
	typename Map::const_iterator iter = map.find(id);

//	if not found
	if(iter == map.end())
	{
		if(bCreate)
		{
		//	try to create the set
			dynamically_create_set(roid, id);

		//	next try to return the set
			return get<TRefElem>(id, false);
		}

		UG_LOG("ERROR in 'LocalShapeFunctionSetProvider::get': "
				"Unknown Trial Space Type "<<id<<" requested for Element"
				" type: "<<roid<<".\n");
		throw(UG_ERROR_LocalShapeFunctionSetNotRegistered(TRefElem::dim, roid, id));
	}

//	return shape function set
	return *(iter->second);
}

template <int dim>
const DimLocalShapeFunctionSet<dim>&
LocalShapeFunctionSetProvider::
get(ReferenceObjectID roid, LFEID id, bool bCreate)
{
//	get type of map
	typedef std::vector<std::map<LFEID, const DimLocalShapeFunctionSet<dim>* > > VecMap;
	typedef std::map<LFEID, const DimLocalShapeFunctionSet<dim>* > Map;


//	init provider and get map
	static VecMap& vMap = inst().get_dim_map<dim>();
	const Map& map = vMap[roid];

//	search for identifier
	typename Map::const_iterator iter = map.find(id);

//	if not found
	if(iter == map.end())
	{
		if(bCreate)
		{
		//	try to create the set
			dynamically_create_set(roid, id);

		//	next try to return the set
			return get<dim>(roid, id, false);
		}

		UG_LOG("ERROR in 'LocalShapeFunctionSetProvider::get': "
				"Unknown Trial Space Type "<<id<<" requested for Element"
				" type: "<<roid<<".\n");
		throw(UG_ERROR_LocalShapeFunctionSetNotRegistered(dim, roid, id));
	}

//	return shape function set
	return *(iter->second);
}

}
#endif /* __H__UG__LIB_DISC__LOCAL_FINITE_ELEMENT__LOCAL_SHAPE_FUNCTION_SET_IMPL__ */
