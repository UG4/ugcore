/*
 * local_shape_function_set_impl.h
 *
 *  Created on: 17.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISCRETIZATION__LOCAL_FINITE_ELEMENT__LOCAL_SHAPE_FUNCTION_SET_IMPL__
#define __H__UG__LIB_DISCRETIZATION__LOCAL_FINITE_ELEMENT__LOCAL_SHAPE_FUNCTION_SET_IMPL__

#include "common/common.h"
#include "local_shape_function_set.h"

namespace ug{

///////////////////////////////////////////
// LocalShapeFunctionSetProvider
///////////////////////////////////////////

template <typename TRefElem>
std::map<LFEID, const ReferenceElemLocalShapeFunctionSet<TRefElem>* >&
LocalShapeFunctionSetProvider::get_map()
{
//	get type of map
	typedef std::map<LFEID, const ReferenceElemLocalShapeFunctionSet<TRefElem>* > Map;

//	create static map
	static Map sShapeFunctionSetMap;

//	return map
	return sShapeFunctionSetMap;
};

template <int dim>
std::map<LFEID, const LocalShapeFunctionSet<dim>* >*
LocalShapeFunctionSetProvider::get_dim_map()
{
//	get type of map
	typedef std::map<LFEID, const LocalShapeFunctionSet<dim>* > Map;

//	create static map
	static Map sShapeFunctionSetMap[NUM_REFERENCE_OBJECTS];

//	return map
	return sShapeFunctionSetMap;
};

template <typename TRefElem>
bool
LocalShapeFunctionSetProvider::
register_set(LFEID type, const ReferenceElemLocalShapeFunctionSet<TRefElem>& set)
{
//	Reference Object type
	static const ReferenceObjectID roid = TRefElem::REFERENCE_OBJECT_ID;

//	get type of map
	typedef std::map<LFEID, const ReferenceElemLocalShapeFunctionSet<TRefElem>* > Map;
	static Map& map = get_map<TRefElem>();
	typedef std::pair<LFEID,const ReferenceElemLocalShapeFunctionSet<TRefElem>*> MapPair;

//	insert into map
	if(map.insert(MapPair(type, &set)).second == false)
	{
		UG_LOG("ERROR in 'LocalShapeFunctionSetProvider::register_set()': "
				"Reference type already registered for trial space: "<<type<<" and "
				" Reference element type "<<roid<<".\n");
		return false;
	}

//	get type of map
	typedef std::map<LFEID, const LocalShapeFunctionSet<TRefElem::dim>* > DimMap;
	static DimMap* vDimMap = get_dim_map<TRefElem::dim>();
	DimMap& dimMap = vDimMap[roid];

	typedef std::pair<LFEID,const LocalShapeFunctionSet<TRefElem::dim>*> DimMapPair;

//	insert into map
	if(dimMap.insert(DimMapPair(type, &set)).second == false)
	{
		UG_LOG("ERROR in 'LocalShapeFunctionSetProvider::register_set()': "
				"Reference type already registered for trial space: "<<type<<" and "
				" Reference element type "<<roid<<".\n");
		return false;
	}

//	all ok
	return true;
}


template <typename TRefElem>
bool
LocalShapeFunctionSetProvider::
unregister_set(LFEID id)
{
//	Reference Object type
	static const ReferenceObjectID roid = TRefElem::REFERENCE_OBJECT_ID;

//	get type of map
	typedef std::map<LFEID, const ReferenceElemLocalShapeFunctionSet<TRefElem>* > Map;

//	init provider and get map
	static Map& map = inst().get_map<TRefElem>();

//	erase element
	if(map.erase(id) != 1) return false;

//	get map
	typedef std::map<LFEID, const LocalShapeFunctionSet<TRefElem::dim>* > DimMap;
	static DimMap* vDimMap = get_dim_map<TRefElem::dim>();
	DimMap& dimMap = vDimMap[roid];

//	erase element
	return (dimMap.erase(id) == 1);
}

template <typename TRefElem>
const ReferenceElemLocalShapeFunctionSet<TRefElem>&
LocalShapeFunctionSetProvider::
get(LFEID id)
{
//	get type of map
	typedef std::map<LFEID, const ReferenceElemLocalShapeFunctionSet<TRefElem>* > Map;
	const static ReferenceObjectID roid = TRefElem::REFERENCE_OBJECT_ID;

//	init provider and get map
	static Map& map = inst().get_map<TRefElem>();

//	search for identifier
	typename Map::const_iterator iter = map.find(id);

//	if not found
	if(iter == map.end())
	{
		UG_LOG("ERROR in 'LocalShapeFunctionSetProvider::get': "
				"Unknown Trial Space Type "<<id<<" requested for Element"
				" type: "<<roid<<".\n");
		throw(UG_ERROR_LocalShapeFunctionSetNotRegistered(TRefElem::dim, roid, id));
	}

//	return shape function set
	return *(iter->second);
}

template <int dim>
const LocalShapeFunctionSet<dim>&
LocalShapeFunctionSetProvider::
get(ReferenceObjectID roid, LFEID id)
{
//	get type of map
	typedef std::map<LFEID, const LocalShapeFunctionSet<dim>* > Map;

//	init provider and get map
	static Map* vMap = inst().get_dim_map<dim>();

//	get map for ref elem type
	const Map& map = vMap[roid];

//	search for identifier
	typename Map::const_iterator iter = map.find(id);

//	if not found
	if(iter == map.end())
	{
		UG_LOG("ERROR in 'LocalShapeFunctionSetProvider::get': "
				"Unknown Trial Space Type "<<id<<" requested for Element"
				" type: "<<roid<<".\n");
		throw(UG_ERROR_LocalShapeFunctionSetNotRegistered(dim, roid, id));
	}

//	return shape function set
	return *(iter->second);
}

}
#endif /* __H__UG__LIB_DISCRETIZATION__LOCAL_FINITE_ELEMENT__LOCAL_SHAPE_FUNCTION_SET_IMPL__ */
