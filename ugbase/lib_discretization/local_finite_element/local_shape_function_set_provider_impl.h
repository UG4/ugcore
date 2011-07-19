/*
 * local_shape_function_set_impl.h
 *
 *  Created on: 17.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISCRETIZATION__LOCAL_SHAPE_FUNCTION_SET_FACTORY_IMPL__
#define __H__UG__LIB_DISCRETIZATION__LOCAL_SHAPE_FUNCTION_SET_FACTORY_IMPL__

#include "common/common.h"
#include "local_shape_function_set_provider.h"

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

template <typename TRefElem>
bool
LocalShapeFunctionSetProvider::
register_set(LFEID type, const LocalShapeFunctionSet<TRefElem>& set)
{
//	get type of map
	typedef std::map<LFEID, const LocalShapeFunctionSet<TRefElem>* > Map;
	static Map& map = get_map<TRefElem>();
	typedef std::pair<LFEID,const LocalShapeFunctionSet<TRefElem>*> MapPair;

//	insert into map
	if(map.insert(MapPair(type, &set)).second == false)
	{
		UG_LOG("ERROR in 'LocalShapeFunctionSetProvider::register_set()': "
				"Reference type already registered for trial space: "<<type<<" and "
				" Reference element type "<<TRefElem::REFERENCE_OBJECT_ID<<".\n");
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
//	get type of map
	typedef std::map<LFEID, const LocalShapeFunctionSet<TRefElem>* > Map;

//	init provider and get map
	static Map& map = inst().get_map<TRefElem>();

//	erase element
	return (map.erase(id) == 1);
}

template <typename TRefElem>
const LocalShapeFunctionSet<TRefElem>&
LocalShapeFunctionSetProvider::
get(LFEID id)
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
		UG_LOG("ERROR in 'LocalShapeFunctionSetProvider::get': "
				"Unknown Trial Space Type "<<id<<" requested for Element"
				" type: "<<roid<<".\n");
		throw(UGFatalError("Trial Space type unknown"));
	}

//	return shape function set
	return *(iter->second);
}

}
#endif /* __H__UG__LIB_DISCRETIZATION__LOCAL_SHAPE_FUNCTION_SET_FACTORY_IMPL__ */
