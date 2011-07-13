/*
 * trialspacefactory_impl.h
 *
 *  Created on: 17.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISCRETIZATION__LOCAL_SHAPE_FUNCTION_SET_FACTORY_IMPL__
#define __H__UG__LIB_DISCRETIZATION__LOCAL_SHAPE_FUNCTION_SET_FACTORY_IMPL__

// include instances
#include "local_shape_function_set_id.h"
#include "lagrange/lagrangep1.h"
#include "lagrange/lagrange.h"


namespace ug{

///////////////////////////////////////////
// LocalShapeFunctionSetProvider
///////////////////////////////////////////

template <typename TRefElem>
bool
LocalShapeFunctionSetProvider::
init_standard_local_shape_function_sets()
{
//	create static Sets
	static LocalShapeFunctionSetWrapper<LagrangeP1<TRefElem, 1> > sSetLagrangeP1;
	static LocalShapeFunctionSetWrapper<LagrangeLSFS<TRefElem, 2> > sSetLagrangeP2;

//	insert into map: P1 Lagrange
	LSFSID type1(LSFSID::LAGRANGE, 1);
	if(!register_local_shape_function_set(type1, sSetLagrangeP1))
		return false;

//	insert into map: P2 Lagrange
	LSFSID type2(LSFSID::LAGRANGE, 2);
	if(!register_local_shape_function_set(type2, sSetLagrangeP2))
		return false;

//	return success
	return true;
}

template <typename TRefElem>
std::map<LSFSID, const LocalShapeFunctionSet<TRefElem>* >&
LocalShapeFunctionSetProvider::get_map()
{
//	get type of map
	typedef std::map<LSFSID, const LocalShapeFunctionSet<TRefElem>* > Map;

//	create static map
	static Map shape_function_set_map;

//	return map
	return shape_function_set_map;
};

template <typename TRefElem>
bool
LocalShapeFunctionSetProvider::
register_local_shape_function_set(LSFSID type, const LocalShapeFunctionSet<TRefElem>& set)
{
//	reference object id
	const ReferenceObjectID id = TRefElem::REFERENCE_OBJECT_ID;

//	get vector of types
	std::vector<const LocalShapeFunctionSetBase*>& vBase = m_baseMap[type];

//	resize vector
	vBase.resize(id+1, NULL);

//	check that no space has been previously registered to this place
	if(vBase[id])
	{
		UG_LOG("ERROR in 'LocalShapeFunctionSetProvider::"
				"register_local_shape_function_set()': "
				"Base type already registered for trial space: "<<type<<" and "
				" Reference element type "<<id<<".\n");
		return false;
	}

//	if ok, add
	vBase[id] = &set;

//	get type of map
	typedef std::map<LSFSID, const LocalShapeFunctionSet<TRefElem>* > Map;
	static Map& map = get_map<TRefElem>();
	typedef std::pair<LSFSID,const LocalShapeFunctionSet<TRefElem>*> MapPair;

//	insert into map
	if(map.insert(MapPair(type, &set)).second == false)
	{
		UG_LOG("ERROR in 'LocalShapeFunctionSetProvider::"
				"register_local_shape_function_set()': "
				"Reference type already registered for trial space: "<<type<<" and "
				" Reference element type "<<id<<".\n");
		return false;
	}

//	all ok
	return true;
}


template <typename TRefElem>
bool
LocalShapeFunctionSetProvider::
unregister_local_shape_function_set(LSFSID id)
{
//	get type of map
	typedef std::map<LSFSID, const LocalShapeFunctionSet<TRefElem>* > Map;

//	init provider and get map
	static Map& map = inst().get_map<TRefElem>();

//	erase element
	bool bRet = true;
	bRet &= (map.erase(id) == 1);
	bRet &= (m_baseMap.erase(id) == 1);
	return bRet;
}

template <typename TRefElem>
const LocalShapeFunctionSet<TRefElem>&
LocalShapeFunctionSetProvider::
get(LSFSID id)
{
//	get type of map
	typedef std::map<LSFSID, const LocalShapeFunctionSet<TRefElem>* > Map;

//	init provider and get map
	static Map& map = inst().get_map<TRefElem>();

//	search for identifier
	typename Map::const_iterator iter = map.find(id);

//	if not found
	if(iter == map.end())
	{
		UG_LOG("ERROR in 'LocalShapeFunctionSetProvider::get': "
				"Unknown Trial Space Type "<<id<<" requested.\n");
		throw(UGFatalError("Trial Space type unknown"));
	}

//	return shape function set
	return *(iter->second);
}

inline
const LocalShapeFunctionSetBase&
LocalShapeFunctionSetProvider::get(LSFSID id, ReferenceObjectID type)
{
//	init provider and search for identifier
	BaseMap::const_iterator iter = inst().m_baseMap.find(id);

//	if not found
	if(iter == m_baseMap.end())
	{
		UG_LOG("ERROR in 'LocalShapeFunctionSetProvider::get': "
				"Unknown Base Trial Space Type "<<id<<" requested.\n");
		throw(UGFatalError("Trial Space type unknown"));
	}

//	get vector
	const std::vector<const LocalShapeFunctionSetBase*>& vBase = iter->second;

//	check that space registered
	if(vBase[type] == NULL)
	{
		UG_LOG("ERROR in 'LocalShapeFunctionSetProvider::get': "
				"Unknown Base Trial Space  for Type "<<id<<" and Reference"
				" element "<<type<<" requested.\n");
		throw(UGFatalError("Trial Space type unknown"));
	}

//	return shape function set
	return *(vBase[type]);
}


}
#endif /* __H__UG__LIB_DISCRETIZATION__LOCAL_SHAPE_FUNCTION_SET_FACTORY_IMPL__ */
