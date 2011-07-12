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
	static bool init = false;
	bool success = true;

//	create static Sets
	static LocalShapeFunctionSetWrapper<LagrangeP1<TRefElem, 1> > sSetLagrangeP1;
	static LocalShapeFunctionSetWrapper<LagrangeLSFS<TRefElem, 2> > sSetLagrangeP2;

	if(!init)
	{
	//	get type of map
		typedef std::map<LSFSID, const LocalShapeFunctionSet<TRefElem>* > Map;
		typedef std::map<LSFSID, const LocalShapeFunctionSetBase* > BaseMap;

	//	get map
		Map& map = get_map<TRefElem>();
		BaseMap& baseMap = get_base_map();

	//	insert into map: P1 Lagrange
		LSFSID type1(LSFSID::LAGRANGE, 1);
		success &= map.insert(std::pair<LSFSID,const LocalShapeFunctionSet<TRefElem>*>
													(type1, &sSetLagrangeP1)).second;
		success &= baseMap.insert(std::pair<LSFSID,const LocalShapeFunctionSetBase*>
													(type1, &sSetLagrangeP1)).second;

	//	insert into map: P2 Lagrange
		LSFSID type2(LSFSID::LAGRANGE, 2);
		success &= map.insert(std::pair<LSFSID, const LocalShapeFunctionSet<TRefElem>*>
													(type2, &sSetLagrangeP2)).second;
		success &= baseMap.insert(std::pair<LSFSID, const LocalShapeFunctionSetBase*>
													(type2, &sSetLagrangeP2)).second;

		init = true;
	}

//	return success
	return success;
}

template <typename TRefElem>
std::map<LSFSID, const LocalShapeFunctionSet<TRefElem>* >&
LocalShapeFunctionSetProvider::get_map()
{
//	get type of map
	typedef std::map<LSFSID, const LocalShapeFunctionSet<TRefElem>* > Map;

//	create static map
	static Map m_shape_function_set_map;

//	return map
	return m_shape_function_set_map;
};

inline
std::map<LSFSID,const LocalShapeFunctionSetBase*>&
LocalShapeFunctionSetProvider::get_base_map()
{
//	get type of map
	typedef std::map<LSFSID, const LocalShapeFunctionSetBase*> Map;

//	create static map
	static Map m_shape_function_set_base_map;

//	return map
	return m_shape_function_set_base_map;
};

template <typename TRefElem>
bool
LocalShapeFunctionSetProvider::
register_local_shape_function_set(	LSFSID id,
									const LocalShapeFunctionSet<TRefElem>& set,
									const LocalShapeFunctionSetBase& baseSet)
{
//	get type of map
	typedef std::map<LSFSID, const LocalShapeFunctionSet<TRefElem>* > Map;
	typedef std::map<LSFSID, const LocalShapeFunctionSetBase*> BaseMap;

//	get map
	static Map& map = inst().get_map<TRefElem>();
	static BaseMap& baseMap = inst().get_base_map();

	bool bRet = true;

//	insert into map
	bRet &= map.insert(std::pair<LSFSID,const LocalShapeFunctionSet<TRefElem>*>(id, &set)).second;
	bRet &= baseMap.insert(std::pair<LSFSID,const LocalShapeFunctionSetBase*>(id, &baseSet)).second;

	return bRet;
}


template <typename TRefElem>
bool
LocalShapeFunctionSetProvider::
unregister_local_shape_function_set(LSFSID id)
{
//	get type of map
	typedef std::map<LSFSID, const LocalShapeFunctionSet<TRefElem>* > Map;
	typedef std::map<LSFSID, const LocalShapeFunctionSetBase*> BaseMap;

//	get map
	static Map& map = inst().get_map<TRefElem>();
	static BaseMap& baseMap = inst().get_base_map();

	bool bRet = true;

//	erase element
	bRet &= (map.erase(id) == 1);
	bRet &= (baseMap.erase(id) == 1);

	return bRet;
}

template <typename TRefElem>
const LocalShapeFunctionSet<TRefElem>&
LocalShapeFunctionSetProvider::
get(LSFSID id)
{
//	get type of map
	typedef std::map<	LSFSID,
						const LocalShapeFunctionSet<TRefElem>* > Map;

//	get map
	static Map& map = inst().get_map<TRefElem>();

//	search for identifier
	typename Map::const_iterator iter = map.find(id);

//	if not found
	if(iter == map.end())
	{
		UG_LOG("ERROR in 'LocalShapeFunctionSetProvider::get': "
				"Unknown Trial Space Type id = "<<id<<" requested.\n");
		throw(UGFatalError("Trial Space type unknown"));
	}

//	return shape function set
	return *(iter->second);
}

inline
const LocalShapeFunctionSetBase&
LocalShapeFunctionSetProvider::get(LSFSID id)
{
//	get type of map
	typedef std::map<LSFSID, const LocalShapeFunctionSetBase* > Map;

//	get map
	static Map& map = inst().get_base_map();

//	search for identifier
	Map::const_iterator iter = map.find(id);

//	if not found
	if(iter == map.end())
	{
		UG_LOG("ERROR in 'LocalShapeFunctionSetProvider::get': "
				"Unknown Trial Space Type id = "<<id<<" requested.\n");
		throw(UGFatalError("Trial Space type unknown"));
	}

//	return shape function set
	return *(iter->second);
}


}
#endif /* __H__UG__LIB_DISCRETIZATION__LOCAL_SHAPE_FUNCTION_SET_FACTORY_IMPL__ */
