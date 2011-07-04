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
		typedef std::map<	LSFSID,
							const LocalShapeFunctionSet<TRefElem>* > Map;

	//	get map
		Map& map = get_map<TRefElem>();

	//	insert into map: P1 Lagrange
		LSFSID type1(LSFSID::LAGRANGE, 1);
		success &= map.insert(
					std::pair<LSFSID,
						  const LocalShapeFunctionSet<TRefElem>*>
							(type1, &sSetLagrangeP1)).second;

	//	insert into map: P2 Lagrange
		LSFSID type2(LSFSID::LAGRANGE, 2);
		success &= map.insert(
					std::pair<LSFSID,
						  const LocalShapeFunctionSet<TRefElem>*>
							(type2, &sSetLagrangeP2)).second;

		init = true;
	}

//	return success
	return success;
}

template <typename TRefElem>
std::map<	LSFSID,
			const LocalShapeFunctionSet<TRefElem>* >&
LocalShapeFunctionSetProvider::
get_map()
{
//	get type of map
	typedef std::map<	LSFSID,
						const LocalShapeFunctionSet<TRefElem>* > Map;

//	create static map
	static Map m_shape_function_set_map;

//	return map
	return m_shape_function_set_map;
};

template <typename TRefElem>
bool
LocalShapeFunctionSetProvider::
register_local_shape_function_set(	LSFSID id,
									const LocalShapeFunctionSet<TRefElem>& set)
{
//	get type of map
	typedef std::map<	LSFSID,
						const LocalShapeFunctionSet<TRefElem>* > Map;

//	get map
	static Map& map = inst().get_map<TRefElem>();

	UG_LOG("Inserting id = " << id << " \n");

//	insert into map
	return map.insert(
			std::pair<LSFSID,
					  const LocalShapeFunctionSet<TRefElem>*>(id, &set)).second;
}


template <typename TRefElem>
bool
LocalShapeFunctionSetProvider::
unregister_local_shape_function_set(LSFSID id)
{
//	get type of map
	typedef std::map<	LSFSID,
						const LocalShapeFunctionSet<TRefElem>* > Map;

//	get map
	static Map& map = inst().get_map<TRefElem>();

//	erase element
	return map.erase(id) == 1;
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
		UG_ASSERT(0, "Unknown Trial Space Type.");
		throw(UG_ERROR_TrialSpaceNotRegistered());
	}

//	return shape function set
	return *(iter->second);
}


}
#endif /* __H__UG__LIB_DISCRETIZATION__LOCAL_SHAPE_FUNCTION_SET_FACTORY_IMPL__ */
