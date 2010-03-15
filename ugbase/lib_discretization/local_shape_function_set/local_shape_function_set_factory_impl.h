/*
 * trialspacefactory_impl.h
 *
 *  Created on: 17.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__LOCAL_SHAPE_FUNCTION_SET_FACTORY_IMPL__
#define __H__LIBDISCRETIZATION__LOCAL_SHAPE_FUNCTION_SET_FACTORY_IMPL__

// include instances
#include "LagrangeP1/p1conform.h"


namespace ug{

///////////////////////////////////////////
// LocalShapeFunctionSetFactory

template <typename TRefElem>
bool
LocalShapeFunctionSetFactory::
init_standard_local_shape_function_sets()
{
	bool init = true;

	init = init && register_local_shape_function_set<TRefElem>(LSFS_LAGRANGEP1, P1conform<TRefElem>::inst());

	return init;
}

template <typename TRefElem>
std::map<LocalShapeFunctionSetID, const LocalShapeFunctionSet<TRefElem>* >&
LocalShapeFunctionSetFactory::
get_local_shape_function_set_map()
{
	typedef std::map<LocalShapeFunctionSetID, const LocalShapeFunctionSet<TRefElem>* > Map;
	static Map m_shape_function_set_map;

	return m_shape_function_set_map;
};

template <typename TRefElem>
bool
LocalShapeFunctionSetFactory::
register_local_shape_function_set(LocalShapeFunctionSetID id, const LocalShapeFunctionSet<TRefElem>& set)
{
	typedef std::map<LocalShapeFunctionSetID, const LocalShapeFunctionSet<TRefElem>* > Map;
	Map& map = get_local_shape_function_set_map<TRefElem>();
	return map.insert(std::pair<LocalShapeFunctionSetID, const LocalShapeFunctionSet<TRefElem>*>(id, &set)).second;
}


template <typename TRefElem>
bool
LocalShapeFunctionSetFactory::
unregister_local_shape_function_set(LocalShapeFunctionSetID id)
{
	typedef std::map<LocalShapeFunctionSetID, const LocalShapeFunctionSet<TRefElem>* > Map;
	Map& map = get_local_shape_function_set_map<TRefElem>();

	return map.erase(id) == 1;
}

template <typename TRefElem>
const LocalShapeFunctionSet<TRefElem>&
LocalShapeFunctionSetFactory::
get_local_shape_function_set(LocalShapeFunctionSetID id)
{
	typedef std::map<LocalShapeFunctionSetID, const LocalShapeFunctionSet<TRefElem>* > Map;
	Map& map = get_local_shape_function_set_map<TRefElem>();

	typename Map::const_iterator iter = map.find(id);

		if(iter == map.end())
		{
			assert(0 && "Unknown Trial Space Type.\n");
		}

		return *(iter->second);
}


}
#endif /* __H__LIBDISCRETIZATION__LOCAL_SHAPE_FUNCTION_SET_FACTORY_IMPL__ */
