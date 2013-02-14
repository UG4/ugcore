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

template <int dim, typename t_shape, typename t_grad>
std::map<LFEID, const LocalShapeFunctionSet<dim, t_shape, t_grad>* >*
LocalShapeFunctionSetProvider::get_map()
{
//	get type of map
	typedef std::map<LFEID, const LocalShapeFunctionSet<dim, t_shape, t_grad>*> Map;

//	create static map
	static Map sShapeFunctionSetMap[NUM_REFERENCE_OBJECTS];

//	return map
	return sShapeFunctionSetMap;
};

template <int dim, typename t_shape, typename t_grad>
std::vector<LocalShapeFunctionSet<dim, t_shape, t_grad>*>&
LocalShapeFunctionSetProvider::get_dynamic_allocated_vector()
{
//	get type of map
	typedef std::vector<LocalShapeFunctionSet<dim, t_shape, t_grad>*> Vec;

//	create static map
	static Vec sShapeFunctionSetMap;

//	return map
	return sShapeFunctionSetMap;
};

template <int dim, typename t_shape, typename t_grad>
void LocalShapeFunctionSetProvider::
register_set(LFEID type,
             const ReferenceObjectID roid,
             const LocalShapeFunctionSet<dim, t_shape, t_grad>& set)
{
//	get type of map
	typedef std::map<LFEID, const LocalShapeFunctionSet<dim, t_shape, t_grad>* > Map;
	Map& map = get_map<dim, t_shape, t_grad>()[roid];

//	insert into map
	typedef typename Map::value_type DimMapPair;
	if(map.insert(DimMapPair(type, &set)).second == false)
		UG_THROW("LocalShapeFunctionSetProvider::register_set(): "
				"Reference type already registered for trial space: "<<type<<" and "
				" Reference element type "<<roid<<".");

	std::map<LFEID, bool>& contMap = get_continuous_map();
	if(contMap.find(type) == contMap.end()){
		contMap.insert(std::map<LFEID, bool>::value_type(type, set.continuous()));
	}else{
		if(contMap[type] != set.continuous())
			UG_THROW("LocalShapeFunctionSetProvider::register_set(): "
					"Reference type says continuous:"<<set.continuous()<<", but "
					" other Reference element say"<<contMap[type]<<".");
	}
}


template <int dim, typename t_shape, typename t_grad>
bool LocalShapeFunctionSetProvider::
unregister_set(LFEID id)
{
	typedef std::map<LFEID, const LocalShapeFunctionSet<dim, t_shape, t_grad>* > Map;
	Map* vMap = get_map<dim, t_shape, t_grad>();

	for(int r = 0; r < NUM_REFERENCE_OBJECTS; ++r)
	{
		const ReferenceObjectID roid = (ReferenceObjectID) r;

	//	get map
		Map& map = vMap[roid];

	//	erase element
		if(!(map.erase(id) == 1)) return false;
	}

	std::map<LFEID, bool>& contMap = get_continuous_map();
	if(contMap.erase(id) == 1) return false;

	return true;
}

template <int dim, typename t_shape, typename t_grad>
const LocalShapeFunctionSet<dim, t_shape, t_grad>&
LocalShapeFunctionSetProvider::
get(ReferenceObjectID roid, LFEID id, bool bCreate)
{
//	init provider and get map
	typedef std::map<LFEID, const LocalShapeFunctionSet<dim, t_shape, t_grad>* > Map;
	const Map& map = inst().get_map<dim, t_shape, t_grad>()[roid];

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
			return get<dim, t_shape, t_grad>(roid, id, false);
		}

		UG_THROW("LocalShapeFunctionSetProvider: Local Shape Function Set not "
				 "found for "<<roid<<" (dim="<<dim<<") and type = "<<id);
	}

//	return shape function set
	return *(iter->second);
}

}
#endif /* __H__UG__LIB_DISC__LOCAL_FINITE_ELEMENT__LOCAL_SHAPE_FUNCTION_SET_IMPL__ */
