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

template <int dim, typename TShape, typename TGrad>
std::map<LFEID, ConstSmartPtr<LocalShapeFunctionSet<dim, TShape, TGrad> > >*
LocalShapeFunctionSetProvider::get_maps()
{
//	get type of map
	typedef std::map<LFEID, ConstSmartPtr<LocalShapeFunctionSet<dim, TShape, TGrad> > > Map;

//	create static map
	static Map sShapeFunctionSetMap[NUM_REFERENCE_OBJECTS];

//	return map
	return sShapeFunctionSetMap;
};

template <int dim, typename TShape, typename TGrad>
void LocalShapeFunctionSetProvider::
register_set(LFEID type,
             const ReferenceObjectID roid,
             ConstSmartPtr<LocalShapeFunctionSet<dim, TShape, TGrad> > set)
{
//	get type of map
	typedef std::map<LFEID, ConstSmartPtr<LocalShapeFunctionSet<dim, TShape, TGrad> > > Map;
	Map& map = get_maps<dim, TShape, TGrad>()[roid];

//	insert into map
	typedef typename Map::value_type DimMapPair;
	if(map.insert(DimMapPair(type, set)).second == false)
		UG_THROW("LocalShapeFunctionSetProvider::register_set(): "
				"Reference type already registered for trial space: "<<type<<" and "
				" Reference element type "<<roid<<".");

	if(m_mContSpace.find(type) == m_mContSpace.end()){
		m_mContSpace.insert(std::map<LFEID, bool>::value_type(type, set->continuous()));
	}else{
		if(m_mContSpace[type] != set->continuous())
			UG_THROW("LocalShapeFunctionSetProvider::register_set(): "
					"Reference type says continuous:"<<set->continuous()<<", but "
					" other Reference element says "<<m_mContSpace[type]<<".");
	}
}


template <int dim, typename TShape, typename TGrad>
bool LocalShapeFunctionSetProvider::
unregister_set(LFEID id)
{
	typedef std::map<LFEID, ConstSmartPtr<LocalShapeFunctionSet<dim, TShape, TGrad> > > Map;
	Map* vMap = get_maps<dim, TShape, TGrad>();

	for(int r = 0; r < NUM_REFERENCE_OBJECTS; ++r)
	{
		const ReferenceObjectID roid = (ReferenceObjectID) r;

	//	get map
		Map& map = vMap[roid];

	//	erase element
		if(!(map.erase(id) == 1)) return false;
	}

	if(m_mContSpace.erase(id) == 1) return false;

	return true;
}

template <int dim, typename TShape, typename TGrad>
ConstSmartPtr<LocalShapeFunctionSet<dim, TShape, TGrad> >
LocalShapeFunctionSetProvider::
getptr(ReferenceObjectID roid, LFEID id, bool bCreate)
{
//	init provider and get map
	typedef std::map<LFEID, ConstSmartPtr<LocalShapeFunctionSet<dim, TShape, TGrad> > > Map;
	const Map& map = inst().get_maps<dim, TShape, TGrad>()[roid];

//	search for identifier
	typename Map::const_iterator iter = map.find(id);

//	if not found
	if(iter == map.end())
	{
		if(bCreate)
		{
		//	try to create the set
			create_set(roid, id);

		//	next try to return the set
			return getptr<dim, TShape, TGrad>(roid, id, false);
		}

		UG_THROW("LocalShapeFunctionSetProvider: Local Shape Function Set not "
				 "found for "<<roid<<" (world dim: "<<dim<<") and type = "<<id<<
				 ". (This is usually due to: a) The function set is not implemented at "
				 " all, or b) The finite element space is discontinuous but the "
				 "evaluation is requested on a subelement, i.e. a grid object "
				 "with dimension less than the dimension where the finite element"
				 " space is defined.)");
	}

//	return shape function set
	return (iter->second);
}

template <int dim, typename TShape, typename TGrad>
const LocalShapeFunctionSet<dim, TShape, TGrad>&
LocalShapeFunctionSetProvider::
get(ReferenceObjectID roid, LFEID id, bool bCreate)
{
	return *getptr<dim,TShape,TGrad>(roid, id, bCreate);
}

}
#endif /* __H__UG__LIB_DISC__LOCAL_FINITE_ELEMENT__LOCAL_SHAPE_FUNCTION_SET_IMPL__ */
