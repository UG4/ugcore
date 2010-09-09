/*
 * trialspacefactory.h
 *
 *  Created on: 17.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__LOCAL_SHAPE_FUNCTION_SET_FACTORY__
#define __H__LIBDISCRETIZATION__LOCAL_SHAPE_FUNCTION_SET_FACTORY__

#include "local_shape_function_set.h"

namespace ug {

enum LocalShapeFunctionSetID {
	LSFS_INVALID = 0,
	LSFS_LAGRANGEP1
};


// Singleton, holding all Trial Spaces available
class LocalShapeFunctionSetFactory {

		// private constructor
		LocalShapeFunctionSetFactory();

		// disallow copy and assignment (intentionally left unimplemented)
		LocalShapeFunctionSetFactory(const LocalShapeFunctionSetFactory&);
		LocalShapeFunctionSetFactory& operator=(const LocalShapeFunctionSetFactory&);

		// private destructor
		~LocalShapeFunctionSetFactory(){};

		// initialize the standard trialspaces (called during construction)
		template <typename TRefElem>
		bool init_standard_local_shape_function_sets();

		// return a map of element_trial_spaces
		template <typename TRefElem>
		std::map<LocalShapeFunctionSetID, const LocalShapeFunctionSet<TRefElem>* >& get_local_shape_function_set_map();

	public:
		// singleton provider
		static LocalShapeFunctionSetFactory& inst()
		{
			static LocalShapeFunctionSetFactory myInst;
			return myInst;
		};

		// register a local shape function set for a given reference element type
		template <typename TRefElem>
		bool register_local_shape_function_set(LocalShapeFunctionSetID id, const LocalShapeFunctionSet<TRefElem>& set);

		// unregister
		template <typename TRefElem>
		bool unregister_local_shape_function_set(LocalShapeFunctionSetID id);

		// get the local shape function set for a given reference element and id
		template <typename TRefElem>
		const LocalShapeFunctionSet<TRefElem>& get_local_shape_function_set(LocalShapeFunctionSetID id);
};

}

#include "local_shape_function_set_factory_impl.h"

#endif /* __H__LIBDISCRETIZATION__LOCAL_SHAPE_FUNCTION_SET_FACTORY__ */
