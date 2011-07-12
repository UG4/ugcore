/*
 * local_shape_function_set
 *
 *  Created on: 17.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISCRETIZATION__LOCAL_SHAPE_FUNCTION_SET_FACTORY__
#define __H__UG__LIB_DISCRETIZATION__LOCAL_SHAPE_FUNCTION_SET_FACTORY__

#include "local_shape_function_set.h"
#include "local_shape_function_set_id.h"

namespace ug {

// LocalShapeFunctionSetProvider
/** class to provide local shape function sets
 *
 *	This class provides references to Local Shape functions sets.
 *	It is implemented as a Singleton.
 */
class LocalShapeFunctionSetProvider {
	private:
	// 	disallow private constructor
		LocalShapeFunctionSetProvider();

	// disallow copy and assignment (intentionally left unimplemented)
		LocalShapeFunctionSetProvider(const LocalShapeFunctionSetProvider&);
		LocalShapeFunctionSetProvider& operator=(const LocalShapeFunctionSetProvider&);

	// 	private destructor
		~LocalShapeFunctionSetProvider(){};

	// 	Singleton provider
		static LocalShapeFunctionSetProvider& inst()
		{
			static LocalShapeFunctionSetProvider myInst;
			return myInst;
		};

	private:
	// 	initialize the standard trialspaces (called during construction)
		template <typename TRefElem>
		bool init_standard_local_shape_function_sets();

	// 	return a map of element_trial_spaces
		template <typename TRefElem>
		static std::map<LSFSID, const LocalShapeFunctionSet<TRefElem>* >&
		get_map();

		inline
		std::map<LSFSID,const LocalShapeFunctionSetBase*>&
		get_base_map();

	public:
	/** register a local shape function set for a given reference element type
	 * This function is used to register a Local Shape Function set for an element
	 * type and the corresponding local shape function set id.
	 *
	 * \param[in]		id 		Identifier for local shape function set
	 * \param[in]		set		Local Shape Function Set to register
	 * \return			bool	true iff registration successful
	 */
		template <typename TRefElem>
		static bool
		register_local_shape_function_set(LSFSID id,
		                                  const LocalShapeFunctionSet<TRefElem>& set,
		                                  const LocalShapeFunctionSetBase& baseSet);

	/** unregister a local shape function set for a given reference element type
	 * This function is used to unregister a Local Shape Function set for an element
	 * type and the corresponding local shape function set id from this Provider.
	 *
	 * \param[in]		id 		Identifier for local shape function set
	 * \return			bool	true iff removal successful
	 */
		template <typename TRefElem>
		static bool unregister_local_shape_function_set(LSFSID id);

	/**	returns the Local Shape Function Set
	 * This function returns the Local Shape Function Set for a reference element
	 * type and an Identifier if a set has been registered for the identifer. Else
	 * an exception is thrown.
	 *
	 * \param[in]	id		Identifier for local shape function set
	 * \return 		set		A const reference to the shape function set
	 */
		// get the local shape function set for a given reference element and id
		template <typename TRefElem>
		static const LocalShapeFunctionSet<TRefElem>&
		get(LSFSID id);

	///	returns the local shape function set base for an id
		inline static const LocalShapeFunctionSetBase& get(LSFSID id);
};


/// Singleton, holding a single local shape function set
/**
 * This class is used to wrap local shape function set into a singleton, such
 * that construction computations is avoided, if the rule is used several times.
 */
class LSFSProvider {

	// 	private constructor
		LSFSProvider();

	// 	disallow copy and assignment (intentionally left unimplemented)
		LSFSProvider(const LSFSProvider&);
		LSFSProvider& operator=(const LSFSProvider&);

	// 	private destructor
		~LSFSProvider(){};

	// 	geometry provider, holding the instance
		template <typename TLSFS>
		inline static TLSFS& inst()
		{
			static TLSFS myInst;
			return myInst;
		};

	public:
	///	returns access to the singleton
		template <typename TLSFS>
		inline static TLSFS& get() {return inst<TLSFS>();}
};


}

#include "local_shape_function_set_provider_impl.h"

#endif /* __H__UG__LIB_DISCRETIZATION__LOCAL_SHAPE_FUNCTION_SET_FACTORY__ */
