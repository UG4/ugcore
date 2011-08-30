/*
 * provider.h
 *
 *  Created on: 20.07.2011
 *      Author: andreasvogel
 */

#ifndef __H__COMMON__UTIL__PROVIDER__
#define __H__COMMON__UTIL__PROVIDER__

namespace ug{

/// Provider, holding a single instance of an object
/**
 * This class is used to wrap an object into a singleton-like provider, such
 * that construction computations is avoided, if the object is used several times.
 *
 * In addition, the object can be shared between unrelated code parts, if the
 * same object is intended to be used, but no passing is possible or wanted.
 */
template <typename TClass>
struct Provider
{
	///	type of provided object
		typedef TClass Type;

	///	returns access to the singleton
		inline static TClass& get()
		{
			static TClass myInst;
			return myInst;
		}
};

} // end namespace ug

#endif /* __H__COMMON__UTIL__PROVIDER__ */
