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
 */
template <typename TClass>
class Provider
{
	public:
		///	type of provided object
		typedef TClass Type;

		///	returns a singleton based on the identifier
		static inline TClass& get(){
			static TClass inst;
			return inst;
		}
};

} // end namespace ug

#endif /* __H__COMMON__UTIL__PROVIDER__ */
