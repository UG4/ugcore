/*
 * provider.h
 *
 *  Created on: 20.07.2011
 *      Author: andreasvogel
 */

#ifndef __H__COMMON__UTIL__PROVIDER__
#define __H__COMMON__UTIL__PROVIDER__

namespace ug{

/// Singleton, holding a single instance of an object
/**
 * This class is used to wrap an object set into a singleton-like provider, such
 * that construction computations is avoided, if the object is used several times.
 *
 * In addition, the object can be shared between unrelated code parts, if the
 * same object is intended to be used, but no passing is possible or wanted.
 */
class Provider
{
	// 	private constructor
		Provider();

	// 	disallow copy and assignment (intentionally left unimplemented)
		Provider(const Provider&);
		Provider& operator=(const Provider&);

	// 	private destructor
		~Provider(){};

	// 	holding the instance
		template <typename TClass>
		inline static TClass& inst()
		{
			static TClass myInst;
			return myInst;
		};

	public:
	///	returns access to the singleton
		template <typename TClass>
		inline static TClass& get() {return inst<TClass>();}
};

} // end namespace ug

#endif /* __H__COMMON__UTIL__PROVIDER__ */
