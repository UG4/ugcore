/*
 * provider.h
 *
 *  Created on: 20.07.2011
 *      Author: andreasvogel
 */

#ifndef __H__COMMON__UTIL__PROVIDER__
#define __H__COMMON__UTIL__PROVIDER__

#include <vector>

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
class Provider
{
	protected:
		/// private constructor
		Provider() {}

		/// singleton provider
		static Provider<TClass>& inst() {
			static Provider<TClass> inst;
			return inst;
		}

		/// vector holding instances
		static std::vector<TClass*> m_vClass;

		/// returns class based on identifier
		static TClass& get_class(size_t p) {
			if(p >= m_vClass.size()){
				m_vClass.resize(p+1, NULL);
			}
			if(m_vClass[p] == NULL){
				m_vClass[p] = new TClass();
			}
			return *m_vClass[p];
		}

		/// clears all instances
		static void clear_classes(){
			for(size_t i = 0; m_vClass.size(); ++i)
				if(m_vClass[i])
					delete m_vClass[i];

			m_vClass.clear();
		}

	public:
		///	type of provided object
		typedef TClass Type;

		///	returns a singleton based on the identifier
		static inline TClass& get(size_t p){
			return inst().get_class(p);
		}

		///	returns a singleton based on the identifier
		static inline TClass& get(){
			static TClass inst;
			return inst;
		}

		///	clears all singletons
		static inline void clear(){
			return inst().clear_classes();
		}
};

template <typename TClass>
std::vector<TClass*> Provider<TClass>::m_vClass = std::vector<TClass*>();


} // end namespace ug

#endif /* __H__COMMON__UTIL__PROVIDER__ */
