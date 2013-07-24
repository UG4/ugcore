/*
 * geom_provider.h
 *
 *  Created on: 20.07.2011
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__DISC_UTIL__GEOM_PROVIDER__
#define __H__UG__LIB_DISC__SPATIAL_DISC__DISC_UTIL__GEOM_PROVIDER__

#include <map>
#include "lib_disc/local_finite_element/local_finite_element_id.h"

namespace ug{


/// Geom Provider, holding a single instance of a geometry
/**
 * This class is used to wrap an object into a singleton-like provider, such
 * that construction computations is avoided, if the object is used several times.
 *
 * In addition, the object can be shared between unrelated code parts, if the
 * same object is intended to be used, but no passing is possible or wanted.
 */
template <typename TGeom>
class GeomProvider
{
	public:
		/// flag indicating if local data may change
		/**
		 * If the local data may change due to run-time setting of trial space
		 * or quadrature order, we must store one object for each type.
		 * In such a case we must store each combination in some data structure
		 * and return the correct one based on the identifier.
		 */
		static const bool staticLocalData = TGeom::staticLocalData;

	protected:
		/// private constructor
		GeomProvider() {m_mLFEIDandOrder.clear();}

		/// singleton provider
		static GeomProvider<TGeom>& inst() {
			static GeomProvider<TGeom> inst;
			return inst;
		}

		/// struct to sort keys
		struct LFEIDandQuadOrder{
				LFEIDandQuadOrder(const LFEID lfeID, const int order)
				: m_lfeID(lfeID), m_order(order) {}

			///	operator <
				bool operator<(const LFEIDandQuadOrder& v) const
				{
					if(m_lfeID != v.m_lfeID) return m_lfeID < v.m_lfeID;
					else return m_order < v.m_order;
				}

			const LFEID m_lfeID;
			const int m_order;
		};

		/// vector holding instances
		typedef std::map<LFEIDandQuadOrder, TGeom*> MapType;
		static MapType m_mLFEIDandOrder;

		/// returns class based on identifier
		static TGeom& get_class(const LFEID lfeID, const int quadOrder) {

			LFEIDandQuadOrder key(lfeID, quadOrder);

			typedef std::pair<typename MapType::iterator,bool> ret_type;
			ret_type ret = m_mLFEIDandOrder.insert(std::pair<LFEIDandQuadOrder,TGeom*>(key,NULL));

			// newly inserted, need construction of data
			if(ret.second == true){
				if(ret.first->second != NULL){
					UG_THROW("Newly inserted element must have Pointer = NULL");
				}
				else{
					ret.first->second = new TGeom();
				}
			}

			return *ret.first->second;
		}

		/// clears all instances
		static void clear_geoms(){
			typedef typename std::map<LFEIDandQuadOrder, TGeom*>::iterator MapIter;
			for(MapIter iter = m_mLFEIDandOrder.begin(); iter != m_mLFEIDandOrder.end(); ++iter)
				if(iter->second)
					delete iter->second;

			m_mLFEIDandOrder.clear();
		}

	public:
		///	type of provided object
		typedef TGeom Type;

		///	returns a singleton based on the identifier
		static inline TGeom& get(const LFEID lfeID, const int quadOrder){
			// in case of static data, use only one object
			if(staticLocalData) return get();

			// return the object based on identifier
			return inst().get_class(lfeID, quadOrder);
		}

		///	returns a singleton based on the identifier
		static inline TGeom& get(){
			static TGeom inst;
			if(!staticLocalData)
				UG_THROW("GeomProvider: accessing geometry without keys, but"
						 " geometry may change local data. Use access by keys instead.");
			return inst;
		}

		///	clears all singletons
		static inline void clear(){
			return inst().clear_geoms();
		}
};

template <typename TGeom>
std::map<typename GeomProvider<TGeom>::LFEIDandQuadOrder, TGeom*> GeomProvider<TGeom>::m_mLFEIDandOrder
						= std::map<typename GeomProvider<TGeom>::LFEIDandQuadOrder, TGeom*>();


} // end namespace ug

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__DISC_UTIL__GEOM_PROVIDER__ */
