/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
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
		static constexpr bool staticLocalData = TGeom::staticLocalData;

	protected:
		/// private constructor
		GeomProvider() {m_mLFEIDandOrder.clear();}

		/// destructor
		~GeomProvider() {clear_geoms();}

		/// singleton provider
		static GeomProvider& inst() {
			static GeomProvider inst;
			return inst;
		}

		/// struct to sort keys
		struct LFEIDandQuadOrder {
				LFEIDandQuadOrder(const LFEID lfeID, const int order)
				: m_lfeID(lfeID), m_order(order) {}

			///	operator <
				bool operator < (const LFEIDandQuadOrder& v) const
				{
					if(m_lfeID != v.m_lfeID) return m_lfeID < v.m_lfeID;
					else return m_order < v.m_order;
				}

			const LFEID m_lfeID;
			const int m_order;
		};

		/// vector holding instances
		using MapType = std::map<LFEIDandQuadOrder, TGeom*>;
		static MapType m_mLFEIDandOrder;

		/// returns class based on identifier
		static TGeom& get_class(const LFEID lfeID, const int quadOrder) {

			LFEIDandQuadOrder key(lfeID, quadOrder);

			using ret_type = std::pair<typename MapType::iterator,bool>;
			ret_type ret = m_mLFEIDandOrder.insert(std::pair<LFEIDandQuadOrder,TGeom*>(key,nullptr));

			// newly inserted, need construction of data
			if(ret.second == true){
				if(ret.first->second != nullptr){
					UG_THROW("Newly inserted element must have Pointer = nullptr");
				}
				else{
					ret.first->second = new TGeom();
				}
			}

			return *ret.first->second;
		}

		/// clears all instances
		static void clear_geoms(){
			using MapIter = typename std::map<LFEIDandQuadOrder, TGeom*>::iterator;
			for(MapIter iter = m_mLFEIDandOrder.begin(); iter != m_mLFEIDandOrder.end(); ++iter)
				if(iter->second)
					delete iter->second;

			m_mLFEIDandOrder.clear();
		}

	public:
		///	type of provided object
		using Type = TGeom;

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
			inst().clear_geoms();
		}
};

template <typename TGeom>
std::map<typename GeomProvider<TGeom>::LFEIDandQuadOrder, TGeom*> GeomProvider<TGeom>::m_mLFEIDandOrder
						= std::map<LFEIDandQuadOrder, TGeom*>();


} // end namespace ug

#endif