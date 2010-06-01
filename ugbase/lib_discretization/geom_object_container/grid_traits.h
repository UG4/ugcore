//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m03 d24

#ifndef __H__UG__GRID_TRAITS__
#define __H__UG__GRID_TRAITS__

namespace ug
{
template <class TGeomObjHandler>
class grid_traits
{
	public:
		/// returns the number of GeomObjects of type 'TGeomObj' for a Subset and a level
		template <class TGeomObj>
		static size_t num(TGeomObjHandler& handler, int subsetInd, int level);

		/// returns iterator to first GeomObject of type 'TGeomObj' for Subset and level
		template <class TGeomObj>
		static typename geometry_traits<TGeomObj>::iterator
		begin(TGeomObjHandler& handler, int subsetInd, int level);

		/// returns iterator to last GeomObject of type 'TGeomObj' for Subset and level
		template <class TGeomObj>
		static typename geometry_traits<TGeomObj>::iterator
		end(TGeomObjHandler& handler, int subsetInd, int level);

		/// returns the number of levels in GeomObjHandler
		static size_t num_levels(TGeomObjHandler& handler);

		/// returns the level of the GeomObj
		template <class TGeomObj>
		static size_t get_level(TGeomObjHandler& handler, TGeomObj* obj);

		/// returns the number of subsets
		static size_t num_subsets(TGeomObjHandler& handler);

		/// returns the subset index of the GeomObj
		template <class TGeomObj>
		static int get_subset_index(TGeomObjHandler& handler, TGeomObj* obj);

		/// returns if the GeomObj is a shadow element
		/**
		 * Returns if a GeomObj is a shadow element. A shadow element is an element
		 * of a surface view of a multigrid, that has a child. Those elements of
		 * lower dimension (e.g. vertices) must be introduced to create the surface
		 * grid on lower levels.
		 */
		template <class TGeomObj>
		static bool is_shadow(TGeomObjHandler& handler, TGeomObj* obj);

		static void enable_subset_attachments(TGeomObjHandler& handler, bool bEnable);

		template <class TGeomObj>
		static inline void attach_to(TGeomObjHandler& handler, IAttachment& attachment, int subsetIndex);

		template <class TGeomObj, class TAttachment>
		static void attach_to_dv(TGeomObjHandler& handler, TAttachment& attachment, int subsetIndex,
						const typename TAttachment::ValueType& defaultValue);

		template <class TGeomObj>
		static void detach_from(TGeomObjHandler& handler, IAttachment& attachment, int subsetIndex);


};


}//	end of namespace

////////////////////////////////
//	include implementation
#include "grid_traits_impl.hpp"

#endif	// __H__UG__GRID_TRAITS__
