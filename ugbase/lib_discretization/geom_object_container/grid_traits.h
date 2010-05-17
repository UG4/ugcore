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
		template <class TGeomObj>
		static typename geometry_traits<TGeomObj>::iterator
		begin(TGeomObjHandler& handler, int subsetInd, int level);

		template <class TGeomObj>
		static typename geometry_traits<TGeomObj>::iterator
		end(TGeomObjHandler& handler, int subsetInd, int level);
		
		
		static size_t num_levels(TGeomObjHandler& handler);

		template <class TGeomObj>
		static size_t get_level(TGeomObjHandler& handler, TGeomObj* obj);
		
		
		static size_t num_subsets(TGeomObjHandler& handler);
		
		template <class TGeomObj>
		static int get_subset_index(TGeomObjHandler& handler, TGeomObj* obj);


		template <class TGeomObj>
		static bool is_shadow(TGeomObjHandler& handler, TGeomObj* obj);
		
};


}//	end of namespace

////////////////////////////////
//	include implementation
#include "grid_traits_impl.hpp"

#endif	// __H__UG__GRID_TRAITS__