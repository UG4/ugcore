//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m03 d24

#ifndef __H__UG__GRID_TRAITS_IMPL__
#define __H__UG__GRID_TRAITS_IMPL__

#include "surface_view.h"

namespace ug
{

template <class THandler> class Surface
{
	public:
		Surface(THandler& handler) : m_handler(handler) {}
		
		inline THandler& handler()	{return m_handler;}
				
	protected:
		THandler& m_handler;
};

////////////////////////////////////////////////////////////////////////
///	specialization of grid_traits for ug::SurfaceView
template <>
class grid_traits<SurfaceView>
{
	public:
		typedef SurfaceView Handler;
		
		template <class TGeomObj>
		static inline typename geometry_traits<TGeomObj>::iterator
		begin(Handler& handler, int subsetInd, int level)
		{
			return handler.begin<TGeomObj>(subsetInd);
		}

		template <class TGeomObj>
		static inline typename geometry_traits<TGeomObj>::iterator
		end(Handler& handler, int subsetInd, int level)
		{
			return handler.end<TGeomObj>(subsetInd);
		}
		
		static inline size_t
		num_levels(Handler& handler)
		{
			return 1;
		}

		template <class TGeomObj>
		static inline size_t
		get_level(Handler& handler, TGeomObj* obj)
		{
			return 0;
		}
		
		
		static inline size_t
		num_subsets(Handler& handler)
		{
			return handler.num_subsets();
		}
		
		template <class TGeomObj>
		static inline int
		get_subset_index(Handler& handler, TGeomObj* obj)
		{
			return handler.get_subset_index(obj);
		}


		template <class TGeomObj>
		static inline bool
		is_shadow(Handler& handler, TGeomObj* obj)
		{
			handler.is_shadow(obj);
		}
};

////////////////////////////////////////////////////////////////////////
///	specialization of grid_traits for ug::MGSubsetHandler
template <>
class grid_traits<MGSubsetHandler>
{
	public:
		typedef MGSubsetHandler Handler;
		
		template <class TGeomObj>
		static inline typename geometry_traits<TGeomObj>::iterator
		begin(Handler& handler, int subsetInd, int level)
		{
			return handler.begin<TGeomObj>(subsetInd, level);
		}

		template <class TGeomObj>
		static inline typename geometry_traits<TGeomObj>::iterator
		end(Handler& handler, int subsetInd, int level)
		{
			return handler.end<TGeomObj>(subsetInd, level);
		}
		
		static inline size_t
		num_levels(Handler& handler)
		{
			return handler.num_levels();
		}

		template <class TGeomObj>
		static inline size_t
		get_level(Handler& handler, TGeomObj* obj)
		{
			return handler.get_level(obj);
		}
		
		
		static inline size_t
		num_subsets(Handler& handler)
		{
			return handler.num_subsets();
		}
		
		template <class TGeomObj>
		static inline int
		get_subset_index(Handler& handler, TGeomObj* obj)
		{
			return handler.get_subset_index(obj);
		}


		template <class TGeomObj>
		static inline bool
		is_shadow(Handler& handler, TGeomObj* obj)
		{
			return false;
		}
};

////////////////////////////////////////////////////////////////////////
///	specialization of grid_traits for ug::MGTopLevelSubsetHandler
template <>
class grid_traits<Surface<MGSubsetHandler> >
{
	public:
		typedef Surface<MGSubsetHandler> Handler;
		
		template <class TGeomObj>
		static inline typename geometry_traits<TGeomObj>::iterator
		begin(Handler& handler, int subsetInd, int level)
		{
			return handler.handler().begin<TGeomObj>(subsetInd,
										handler.handler().num_levels() - 1);
		}

		template <class TGeomObj>
		static inline typename geometry_traits<TGeomObj>::iterator
		end(Handler& handler, int subsetInd, int level)
		{
			return handler.handler().end<TGeomObj>(subsetInd,
										handler.handler().num_levels() - 1);
		}
		
		static inline size_t
		num_levels(Handler& handler)
		{
			return 1;
		}

		template <class TGeomObj>
		static inline size_t
		get_level(Handler& handler, TGeomObj* obj)
		{
			return 0;
		}
		
		
		static inline size_t
		num_subsets(Handler& handler)
		{
			return handler.handler().num_subsets();
		}
		
		template <class TGeomObj>
		static inline int
		get_subset_index(Handler& handler, TGeomObj* obj)
		{
			return handler.handler().get_subset_index(obj);
		}


		template <class TGeomObj>
		static inline bool
		is_shadow(Handler& handler, TGeomObj* obj)
		{
			return false;
		}
};

////////////////////////////////////////////////////////////////////////
/**
 * This method calls all methods declared in the GeomObjContainer-
 * implementation that is associated with the given handler.
 *
 * This helps to find errors during compile-time.
 */
template <class TGeomObjHandler>
void TestGeomObjContainerImplementation(TGeomObjHandler& handler)
{
	typedef grid_traits<TGeomObjHandler> ObjCon;
	
	UG_LOG("  num levels: " << ObjCon::num_levels(handler) << "\n");
	UG_LOG("  num subsets: " << ObjCon::num_subsets(handler) << "\n");
	
	bool bEmpty = ObjCon::template begin<VertexBase>(handler, 0, 0)
					== ObjCon::template end<VertexBase>(handler, 0, 0);
	if(bEmpty){
		UG_LOG("  handler has no vertices on level 0 subset 0.\n");
	}
	else{
		UG_LOG("  handler has vertices on level 0 subset 0.\n");
	}
	
	if(!bEmpty){
		VertexBase* v = *ObjCon::template begin<VertexBase>(handler, 0, 0);
		UG_LOG("  infos on first vertex on level 0, subset 0:\n");
		UG_LOG("    level: " << ObjCon::get_level(handler, v) << "\n");
		UG_LOG("    subset: " << ObjCon::get_subset_index(handler, v) << "\n");
		if(ObjCon::is_shadow(handler, v)){
			UG_LOG("    is shadow\n");
		}
		else{
			UG_LOG("    is no shadow\n");
		}
	}
}

}//	end of namespace

#endif	// __H__UG__GRID_TRAITS_IMPL__