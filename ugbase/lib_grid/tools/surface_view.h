//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m03 d24

#ifndef __H__UG__SURFACE_VIEW__
#define __H__UG__SURFACE_VIEW__

#include <vector>
#include "lib_grid/lg_base.h"

namespace ug
{

/**
 * The SurfaceView is a special form of a SubsetHandler.
 * It allows registration of GridObservers and notifies them
 * when elements have been created (assigned to a subset != -1)
 * or are to be erased (assigned to subset -1).
 * The SurfaceView passes the associated grid to those methods.
 *
 * \todo: instead of (ab)using GridObservers it would probably make
 *		sense to introduce a SurfaceViewObserver...
 */
class SurfaceView : public SubsetHandler
{
	public:
		SurfaceView();
		SurfaceView(MultiGrid& mg);
		virtual ~SurfaceView();
		
	///	warning: this method is not virtual!
		void assign_grid(MultiGrid& mg);
		
		virtual void assign_subset(VertexBase* elem, int subsetIndex);
		virtual void assign_subset(EdgeBase* elem, int subsetIndex);
		virtual void assign_subset(Face* elem, int subsetIndex);
		virtual void assign_subset(Volume* elem, int subsetIndex);
		
		void register_observer(GridObserver* observer, uint observerType = OT_FULL_OBSERVER);
		void unregister_observer(GridObserver* observer);
		
	///	returns if the element is a shadow, i.e. has children
		template <class TGeomObj>
		inline bool is_shadow(TGeomObj* obj) const	{return m_pMG->has_children(obj);}

	///	returns if the element shadows another
		template <class TGeomObj>
		inline bool shadows(TGeomObj* obj) const
		{
			GeometricObject* parent = m_pMG->get_parent(obj);
			if(parent == NULL) return false;
			else return (get_subset_index(parent) != -1);
		}

	///	returns father of a shadowing element
		template <class TGeomObj>
		inline GeometricObject* get_parent(TGeomObj* obj) const {return m_pMG->get_parent(obj);}

	///	returns child for shadow
		inline VertexBase* get_child(VertexBase* vrt) const {return m_pMG->get_child_vertex<VertexBase>(vrt);}
		
	///	returns the level in grid hierarchy of an element in the surface
		template <class TGeomObj>
		inline int get_level(TGeomObj* obj) const	{return m_pMG->get_level(obj);}

		virtual void grid_to_be_destroyed(Grid* grid);
		
	protected:
		typedef std::vector<GridObserver*>	ObserverContainer;
		
	protected:
		MultiGrid*			m_pMG;
		
	//	observer handling
		ObserverContainer	m_gridObservers;
		ObserverContainer	m_vertexObservers;
		ObserverContainer	m_edgeObservers;
		ObserverContainer	m_faceObservers;
		ObserverContainer	m_volumeObservers;
};

}//	end of namespace

#endif // __H__UG__SURFACE_VIEW__
