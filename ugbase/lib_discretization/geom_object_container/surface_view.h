//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m03 d24

#ifndef __H__UG__SURFACE_VIEW__
#define __H__UG__SURFACE_VIEW__

#include <vector>
#include "lib_grid/lib_grid.h"

namespace ug
{

/**
 * The SurfaceView is a special form of a SubsetHandler.
 * It allows registration of GridObservers and notifies them
 * when elements have been created (assigned to a subset != -1)
 * or are to be erased (assigned to subset -1).
 * The SurfaceView passes the associated grid to those methods.
 *
 * Please note that registered_at_grid will not be called.
 * unregistered_from_grid will be called with a NULL pointer.
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

		virtual void registered_at_grid(Grid* grid);

		template <class TGeomObj>
		inline bool is_shadow(TGeomObj* obj)	{return m_pMG->has_children(obj);}

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
