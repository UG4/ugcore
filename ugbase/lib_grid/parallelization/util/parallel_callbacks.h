#ifndef __H__UG__parallel_callbacks__
#define __H__UG__parallel_callbacks__

#include "../distributed_grid.h"
#include "common/assert.h"

namespace ug
{

///	Returns true if an element is a regular surface element.
/**	Regular surface elements are elements which lie on the
 * surface and are not a ghost element.
 * (Ghost are vertical masters, which do not lie in any other interface).
 */
class IsRegularSurfaceElem
{
	public:
		IsRegularSurfaceElem(const DistributedGridManager& dgm) :
			m_dgm(dgm), m_mg(dgm.get_assigned_grid())
		{
			UG_ASSERT(m_mg, "A grid has to be assigned to the distributed grid manager.");
		}

		bool operator() (Vertex* v)	{return is_ok(v);}
		bool operator() (Edge* e)	{return is_ok(e);}
		bool operator() (Face* f)		{return is_ok(f);}
		bool operator() (Volume* v)		{return is_ok(v);}

	private:
		template <class TElem>
		inline bool is_ok(TElem* e)
		{
			return !(m_mg->has_children(e) || m_dgm.is_ghost(e));
		}

		const DistributedGridManager& 	m_dgm;
		const MultiGrid*				m_mg;
};

}//	end of namespace

#endif
