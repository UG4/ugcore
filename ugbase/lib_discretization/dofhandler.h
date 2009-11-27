/*
 * dofhandler.h
 *
 *  Created on: 23.04.2009
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__DOFHANDLER__
#define __H__LIBDISCRETIZATION__DOFHANDLER__

#include "common/types.h"
#include "lib_grid/lib_grid.h"

namespace ug {

typedef uint DOFIndex;

////////////////////////////////////////////////////////////////////////
//	DOFHandler
///	takes care of the degrees of freedom of a grid.
/**
 * The DOFHandler manages the degrees of freedom of a grid. He takes care
 * of the creation and destruction of DOFs and for a geometric object
 * he knows the associated degrees of freedom. The class DOFHandler is
 * just a general structure from which specific DOFHandlers can be
 * derived.
 */
class DoFHandler : public GridObserver
{
	public:
		DoFHandler(){};
		DoFHandler(Grid& grid){assert(0);};
		DoFHandler(const DoFHandler& sh){assert(0);};

		virtual DOFIndex get_dofindex_at_vertex(VertexBase* vrt) = 0;
		virtual void reset_dofindeces() = 0;
		virtual uint total_number_of_dof() const = 0;

		virtual ~DoFHandler(){};

	protected:
		Grid* m_pGrid; // associated Grid


};

class NodalDoFHandler : public DoFHandler
{
	public:
		NodalDoFHandler(){};
		NodalDoFHandler(Grid& grid) {m_pGrid = NULL; assign_grid(grid);}
		NodalDoFHandler(const NodalDoFHandler& dofhandler){assert(0);};

		virtual ~NodalDoFHandler()
		{
			m_pGrid->unregister_observer(this);
		};

		void assign_grid(Grid& grid)
		{
			grid.register_observer(this, OT_GRID_OBSERVER | OT_VERTEX_OBSERVER );
		}

		virtual void registered_at_grid(Grid* grid)
		{
			if(m_pGrid)
				m_pGrid->unregister_observer(this);
			m_pGrid = grid;

			m_pGrid->attach_to_vertices(m_aDOFIndex);
			m_aaDOFIndexVRT.access(*m_pGrid, m_aDOFIndex);
		}

		virtual void unregistered_from_grid(Grid* grid)
		{
			m_pGrid->detach_from_vertices(m_aDOFIndex);
			m_aaDOFIndexVRT.invalidate();
		}


		virtual DOFIndex get_dofindex_at_vertex(VertexBase* vrt)
		{
			assert(m_aaDOFIndexVRT.valid());
			return m_aaDOFIndexVRT[vrt];
		}

		virtual void reset_dofindeces()
		{
			assert(m_aaDOFIndexVRT.valid());

			int n = 0;
			VertexIterator iterBegin, iterEnd, iter;

			iterBegin = m_pGrid->begin<Vertex>();
			iterEnd = m_pGrid->end<Vertex>();

			for(iter = iterBegin; iter != iterEnd; iter++)
			{
				Vertex* vrt = *iter;
				m_aaDOFIndexVRT[vrt] = n++;
			}

			std::cout << n << " DoFs assigned to grid" << std::endl;
			m_numberOfDOF = n;
		}

		virtual uint total_number_of_dof() const
		{
			return m_numberOfDOF;
		}

	protected:
		typedef ug::Attachment<DOFIndex> ADOFIndex;

	protected:
		ADOFIndex m_aDOFIndex;
		Grid::VertexAttachmentAccessor<ADOFIndex> m_aaDOFIndexVRT;
		uint m_numberOfDOF;

};

} // end namespace

#endif /* __H__LIBDISCRETIZATION__DOFHANDLER__ */


