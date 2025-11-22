/*
 * Copyright (c) 2015:  G-CSC, Goethe University Frankfurt
 * Author: Markus Breit
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

#ifndef LIB_GRID__COPY_ATTACHMENT_HANDLER_H_
#define LIB_GRID__COPY_ATTACHMENT_HANDLER_H_

#include "common/util/smart_pointer.h"  // for SmartPtr
#include "lib_grid/grid/grid_observer.h"  // for GridObserver
#include "lib_grid/multi_grid.h"  // for MultiGrid

#ifdef UG_PARALLEL
    #include "lib_grid/parallelization/util/compol_copy_attachment.h" // ComPol_CopyAttachment
    #include "lib_grid/parallelization/distributed_grid.h"        // DistributedGridManager
#endif

namespace ug{


/**
 * @brief handler for attachments in a multi-grid
 *
 * When using attachments in a multi-grid hierarchy, one
 * EITHER needs to create the attachment after the geometry is refined (and distributed)
 * up to the state in which the attachments are used - in this case, every geometric
 * object which the data is attached to will have its attachment and there is no problem
 * OR one creates the attachment before the geometry is refined (and distributed) - in
 * this case, if the attachment is required on elements that are created after creation
 * of the attachment, one needs to take care to supply those elements with attachment
 * values. This class will do exactly this, but only to a very limited extent:
 *
 * When this class is assigned a multi-grid, it will propagate attachment values from the
 * base level upwards to the children <b>of the same type<\b> by copying.
 * Vertical communication is performed on each level to make sure that all those children
 * will get their values.
 *
 * Moreover, making use of the GridObserver interface, this class will copy attachment
 * values to each child element (of equal type as their parent) created after assignment
 * of the grid.
 *
 *
 * This will not, in general, guarantee that every element of the required type will
 * get an attachment value. It will however suffice if the element type the attachment is
 * attached to is full-dimensional (in its domain); it will obviously also suffice if the
 * attachment values are only required on direct children (of equal type) of coarse grid
 * elements, e.g. on d-dimensional elements on a d-dimensional manifold.
 *
 * Example:
 *
 * If, in the base level, attachments (a1, a2) are attached to their resp. nodes like this:
 *
 *  a1 --------- a2    <-- lvl 0
 *
 * then, on the next level, the node marked by "o" will not have an attachment value
 * propagated from below (since it does not have a parent of the same type):
 *
 *  a1 --- o --- a2    <-- lvl 1
 *
 *  a1 --------- a2    <-- lvl 0
 *
 * You may, however, choose to derive another handler class from this implementation
 * and treat the cases where an upper-level element is child of an element of different
 * type there by overloading copy_from_other_elem_type().
 * In the same manner, you might use derivation and overloading (of the copy() method)
 * to perform a more specialized way of copying.
 *
 * @note Implementation relies on a template trick: One cannot - as decreed by the
 * standard - specialize template methods of a not fully specialized template class;
 * it is possible, however, to @e partially specialize nested template classes of a
 * template class. This is what is used here.
 *
 */
template <typename TElem, typename TAttachment>
class CopyAttachmentHandler : public GridObserver
{
	public:
		/// constructor
		CopyAttachmentHandler()
		: m_bEnableVertComm(true) {}

		/// destructor
		~CopyAttachmentHandler() override {
			if (m_spMG.valid())
				m_spMG->unregister_observer(this);

			m_aa.invalidate();
		};

		void set_attachment(const TAttachment& attch)
		{
			m_a = attch;

			// if we already have a grid, then this must be a reset
			if (m_spMG.valid())
			{
				m_aa.invalidate();
				m_aa.access(*m_spMG, m_a);
			}
		}

		void set_grid(SmartPtr<MultiGrid> mg)
		{
			// do nothing if given grid is already set
			if (m_spMG == mg) return;

			// if another grid is already given, unregister from old grid and invalidate accessor
			if (m_spMG.valid())
			{
				m_spMG->unregister_observer(this);
				m_aa.invalidate();
			}

			// set new grid
			m_spMG = mg;
			UG_COND_THROW(!m_spMG.valid(), "No valid multigrid given!");

			// register as appropriate observer
			register_as_observer<TElem, void>(m_spMG, this);

			// check that attachment is attached to grid
			UG_COND_THROW(!m_spMG->has_attachment<TElem>(m_a),
						  "Given grid does not have given attachment attached to it.");

			// init accessor
			m_aa.access(*m_spMG, m_a);

			// copy to all levels that are already present
			propagate_to_levels();
		}

		void enable_vertical_communication(bool b)
		{
			m_bEnableVertComm = b;
		}


		// GridObserver implementations
		void vertex_created
		(
			Grid* grid,
			Vertex* vrt,
			GridObject* pParent = nullptr,
			bool replacesParent = false
		) override {
			propagate<Vertex, void>(vrt, pParent, this);
		}

		void edge_created
		(
			Grid* grid,
			Edge* e,
			GridObject* pParent = nullptr,
			bool replacesParent = false
		) override {
			propagate<Edge, void>(e, pParent, this);
		}

		void face_created
		(
			Grid* grid,
			Face* f,
			GridObject* pParent = nullptr,
			bool replacesParent = false
		) override {
			propagate<Face, void>(f, pParent, this);
		}

		void volume_created
		(
			Grid* grid,
			Volume* vol,
			GridObject* pParent = nullptr,
			bool replacesParent = false
		) override {
			propagate<Volume, void>(vol, pParent, this);
		}

	protected:

		// this template will be used if created is called with anything BUT TElem (in that case: do nothing)
		template <typename TCreatedElem, typename Dummy>
		struct propagate
		{
			propagate(TCreatedElem* elem, GridObject* pParent, CopyAttachmentHandler* cah) {};
		};

		// this template will be used if created is called with an elem of type TElem
		template <typename Dummy>
		struct propagate<TElem, Dummy>
		{
			propagate(TElem* elem, GridObject* pParent, CopyAttachmentHandler* cah)
			{
				// check that parent is given
				if (!pParent) return;

				// if parent is of different elem type than child
				TElem* par = dynamic_cast<TElem*>(pParent);
				if (!par)
				{
					cah->copy_from_other_elem_type(pParent, elem);
					return;
				}

				// copy attachment value from parent to child
				cah->copy(par, elem);
			}
		};

		friend struct propagate<TElem, void>;

		template <typename TObserverElem, typename Dummy>
		struct register_as_observer
		{
			register_as_observer(SmartPtr<MultiGrid> mg, CopyAttachmentHandler<TObserverElem, TAttachment>* cah) {};
		};

		template <typename Dummy>
		struct register_as_observer<Vertex, Dummy>
		{
			register_as_observer(SmartPtr<MultiGrid> mg, CopyAttachmentHandler<Vertex, TAttachment>* cah)
			{
				mg->register_observer(cah, OT_VERTEX_OBSERVER);
			}
		};
		template <typename Dummy>
		struct register_as_observer<Edge, Dummy>
		{
			register_as_observer(SmartPtr<MultiGrid> mg, CopyAttachmentHandler<Edge, TAttachment>* cah)
			{
				mg->register_observer(cah, OT_EDGE_OBSERVER);
			}
		};
		template <typename Dummy>
		struct register_as_observer<Face, Dummy>
		{
			register_as_observer(SmartPtr<MultiGrid> mg, CopyAttachmentHandler<Face, TAttachment>* cah)
			{
				mg->register_observer(this, OT_FACE_OBSERVER);
			}
		};
		template <typename Dummy>
		struct register_as_observer<Volume, Dummy>
		{
			register_as_observer(SmartPtr<MultiGrid> mg, CopyAttachmentHandler<Volume, TAttachment>* cah)
			{
				mg->register_observer(this, OT_VOLUME_OBSERVER);
			}
		};


		void propagate_to_levels()
		{
			// ensure a final vertical communication on the top level
			const size_t nProp = m_bEnableVertComm ? m_spMG->num_levels() + 1 : m_spMG->num_levels();
			for (size_t i = 1; i < nProp; ++i)
				propagate_to_level(i, m_bEnableVertComm);
		}

		void propagate_to_level(size_t lvl, bool vertComm = true)
		{
#ifdef UG_PARALLEL
			// copy from vMasters to vSlaves
			if (pcl::NumProcs() > 1 && vertComm)
			{
				using layout_type = typename GridLayoutMap::Types<TElem>::Layout;
				DistributedGridManager& dgm = *m_spMG->distributed_grid_manager();
				GridLayoutMap& glm = dgm.grid_layout_map();
				pcl::InterfaceCommunicator<layout_type> icom;

				ComPol_CopyAttachment<layout_type, TAttachment> compolCopy(*m_spMG, m_a);
				icom.exchange_data(glm, INT_V_MASTER, INT_V_SLAVE, compolCopy);
				icom.communicate();
			}
#endif
			// iterate over all TElems of level
			using iter_type = typename geometry_traits<TElem>::const_iterator;
			iter_type iter = m_spMG->begin<TElem>(lvl);
			iter_type iter_end = m_spMG->end<TElem>(lvl);
			for (; iter != iter_end; ++iter)
			{
				TElem* child = *iter;
				GridObject* parent = m_spMG->get_parent(*iter);

				propagate<TElem, void>(child, parent, this);
			}
		}

		virtual void copy(TElem* parent, TElem* child)
		{
			m_aa[child] = m_aa[parent];
		}

		virtual void copy_from_other_elem_type(GridObject* parent, TElem* child) {};

	protected:
		/// multigrid
		SmartPtr<MultiGrid> m_spMG;

		bool m_bEnableVertComm;

		/// attachment to be handled
		TAttachment m_a;

		/// accessor for handled attachment
		Grid::AttachmentAccessor<TElem, TAttachment> m_aa;

};

} // end namespace ug

#endif