/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
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

#ifndef __H__UG__refmark_adjuster_interface__
#define __H__UG__refmark_adjuster_interface__

#include "refiner_interface.h"
#include "common/util/smart_pointer.h"

namespace ug{

///	\addtogroup lib_grid_algorithms_refinement
///	@{

/**	Used to steer the mark-adjustment process in classes derived from IRefiner.
 * Instances can be registered at IRefiner and are then used to mark new elements
 * based on already marked elements. Different specializations can thereby implement
 * different rules (e.g. for parallel refinement or anisotropic refinement).
 * The different specialization should be constructed orthogonal to each other,
 * so that they can be chained.
 *
 * The callback-methods 'ref_marks_changed' and 'coarsen_marks_changed' are called
 * with lists of elements which have been marked since the last adjustment-round.
 * The callbacks can then use IRefiner::mark to adjust marks as required.
 * If a callback changes a ref-mark on an element,
 * this will automatically trigger a new round of adjustments
 * (then only considering elements which have been marked in the current round).
 *
 * Adjustment normally stops as soon as no adjuster marked any new elements.
 */
class IRefMarkAdjuster
{
	public:
		IRefMarkAdjuster() :
			m_enabled(true),
			m_nodeDependencyOrder1(true)
		{}

		virtual ~IRefMarkAdjuster()	= default;

		virtual void ref_marks_changed(IRefiner& ref,
										const std::vector<Vertex*>& vrts,
										const std::vector<Edge*>& edges,
										const std::vector<Face*>& faces,
										const std::vector<Volume*>& vols)
		{}

		virtual void coarsen_marks_changed(IRefiner& ref,
										const std::vector<Vertex*>& vrts,
										const std::vector<Edge*>& edges,
										const std::vector<Face*>& faces,
										const std::vector<Volume*>& vols)
		{}

		virtual void enable(bool enable)	{m_enabled = enable;}
		virtual bool enabled() const		{return m_enabled;}

	///	enables or disables node-dependency-order-1.
	/**	\{
	 * If enabled, hanging nodes may only depend on non-hanging nodes.
	 * An edge containing a hanging node thus will not have a hanging-node
	 * as a corner vertex.
	 *
	 * Enabled by default.*/
		void enable_node_dependency_order_1(bool bEnable)	{m_nodeDependencyOrder1 = bEnable;}
		bool node_dependency_order_1_enabled()				{return m_nodeDependencyOrder1;}
	/**	\} */

	private:
		bool m_enabled;
		bool m_nodeDependencyOrder1;
};

using SPIRefMarkAdjuster = SmartPtr<IRefMarkAdjuster>;

/// @}

}// end of namespace

#endif
