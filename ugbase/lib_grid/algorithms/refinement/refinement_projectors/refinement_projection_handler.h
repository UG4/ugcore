/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__REFINEMENT_PROJECTION_HANDLER__
#define __H__UG__REFINEMENT_PROJECTION_HANDLER__

#include "../refinement_callbacks.h"
#include "lib_grid/tools/subset_handler_interface.h"

namespace ug{

///	Register different refinement callbacks for different subsets
/**	Using set_callback you may specify a callback which will be executed,
 * if the parent element is contained in the specified subset.
 * You may also specify a default callback which is executed if no parent is
 * present, if the parent isn't assigned to a subset or if no callback has
 * been assigned to the parent's subset. By default the default callback is
 * RefinementCallbackLinear.
 *
 * Use IRefiner::set_refinement_callback to set an instance of this class as
 * the refinement callback of a refiner.
 */
template <class TAPosition>
class RefinementProjectionHandler : public IRefinementCallback {
	public:
		RefinementProjectionHandler(SmartPtr<ISubsetHandler> sh, TAPosition aPos);

		virtual ~RefinementProjectionHandler();

		void set_default_callback(SmartPtr<IRefinementCallback> callback);
		void set_callback(int subsetIndex, SmartPtr<IRefinementCallback> callback);
		void set_callback(std::string subsetName, SmartPtr<IRefinementCallback> callback);

	////////////////////////////////////////
	//	IMPLEMENTATION OF IRefinementCallback
	///	called when a new vertex was created from an old vertex.
		virtual void new_vertex(Vertex* vrt, Vertex* parent);
	///	called when a new vertex was created from an old edge.
		virtual void new_vertex(Vertex* vrt, Edge* parent);
	///	called when a new vertex was created from an old face.
		virtual void new_vertex(Vertex* vrt, Face* parent);
	///	called when a new vertex was created from an old volume.
		virtual void new_vertex(Vertex* vrt, Volume* parent);

	///	callback for vertices in flat grids.
		virtual void flat_grid_vertex_encountered(Vertex* vrt);

	///	returns the position of the given vertex.
		virtual int current_pos(number* coordsOut, Vertex* vrt, int maxCoords);

	private:
		template <class TParent>
		void handle_new_vertex(Vertex* vrt, TParent* parent);

		SmartPtr<ISubsetHandler>					m_sh;
		Grid::VertexAttachmentAccessor<TAPosition>	m_aaPos;
		std::vector<SmartPtr<IRefinementCallback> >	m_callbacks;
		SmartPtr<IRefinementCallback>				m_defaultCallback;
};

}	//	end of namespace

////////////////////
#include "refinement_projection_handler_impl.hpp"

#endif /* REFINEMENT_PROJECTION_HANDLER_H_ */
