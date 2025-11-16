/*
 * Copyright (c) 2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG_horizontal_anisotropy_adjuster
#define __H__UG_horizontal_anisotropy_adjuster

#include "../ref_mark_adjuster_interface.h"

namespace ug{

///	Selects additional edges to preserve layers in a grid with horizontal anisotropies
/** Selects additional edges to preserve layers in a grid with horizontal anisotropies
 * and builds a closure so that no hanging nodes are generated.
 * 
 * \note	This adjuster regards the grid as a serial grid. If the grid represents a part
 * 			of a distributed grid, then the additional use of a parallel adjuster is required.
 */
template <class TAPos>
class HorizontalAnisotropyAdjuster : public IRefMarkAdjuster
{
	public:
		using position_attachment_t = TAPos;
		using position_accessor_t = Grid::VertexAttachmentAccessor<TAPos>;

		HorizontalAnisotropyAdjuster(TAPos aPos) : m_aPos(aPos) {}

		virtual ~HorizontalAnisotropyAdjuster()	{}

		virtual void ref_marks_changed(IRefiner& ref,
										const std::vector<Vertex*>& vrts,
										const std::vector<Edge*>& edges,
										const std::vector<Face*>& faces,
										const std::vector<Volume*>& vols);
	private:
		TAPos	m_aPos;
};

}//	end of namespace

#endif	//__H__UG_horizontal_anisotropy_adjuster
