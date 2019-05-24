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

#include "shadow_copy_adjuster.h"


namespace ug {


ShadowCopyAdjuster::~ShadowCopyAdjuster()
{}


void ShadowCopyAdjuster:: ref_marks_changed
(
	IRefiner& ref,
	const std::vector<Vertex*>& vrts,
	const std::vector<Edge*>& edges,
	const std::vector<Face*>& faces,
	const std::vector<Volume*>& vols
)
{
	Grid::face_traits::secure_container fl;
	const size_t nv = vols.size();
	for (size_t i = 0; i < nv; ++i)
	{
		Volume* vol = vols[i];

		if (!ref.marked_full(vol))
			continue;

		bool allSidesRefineable = true;
		ref.grid()->associated_elements(fl, vol);
		const size_t flSz = fl.size();
		for (size_t j = 0; j < flSz; ++j)
		{
			Face* f = fl[j];
			RefinementMark curMark = ref.get_mark(f);
			if (!ref.mark(f, RM_FULL))
			{
				ref.mark(f, curMark);
				allSidesRefineable = false;
				break;
			}
			ref.mark(f, curMark);
		}

		if (!allSidesRefineable)
			ref.mark(vol, RM_CLOSURE);
	}

	Grid::edge_traits::secure_container el;
	const size_t nf = faces.size();
	for (size_t i = 0; i < nf; ++i)
	{
		Face* f = faces[i];

		if (!ref.marked_full(f))
			continue;

		bool allSidesRefineable = true;
		ref.grid()->associated_elements(el, f);
		const size_t elSz = el.size();
		for (size_t j = 0; j < elSz; ++j)
		{
			Edge* e = el[j];
			RefinementMark curMark = ref.get_mark(e);
			if (!ref.mark(e, RM_FULL))
			{
				ref.mark(e, curMark);
				allSidesRefineable = false;
				break;
			}
			ref.mark(e, curMark);
		}

		if (!allSidesRefineable)
			ref.mark(f, RM_CLOSURE);
	}
}


}  // namespace ug


