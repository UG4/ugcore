/*
 * Copyright (c) 2009-2019:  G-CSC, Goethe University Frankfurt
 * Author: Stephan Grein
 * Creation date: 2019-07-05
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

#ifndef __H__LIB_GRID__FILE_IO_IMPL__
#define __H__LIB_GRID__FILE_IO_IMPL__

#include "file_io.h"
#include "lib_grid/algorithms/attachment_util.h"
#include "lib_grid/algorithms/subset_util.h"

namespace ug {
	template<class TElem>
	void CopyGridElements(Grid& srcGrid, Grid& destGrid,
					      ISubsetHandler& srcSH, ISubsetHandler& destSH,
						  Attachment<Vertex*>& aNewVrt)
	{
		Grid::VertexAttachmentAccessor<Attachment<Vertex*> > aaNewVrt(srcGrid, aNewVrt);
		GridObjectCollection goc = srcGrid.get_grid_objects();
		CustomVertexGroup vrts;

		using iter_t = typename Grid::traits<TElem>::iterator;

		for (iter_t eIter = goc.begin<TElem>(); eIter != goc.end<TElem>(); ++eIter)
		{
			TElem* e = *eIter;
			vrts.resize(e->num_vertices());

			for (size_t iv = 0; iv < e->num_vertices(); ++iv)
			{
				vrts.set_vertex(iv, aaNewVrt[e->vertex(iv)]);
			}

			TElem* ne = *destGrid.create_by_cloning(e, vrts);
			destSH.assign_subset(ne, srcSH.get_subset_index(e));
		}
	}
}

#endif /* __H__LIB_GRID__FILE_IO_IMPL__ */
