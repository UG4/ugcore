/*
 * Copyright (c) 2017:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG_grid_util
#define __H__UG_grid_util

#include "../grid/grid.h"
#include "../tools/selector_grid.h"

namespace ug{

/** Make sure that aNewVrt is attached to srcSel.grid() and contains a
 * pointer to a valid vertex in destGrid for each selected vertex in srcGrid.
 * Also make sure that all vertices belonging to a selected element have been
 * selected, too.*/
template <class TElem>
inline void CopySelectedElements(Selector& srcSel, Grid& destGrid, AVertex aNewVrt)
{
	Grid& srcGrid						= *srcSel.grid();

	Grid::VertexAttachmentAccessor<AVertex> aaNewVrt(srcGrid, aNewVrt);

	CustomVertexGroup vrts;
	typedef typename Grid::traits<TElem>::iterator iter_t;

	for(iter_t eiter = srcSel.begin<TElem>();
		eiter != srcSel.end<TElem>(); ++eiter)
	{
		TElem* e = *eiter;
		vrts.resize(e->num_vertices());
		for(size_t iv = 0; iv < e->num_vertices(); ++iv)
			vrts.set_vertex(iv, aaNewVrt[e->vertex(iv)]);

		destGrid.create_by_cloning(e, vrts);
	}
}

template <class TAAPosSrc, class TAAPosDest>
inline void CopySelection(Selector& srcSel, Grid& destGrid,
						  TAAPosSrc aaPosSrc, TAAPosDest aaPosDest)
{
	UG_COND_THROW(!srcSel.grid(), "The specified selector has to operate on a grid.");

	Grid& srcGrid						= *srcSel.grid();

	AVertex aNewVrt;
	srcGrid.attach_to_vertices(aNewVrt);
	Grid::VertexAttachmentAccessor<AVertex> aaNewVrt(srcGrid, aNewVrt);

	SelectAssociatedGridObjects(srcSel);

//	create new vertices in destGrid
	for(VertexIterator viter = srcSel.begin<Vertex>();
		viter != srcSel.end<Vertex>(); ++viter)
	{
		Vertex* v = *viter;
		Vertex* nv = *destGrid.create_by_cloning(v);
		aaNewVrt[v] = nv;
		aaPosDest[nv] = aaPosSrc[v];
	}

	CopySelectedElements<Edge>(srcSel, destGrid, aNewVrt);
	CopySelectedElements<Face>(srcSel, destGrid, aNewVrt);
	CopySelectedElements<Volume>(srcSel, destGrid, aNewVrt);

	srcGrid.detach_from_vertices(aNewVrt);
}

}//	end of namespace

#endif	//__H__UG_grid_util
