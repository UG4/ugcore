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

#ifndef __H__UG__debug_util_impl__
#define __H__UG__debug_util_impl__

#include <fstream>
#include <sstream>
#include "lib_grid/lg_base.h"
#include "lib_grid/algorithms/geom_obj_util/geom_obj_util.h"
#include "common/util/table.h"
#include "lib_grid/tools/periodic_boundary_manager.h"

namespace ug
{

template <typename TElem>
vector3 GetGridObjectCenter(Grid& g, TElem* elem)
{
	if(g.has_vertex_attachment(aPosition)){
		Grid::VertexAttachmentAccessor aaPos(g, aPosition);
		return CalculateCenter(elem, aaPos);
	}
	else if(g.has_vertex_attachment(aPosition2)){
		Grid::VertexAttachmentAccessor aaPos(g, aPosition2);
		vector2 v = CalculateCenter(elem, aaPos);
		return vector3(v.x(), v.y(), 0);
	}
	if(g.has_vertex_attachment(aPosition1)){
		Grid::VertexAttachmentAccessor aaPos(g, aPosition1);
		vector1 v = CalculateCenter(elem, aaPos);
		return vector3(v.x(), 0, 0);
	}

	UG_LOG("GetGridObjectCenter failed! No standard position attachment found.\n");
	return vector3(0, 0, 0);
}


inline vector3 GetGridObjectCenter(Grid& g, GridObject* elem)
{
	switch(elem->base_object_id()){
		case VERTEX:	return GetGridObjectCenter(g, static_cast<Vertex*>(elem));
		case EDGE:		return GetGridObjectCenter(g, static_cast<Edge*>(elem));
		case FACE:		return GetGridObjectCenter(g, static_cast<Face*>(elem));
		case VOLUME:	return GetGridObjectCenter(g, static_cast<Volume*>(elem));
		default:		UG_THROW("Unknown base object type."); 
	}
	return vector3(0, 0, 0);
}


template <typename TElem>
int GetGridObjectIndex(Grid& g, TElem* elem)
{
	using TBase = typename Grid::traits<TElem>::base_object;

	int counter = 0;
	for(typename Grid::traits<TBase>::iterator iter = g.begin<TBase>();
		iter != g.end<TBase>(); ++iter, ++counter)
	{
		if(*iter == elem)
			return counter;
	}
	return -1;
}

template <typename TElem, typename TAValue>
void WriteDebugValuesToFile(const char* filename, Grid& grid,
							TAValue& aVal, bool levelWise)
{
	using TBase = typename Grid::traits<TElem>::base_object;

	std::ofstream out(filename);
	if(!out)
		return;

	Grid::AttachmentAccessor<TBase, TAValue> aaVal(grid, aVal);

	Table<std::stringstream> table(grid.num<TElem>() + 1, 3);
	table(0, 0) << "lvl";	table(0, 1) << "center";	table(0, 2) << "value";

	size_t row = 1;
	if(levelWise){
		GridObjectCollection goc = grid.get_grid_objects();
		for(size_t lvl = 0; lvl < goc.num_levels(); ++lvl){
			for(typename Grid::traits<TElem>::iterator iter = goc.begin<TElem>(lvl);
				iter != goc.end<TElem>(lvl); ++iter, ++row)
			{
				table(row, 0) << lvl;
				table(row, 1) << GetGridObjectCenter(grid, *iter);
				table(row, 2) << aaVal[*iter];
			}
		}
	}
	else if(MultiGrid* pmg = dynamic_cast<MultiGrid*>(&grid)){
		MultiGrid& mg = *pmg;
		for(typename Grid::traits<TElem>::iterator iter = grid.begin<TElem>();
			iter != grid.end<TElem>(); ++iter, ++row)
		{
			table(row, 0) << mg.get_level(*iter);
			table(row, 1) << GetGridObjectCenter(grid, *iter);
			table(row, 2) << aaVal[*iter];
		}
	}
	else{
		for(typename Grid::traits<TElem>::iterator iter = grid.begin<TElem>();
			iter != grid.end<TElem>(); ++iter, ++row)
		{
			table(row, 0) << 0;
			table(row, 1) << GetGridObjectCenter(grid, *iter);
			table(row, 2) << aaVal[*iter];
		}
	}

	out << table;
	out.close();
}

}//	end of namespace

#endif
