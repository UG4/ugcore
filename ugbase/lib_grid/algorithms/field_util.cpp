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

#include "field_util.h"
#include "lib_grid/algorithms/geom_obj_util/vertex_util.h"

using namespace std;

namespace ug{

////////////////////////////////////////////////////////////////////////////////
void CreateGridFromField(Grid& grid,
						 const Field<number>& field,
						 const vector2& cellSize,
						 const vector2& offset,
						 number noDataValue,
						 Grid::VertexAttachmentAccessor<APosition> aaPos)
{
	Field<Vertex*> vrtField(field.width(), field.height());
	vrtField.fill_all(0);

//	iterate over cells and create triangles or quadrilaterals between
//	adjacent vertices
	const int fieldWidth = (int)field.width();
	const int fieldHeight = (int)field.height();

	for(int irow = 0; irow + 1 < fieldHeight; ++irow){
		for(int icol = 0; icol + 1 < fieldWidth; ++icol){
		//	corner indices in ccw order
		//	index order:
		//	0--3
		//	|  |
		//	1--2
			const int ci[4][2] = {{icol, irow},
								  {icol, irow + 1},
								  {icol + 1, irow + 1},
								  {icol + 1, irow}};

		//	cell corner values in ccw order
			number val[4];
			size_t numVals = 0;
			int firstInvalid = -1;
			for(size_t i = 0; i < 4; ++i){
				val[i] = field.at(ci[i][0], ci[i][1]);
				if(val[i] != noDataValue)
					++numVals;
				else if(firstInvalid == -1)
					firstInvalid = (int)i;
			}

			if(numVals < 3)
				continue;
			
		//	create vertices as necessary and store corner vertices in vrts[]
			Vertex* vrts[4];
			for(size_t i = 0; i < 4; ++i){
				const int ic = ci[i][0];
				const int ir = ci[i][1];
				Vertex* vrt = vrtField.at(ic, ir);
				if((val[i] != noDataValue) && (!vrt)){
					vrt = *grid.create<RegularVertex>();
					aaPos[vrt] =
						vector3(offset.x() + cellSize.x() * (number)ic,
								offset.y() + cellSize.y() * (number)ir,
								field.at(ic, ir));
					vrtField.at(ic, ir) = vrt;
				}
				vrts[i] = vrt;
			}

			if(numVals == 4){
			//	create a quad
				grid.create<Quadrilateral>(
					QuadrilateralDescriptor(vrts[3], vrts[2], vrts[1], vrts[0]));
			}
			else if(numVals == 3){
			//	create a tri
				grid.create<Triangle>(
					TriangleDescriptor(vrts[(firstInvalid + 3) % 4],
									   vrts[(firstInvalid + 2) % 4],
									   vrts[(firstInvalid + 1) % 4]));
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////////////
void CreateGridFromFieldBoundary(Grid& grid,
								 const Field<number>& field,
								 const vector2& cellSize,
								 const vector2& offset,
								 number noDataValue,
								 Grid::VertexAttachmentAccessor<APosition> aaPos)
{
	Field<Vertex*> vrtField(field.width(), field.height());
	vrtField.fill_all(0);

//	iterate over cells and create triangles or quadrilaterals between
//	adjacent vertices
	const int fieldWidth = (int)field.width();
	const int fieldHeight = (int)field.height();

	for(int irow = 0; irow + 1 < fieldHeight; ++irow){
		for(int icol = 0; icol + 1 < fieldWidth; ++icol){
		//	corner indices in ccw order
		//	index order:
		//	0--3
		//	|  |
		//	1--2
			const int ci[4][2] = {{icol, irow},
								  {icol, irow + 1},
								  {icol + 1, irow + 1},
								  {icol + 1, irow}};

		//	cell corner values in ccw order
			number val[4];
			size_t numVals = 0;
			int firstInvalid = -1;
			for(size_t i = 0; i < 4; ++i){
				val[i] = field.at(ci[i][0], ci[i][1]);
				if(val[i] != noDataValue)
					++numVals;
				else if(firstInvalid == -1)
					firstInvalid = (int)i;
			}

			if(numVals < 3)
				continue;
			
		//	create vertices as necessary and store corner vertices in vrts[]
			Vertex* vrts[4];

			if(numVals == 3){
			//	create a diagonal which doesn't include the firstInvalid corner
				int ind[2]  = {(firstInvalid + 1) % 4, (firstInvalid + 3) % 4};

				for(int i = 0; i < 2; ++i){
					int ic = ci[ind[i]][0];
					int ir = ci[ind[i]][1];
					Vertex* vrt = vrtField.at(ic, ir);
					if((val[ind[i]] != noDataValue) && (!vrt)){
						vrt = *grid.create<RegularVertex>();
						aaPos[vrt] =
							vector3(offset.x() + cellSize.x() * (number)ic,
									offset.y() + cellSize.y() * (number)ir,
									val[ind[i]]);
						vrtField.at(ic, ir) = vrt;
					}
					vrts[ind[i]] = vrt;
				}

				grid.create<RegularEdge>(EdgeDescriptor(vrts[ind[0]], vrts[ind[1]]));
			}

		//	now check for each inner side whether it's a neighbor to an empty cell
		constexpr int colAdd[4] = {-1, 0, 1, 0};
		constexpr int rowAdd[4] = {0, 1, 0, -1};
			for(int iside = 0; iside < 4; ++iside){
				const int ind[2] = {iside, (iside + 1) % 4};
				if(val[ind[0]] != noDataValue && val[ind[1]] != noDataValue){
				//	the edge is either inner or boundary
					bool gotOne = false;
					for(int i = 0; i < 2; ++i){
						const int ic = ci[ind[i]][0] + colAdd[iside];
						const int ir = ci[ind[i]][1] + rowAdd[iside];
						if(ic >= 0 && ir >= 0 && ic < fieldWidth && ir < fieldHeight){
							if(field.at(ic, ir) != noDataValue){
								gotOne = true;
								break;
							}
						}
					}

					if(!gotOne){
					//	the current side has to be represented by an edge
						for(int i = 0; i < 2; ++i){
							const int ic = ci[ind[i]][0];
							const int ir = ci[ind[i]][1];
							Vertex* vrt = vrtField.at(ic, ir);
							if((val[ind[i]] != noDataValue) && (!vrt)){
								vrt = *grid.create<RegularVertex>();
								aaPos[vrt] =
									vector3(offset.x() + cellSize.x() * (number)ic,
											offset.y() + cellSize.y() * (number)ir,
											val[ind[i]]);
								vrtField.at(ic, ir) = vrt;
							}
							vrts[ind[i]] = vrt;
						}
						grid.create<RegularEdge>(EdgeDescriptor(vrts[ind[0]], vrts[ind[1]]));
					}
				}
			}
		}
	}
}

}//	end of namespace
