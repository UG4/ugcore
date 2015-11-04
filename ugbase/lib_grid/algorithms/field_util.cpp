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
			const int colAdd[4] = {-1, 0, 1, 0};
			const int rowAdd[4] = {0, 1, 0, -1};
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
