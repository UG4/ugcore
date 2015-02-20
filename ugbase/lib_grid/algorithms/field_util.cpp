// created by Sebastian Reiter
// s.b.reiter@gmail.com

#include "field_util.h"

namespace ug{
	
void CreateGridFromField(Grid& grid,
						 const Field<number>& field,
						 const vector2& cellSize,
						 const vector2& offset,
						 number noDataValue,
						 Grid::VertexAttachmentAccessor<APosition> aaPos)
{
	Field<Vertex*> vrtField(field.width(), field.height());

//	iterate over points. We create a vertex for each point which will be
//	part of a triangle or quadrilateral
	for(size_t irow = 0; irow < field.height(); ++irow){
		for(size_t icol = 0; icol < field.width(); ++icol){
			if(field.at(icol, irow) != noDataValue){
				bool left = false;
				bool right = false;
				bool up = false;
				bool down = false;

				if(icol > 0)
					left = (field.at(icol - 1, irow) != noDataValue);
				if(icol + 1 < field.width())
					right = (field.at(icol + 1, irow) != noDataValue);
				if(irow > 0)
					up = (field.at(icol, irow - 1) != noDataValue);
				if(irow + 1 < field.height())
					down = (field.at(icol, irow + 1) != noDataValue);

				if((left && up) || (up && right) || (right && down) || (down && left)){
					vrtField.at(icol, irow) = *grid.create<RegularVertex>();
					aaPos[vrtField.at(icol, irow)] =
						vector3(offset.x() + cellSize.x() * (number)icol,
								offset.y() + cellSize.y() * (number)irow,
								field.at(icol, irow));
				}
				else
					vrtField.at(icol, irow) = NULL;
			}
			else
				vrtField.at(icol, irow) = NULL;
		}
	}

//	iterate over cells and create triangles or quadrilaterals between
//	adjacent vertices
	for(size_t irow = 0; irow + 1 < field.height(); ++irow){
		for(size_t icol = 0; icol + 1 < field.width(); ++icol){
		//	cell vertices in ccw order
			Vertex* vrts[4];
			vrts[0] = vrtField.at(icol, irow);
			vrts[1] = vrtField.at(icol, irow + 1);
			vrts[2] = vrtField.at(icol + 1, irow + 1);
			vrts[3] = vrtField.at(icol + 1, irow);

			size_t numVrts = 0;
			int firstZero = -1;
			for(size_t i = 0; i < 4; ++i){
				if(vrts[i])
					++numVrts;
				else if(firstZero == -1)
					firstZero = (int)i;
			}

			if(numVrts == 4){
			//	create a quad
				grid.create<Quadrilateral>(
					QuadrilateralDescriptor(vrts[0], vrts[1], vrts[2], vrts[3]));
			}
			else if(numVrts == 3){
			//	create a tri
				grid.create<Triangle>(
					TriangleDescriptor(vrts[(firstZero + 1) % 4],
									   vrts[(firstZero + 2) % 4],
									   vrts[(firstZero + 3) % 4]));
			}
		}
	}
}

}//	end of namespace
