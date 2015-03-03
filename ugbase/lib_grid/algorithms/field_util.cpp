// created by Sebastian Reiter
// s.b.reiter@gmail.com

#include "field_util.h"
#include "lib_grid/algorithms/geom_obj_util/vertex_util.h"

namespace ug{

using namespace std;

////////////////////////////////////////////////////////////////////////////////
Heightfield::Heightfield() :
	m_cellSize(1, 1),
	m_offset(0, 0),
	m_noDataValue(1e30)
{
}

number Heightfield::
interpolate(number x, number y) const
{
	pair<int, int> ind = coordinate_to_index(x, y);
	if( ind.first >= 0 && ind.first < (int)m_field.width() &&
		ind.second >= 0 && ind.second < (int)m_field.height())
	{
		return m_field.at(ind.first, ind.second);
	}
	return m_noDataValue;
}

std::pair<int, int> Heightfield::
coordinate_to_index(number x, number y) const
{
//	roundOffset 0: value is constant across each cell
//	roundOffset 0.5: value is constant around each node
	const number roundOffset = 0.5;

	pair<int, int> c;
	if(m_cellSize.x() != 0)
		c.first = (int)(roundOffset + (x - m_offset.x()) / m_cellSize.x());
	else
		c.first = 0;
	
	if(m_cellSize.y() != 0)
		c.second = (int)(roundOffset + (y - m_offset.y()) / m_cellSize.y());
	else
		c.second = 0;
	return c;
}

vector2 Heightfield::index_to_coordinate(int ix, int iy) const
{
	return vector2(	m_offset.x() + (number)ix * m_cellSize.x(),
					m_offset.y() + (number)iy * m_cellSize.y());
}

vector2 Heightfield::extent() const
{
	return vector2(m_cellSize.x() * (number)m_field.width(),
				   m_cellSize.y() * (number)m_field.height());
}

////////////////////////////////////////////////////////////////////////////////
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

////////////////////////////////////////////////////////////////////////////////
void CreateGridFromField(Grid& grid,
						 const Heightfield& hfield,
						 Grid::VertexAttachmentAccessor<APosition> aaPos)
{
	CreateGridFromField(grid, hfield.field(), hfield.cell_size(), hfield.offset(),
						hfield.no_data_value(), aaPos);
}


////////////////////////////////////////////////////////////////////////////////
void CreateGridFromFieldBoundary(Grid& grid,
								 const Field<number>& field,
								 const vector2& cellSize,
								 const vector2& offset,
								 number noDataValue,
								 Grid::VertexAttachmentAccessor<APosition> aaPos)
{
//	iterate over cells and create edges between neighbored boundary vertices
//	we'll create duplicate vertices and remove them later on
	VertexIterator firstVrtIter;
	bool gotFirstVrt = false;

	for(int irow = -1; irow < (int)field.height(); ++irow){
		for(int icol = -1; icol < (int)field.width(); ++icol){
			int corner[4][2] = {{icol, irow},
								{icol, irow + 1},
								{icol + 1, irow + 1},
								{icol + 1, irow}};
		//	values of cell corners
			number values[4];
			for(int i = 0; i < 4; ++i){
				if(	(corner[i][0] >= 0) && (corner[i][0] < (int)field.width()) &&
					(corner[i][1] >= 0) && (corner[i][1] < (int)field.height()))
				{
					values[i] = field.at(corner[i][0], corner[i][1]);
				}
				else
					values[i] = noDataValue;
			}

			int validCorners[4];
			size_t numValid = 0;
			int firstValid = -1;
			int firstInvalid = -1;
			for(size_t i = 0; i < 4; ++i){
				if(values[i] != noDataValue){
					validCorners[numValid] = i;
					++numValid;
					if(firstValid == -1)
						firstValid = (int)i;
				}
				else if(firstInvalid == -1)
					firstInvalid = (int)i;
			}

			bool createEdge = false;
			int ci[2];
			if(numValid == 3){
			//	we'll create an edge which separates the invalid corner
				ci[0] = (firstInvalid + 1) % 4;
				ci[1] = (firstInvalid + 3) % 4;
				createEdge = true;
			}
			else if((numValid == 2) &&
				((validCorners[1] == (validCorners[0] + 1)) ||
				 (validCorners[0] == (validCorners[1] + 1) % 4)))
			{
			//	we only create an edge if the two valid corners are direct neighbors
				ci[0] = validCorners[0];
				ci[1] = validCorners[1];
				createEdge = true;
			}

			if(createEdge){
				Vertex* vrts[2];
				for(int i = 0; i < 2; ++i){
					VertexIterator vi = grid.create<RegularVertex>();
					aaPos[*vi] = vector3(offset.x() + cellSize.x() * (number)corner[ci[i]][0],
										 offset.y() + cellSize.y() * (number)corner[ci[i]][1],
										 field.at(corner[ci[i]][0], corner[ci[i]][1]));
					vrts[i] = *vi;
					if(!gotFirstVrt){
						gotFirstVrt = true;
						firstVrtIter = vi;
					}
				}

				grid.create<RegularEdge>(EdgeDescriptor(vrts[0], vrts[1]));
			}
		}
	}

	if(gotFirstVrt){
	//	we have to remove double vertices
		RemoveDoubles<3>(grid, firstVrtIter, grid.end<Vertex>(), aaPos, SMALL);
	}
}


void CreateGridFromFieldBoundary(Grid& grid,
								 const Heightfield& hfield,
								 Grid::VertexAttachmentAccessor<APosition> aaPos)
{
	CreateGridFromFieldBoundary(grid, hfield.field(), hfield.cell_size(), hfield.offset(),
								hfield.no_data_value(), aaPos);
}

}//	end of namespace
