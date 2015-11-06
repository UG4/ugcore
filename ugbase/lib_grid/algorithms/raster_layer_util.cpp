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

#include "raster_layer_util.h"
#include "field_util.h"
#include "lib_grid/file_io/file_io_asc.h"

using namespace std;

namespace ug{

void RasterLayers::
load_from_files(const std::vector<LayerDesc>& layerDescs)
{
	resize(layerDescs.size());
	for(size_t i = 0; i < layerDescs.size(); ++i){
		LoadHeightfieldFromASC(heightfield(i), layerDescs[i].filename().c_str());
		set_min_height(i, layerDescs[i].min_height());
	}
}

void RasterLayers::
load_from_files(const std::vector<SPLayerDesc>& layerDescs)
{
	vector<LayerDesc>	descs;
	descs.reserve(layerDescs.size());
	for_each_in_vec(const SPLayerDesc& d, layerDescs){
		descs.push_back(*d);	
	}end_for;
	load_from_files(descs);
}

void RasterLayers::
load_from_files(const std::vector<std::string>& filenames, number minHeight)
{
	vector<LayerDesc>	descs;
	descs.reserve(filenames.size());
	for_each_in_vec(const std::string& fname, filenames){
		descs.push_back(LayerDesc(fname, minHeight));	
	}end_for;
	load_from_files(descs);
}

void RasterLayers::
resize(size_t newSize)
{
	size_t oldSize = size();
	m_layers.resize(newSize);
	for(size_t i = oldSize; i < newSize; ++i){
		m_layers[i] = make_sp(new layer_t);
	}
}

void RasterLayers::
invalidate_flat_cells()
{
//	since the top-layer is considered to be the terrain surface, we'll start
//	one layer below the top.
	if(size() <= 1)
		return;
	for(int lvl = (int)size() - 2; lvl >= 0; --lvl){
		Heightfield& curHF = heightfield(lvl);
		Field<number>& innerField = curHF.field();
		number minHeight = min_height(lvl);

	//	iterate over the cells of each heightfield and invalidate values
		for(int iy = 0; iy < (int)innerField.height(); ++iy){
			for(int ix = 0; ix < (int)innerField.width(); ++ix){
				number val = innerField.at(ix, iy);
				if(val != curHF.no_data_value()){
					bool gotDataVal = false;
					vector2 c = curHF.index_to_coordinate(ix, iy);
				//	compare the value against values from higher layers
					for(int ulvl = lvl + 1; ulvl < (int)size(); ++ulvl){
						Heightfield& uHF = heightfield(ulvl);
						number uval = uHF.interpolate(c);
						if(uval != uHF.no_data_value()){
							gotDataVal = true;
							if(val + minHeight > uval){
							//	the height in the current cell is too small
								innerField.at(ix, iy) = curHF.no_data_value();
							}
							break;
						}
					}

					if(!gotDataVal){
					//	The cell lies outside of the domain which is specified
					//	by the upmost field. The cell thus has to be invalidated
						innerField.at(ix, iy) = curHF.no_data_value();
					}
				}
			}
		}
	}
}


void  RasterLayers::
invalidate_small_lenses(number minArea)
{
	for(size_t lvl = 0; lvl < size(); ++lvl){
		Heightfield& hf = heightfield(lvl);
		number cellSize = hf.cell_size().x() * hf.cell_size().y();
		if(cellSize > 0){
			InvalidateSmallLenses(hf.field(), minArea / cellSize,
								  hf. no_data_value());
		}
	}
}

void RasterLayers::
remove_small_holes(number maxArea)
{
	using namespace std;

	const size_t numNbrs = 4;
	const int xadd[numNbrs] = {0, -1, 1, 0};
	const int yadd[numNbrs] = {-1, 0, 0, 1};

	Field<bool>	visited;
	vector<pair<int, int> > cells;

	for(int lvl = (int)size() - 2; lvl >= 0; --lvl){
		Heightfield& hf = heightfield(lvl);
		Field<number>& field = hf.field();
		number noDataValue = hf.no_data_value();
		number minHeight = min_height(lvl);

		number cellSize = hf.cell_size().x() * hf.cell_size().y();
		if(cellSize <= 0)
			continue;

		size_t thresholdCellCount(maxArea / cellSize);


		visited.resize_no_copy(field.width(), field.height());
		visited.fill_all(false);

		const int fwidth = (int)field.width();
		const int fheight = (int)field.height();

		for(int outerIy = 0; outerIy < fheight; ++outerIy){
			for(int outerIx = 0; outerIx < fwidth; ++outerIx){
				if(visited.at(outerIx, outerIy)
				   || (field.at(outerIx, outerIy) != noDataValue))
				{
					continue;
				}

				cells.clear();
				cells.push_back(make_pair(outerIx, outerIy));
				size_t curCell = 0;
				while(curCell < cells.size()){
					int ix = cells[curCell].first;
					int iy = cells[curCell].second;

					for(size_t inbr = 0; inbr < numNbrs; ++inbr){
						int nx = ix + xadd[inbr];
						int ny = iy + yadd[inbr];
						if((nx >= 0 && nx < fwidth && ny >= 0 && ny < fheight)
						   && !visited.at(nx, ny))
						{
							visited.at(nx, ny) = true;
							if(field.at(nx, ny) == noDataValue){
								cells.push_back(make_pair(nx, ny));
							}
						}
					}
					++curCell;
				}
				if(cells.size() < thresholdCellCount){
					for(size_t i = 0; i < cells.size(); ++i){
						int ix = cells[i].first;
						int iy = cells[i].second;
						vector2 c = hf.index_to_coordinate(ix, iy);
						pair<int, number> result = trace_line_up(c, (size_t)lvl);
						if(result.first > lvl){
						//	we have to adjust not only the value in the current layer
						//	but also possibly values in lower levels, to avoid
						//	creating new thin layers
							number curVal = result.second - minHeight;
							field.at(ix, iy) = curVal;
							for(int lowerLvl = lvl - 1; lowerLvl >= 0 ; --lowerLvl){
								Heightfield& lhf = heightfield(lowerLvl);
								Field<number>& lfield = lhf.field();
								pair<int, int> index = lhf.coordinate_to_index(c.x(), c.y());
								int lx = index.first;
								int ly = index.second;
								if(lx >= 0 && lx < (int)lfield.width()
									&& ly >= 0 && ly < (int)lfield.height())
								{
									number lval = lfield.at(lx, ly);
									if(lval != lhf.no_data_value()){
										if(curVal - lval < minHeight){
											curVal -= minHeight;
											lfield.at(lx, ly) = curVal;
										}
										else
											break;
									}
								}
							}
						}
					}
				}
			}
		}
	}
}


void RasterLayers::
snap_cells_to_higher_layers()
{
	if(size() <= 1)
		return;

	for(int lvl = (int)size() - 2; lvl >= 0; --lvl){
		Heightfield& curHF = heightfield(lvl);
		Field<number>& curField = curHF.field();
		Heightfield& upperHF = heightfield(lvl + 1);
		number minHeight = min_height(lvl);

		for(size_t iy = 0; iy < curField.height(); ++iy){
			for(size_t ix = 0; ix < curField.width(); ++ix){
				number curVal = curField.at(ix, iy);
				vector2 c = curHF.index_to_coordinate(ix, iy);
				number upperVal = upperHF.interpolate(c);
				if(upperVal == upperHF.no_data_value())
					continue;

				if((curVal == curHF.no_data_value()) || (curVal > upperVal)
					|| (upperVal - curVal < minHeight))
				{
					curField.at(ix, iy) = upperVal;
				}
			}
		}
	}
}


void RasterLayers::
eliminate_invalid_cells()
{
	for(size_t lvl = 0; lvl < size(); ++lvl){
		Heightfield& hf = heightfield(lvl);
		EliminateInvalidCells(hf.field(), hf.no_data_value());
	}
}


void RasterLayers::
blur_layers(number alpha, size_t numIterations)
{
	for(size_t i = 0; i < size(); ++i){
		BlurField(heightfield(i).field(), alpha, numIterations, heightfield(i).no_data_value());
	}
}

std::pair<int, number> RasterLayers::
trace_line_down(const vector2& c, size_t firstLayer) const
{
	for(int i = (int)firstLayer; i >= 0; --i){
		number val = heightfield(i).interpolate(c);
		if(val != heightfield(i).no_data_value()){
			return make_pair(i, val);
		}
	}

	return make_pair<int, number>(-1, 0);
}

std::pair<int, number> RasterLayers::
trace_line_up(const vector2& c, size_t firstLayer) const
{
	for(size_t i = firstLayer; i < size(); ++i){
		number val = heightfield(i).interpolate(c);
		if(val != heightfield(i).no_data_value()){
			return make_pair((int)i, val);
		}
	}

	return make_pair<int, number>(-1, 0);
}

}//	end of namespace
