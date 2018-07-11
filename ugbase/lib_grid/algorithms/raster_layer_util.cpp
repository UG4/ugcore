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

#include <vector>
#include "raster_layer_util.h"
#include "field_util.h"
#include "lib_grid/file_io/file_io_asc.h"

using namespace std;

namespace ug{

//	initialize all values with simple-values and record a list of all invalid ones
struct CellIdx {
	CellIdx()	{}
	CellIdx(size_t x, size_t y) : ix(x), iy(y)	{}
	size_t ix;
	size_t iy;
};

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

std::pair<int, int> RasterLayers::
get_layer_indices(const vector3& c) const
{
	vector2 xy(c.x(), c.y());
	
	int upperLayer = -1;
	int lowerLayer = -1;
	while(1){
		std::pair<int, number> cut = trace_line_down(xy, upperLayer + 1);
		if(cut.first == -1)
			break;
		else if(cut.second > c.z()){
			lowerLayer = cut.first;
			break;
		}
		else
			upperLayer = cut.first;
	}

	return make_pair(upperLayer, lowerLayer);
}

number RasterLayers::
relative_to_global_height(const vector2& c, number relHeight) const
{
	static const int order = 1;

	if(m_relativeToGlobalHeights.empty())
		return relative_to_global_height_simple(c, relHeight);
	else{
		const int relHeightLow = (int)(relHeight + SMALL);
		const int topLvl = (int)m_relativeToGlobalHeights.size() - 1;
		if((relHeightLow >= 0)
			&& (relHeightLow + 1 <= topLvl))
		{
			number rel = relHeight - (number)relHeightLow;
			return (1. - rel) * m_relativeToGlobalHeights[relHeightLow]
									->interpolate(c, order)
					+ rel * m_relativeToGlobalHeights[relHeightLow + 1]
									->interpolate(c, order);
		}
		else if(relHeightLow >= topLvl)
			return m_relativeToGlobalHeights[topLvl]->interpolate(c, order);
		else
			return m_relativeToGlobalHeights[0]->interpolate(c, order);
	}
}

number RasterLayers::
relative_to_global_height_simple(const vector2& c, number relHeight) const
{
	const int relHeightLower = max((int)(relHeight + SMALL), 0);
	const int relHeightUpper = min(relHeightLower + 1, (int)num_layers() - 1);
	const std::pair<int, number> lower = trace_line_down(c, relHeightLower);
	const std::pair<int, number> upper = trace_line_up(c, relHeightUpper);
	
	UG_COND_THROW(lower.first < 0, "Invalid lower layer for coordinate " << c
				  << " and relative height " << relHeight);
	UG_COND_THROW(upper.first < 0, "Invalid upper layer for coordinate " << c
				  << " and relative height " << relHeight);

	if(upper.second < lower.second)
		return upper.second;
	
	number layerDiff = max(1, upper.first - lower.first);
	number rel = (clip<number>(relHeight - (number)relHeightLower, 0, 1)
				  + relHeightLower - lower.first)
				/ layerDiff;
	return (1. - rel) * lower.second + rel * upper.second;
}

void RasterLayers::
construct_relative_to_global_height_table(size_t iterations, number alpha)
{
	m_relativeToGlobalHeights.clear();
	for(size_t i = 0; i < m_layers.size(); ++i){
		m_relativeToGlobalHeights.push_back(
				make_sp(new Heightfield(m_layers[i]->heightfield)));
	}


//	make sure that the bottom-layer has no holes (compared to the top layer).
//	we'll set it to the same height as its local upper layer.
	{
		Heightfield& baseHF = *m_relativeToGlobalHeights[0];
		Field<number>& baseField = baseHF.field();

		for(size_t iy = 0; iy < baseField.height(); ++iy){
			for(size_t ix = 0; ix < baseField.width(); ++ix){
				if(baseField.at(ix, iy) == baseHF.no_data_value()){
					vector2 c = baseHF.index_to_coordinate(ix, iy);
					std::pair<int, number> upper = trace_line_up(c, 0);
					if(upper.first != -1)
						baseField.at(ix, iy) = upper.second;
				}
			}
		}
	}


	std::vector<std::vector<CellIdx> > allCells;
	allCells.resize(m_relativeToGlobalHeights.size());

	for(int lvl = 1; lvl + 1 < (int)m_relativeToGlobalHeights.size(); ++lvl){
		Heightfield& curHF = *m_relativeToGlobalHeights[lvl];
		Field<number>& curField = curHF.field();

		std::vector<CellIdx>& cells = allCells[lvl];

		for(size_t iy = 0; iy < curField.height(); ++iy){
			for(size_t ix = 0; ix < curField.width(); ++ix){
				if(curField.at(ix, iy) == curHF.no_data_value()){
					vector2 c = curHF.index_to_coordinate(ix, iy);
					if(trace_line_up(c, lvl).first != -1){
						curField.at(ix, iy) = relative_to_global_height_simple(c, lvl);
						cells.push_back(CellIdx(ix, iy));
					}
				}
			}
		}
	}

//	now iterate over all recorded cells and perform relaxation
	for(size_t iteration = 0; iteration < iterations; ++iteration){
		for(int lvl = 1; lvl + 1 < (int)m_relativeToGlobalHeights.size(); ++lvl){
			Field<number>& cur = m_relativeToGlobalHeights[lvl]->field();
			Field<number>& lower = m_relativeToGlobalHeights[lvl-1]->field();
			Field<number>& upper = m_relativeToGlobalHeights[lvl+1]->field();

			std::vector<CellIdx>& cells = allCells[lvl];
			for(size_t icell = 0; icell < cells.size(); ++icell){
				const CellIdx ci = cells[icell];
			//	average dist-relations of neighbors
				number nbrVal = 0;
				number numNbrs = 0;

				if(ci.ix > 0){
					nbrVal += upper_lower_dist_relation(
									lower, cur,
									upper, ci.ix - 1, ci.iy);
					++numNbrs;
				}
				if(ci.ix + 1 < cur.width()){

					nbrVal += upper_lower_dist_relation(
									lower, cur,
									upper, ci.ix + 1, ci.iy);
					++numNbrs;
				}
				if(ci.iy > 0){
					nbrVal += upper_lower_dist_relation(
									lower, cur,
									upper, ci.ix, ci.iy - 1);
					++numNbrs;
				}
				if(ci.iy + 1 < cur.height()){
					nbrVal += upper_lower_dist_relation(
									lower, cur,
									upper, ci.ix, ci.iy + 1);
					++numNbrs;
				}

				if(numNbrs > 0){
					number locVal = upper_lower_dist_relation(
											lower, cur,
											upper, ci.ix, ci.iy);

					locVal = (1. - alpha) * locVal + alpha * nbrVal / numNbrs;

				//	reconstruct height value from dist-relation
					const number upperVal = upper.at(ci.ix, ci.iy);
					const number lowerVal = lower.at(ci.ix, ci.iy);

					number dirTotal = upperVal - lowerVal;
					cur.at(ci.ix, ci.iy) = upperVal - locVal * dirTotal;

				//	enforce minWidth constraints
					const number upperDist = upperVal - cur.at(ci.ix, ci.iy);
					const number lowerDist = cur.at(ci.ix, ci.iy) - lowerVal;
					if(upperDist < layer(lvl).minHeight)
						if(lowerDist < layer(lvl-1).minHeight)
							cur.at(ci.ix, ci.iy)
								= 0.5 * (upperVal - upperDist + lowerVal + lowerDist);
						else
							cur.at(ci.ix, ci.iy) = upperVal - upperDist;
					else if(lowerDist < layer(lvl-1).minHeight)
						cur.at(ci.ix, ci.iy) = lowerVal + lowerDist;
				}
			}
		}
	}
}


void RasterLayers::
invalidate_relative_to_global_height_table ()
{
	m_relativeToGlobalHeights.clear();
}

number RasterLayers::
upper_lower_dist_relation (	Field<number>&lower,
							Field<number>& middle,
							Field<number>& upper,
							size_t ix,
							size_t iy)
{
	number totalDist = fabs(upper.at(ix, iy) - lower.at(ix, iy));
	if(totalDist > 0)
		return (upper.at(ix, iy) - middle.at(ix, iy)) / totalDist;
	return 0;
}

// void RasterLayers::
// save_to_files(const char* filenamePrefix_cstr)
// {
// 	std::string filenamePrefix = filenamePrefix_cstr;
// 	for(size_t i = 0; i < m_layers.size(); ++i){
		
// 	}
// }

// void RasterLayers::
// save_rel_to_glob_table_to_files(const char* filenamePrefix_cstr)
// {

// }



}//	end of namespace
