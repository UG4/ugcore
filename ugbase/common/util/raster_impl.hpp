/*
 * Copyright (c) 2016:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG_raster_impl
#define __H__UG_raster_impl

#include <limits>
#include <cstring>
#include <algorithm>
#include <fstream>
#include "common/error.h"
#include "common/util/file_util.h"
#include "common/util/string_util.h"
#include "raster_kernels.h"

namespace ug{

////////////////////////////////////////////////////////////////////////////////
//	Raster::MultiIndex

template <class T, int TDIM>
Raster<T, TDIM>::MultiIndex::
MultiIndex()
{}

template <class T, int TDIM>
Raster<T, TDIM>::MultiIndex::
MultiIndex(size_t i)
{
	set(i);
}

template <class T, int TDIM>
int Raster<T, TDIM>::MultiIndex::
dim () const				
{
	return TDIM;
}

template <class T, int TDIM>
void Raster<T, TDIM>::MultiIndex::
set (size_t i)				
{
	for(int d = 0; d < TDIM; ++d)
		m_ind[d] = i;
}

template <class T, int TDIM>
size_t& Raster<T, TDIM>::MultiIndex::
operator[] (int d)			
{
	return m_ind[d];
}

template <class T, int TDIM>
size_t Raster<T, TDIM>::MultiIndex::
operator[] (int d) const	
{
	return m_ind[d];
}


////////////////////////////////////////////////////////////////////////////////
//	Raster::Coordinate

template <class T, int TDIM>
Raster<T, TDIM>::Coordinate::
Coordinate()
{}

template <class T, int TDIM>
Raster<T, TDIM>::Coordinate::
Coordinate(number c)
{
	set(c);
}

template <class T, int TDIM>
Raster<T, TDIM>::Coordinate::
Coordinate(const MathVector<TDIM, number>& v)
{
	for(int d = 0; d < TDIM; ++d)
		m_coord[d] = v[d];
}

template <class T, int TDIM>
int Raster<T, TDIM>::Coordinate::
dim () const
{
	return TDIM;
}

template <class T, int TDIM>
void Raster<T, TDIM>::Coordinate::
set (number c)
{
	for(int d = 0; d < TDIM; ++d)
		m_coord[d] = c;
}


template <class T, int TDIM>
number& Raster<T, TDIM>::Coordinate::
operator[] (int d)
{
	return m_coord[d];
}

template <class T, int TDIM>
number Raster<T, TDIM>::Coordinate::
operator[] (int d) const
{
	return m_coord[d];
}

template <class T, int TDIM>
typename Raster<T, TDIM>::Coordinate& Raster<T, TDIM>::Coordinate::
operator+= (const Coordinate& c)
{
	for(int d = 0; d < TDIM; ++d)
		m_coord[d] += c[d];
}

template <class T, int TDIM>
typename Raster<T, TDIM>::Coordinate& Raster<T, TDIM>::Coordinate::
operator-= (const Coordinate& c)
{
	for(int d = 0; d < TDIM; ++d)
		m_coord[d] -= c[d];
}

template <class T, int TDIM>
typename Raster<T, TDIM>::Coordinate& Raster<T, TDIM>::Coordinate::
operator*= (number s)
{
	for(int d = 0; d < TDIM; ++d)
		m_coord[d] *= s;
}



////////////////////////////////////////////////////////////////////////////////
//	Raster - public

template <class T, int TDIM>
Raster<T, TDIM>::
Raster () :
	m_data(NULL),
	m_numNodes(0),
	m_selNode(0),
	m_minCorner(0),
	m_extension(1),
	m_cellExtension(1),
	m_cursor(0),
	m_numNodesTotal(0),
	m_noDataValue(std::numeric_limits<T>::max())
{
	update_num_nodes_total();
	update_cell_extension();
}


template <class T, int TDIM>
Raster<T, TDIM>::
Raster (const Raster<T, TDIM>& raster) :
	m_data(NULL),
	m_numNodes(raster.m_numNodes),
	m_selNode(raster.m_selNode),
	m_minCorner(raster.m_minCorner),
	m_extension(raster.m_extension),
	m_cellExtension(raster.m_cellExtension),
	m_cursor(raster.m_cursor),
	m_numNodesTotal(raster.m_numNodesTotal),
	m_noDataValue(raster.m_noDataValue)
{
	update_num_nodes_total();
	update_cell_extension();

	create();

	memcpy(m_data, raster.m_data, num_nodes_total() * sizeof(T));
}

template <class T, int TDIM>
Raster<T, TDIM>::
Raster (const MultiIndex& numNodes) :
	m_data(NULL),
	m_numNodes(numNodes),
	m_selNode(0),
	m_minCorner(0),
	m_extension(1),
	m_cellExtension(1),
	m_cursor(0),
	m_numNodesTotal(0),
	m_noDataValue(std::numeric_limits<T>::max())
{
	update_num_nodes_total();

	for(int d = 0; d < TDIM; ++d)
		m_extension[d] = numNodes[d] - 1;

	update_cell_extension();

	create();
}

template <class T, int TDIM>
Raster<T, TDIM>::
Raster (const MultiIndex& numNodes,
		const Coordinate& extension,
		const Coordinate& minCorner) :
	m_data(NULL),
	m_numNodes(numNodes),
	m_selNode(0),
	m_minCorner(minCorner),
	m_extension(extension),
	m_cellExtension(1),
	m_cursor(0),
	m_numNodesTotal(0),
	m_noDataValue(std::numeric_limits<T>::max())
{
	update_num_nodes_total();
	update_cell_extension();
	create();
}


template <class T, int TDIM>
Raster<T, TDIM>::
~Raster ()
{
	if(m_data)
		delete[] m_data;
}

template <class T, int TDIM>
Raster<T, TDIM>& Raster<T, TDIM>::
operator= (const Raster& raster)
{
	m_numNodes		= raster.m_numNodes;
	m_selNode		= raster.m_selNode;
	m_minCorner		= raster.m_minCorner;
	m_extension		= raster.m_extension;
	m_cursor		= raster.m_cursor;
	m_noDataValue	= raster.m_noDataValue;

	update_num_nodes_total();
	update_cell_extension();

	create();

	memcpy(m_data, raster.m_data, num_nodes_total() * sizeof(T));

	return *this;
}

template <class T, int TDIM>
int Raster<T, TDIM>::
dim () const
{
	return TDIM;
}

template <class T, int TDIM>
void Raster<T, TDIM>::
set_num_nodes (int dim, size_t num)
{
	m_numNodes[dim] = num;
	update_num_nodes_total();
	update_cell_extension(dim);
}

template <class T, int TDIM>
void Raster<T, TDIM>::
set_num_nodes (const typename Raster<T, TDIM>::MultiIndex& mi)
{
	m_numNodes = mi;
	update_num_nodes_total();
	update_cell_extension();
}

template <class T, int TDIM>
size_t Raster<T, TDIM>::
num_nodes_total () const
{
	return m_numNodesTotal;
}

template <class T, int TDIM>
size_t Raster<T, TDIM>::
num_nodes (int dim) const
{
	return m_numNodes[dim];
}

template <class T, int TDIM>
const typename Raster<T, TDIM>::MultiIndex& Raster<T, TDIM>::
num_nodes () const
{
	return m_numNodes;
}


template <class T, int TDIM>
void Raster<T, TDIM>::
create ()
{
	if(m_data){
		delete[] m_data;
		m_data = NULL;
	}

	update_num_nodes_total(); // this isn't strictly necessary if everything works right.

	const size_t num = num_nodes_total();
	if(num){
		m_data = new T[num];
	}
}


template <class T, int TDIM>
T& Raster<T, TDIM>::
node_value (const MultiIndex& mi)
{
	return m_data[data_index(mi)];
}

template <class T, int TDIM>
T Raster<T, TDIM>::
node_value (const MultiIndex& mi) const
{
	return m_data[data_index(mi)];
}


template <class T, int TDIM>
void Raster<T, TDIM>::
set_min_corner (int dim, number coord)
{
	m_minCorner[dim] = coord;
}

template <class T, int TDIM>
void Raster<T, TDIM>::
set_min_corner (const Coordinate& coord)
{
	m_minCorner = coord;
}

template <class T, int TDIM>
const typename Raster<T, TDIM>::Coordinate& Raster<T, TDIM>::
min_corner () const
{
	return m_minCorner;
}

template <class T, int TDIM>
number Raster<T, TDIM>::
min_corner (int dim) const
{
	return m_minCorner[dim];
}

template <class T, int TDIM>
void Raster<T, TDIM>::
set_extension (int dim, number ext)
{
	m_extension[dim] = ext;
	update_cell_extension(dim);
}

template <class T, int TDIM>
void Raster<T, TDIM>::
set_extension (const Coordinate& ext)
{
	m_extension = ext;
	update_cell_extension();
}

template <class T, int TDIM>
const typename Raster<T, TDIM>::Coordinate& Raster<T, TDIM>::
extension () const
{
	return m_extension;
}

template <class T, int TDIM>
number Raster<T, TDIM>::
extension (int dim) const
{
	return m_extension[dim];
}

template <class T, int TDIM>
const typename Raster<T, TDIM>::MultiIndex&  Raster<T,TDIM>::
node_index(const Coordinate& coord, int order) const
{
	MultiIndex mi(-1);
	for(size_t d = 0; d < TDIM; ++d)
	{ 
		switch(order){
		case 0: {
			mi[d] = static_cast<int>(0.5 + (coord[d] - m_minCorner[d]) / m_cellExtension[d]);
			if(mi[d] < 0)					mi[d] = 0;
			else if(mi[d] >= num_nodes(d))	mi[d] = num_nodes(d) - 1;
		}break;

		case 1:{
			mi[d] = static_cast<int>((coord[d] - m_minCorner[d]) / m_cellExtension[d]);
			if(mi[d] < 0)					mi[d] = 0;
			else if(mi[d]+1 >= num_nodes(d))	mi[d] = num_nodes(d) - 2;
		}break;

		default:
			UG_THROW("Raster::interpolate(): Unsupported interpolation order: " << order);
		}
	}
	return mi;
}

template<class T, int TDIM>
const AABox<number> Raster<T,TDIM>::
bounding_box(const MultiIndex& mi) const{
	AABox<MathVector<TDIM, number> > box;
	for(size_t d = 0; d < TDIM; ++d)
	{
		box.min[d]=(number)mi[d] * m_cellExtension[d]+ m_minCorner[d];
		box.max[d]=box.min[d]+m_cellExtension[d];
	}	
	return box;
}

template <class T, int TDIM>
T Raster<T, TDIM>::
interpolate (const Coordinate& coord, int order) const
{
	switch(order){
		case 0: {
			MultiIndex mi;
			for(size_t d = 0; d < TDIM; ++d)
			{
				mi[d] = static_cast<int>(0.5 + (coord[d] - m_minCorner[d]) / m_cellExtension[d]);
				if(mi[d] < 0)					mi[d] = 0;
				else if(mi[d] >= num_nodes(d))	mi[d] = num_nodes(d) - 1;
			}
			return node_value(mi);
		} break;

		case 1: {
			MultiIndex mi;
			Coordinate lc;
			for(size_t d = 0; d < TDIM; ++d)
			{
				mi[d] = static_cast<int>((coord[d] - m_minCorner[d]) / m_cellExtension[d]);
				if(mi[d] < 0){
					mi[d] = 0;
					lc[d] = 0;
				}
				else if(mi[d] + 1 >= num_nodes(d)){
					mi[d] = num_nodes(d) - 2;
					lc[d] = 1;
				}
				else{
					lc[d] = ( coord[d]
							  - ((number)mi[d] * m_cellExtension[d]
							  	 + m_minCorner[d]))
							/ m_cellExtension[d];
				}
			}

			return interpolate_linear(mi, lc);
		} break;

		default:
			UG_THROW("Raster::interpolate(): Unsupported interpolation order: " << order);
	}
	return m_noDataValue;
}


template <class T, int TDIM>
void Raster<T, TDIM>::
set_no_data_value(T val)
{
	m_noDataValue = val;
}

template <class T, int TDIM>
T Raster<T, TDIM>::
no_data_value() const
{
	return m_noDataValue;
}


template <class T, int TDIM>
void Raster<T, TDIM>::
blur(T alpha, size_t iterations)
{
	raster_kernels::Blur<T, TDIM> blurKernel (alpha);
	const MultiIndex start(0);

	for(size_t iiter = 0; iiter < iterations; ++iiter)
		run_on_all (blurKernel);
}


template <class T, int TDIM>
template <class TKernel>
typename TKernel::result_t Raster<T, TDIM>::
run_on_all()
{
	TKernel kernel;
	run_on_all (MultiIndex(0), kernel, TDIM - 1);
	return kernel.result();
}

template <class T, int TDIM>
template <class TKernel>
void Raster<T, TDIM>::
run_on_all(TKernel& kernel)
{
	run_on_all (MultiIndex(0), kernel, TDIM - 1);
}


template <class T, int TDIM>
template <class TKernel>
void Raster<T, TDIM>::
run_on_all(const MultiIndex& start, TKernel& kernel, int curDim)
{
	if(curDim > 0) {
		const size_t numNodes = num_nodes(curDim);
		for(MultiIndex cur = start; cur[curDim] < numNodes; ++cur[curDim]){
			run_on_all(cur, kernel, curDim - 1);
		}
	}
	else {
		const size_t numNodes = num_nodes(0);
		for(MultiIndex cur = start; cur[0] < numNodes; ++cur[0]){
			kernel (*this, cur);
		}
	}
}


template <class T, int TDIM>
template <class TKernel>
typename TKernel::result_t Raster<T, TDIM>::
run_on_nbrs(const MultiIndex& center)
{
	TKernel kernel;
	run_on_nbrs(center, kernel, TDIM - 1);
	return kernel.result();
}


template <class T, int TDIM>
template <class TKernel>
void Raster<T, TDIM>::
run_on_nbrs(const MultiIndex& center, TKernel& kernel)
{
	run_on_nbrs(center, kernel, TDIM - 1);
}


template <class T, int TDIM>
template <class TKernel>
void Raster<T, TDIM>::
run_on_nbrs(const MultiIndex& center, TKernel& kernel, int curDim)
{
	if(curDim > 0)
		run_on_nbrs(center, kernel, curDim - 1);

	if(center[curDim] > 0){
		MultiIndex c = center;
		--c[curDim];
		kernel(*this, c);
	}

	if(center[curDim] + 1 < num_nodes(curDim)){
		MultiIndex c = center;
		++c[curDim];
		kernel(*this, c);
	}
}


template <class T, int TDIM>
void Raster<T, TDIM>::
select_node (int dim, size_t index)
{
	m_selNode[dim] = index;
}

template <class T, int TDIM>
void Raster<T, TDIM>::
select_node (const MultiIndex& mi)
{
	m_selNode = mi;
}


template <class T, int TDIM>
void Raster<T, TDIM>::
set_selected_node_value (T val)
{
	node_value(m_selNode) = val;
}

template <class T, int TDIM>
T Raster<T, TDIM>::
selected_node_value () const
{
	return node_value(m_selNode);
}


template <class T, int TDIM>
void Raster<T, TDIM>::
set_cursor (int dim, number coord)
{
	m_cursor[dim] = coord;
}

template <class T, int TDIM>
void Raster<T, TDIM>::
set_cursor (const Coordinate& coord)
{
	m_cursor = coord;
}


template <class T, int TDIM>
T Raster<T, TDIM>::
interpolate_at_cursor (int order) const
{
	return interpolate(m_cursor, order);
}


template <class T, int TDIM>
void Raster<T, TDIM>::
load_from_asc (const char* filename)
{
	using namespace std;

	#define LFA_ERR_WHERE "Error in Raster::load_from_asc('" << filename << "'): "
//	this macro helps with error-checks for bad dimensions
	#define LFA_CHECK_DIM(d, line)\
				UG_COND_THROW(d >= TDIM, LFA_ERR_WHERE << "Bad dimension '" << d <<\
					"' in line " << line << " of file " << filename << "," <<\
					"while trying to read a " << TDIM << "d raster.");

	std::string fullFileName = FindFileInStandardPaths(filename);
	UG_COND_THROW(fullFileName.empty(),
				  LFA_ERR_WHERE << "Couldn't find the specified file in any of the standard paths.");

	ifstream in(fullFileName.c_str());
	UG_COND_THROW(!in, LFA_ERR_WHERE << "Couldn't access file.");

	MultiIndex numNodes(0);

//	indicate whether minCoord was specified as cell-center
	MultiIndex minCoordIsCenter(0);

	Coordinate minCoord(0);
	Coordinate cellSize(0);

	T noDataValue = T();

//	parse header
//	the header lenght varies between different asc files. This depends on the
//	dimension and whether equlateral cells are specified or not, i.e., whether
//	'cellsize' or whether 'xcellsize', 'ycellsize', ... was specified.
//	We're trying to guess the correct length here.
	int headerLen = 0;
	if(TDIM == 1)		headerLen = 4;
	else if(TDIM == 2)	headerLen = 6;
	else if(TDIM == 3)	headerLen = 8;
	else{
		UG_THROW("Raster::load_from_asc only supports 1, 2, and 3 dimensions\n");
	}
	
	for(int i = 0; i < headerLen; ++i){
		string name;
		double value;
		in >> name >> value;
		UG_COND_THROW(!in, LFA_ERR_WHERE <<
					  "Couldn't parse expected name-value pair in row " << i);

		name = ToLower(name);

		if(name.compare("ncols") == 0){
			LFA_CHECK_DIM(0, i);
			numNodes[0] = (int)value;
		}
		else if(name.compare("nrows") == 0){
			LFA_CHECK_DIM(1, i);
			numNodes[1] = (int)value;
		}
		else if(name.compare("nstacks") == 0){
			LFA_CHECK_DIM(2, i);
			numNodes[2] = (int)value;
		}

		else if(name.compare("xllcenter") == 0){
			LFA_CHECK_DIM(0, i);
			minCoord[0] = value;
			minCoordIsCenter[0] = 1;
		}
		else if(name.compare("yllcenter") == 0){
			LFA_CHECK_DIM(1, i);
			minCoord[1] = value;
			minCoordIsCenter[1] = 1;
		}
		else if(name.compare("zllcenter") == 0){
			LFA_CHECK_DIM(2, i);
			minCoord[2] = value;
			minCoordIsCenter[2] = 1;
		}
		else if(name.compare("xllcorner") == 0){
			LFA_CHECK_DIM(0, i);
			minCoord[0] = value;
		}
		else if(name.compare("yllcorner") == 0){
			LFA_CHECK_DIM(1, i);
			minCoord[1] = value;
		}
		else if(name.compare("zllcorner") == 0){
			LFA_CHECK_DIM(2, i);
			minCoord[2] = value;
		}

		else if(name.compare("cellsize") == 0){
			for(int d = 0; d < TDIM; ++d)
				cellSize[d] = value;
		}

		else if(name.compare("xcellsize") == 0){
			LFA_CHECK_DIM(0, i);
		//	we have to read additional cell-sizes for the other dimensions
			headerLen += (TDIM - 1);
			cellSize[0] = value;
		}
		else if(name.compare("ycellsize") == 0){
			LFA_CHECK_DIM(1, i);
			cellSize[1] = value;
		}
		else if(name.compare("zcellsize") == 0){
			LFA_CHECK_DIM(2, i);
			cellSize[2] = value;
		}

		else if(name.compare("nodata_value") == 0){
			noDataValue = value;
		}

		else{
			UG_THROW(LFA_ERR_WHERE << "unknown identifier in header: " << name);
		}
	}
	
	for(int d = 0; d < TDIM; ++d){
		if(minCoordIsCenter[d])
			minCoord[d] -= 0.5 * cellSize[d];
	}

//	check validity
	for(int d = 0; d < TDIM; ++d){
		UG_COND_THROW(numNodes[d] == 0, LFA_ERR_WHERE << "Num nodes may not be 0 for dim " << d);
		UG_COND_THROW(cellSize[d] <= 0, LFA_ERR_WHERE << "cell-size must be bigger than 0 for dim " << d);
	}

	set_num_nodes(numNodes);
	set_min_corner(minCoord);
	Coordinate extension = cellSize;
	for(int d = 0; d < TDIM; ++d)
		extension[d] *= (number)(numNodes[d] - 1);
	set_extension(extension);
	set_no_data_value(noDataValue);

	create();

//	parse values
//	y and z are inverted
	size_t num[3] = {0, 1, 1};
	for(size_t i = 0; i < TDIM; ++i)
		num[i] = m_numNodes[i];
		
	for(size_t iz = 0; iz < num[2]; ++iz){
		for(size_t iy = 0; iy < num[1]; ++iy){
			for(size_t ix = 0; ix < num[0]; ++ix)
			{
				const size_t ty = num[1] - 1 - iy;
				const size_t tz = num[2] - 1 - iz;
				in >> m_data[ix + num[0] * (ty + num[1] * tz)];
				UG_COND_THROW(!in, LFA_ERR_WHERE << "Couldn't read value for at ("
							  << ix << ", " << iy << ", " << iz << ")");
			}
		}
	}
}

template <class T, int TDIM>
void Raster<T, TDIM>::
save_to_asc (const char* filename) const
{
	using namespace std;
	#define STA_ERR_WHERE "Error in Raster::save_to_asc('" << filename << "'): "

	UG_COND_THROW(!m_data, STA_ERR_WHERE << "Can't write an unitinialized raster."
				  "Please call 'create' or 'load_from_asc' first.");

	ofstream out(filename);
	UG_COND_THROW(!out, STA_ERR_WHERE << "Couldn't open file for writing.");

	out << "ncols         " << num_nodes(0) << endl;
	if(TDIM > 1)
		out << "nrows         " << num_nodes(1) << endl;
	if(TDIM > 2)
		out << "nstacks       " << num_nodes(2) << endl;

	out << "xllcorner     " << setprecision(16) << min_corner(0) << endl;
	if(TDIM > 1)
		out << "yllcorner     " << setprecision(16) << min_corner(1) << endl;
	if(TDIM > 2)
		out << "zllcorner     " << setprecision(16) << min_corner(2) << endl;

	bool equlateralCells = true;
	for(int d = 1; d < TDIM; ++d){
		if(m_cellExtension[d] != m_cellExtension[0]){
			equlateralCells = false;
			break;
		}
	}

	if(equlateralCells)
		out << "cellsize      " << m_cellExtension[0] << endl;
	else{
		out << "xcellsize     " << m_cellExtension[0] << endl;
		if(TDIM > 1)
			out << "ycellsize     " << m_cellExtension[1] << endl;
		if(TDIM > 2)
			out << "zcellsize     " << m_cellExtension[2] << endl;
	}

	out << "NODATA_value  " << no_data_value() << endl;

//	write values
//	y and z are inverted
	size_t num[3] = {0, 1, 1};
	for(size_t i = 0; i < TDIM; ++i)
		num[i] = num_nodes(i);
		
	for(size_t iz = 0; iz < num[2]; ++iz){
		for(size_t iy = 0; iy < num[1]; ++iy){
			for(size_t ix = 0; ix < num[0]; ++ix)
			{
				if(ix > 0)
					out << " ";
				const size_t ty = num[1] - 1 - iy;
				const size_t tz = num[2] - 1 - iz;
				out << m_data[ix + num[0] * (ty + num[1] * tz)];
			}
			out << endl;
		}
		out << endl;
	}
	out << endl;
}

////////////////////////////////////////////////////////////////////////////////
//	Raster - private
template <class T, int TDIM>
size_t Raster<T, TDIM>::
data_index (const MultiIndex& mi, int curDim, size_t curVal) const
{
	if(curDim == 0)
		return curVal + mi[0];
	else{
		return data_index(mi, curDim - 1, m_numNodes[curDim - 1] * (mi[curDim] + curVal));
	}
}

template <class T, int TDIM>
void Raster<T, TDIM>::
update_num_nodes_total()
{
	m_numNodesTotal = 1;
	for(int d = 0; d < TDIM; ++d)
		m_numNodesTotal *= num_nodes(d);
}

template <class T, int TDIM>
void Raster<T, TDIM>::
update_cell_extension()
{
	for(int d = 0; d < TDIM; ++d)
		update_cell_extension(d);
}

template <class T, int TDIM>
void Raster<T, TDIM>::
update_cell_extension(int dim)
{
	if(m_numNodes[dim] > 1 && m_extension[dim] > 0)
		m_cellExtension[dim] = m_extension[dim] / (m_numNodes[dim] - 1);
	else
		m_cellExtension[dim] = 1;
}


template <class T, int TDIM>
T Raster<T, TDIM>::
interpolate_linear (
		const MultiIndex& minNodeInd,
		Coordinate& localCoord,
		int curDim) const
{
	if(curDim == 0)
		return node_value(minNodeInd);

	MultiIndex miMax = minNodeInd;
	miMax[curDim - 1] += 1;

	T val0 = interpolate_linear(minNodeInd, localCoord, curDim - 1);
	T val1 = interpolate_linear(miMax, localCoord, curDim - 1);

//	perform linear interpolation
	val0 *= (1. - localCoord[curDim - 1]);
	val1 *= localCoord[curDim - 1];
	val0 += val1;
	return val0;
}

}//	end of namespace

#endif	//__H__UG_raster_impl
