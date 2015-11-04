/**
 * 	connection_viewer_input.h
 *
 *  Created by Martin Rupp on 26.09.12.
 *  Copyright 2012 . All rights reserved.
 *
 */

#ifndef CONNECTION_VIEWER_IO_H
#define CONNECTION_VIEWER_IO_H

#include <fstream>
#include <vector>
#include <string>
#include <cassert>
#include "common/progress.h"
#include "common/math/misc/math_util.h"


namespace ug{

namespace ConnectionViewer{

template<typename matrix_type>
bool ReadMatrix(std::string filename, matrix_type &matrix, std::vector<MathVector<3> > &grid, int &dimension)
{
	std::cout << " Reading matrix from " <<  filename << "... ";
	std::fstream matfile(filename.c_str(), std::ios::in);
	if(matfile.is_open() == false) { std::cout << "failed.\n"; return false; }

	int version=-1, gridsize=-1;
	dimension=-1;

	matfile >> version;
	matfile >> dimension;
	matfile >> gridsize;

	assert(version == 1);
	assert(dimension == 2 || dimension == 3);
	assert(gridsize != -1);

	std::cout << "v=" << version << " dim=" << dimension << " size=" << gridsize;
	PROGRESS_START(prog, gridsize*2, "ReadMatrix " << dimension << "d from " << filename << " , " << gridsize << " x " << gridsize);
	grid.resize(gridsize);
	for(int i=0; i<gridsize; i++)
	{
		if(i%100) { PROGRESS_UPDATE(prog, i); }
		if(matfile.eof())
		{
			assert(0);
			std::cout << " failed.\n";
			return false;
		}

		matfile >> grid[i].x() >> grid[i].y();
		if(dimension==3) matfile >> grid[i].z();
		else grid[i].z() = 0;
	}

	int printStringsInWindow;
	matfile >> printStringsInWindow;

	matrix.resize_and_clear(gridsize, gridsize);
	bool bEOF = matfile.eof();

	while(!bEOF)
	{
		int from, to; double value;
        char c = matfile.peek();
		if(c == -1 || c == 'c' || c == 'v' || matfile.eof())
			break;

		matfile >> from >> to >> value;
		if(value != 0.0)
			matrix(from, to) = value;

		if(from%100) { PROGRESS_UPDATE(prog, gridsize+from); }
		bEOF = matfile.eof();
	}
	matrix.defragment();

	return true;
}

template<typename vector_type>
bool ReadVector(std::string filename, vector_type &vec)
{
    Progress p;
	std::cout << " Reading std::vector from " <<  filename << "... ";
	std::fstream matfile(filename.c_str(), std::ios::in);
	if(matfile.is_open() == false) { std::cout << "failed.\n"; return false; }

	int version=-1, gridsize=-1, dimension=-1;

	matfile >> version;
	matfile >> dimension;
	matfile >> gridsize;

	assert(version == 1);
	assert(dimension == 2 || dimension == 3);
	assert(gridsize != -1);

	PROGRESS_START(prog, gridsize*2, "ReadVector " << dimension << "d from " << filename << " , " << gridsize << " x " << gridsize);
	for(int i=0; i<gridsize; i++)
	{
		if(i%100) { PROGRESS_UPDATE(prog, i); }
		if(matfile.eof())
		{
			std::cout << " failed.\n";
			assert(0);
			return false;
		}
		double x, y, z;
		matfile >> x >> y;
		if(dimension==3) matfile >> z;
	}

	int printStringsInWindow;
	matfile >> printStringsInWindow;

	vec.resize(gridsize);
	bool bEOF = matfile.eof();
	while(!bEOF)
	{
		int from, to; double value;
		char c = matfile.peek();
		if(c == -1 || c == 'c' || c == 'v' || matfile.eof())
			break;

		matfile >> from >> to >> value;
		assert(from == to);
		vec[from] = value;
		if(from%100) { PROGRESS_UPDATE(prog, from); }
		bEOF = matfile.eof();
	}
	return true;
}

/*
template<typename vector_type>
bool WriteVector(std::string filename, const vector_type &vec, const std::vector<MathVector<3> > &grid, int dimension)
{
	std::cout << " Writing std::vector (size " << vec.size() << " to " << filename << "\n";
	assert(vec.size() == grid.size());
	std::fstream f(filename.c_str(), std::ios::out);
	WriteGridHeader(f, grid, dimension);
	for(size_t i=0; i<vec.size(); i++)
		f << i << " " << i << " " << vec[i] << "\n";
	return true;
}


template<typename matrix_type>
bool WriteMatrix(std::string filename, const matrix_type &mat, const std::vector<MathVector<3> > &grid, int dimension)
{
	std::cout << " Writing matrix (size " << mat.num_rows() << " x " << mat.num_cols() << ") to " << filename << "\n";
	assert(mat.num_rows() == mat.num_cols() && mat.num_cols() == grid.size());
	std::fstream f(filename.c_str(), std::ios::out);
	WriteGridHeader(f, grid, dimension);
	for(size_t i=0; i<mat.num_rows(); i++)
	{
		for(typename matrix_type::const_row_iterator it = mat.begin_row(i); it != mat.end_row(i); ++it)
			f << i << " " << index(it) << " " << value(it) << "\n";
	}
	return true;
}

template<typename matrix_type>
bool WriteProlongation(std::string filename, const matrix_type &mat, const std::vector<MathVector<3> > &finegrid, int dimension, std::vector<int> newIndex, int nCoarse)
{
	std::cout << " Writing prolongation matrix (size " << mat.num_rows() << " x " << mat.num_cols() << ") to " << filename << "\n";
	assert(mat.num_rows() == finegrid.size() && mat.num_cols() == nCoarse);
	std::fstream f(filename.c_str(), std::ios::out);
	WriteGridHeader(f, finegrid, dimension);

	std::vector<int> parentIndex(nCoarse, -1);
	for(size_t i=0; i<newIndex.size(); i++)
		if(newIndex[i] != -1)
			parentIndex[newIndex[i]] = i;

	for(size_t i=0; i<mat.num_rows(); i++)
	{
		for(typename matrix_type::const_row_iterator it = mat.begin_row(i); it != mat.end_row(i); ++it)
			f << i << " " << parentIndex[index(it)] << " " << value(it) << "\n";
	}
	return true;
}

template<typename matrix_type>
bool WriteRestriction(std::string filename, const matrix_type &mat, const std::vector<MathVector<3> > &finegrid, int dimension, std::vector<int> newIndex, int nCoarse)
{
	std::cout << " Writing restriction matrix (size " << mat.num_rows() << " x " << mat.num_cols() << ") to " << filename << "\n";
	assert(mat.num_rows() == nCoarse);
	assert(mat.num_cols() == finegrid.size());
	std::fstream f(filename.c_str(), std::ios::out);
	WriteGridHeader(f, finegrid, dimension);

	std::vector<int> parentIndex(nCoarse, -1);
	for(size_t i=0; i<newIndex.size(); i++)
		if(newIndex[i] != -1)
			parentIndex[newIndex[i]] = i;

	for(size_t i=0; i<mat.num_rows(); i++)
	{
		for(typename matrix_type::const_row_iterator it = mat.begin_row(i); it != mat.end_row(i); ++it)
			f << parentIndex[i] << " " << index(it) << " " << value(it) << "\n";
	}
	return true;
}

template<typename matrix_type>
bool WriteGraph(std::string filename, const matrix_type &g, const std::vector<MathVector<3> > &grid, int dimension)
{
	std::cout << " Writing graph (size " << g.size() << ") to " << filename << "\n";
	assert(g.size() == grid.size());
	std::fstream f(filename.c_str(), std::ios::out);
	WriteGridHeader(f, grid, dimension);
	for(size_t i=0; i<g.size(); i++)
	{
		for(typename matrix_type::const_row_iterator it = g.begin_row(i); it != g.end_row(i); ++it)
			f << i << " " << index(it) << " 1\n";
	}
	return true;
}
*/

}
} // namespace ug

#endif // CONNECTION_VIEWER_IO_H
