/*
 * connection_viewer_output.h
 *
 *  Created on: 03.08.2011
 *      Author: mrupp
 */

#ifndef CONNECTIONVIEWEROUTPUT_H_
#define CONNECTIONVIEWEROUTPUT_H_

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "common/progress.h"
#ifdef UG_PARALLEL
#include "pcl/pcl.h"
#include "pcl/parallel_archive.h"
#endif


#define CONNECTION_VIEWER_VERSION 1

namespace ug
{
namespace ConnectionViewer
{

#ifdef UG_PARALLEL
/**
 * extends the filename (add p000X extension in parallel) and writes a parallel pvec/pmat "header" file
 */
std::string GetParallelName(std::string name, const pcl::ProcessCommunicator &pc, bool bWriteHeader=true);

inline std::string GetParallelName(std::string name)
{
	pcl::ProcessCommunicator pc;
	return GetParallelName(name, pc);

}
template<typename T>
inline std::string GetParallelName(T &t, std::string name)
{
	return GetParallelName(name, t.layouts()->proc_comm());
}
#else
template<typename T>
inline std::string GetParallelName(T &t, std::string name)
{
	return name;
}
#endif


bool AddMarkers(std::string filename, std::string markfilename);
bool WriteMarkers(std::string markfilename, std::vector<bool> markers, float r, float g, float b, float alpha, int size);

template<typename TPosition>
bool WriteGridHeader(std::ostream &f, const TPosition &positions, size_t N, int dimension)
{
	assert(dimension == 2 || dimension == 3);
	f << 1 << "\n";
	f << dimension << "\n";
	f << N << "\n";
	if(dimension == 1)
		for(size_t i=0; i < N; i++)
			f << positions[i][0] << " 0.0\n";
	else if(dimension == 2)
		for(size_t i=0; i < N; i++)
			f << positions[i][0] << " " << positions[i][1] << "\n";
	else
		for(size_t i=0; i < N; i++)
			f << positions[i][0] << " " << positions[i][1] << " " << positions[i][2] << "\n";

	f << 1 << "\n"; // stringsInWindow
	return true;
}

template<typename TPosition>
bool WriteGridHeader(std::ostream &f, const TPosition &positions, int dimension)
{
	return WriteGridHeader(f, positions, positions.size(), dimension);
}

// WriteMatrixToConnectionViewer
//--------------------------------------------------
/**
 * \brief writes to a file in somewhat SparseMatrix-market format (for connection viewer)
 * \param filename Filename to write matrix to
 * \param A SparseMatrix A.
 * \param positions Positions, there has to be one position for each i in (0, ..., max(A.num_rows(), A.num_cols())).
 * \param dimensions Dimension of positions
 */
template<typename Matrix_type, typename postype>
void WriteMatrix(std::ostream &file, const Matrix_type &A, postype *positions, int dimensions)
{
	PROFILE_FUNC_GROUP("debug");
	size_t rows = A.num_rows();

	WriteGridHeader(file, positions, rows, dimensions);
	PROGRESS_START(prog, rows, "WriteMatrixToConnectionViewer " << dimensions << "d, " << rows << "x" << rows);

	// write connections
	for(size_t i=0; i < rows; i++)
	{
		PROGRESS_UPDATE(prog, i);
		for(typename Matrix_type::const_row_iterator conn = A.begin_row(i); conn != A.end_row(i); ++conn)
			if(conn.value() != 0.0)
				file << i << " " << conn.index() << " " << conn.value() <<		std::endl;
			else
				file << i << " " << conn.index() << " 0" <<	std::endl;
	}
}

template<typename Matrix_type, typename postype>
void WriteMatrix(std::string filename, const Matrix_type &A, postype *positions, int dimensions)
{
	std::fstream file(filename.c_str(), std::ios::out);
	WriteMatrix(file, A, positions, dimensions);
}


template<typename Matrix_type, typename postype>
void WriteMatrixPar(std::string name, const Matrix_type &A, const postype *positions, int dimensions)
{
	WriteMatrix(GetParallelName(A, name), A, positions, dimensions);
/*// the new version using parallel archive
#ifndef UG_PARALLEL
	WriteMatrix(GetParallelName(A, filename), A, positions, dimensions);
#else
	const pcl::ProcessCommunicator &pc = A.layouts()->proc_comm();
	if(pcl::GetNumProcesses() == 1)
		WriteMatrix(name, A, positions, dimensions);

	char buf[20];
	int rank = pcl::GetProcRank();


	std::string fname = FilenameWithoutExtension(name);
	std::string ext = GetFilenameExtension(name);

	pcl::ParallelArchive ar(name + ".tar", pc);
	if(rank == pc.get_proc_id(0))
	{
		// create the .pmat "header" file
		std::stringstream &file = ar.create_stringstream_file(fname + ".p" + ext);
		file << pc.size() << "\n";
		for(size_t i=0; i<pc.size(); i++)
		{
			sprintf(buf, "_p%04d.%s", pc.get_proc_id(i), ext.c_str());
			file << fname << buf << "\n";
		}
	}
	// create the .mat file
	sprintf(buf, "_p%04d.%s", rank, ext.c_str());
	WriteMatrix(ar.create_stringstream_file(fname + buf), A, positions, dimensions);
#endif
*/
}


//	NOTE:	The commented version below contains a bug which only occurs in rare situations,
//			resulting in a matrix output, where one entry has connections to
//			many other entries. Most presumably a problem with the mappings.
//			Instead of fixing it, the method was instead replaced by the method below
//			the commented section, which simplifies the output. Since the version with
//			the mapping might still be of interest, it remains in a commented form.
//// WriteMatrix
////--------------------------------------------------
///**
// * this version can handle different from and to spaces
// */
//template <typename Matrix_type, typename postype>
//bool WriteMatrixToConnectionViewerMapped(std::string filename,
//									const Matrix_type &A,
//									std::vector<postype> &positionsFrom, std::vector<postype> &positionsTo, size_t dimensions)
//{
//	PROFILE_FUNC_GROUP("debug");
//#ifdef UG_PARALLEL
//	filename = GetParallelName(A, filename);
//#endif
//
//	/*const char * p = strstr(filename, ".mat");
//	if(p == NULL)
//	{
//		UG_LOG("Currently only '.mat' format supported for domains.\n");
//		return false;
//	}*/
//
//	if(positionsFrom.size() != A.num_cols())
//	{
//		UG_LOG("uFrom.size() != A.num_cols() !\n");
//		return false;
//	}
//	if(positionsTo.size() != A.num_rows())
//	{
//		UG_LOG("uTo.size() != A.num_rows() !\n");
//		return false;
//	}
//
//	std::vector<postype> positions;
//	std::vector<size_t> mapFrom, mapTo;
//	mapFrom.resize(positionsFrom.size());
//	mapTo.resize(positionsTo.size());
//
//	if(positionsFrom.size() > positionsTo.size())
//	{
//		positions.resize(positionsFrom.size());
//		for(size_t i=0; i<positionsFrom.size(); i++)
//		{
//			positions[i] = positionsFrom[i];
//			mapFrom[i] = i;
//		}
//
//
//		for(size_t i=0; i<positionsTo.size(); i++)
//		{
//			size_t j;
//			for(j=0; j<positionsFrom.size(); j++)
//			{
//				if(positionsTo[i] == positionsFrom[j])
//					break;
//			}
//			mapTo[i] = j;
//			if(j == positionsFrom.size())
//				positions.push_back(positionsTo[i]);
//		}
//	}
//	else
//	{
//		positions.resize(positionsTo.size());
//		for(size_t i=0; i<positionsTo.size(); i++)
//		{
//			positions[i] = positionsTo[i];
//			mapTo[i] = i;
//		}
//
//		for(size_t  i=0; i<positionsFrom.size(); i++)
//		{
//			size_t j;
//			for(j=0; j<positionsTo.size(); j++)
//			{
//				if(positionsFrom[i] == positionsTo[j])
//					break;
//			}
//			mapFrom[i] = j;
//			if(j == positionsTo.size())
//				positions.push_back(positionsFrom[i]);
//		}
//	}
//
//
//	std::fstream file(filename.c_str(), std::ios::out);
//	file << CONNECTION_VIEWER_VERSION << std::endl;
//	file << dimensions << std::endl;
//
//	// write positions
//	file << positions.size() << std::endl;
//
//	if(dimensions == 1)
//		for(size_t i=0; i < positions.size(); i++)
//			file << positions[i][0] << " 0.0"  << std::endl;
//	else if(dimensions == 2)
//
//		for(size_t i=0; i < positions.size(); i++)
//			file << positions[i][0] << " " << positions[i][1] << std::endl;
//	else
//		for(size_t i=0; i < positions.size(); i++)
//		  file << positions[i][0] << " " << positions[i][1] << " " << positions[i][2] << std::endl;
//
//	file << 1 << std::endl; // show all cons
//	// write connections
//	for(size_t i=0; i < A.num_rows(); i++)
//	{
//		for(typename Matrix_type::const_row_iterator conn = A.begin_row(i); conn != A.end_row(i); ++conn)
//			if(conn.value() != 0.0)
//				file << mapTo[i] << " " << mapFrom[conn.index()] << " " << conn.value() <<		std::endl;
//			else
//				file << mapTo[i] << " " << mapFrom[conn.index()] << " 0" << std::endl;
//	}
//	return true;
//}

// WriteMatrixToConnectionViewer
//--------------------------------------------------
/**
 * this version can handle different from and to spaces
 */
template <typename Matrix_type, typename postype>
bool WriteMatrix(	std::string filename,
									const Matrix_type &A,
									const std::vector<postype> &positionsFrom,
									const std::vector<postype> &positionsTo, size_t dimensions)
{
	PROFILE_FUNC_GROUP("debug");

	if(positionsFrom.size() != A.num_cols())
	{
		UG_LOG("uFrom.size() != A.num_cols() !\n");
		return false;
	}
	if(positionsTo.size() != A.num_rows())
	{
		UG_LOG("uTo.size() != A.num_rows() !\n");
		return false;
	}

	size_t fromOffset = positionsTo.size();


	std::fstream file(filename.c_str(), std::ios::out);
	file << CONNECTION_VIEWER_VERSION << std::endl;
	file << dimensions << std::endl;

	// write positions
	file << positionsFrom.size() + positionsTo.size() << std::endl;

	if(dimensions == 1){
		for(size_t i=0; i < positionsTo.size(); i++)
			file << positionsTo[i][0] << " 0.0"  << std::endl;
		for(size_t i=0; i < positionsFrom.size(); i++)
			file << positionsFrom[i][0] << " 0.0"  << std::endl;
	}
	else if(dimensions == 2){
		for(size_t i=0; i < positionsTo.size(); i++)
			file << positionsTo[i][0] << " " << positionsTo[i][1] << std::endl;
		for(size_t i=0; i < positionsFrom.size(); i++)
			file << positionsFrom[i][0] << " " << positionsFrom[i][1] << std::endl;
	}
	else{
		for(size_t i=0; i < positionsTo.size(); i++)
		  file << positionsTo[i][0] << " " << positionsTo[i][1] << " " << positionsTo[i][2] << std::endl;
		for(size_t i=0; i < positionsFrom.size(); i++)
		  file << positionsFrom[i][0] << " " << positionsFrom[i][1] << " " << positionsFrom[i][2] << std::endl;
	}

	file << 1 << std::endl; // show all cons
	// write connections

	PROGRESS_START(prog, A.num_rows(), "WriteMatrixToConnectionViewer " << dimensions << "d, " << A.num_rows() << "x" << A.num_rows() );
	for(size_t i=0; i < A.num_rows(); i++)
	{
		PROGRESS_UPDATE(prog, i);
		for(typename Matrix_type::const_row_iterator conn = A.begin_row(i); conn != A.end_row(i); ++conn)
			if(conn.value() != 0.0)
				file << i << " " << conn.index() + fromOffset << " " << conn.value() <<	std::endl;
			else
				file << i << " " << conn.index() + fromOffset << " 0" << std::endl;
	}
	return true;
}


template <typename Matrix_type, typename postype>
bool WriteMatrixPar(std::string filename, const Matrix_type &A,
					const std::vector<postype> &positionsFrom, const std::vector<postype> &positionsTo, size_t dimensions)
{
	return WriteMatrix(GetParallelName(A, filename), A, positionsFrom, positionsTo, dimensions);
}

// WriteVectorToConnectionViewer
//--------------------------------------------------
/**
 * \brief writes to a file in somewhat SparseMatrix-market format (for connection viewer)
 * \param filename Filename to write matrix to
 * \param b Vector
 * \param positions Positions, there has to be one position for each i in (0, ..., max(A.num_rows(), A.num_cols())).
 * \param dimensions	Dimensions of Positions
 */
template<typename Vector_type, typename postype>
void WriteVector(std::string filename, const Vector_type &b, const postype *positions, int dimensions, const Vector_type *compareVec=NULL)
{
	PROFILE_FUNC_GROUP("debug");

	std::fstream file(filename.c_str(), std::ios::out);
	size_t rows = b.size();
	WriteGridHeader(file, positions, rows, dimensions);

	// write connections
	if(compareVec == NULL)
		for(size_t i=0; i < rows; i++)
			file << i << " " << i << " " << b[i] <<		std::endl;
	else
		for(size_t i=0; i < rows; i++)
			file << i << " " << i << " " << b[i]-(*compareVec)[i] << std::endl;
}

template<typename Vector_type, typename postype>
void WriteVectorPar(std::string filename, const Vector_type &b, const postype *positions, int dimensions, const Vector_type *compareVec=NULL)
{
	WriteVector(GetParallelName(b, filename), b, positions, dimensions, compareVec);
}

template<typename Matrix_type, typename Vector_type, typename postype>
void WriteVector(std::string filename, const Matrix_type &A, const Vector_type &v,
		postype *positions, int dimensions, const Vector_type *compareVec=NULL)
{
	PROFILE_FUNC_GROUP("debug");
	if(dimensions != 2)
	{
		UG_LOG(__FILE__ << ":" << __LINE__ << " WriteVectorToConnectionViewer: only dimension=2 supported");
		return;
	}
	filename = GetParallelName(A, filename);

	size_t rows = A.num_rows();
	std::fstream file(filename.c_str(), std::ios::out);
	WriteGridHeader(file, positions, rows, dimensions);

	PROGRESS_START(prog, rows, "WriteVectorToConnectionViewer " << dimensions << "d, " << A.num_rows() << "x" << A.num_rows() );
	// write connections
	for(size_t i=0; i < rows; i++)
	{
		PROGRESS_UPDATE(prog, i);
		for(typename Matrix_type::const_row_iterator conn = A.begin_row(i); conn != A.end_row(i); ++conn)
			if(conn.value() != 0.0)
				file << i << " " << conn.index() << " " << conn.value() <<		std::endl;
			else
				file << i << " " << conn.index() << " 0" <<	std::endl;
	}


	std::string nameValues = filename;
	nameValues.resize(nameValues.find_last_of("."));
	nameValues.append(".values");
	file << "v " << nameValues << "\n";

	std::fstream fileValues(nameValues.c_str(), std::ios::out);
	if(compareVec == NULL)
	{
		for(size_t i=0; i < rows; i++)
			fileValues << i << " " << v[i] <<	"\n";
	}
	else
	{
		typename Vector_type::value_type t;
		for(size_t i=0; i < rows; i++)
		{
			t = v[i];
			t -= (*compareVec)[i];
			fileValues << i << " " << t <<	"\n";
		}
	}

}

template<typename Matrix_type, typename Vector_type, typename postype>
void WriteVectorPar(std::string filename, const Matrix_type &A, const Vector_type &v,
		postype *positions, int dimensions, const Vector_type *compareVec=NULL)
{
	WriteVector(GetParallelName(A, filename), A, v, positions, dimensions, compareVec);
}

#if 0
template<typename Vector_type, typename postype>
void WriteVectorNEW(std::string filename, const Vector_type &b, postype *positions, int dimensions)
{
	PROFILE_FUNC_GROUP("debug");
	filename = GetParallelName(filename);

	std::fstream file(filename.c_str(), std::ios::out);
	file << CONNECTION_VIEWER_VERSION << std::endl;
	file << 3 << std::endl;

	double nmax=0;
	for(size_t i=0; i<b.size(); i++)
		if(nmax < BlockNorm(b[i])) nmax = BlockNorm(b[i]);

	nmax*=4;
	double scale = 1.0/nmax;
	size_t rows = b.size();
	// write positions
	file << rows << std::endl;
	if(dimensions == 1)
		for(size_t i=0; i < rows; i++)
			file << positions[i][0] << " 0.0" << std::endl;
	else if(dimensions == 2)
		for(size_t i=0; i < rows; i++)
			file << positions[i][0] << " " << positions[i][1] << " " << b[i] * scale << std::endl;
	else
		for(size_t i=0; i < rows; i++)
		  file << positions[i][0] << " " << positions[i][1] << " " << positions[i][2] << std::endl;

	file << 1 << std::endl; // show all cons
	// write connections
	for(size_t i=0; i < rows; i++)
	{
		file << i << " " << i << " " << b[i] <<		std::endl;
	}
}
#endif

} // namespace ConnectionViewer
} // namespace ug
#endif /* CONNECTIONVIEWEROUTPUT_H_ */
