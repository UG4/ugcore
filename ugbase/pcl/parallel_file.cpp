/*
 * Copyright (c) 2014-2015:  G-CSC, Goethe University Frankfurt
 * Author: Martin Rupp
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

#include "pcl_process_communicator.h"
#include "common/util/binary_buffer.h"
#include "common/log.h"
#include <map>
#include <string>
#include <mpi.h>

namespace pcl{


void WriteCombinedParallelFile(ug::BinaryBuffer &buffer, std::string strFilename, pcl::ProcessCommunicator pc)
{
		MPI_Status status;
	MPI_Comm m_mpiComm = pc.get_mpi_communicator();
	MPI_File fh;
	bool bFirst = pc.get_proc_id(0) == pcl::ProcRank();

	char filename[1024];
	strcpy(filename, strFilename.c_str());

	if(MPI_File_open(m_mpiComm, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh))
		UG_THROW("could not open "<<filename);

	long long mySize = buffer.write_pos();
	long long myNextOffset = 0;
	MPI_Scan(&mySize, &myNextOffset, 1, MPI_LONG_LONG, MPI_SUM, m_mpiComm);


	std::vector<long long> allNextOffsets;
	allNextOffsets.resize(pc.size(), 0);
	//else allNextOffsets.resize(1);

	myNextOffset += (pc.size())*sizeof(long long) + sizeof(int);
	MPI_Gather(&myNextOffset, 1, MPI_LONG_LONG, &allNextOffsets[0], 1, MPI_LONG_LONG, pc.get_proc_id(0), m_mpiComm);

	if(bFirst)
	{
		int numProcs = pcl::NumProcs();
		MPI_File_write(fh, &numProcs, sizeof(numProcs), MPI_BYTE, &status);
		for(size_t i=0; i<allNextOffsets.size(); i++)
		{
//			UG_LOG("allNextOffsets[" << i << "] = " << allNextOffsets[i] << "\n");
			MPI_File_write(fh, &allNextOffsets[i], sizeof(allNextOffsets[i]), MPI_BYTE, &status);
		}
	}

	long long myOffset = myNextOffset - mySize;
	MPI_File_seek(fh, myOffset, MPI_SEEK_SET);

//	UG_LOG_ALL_PROCS("MySize = " << mySize << "\n" << " myOffset = " << myOffset << "\n");
//	UG_LOG_ALL_PROCS("buffer.write_pos() = " << buffer.write_pos() << "\n" << "(pc.size()+1)*sizeof(size_t) = " << (pc.size()+1)*sizeof(size_t) << "\n");

	MPI_File_write(fh, buffer.buffer(), mySize, MPI_BYTE, &status);

	MPI_File_close(&fh);
}

void ReadCombinedParallelFile(ug::BinaryBuffer &buffer, std::string strFilename, pcl::ProcessCommunicator pc)
{
	MPI_Status status;
	MPI_Comm m_mpiComm = pc.get_mpi_communicator();
	MPI_File fh;

	char filename[1024];
	strcpy(filename, strFilename.c_str());
	if(MPI_File_open(m_mpiComm, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh))
		UG_THROW("could not open "<<filename);

	std::vector<long long> allNextOffsets;
	allNextOffsets.resize(pc.size()+1);

	allNextOffsets[0] = (pc.size())*sizeof(long long) + sizeof(int);
	bool bFirst = (pc.get_proc_id(0) == pcl::ProcRank());
	if (bFirst)
	{
		int numProcs;
		MPI_File_read(fh, &numProcs, sizeof(numProcs), MPI_BYTE, &status);
		UG_COND_THROW(numProcs != pcl::NumProcs(), "checkPoint numProcs = " << numProcs << ", but running on " << pcl::NumProcs());

		for(size_t i=1; i<allNextOffsets.size(); i++)
		{
			MPI_File_read(fh, &allNextOffsets[i], sizeof(allNextOffsets[i]), MPI_BYTE, &status);
//			UG_LOG("allNextOffsets[" << i << "] = " << allNextOffsets[i] << "\n");
		}
	}
	long long myNextOffset, myNextOffset2;
	MPI_Scatter(&allNextOffsets[0], 1, MPI_LONG_LONG, &myNextOffset, 1, MPI_LONG_LONG, pc.get_proc_id(0), m_mpiComm);
	MPI_Scatter(&allNextOffsets[1], 1, MPI_LONG_LONG, &myNextOffset2, 1, MPI_LONG_LONG, pc.get_proc_id(0), m_mpiComm);

	long long mySize = myNextOffset2-myNextOffset;

//	UG_LOG_ALL_PROCS("MySize = " << mySize << "\n" << "myNextOffset = " << myNextOffset << " - " << myNextOffset2 << "\n");

	MPI_File_seek(fh, myNextOffset, MPI_SEEK_SET);

	char *p = new char[mySize];
	MPI_File_read(fh, p, mySize, MPI_BYTE, &status);
	buffer.clear();
	buffer.reserve(mySize);
	buffer.write(p, mySize);
	delete[] p;

	MPI_File_close(&fh);
	//	UG_LOG("File read.\n");
}

}
