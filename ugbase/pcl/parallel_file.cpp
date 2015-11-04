
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

	int mySize = buffer.write_pos();
	int myNextOffset = 0;
	MPI_Scan(&mySize, &myNextOffset, 1, MPI_INT, MPI_SUM, m_mpiComm);


	std::vector<int> allNextOffsets;
	allNextOffsets.resize(pc.size(), 0);
	//else allNextOffsets.resize(1);

	myNextOffset += (pc.size()+1)*sizeof(int);
	MPI_Gather(&myNextOffset, 1, MPI_INT, &allNextOffsets[0], 1, MPI_INT, pc.get_proc_id(0), m_mpiComm);

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

	int myOffset = myNextOffset - buffer.write_pos();
	MPI_File_seek(fh, myOffset, MPI_SEEK_SET);

//	UG_LOG_ALL_PROCS("MySize = " << mySize << "\n" << " myOffset = " << myOffset << "\n");
//	UG_LOG_ALL_PROCS("buffer.write_pos() = " << buffer.write_pos() << "\n" << "(pc.size()+1)*sizeof(size_t) = " << (pc.size()+1)*sizeof(size_t) << "\n");

	MPI_File_write(fh, buffer.buffer(), buffer.write_pos(), MPI_BYTE, &status);

	MPI_File_close(&fh);
}

void ReadCombinedParallelFile(ug::BinaryBuffer &buffer, std::string strFilename, pcl::ProcessCommunicator pc)
{
	MPI_Status status;
	MPI_Comm m_mpiComm = pc.get_mpi_communicator();
	MPI_File fh;

	char filename[1024];
	strcpy(filename, strFilename.c_str());
	if(MPI_File_open(m_mpiComm, filename, MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &fh))
		UG_THROW("could not open "<<filename);

	std::vector<int> allNextOffsets;
	allNextOffsets.resize(pc.size()+1);

	allNextOffsets[0] = (pc.size()+1)*sizeof(int);
	bool bFirst = pc.get_proc_id(0) == pcl::ProcRank();
	if(bFirst)
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
	int myNextOffset, myNextOffset2;
	MPI_Scatter(&allNextOffsets[0], 1, MPI_INT, &myNextOffset, 1, MPI_INT, pc.get_proc_id(0), m_mpiComm);
	MPI_Scatter(&allNextOffsets[1], 1, MPI_INT, &myNextOffset2, 1, MPI_INT, pc.get_proc_id(0), m_mpiComm);

	int mySize = myNextOffset2-myNextOffset;

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
