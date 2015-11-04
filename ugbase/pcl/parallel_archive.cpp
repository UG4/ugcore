
#include "parallel_archive.h"
#include "common/util/tar.h"

namespace pcl{
using namespace ug;

size_t Get512Padding(size_t s)
{
	return s%512 == 0 ? 0 : (512- (s%512));
}

void WriteParallelArchive(ProcessCommunicator &pc, std::string strFilename, const std::vector<FileBufferDescriptor> &files)
{

	char padding[1024];
	memset(padding, 0, 1024);

	MPI_Offset my_current_offset;

	MPI_Status status;
	MPI_Comm m_mpiComm = pc.get_mpi_communicator();
	MPI_File fh;

	bool bLast = pc.get_local_proc_id()+1 == (int)pc.size();
	bool bFirst = pc.get_proc_id(0) == pcl::ProcRank();


	char filename[1024];
	strcpy(filename, strFilename.c_str());
	MPI_File_open(m_mpiComm, filename,
					MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &fh);


	int mySize = 0;
	if(bLast) mySize += 1024;
	size_t ug4tarLookupSize = sizeof(size_t)*pc.size();
	if(bFirst) mySize += sizeof(TarHeader) + ug4tarLookupSize + Get512Padding(ug4tarLookupSize);


	for(size_t i=0; i<files.size(); i++)
	{
		// byte positions need to be 512-aligned
		size_t s = files[i].size;
		mySize += sizeof(TarHeader) + s + Get512Padding(s);
	}

	int myOffset = 0;
	MPI_Scan(&mySize, &myOffset, 1, MPI_INT, MPI_SUM, m_mpiComm);
	myOffset-=mySize;
//	UG_LOG_ALL_PROCS("MySize = " << mySize << "\n" << "MyOffset = " << myOffset << "\n");


	MPI_File_seek(fh, myOffset, MPI_SEEK_SET);

	std::vector<int> allOffsets;
	if(bFirst) allOffsets.resize(pc.size());
	else allOffsets.resize(1);
	MPI_Gather(&myOffset, 1, MPI_INT,
			&allOffsets[0], 1, MPI_INT, pc.get_proc_id(0), m_mpiComm);

	if(bFirst)
	{
//		UG_LOG_ALL_PROCS("I am first, writing .ug4_tar_lookup_table of size " << ug4tarLookupSize <<"\n");
		TarHeader t;
		t.set_filename(".tar_lookup_table");
		t.set_filesize(ug4tarLookupSize);
		t.set_checksum();

		 // write header
		 MPI_File_write(fh, (void*)&t, sizeof(t), MPI_BYTE, &status);
		 // write file
		 MPI_File_write(fh, (void*)&allOffsets[0], ug4tarLookupSize, MPI_BYTE, &status);
		 // write padding
		 MPI_File_write(fh, (void*)padding, Get512Padding(ug4tarLookupSize), MPI_BYTE, &status);
	}



	for(size_t i=0; i<files.size(); i++)
	{
		 MPI_File_get_position(fh, &my_current_offset);
//		 UG_LOG_ALL_PROCS("Writing file " << files[i].name << " of size " << files[i].size
//				 << " " << 512-files[i].size%512 << " " << " at pos " << my_current_offset << "\n");

		 TarHeader t;
		 t.set_filename(files[i].name);
		 t.set_filesize(files[i].size);
		 t.set_checksum();

		 // write header
		 MPI_File_write(fh, (void*)&t, sizeof(t), MPI_BYTE, &status);
		 // write file
		 MPI_File_write(fh, (void*)files[i].buf, files[i].size, MPI_BYTE, &status);
		 // write padding
		 MPI_File_write(fh, (void*)padding, Get512Padding(files[i].size), MPI_BYTE, &status);
	}


	if(bLast)
		MPI_File_write(fh, (void*)padding, 2*512, MPI_BYTE, &status);

	MPI_File_close(&fh);
}


}
