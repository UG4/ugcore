/*
 * parallel_file.h
 *
 *  Created on: 01.07.2014
 *      Author: mrupp
 */

#ifndef PARALLEL_FILE_H_
#define PARALLEL_FILE_H_

#include "pcl_process_communicator.h"
#include "common/util/binary_buffer.h"

namespace pcl{

/**
 * This function writes a binarybuffers from all participating cores into a parallel file.
 *
 * NOTE: you have to use this function to do i/o from a lot of cores (1000+),
 * otherwise you will get big i/o problems.
 *
 * The file format is as follows:
 *
 * size_t 	numProcs
 * int		nextOffset[numProcs]
 * byte data1[...]
 * byte data2[...]
 * ...
 *
 * That means, in the file the first entry is size_t numProcs, then an array of ints with the offset of
 * the next data set (see more description below), and then the actual data.
 * This function is executed in parallel, so if core 0 has 1024 bytes of data, core 1 has 256 bytes of data, and core 3 has 500 bytes,
 * we have sizeof(size_t) + 3*sizeof(int) = 4*8 = 32 bytes of header, so
 * (size_t) 3
 * (int) 32+1024
 * (int) 32+1024+256
 * (int) 32+1024+256+500
 * 32: data1
 * 32+1024: data2
 * 32+1024+256: data2
 *
 * We store the nextOffset to get access to the size of the data written.
 *
 * @param buffer		a Binary buffer with data
 * @param strFilename	the filename
 * @param pc			a processes communicator (default pcl::World)
 */
void ParallelFileWrite(ug::BinaryBuffer &buffer, std::string strFilename, pcl::ProcessCommunicator pc = pcl::ProcessCommunicator(pcl::PCD_WORLD));


/**
 * This function reads a binarybuffers to all participating cores from a parallel file.
 *
 * NOTE: you have to use this function to do i/o from a lot of cores (1000+),
 * otherwise you will get big i/o problems.
 *
 * @param buffer		a Binary buffer to read data to
 * @param strFilename	the filename
 * @param pc			a processes communicator (default pcl::World)
 */
void ParallelFileRead(ug::BinaryBuffer &buffer, std::string strFilename, pcl::ProcessCommunicator pc = pcl::ProcessCommunicator(pcl::PCD_WORLD));

}
#endif /* PARALLEL_ARCHIVE_H_ */
