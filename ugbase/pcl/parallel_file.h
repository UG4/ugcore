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

#ifndef PARALLEL_FILE_H_
#define PARALLEL_FILE_H_


#include "common/util/binary_buffer.h"

#include "pcl_process_communicator.h"

namespace pcl {

/**
 * This function writes a binarybuffers from all participating cores into one combined parallel file.
 * Note that to read these files, you HAVE to use ReadCombinedParallelFile
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
void WriteCombinedParallelFile(ug::BinaryBuffer &buffer,
                               const std::string &strFilename,
                               const ProcessCommunicator &pc = ProcessCommunicator(PCD_WORLD));


/**
 * This function reads a binarybuffers to all participating cores from a combined parallel file,
 * that is one file which contains data for each core.
 * Note that this is not ParallelReadFile, so each core gets DIFFERENT data.
 * It HAS to be used together with WriteCombinedParallelFile
 *
 * NOTE: you have to use this function to do i/o from a lot of cores (1000+),
 * otherwise you will get big i/o problems.
 *
 * @param buffer		a Binary buffer to read data to
 * @param strFilename	the filename
 * @param pc			a processes communicator (default pcl::World)
 */
void ReadCombinedParallelFile(ug::BinaryBuffer &buffer,
                              const std::string& strFilename,
                              const ProcessCommunicator& pc = ProcessCommunicator(PCD_WORLD));

}
#endif