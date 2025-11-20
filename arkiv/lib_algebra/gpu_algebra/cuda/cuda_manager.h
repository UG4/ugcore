/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef CUDAManager_H
#define	CUDAManager_H

#define USE_CUSPARSE

/* Using updated (v2) interfaces to cublas and cusparse */
#include <cuda_runtime.h>

#ifdef USE_CUSPARSE
#include <cusparse_v2.h>
#endif

#include <cublas_v2.h>
#include <vector>

// Utilities and system includes
#include "common/error.h"
#include "common/log.h"

#include <string>

namespace ug{
extern DebugID DID_CUDA;

const char *CUDAError(int err);

template<typename T>
inline void CudaCheckStatus(T status, const char * file, int line)
{
	unsigned int s = static_cast<unsigned int>(status );
	UG_COND_THROW(status != 0, "CUDA error at " << file << ":" << line << " " << s << " = " << ug::CUDAError(status) );
}


#define CUDA_CHECK_STATUS(status ) CudaCheckStatus(status, __FILE__, __LINE__)

#define CUDA_CHECK_SUCCESS(err, desc) \
if(err != cudaSuccess)\
{\
    UG_THROW("Error in " << __FUNCTION__ << ": CUDA ERROR " << err <<":\n" <<\
    		ug::CUDAError(err) << "\n----------------------------\n" << desc << "\n");\
}


template<typename T>
T *MyCudaAlloc(size_t N)
{
	UG_DLOG(DID_CUDA, 2, "CUDA: Allocating " <<  sizeof(T)*N  << " bytes.\n");

	T *p;
	cudaError_t err = cudaMalloc ((void**) &p, sizeof(T)*N);
	if(err != cudaSuccess)
	{
		UG_THROW("Error in " << __FUNCTION__ << "when allocating " << sizeof(T)*N << " bytes. CUDA ERROR " << err <<": " <<
				ug::CUDAError(err));
	}
	return p;
}


class CUDAManager
{
public:
    virtual ~CUDAManager();
    void init();
    static CUDAManager &get_instance();
    
#ifdef USE_CUSPARSE
public:
    static inline cusparseHandle_t get_cusparseHandle() { return get_instance().cusparseHandle; }
private:
    cusparseHandle_t cusparseHandle;
#endif
    
public:
    static inline cublasHandle_t get_cublasHandle() { return get_instance().cublasHandle; }
    size_t m_maxThreadsPerBlock;

    template<typename T>
    T *get_temp_buffer(size_t n)
    {
    	size_t N = n*sizeof(T);
    	if(N < m_tempSize) return (T*)m_tempBuffer;

    	UG_DLOG(DID_CUDA, 2, "CUDA: Allocating Temp Buffer " <<  N  << " bytes.\n");
    	if(m_tempBuffer)
    		cudaFree(m_tempBuffer);

    	m_tempBuffer = MyCudaAlloc<char>(n);

    	return (T*)m_tempBuffer;
    }

    template<typename T>
    T *get_temp_return_buffer()
	{
    	return (T*)m_tempRetBuffer;
	}

    static void get_cuda_devices(std::vector<cudaDeviceProp> &devices);
    static int get_max_multiprocessor_cuda_device();

private:    
    cublasHandle_t cublasHandle;
    void *m_tempBuffer;
    void *m_tempRetBuffer;
    size_t m_tempSize;
};


template<typename T>
inline void CudaCpyToDevice(typename T::value_type *dest, T &vec)
{
	UG_DLOG(DID_CUDA, 2, "Copying " << vec.size() << " to device\n");
	//std::cout << "copy!\n";
	CUDA_CHECK_SUCCESS( cudaMemcpy(dest, &vec[0], vec.size()*sizeof(typename T::value_type), cudaMemcpyHostToDevice),
			"cudaMemcpy vec size " << vec.size());
}

template<typename T>
inline void CudaCpyToHost(T &dest, typename T::value_type *src)
{
	UG_DLOG(DID_CUDA, 2, "Copying " << dest.size() << " to host\n");
	//std::cout << "copy!\n";
	CUDA_CHECK_SUCCESS( cudaMemcpy(&dest[0], src, dest.size()*sizeof(typename T::value_type), cudaMemcpyDeviceToHost),
			"cudaMemcpy dest size " << dest.size())
}


template<typename T>
inline typename T::value_type *CudaCreateAndCopyToDevice(T &vec)
{
	UG_DLOG(DID_CUDA, 2, "Create and Copying " << vec.size() << " to host\n");
	typename T::value_type *dest;
	int N = vec.size()*sizeof(typename T::value_type);
	CUDA_CHECK_SUCCESS( cudaMalloc((void **)&dest, N),
			"Error at cudaMalloc of " << N << " bytes");

	CudaCpyToDevice(dest, vec);
	return dest;
}

template<typename T>
T CUDA_GetElementFromDevice(T *p, size_t i=0)
{
	T t;
	cudaMemcpy(&t, p+i, sizeof(T), cudaMemcpyDeviceToHost);
	return t;
}
}
#endif