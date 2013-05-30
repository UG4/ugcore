/* 
 * File:   CUDAHelper.h
 * Author: mrupp
 *
 * Created on 16. Oktober 2012, 13:07
 */

#ifndef CUDAHELPER_H
#define	CUDAHELPER_H

#define USE_CUSPARSE

/* Using updated (v2) interfaces to cublas and cusparse */
#include <cuda_runtime.h>

#ifdef USE_CUSPARSE
#include <cusparse_v2.h>
#endif

#include <cublas_v2.h>
#include <vector>

// Utilities and system includes
#include "helper_cuda.h"
#include "common/error.h"

#include <string>
namespace ug{
std::string CUDAError(int err);

class CUDAHelper
{
public:
    virtual ~CUDAHelper();
    void init();
    static CUDAHelper &get_instance();
    
#ifdef USE_CUSPARSE
public:
    static inline cusparseHandle_t get_cusparseHandle() { return get_instance().cusparseHandle; }
private:
    cusparseHandle_t cusparseHandle;
#endif
    
public:
    static inline cublasHandle_t get_cublasHandle() { return get_instance().cublasHandle; }
    size_t m_maxThreadsPerBlock;

private:    
    cublasHandle_t cublasHandle;
};


template<typename T>
inline void CudaCpyToDevice(typename T::value_type *dest, T &vec)
{
	UG_LOG("Copying " << vec.size() << " to device\n");
	//std::cout << "copy!\n";
    cudaError_t err = cudaMemcpy(dest, &vec[0], vec.size()*sizeof(typename T::value_type), cudaMemcpyHostToDevice);
    if(err != cudaSuccess)
	{
        UG_THROW("Error in " << __FUNCTION__ << ": CUDA ERROR " << err <<": " <<
                cudaGetErrorString(err) << "\n");
	}
}

template<typename T>
inline void CudaCpyToHost(T &dest, typename T::value_type *src)
{
	UG_LOG("Copying " << dest.size() << " to host\n");
	//std::cout << "copy!\n";
    cudaError_t err = cudaMemcpy(&dest[0], src, dest.size()*sizeof(typename T::value_type), cudaMemcpyDeviceToHost);
    if(err != cudaSuccess)
	{
        UG_THROW("Error in " << __FUNCTION__ << ": CUDA ERROR " << err <<": " <<
                cudaGetErrorString(err) << "\n");
	}
}


template<typename T>
inline typename T::value_type *CudaCreateAndCopyToDevice(T &vec)
{
	UG_LOG("Create and Copying " << vec.size() << " to host\n");
	typename T::value_type *dest;
	int N = vec.size()*sizeof(typename T::value_type);
	cudaError_t err = cudaMalloc((void **)&dest, N);
//	UG_LOG("CudaCreateAndCopyToDevice, cudaMalloc(" << N << ") = " << dest << "\n");
	//UG_LOG("cudaAlloc " << vec.size() << "\n");
	if(err != cudaSuccess)
	{
        UG_THROW("Error at cudaMalloc of " << N << " bytes: CUDA ERROR " << err <<": " <<
                cudaGetErrorString(err) << "\n");
	}

	CudaCpyToDevice(dest, vec);
	return dest;
}


}
#endif	/* CUDAHELPER_H */
