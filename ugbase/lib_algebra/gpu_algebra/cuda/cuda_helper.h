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
#include "common/log.h"

#include <string>
namespace ug{
std::string CUDAError(int err);

template<typename T>
inline void CudaCheckStatus(T status, const char * file, int line)
{
	unsigned int s = static_cast<unsigned int>(status );
	UG_COND_THROW(status != 0, "CUDA error at " << file << ":" << line << " " << s << " = " << _cudaGetErrorEnum(status) );
}

#define CUDA_CHECK_STATUS(status ) CudaCheckStatus(status, __FILE__, __LINE__)

#define CUDA_CHECK_SUCCESS(err, desc) \
if(err != cudaSuccess)\
{\
    UG_THROW("Error in " << __FUNCTION__ << ": CUDA ERROR " << err <<": " <<\
            cudaGetErrorString(err) << "\n" << desc << "\n");\
}


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
	CUDA_CHECK_SUCCESS( cudaMemcpy(dest, &vec[0], vec.size()*sizeof(typename T::value_type), cudaMemcpyHostToDevice),
			"cudaMemcpy vec size " << vec.size());
}

template<typename T>
inline void CudaCpyToHost(T &dest, typename T::value_type *src)
{
	UG_LOG("Copying " << dest.size() << " to host\n");
	//std::cout << "copy!\n";
	CUDA_CHECK_SUCCESS( cudaMemcpy(&dest[0], src, dest.size()*sizeof(typename T::value_type), cudaMemcpyDeviceToHost),
			"cudaMemcpy dest size " << dest.size())
}


template<typename T>
inline typename T::value_type *CudaCreateAndCopyToDevice(T &vec)
{
	UG_LOG("Create and Copying " << vec.size() << " to host\n");
	typename T::value_type *dest;
	int N = vec.size()*sizeof(typename T::value_type);
	CUDA_CHECK_SUCCESS( cudaMalloc((void **)&dest, N),
			"Error at cudaMalloc of " << N << " bytes");

	CudaCpyToDevice(dest, vec);
	return dest;
}


}
#endif	/* CUDAHELPER_H */
