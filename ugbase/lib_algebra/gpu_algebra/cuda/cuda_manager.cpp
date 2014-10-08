/* 
 * File:   CUDAManager.cpp
 * Author: mrupp
 * 
 * Created on 16. Oktober 2012, 13:07
 */
//#ifdef USE_CUDA
#include "cuda_manager.h"
#include <iostream>

using namespace std;

namespace ug{

DebugID DID_CUDA("CUDA");


CUDAManager::~CUDAManager()
{
	if(m_tempBuffer)
		cudaFree(m_tempBuffer);
	cudaFree(m_tempRetBuffer);

#ifdef USE_CUSPARSE
	if(cusparseHandle) cusparseDestroy(cusparseHandle);
#endif
    if(cublasHandle) cublasDestroy(cublasHandle);
	cudaDeviceReset();    
	cout << "Cleaned up CUDA." << endl;
}


void CUDAManager::init()
{
	//cudaDeviceReset();
	
     // This will pick the best possible CUDA capable device
    cudaDeviceProp deviceProp;
    int devID;
	int device_count;
	cudaGetDeviceCount(&device_count);
	cout << device_count << " CUDA devices.\n";
    devID = gpuGetMaxGflopsDeviceId();
    checkCudaErrors(cudaSetDevice(devID));
	cudaSetDevice(devID);

    if (devID < 0)
    {
    	UG_THROW("no CUDA device found.\n");
    }

    CUDA_CHECK_STATUS(cudaGetDeviceProperties(&deviceProp, devID));


    // Statistics about the GPU device
    printf("> GPU device has %d Multi-Processors, SM %d.%d compute capabilities\n\n", deviceProp.multiProcessorCount, deviceProp.major, deviceProp.minor);
    int version = (deviceProp.major * 0x10 + deviceProp.minor);
    if (version < 0x11)
    {

        cudaDeviceReset();
        UG_THROW("Requires a minimum CUDA compute 1.1 capability\n");
    }
    
    m_maxThreadsPerBlock = deviceProp.maxThreadsPerBlock;
    UG_DLOG(DID_CUDA, 0, "CUDA Initialized:"
    		"\n - CUDA Device '" << deviceProp.name << "': " <<
            "\n - Total global Memory: " << deviceProp.totalGlobalMem/(1024*1024*1024.0) << " GB"
            "\n - Shared Mem per Block: " << deviceProp.sharedMemPerBlock/(1024.0) << " KB"
            "\n - Regs per Block: " << deviceProp.regsPerBlock <<
            "\n - Warp Size: " << deviceProp.warpSize <<
            "\n - Maximum Number of Threads per Block: " << deviceProp.maxThreadsPerBlock <<
            "\n - Max Thread Dim: (" << deviceProp.maxThreadsDim[0] << ", " << deviceProp.maxThreadsDim[1] << ", " << deviceProp.maxThreadsDim[2] << ")" <<
            "\n - Max Grid Size: (" << deviceProp.maxGridSize[0] << ", " << deviceProp.maxGridSize[1] << ", " << deviceProp.maxGridSize[2] << ")" <<
            "\n - Clock Rate: " << deviceProp.clockRate/1000.0 << " Mhz"
            "\n - Total Const Mem: " << deviceProp.totalConstMem/(1024.0) << " KB"
            "\n - Compute Capability: " << deviceProp.major << "." << deviceProp.minor <<
            "\n - Number of multiprocessors: " << deviceProp.multiProcessorCount <<
            "\n - Maximum Texture Size 1D: " << deviceProp.maxTexture1D <<
            "\n - Maximum Texture Size 2D: " << deviceProp.maxTexture2D[0] << " x " << deviceProp.maxTexture2D[1] <<
            "\n - Maximum Texture Size 3D: " << deviceProp.maxTexture3D[0] << " x " << deviceProp.maxTexture3D[1] << " x " << deviceProp.maxTexture3D[2] <<
            "\n - Memory Clock Rate: " << deviceProp.memoryClockRate/(1000000.0)<< " Mhz"
            "\n - Memory Bus Width: " << deviceProp.memoryBusWidth <<
            "\n - L2 Cache Size: " << deviceProp.l2CacheSize <<
            "\n - Max Threads per multiprocessor: " << deviceProp.maxThreadsPerMultiProcessor <<
            "\n");
            
            
      /* Get handle to the CUBLAS context */
    cublasHandle = 0;
    cublasStatus_t cublasStatus = cublasCreate(&cublasHandle);
    CUDA_CHECK_STATUS(cublasStatus);

#ifdef USE_CUSPARSE
    /* Get handle to the CUSPARSE context */
    cusparseHandle = 0;
    cusparseStatus_t cusparseStatus = cusparseCreate(&cusparseHandle);
    CUDA_CHECK_STATUS(cusparseStatus);

    cusparseMatDescr_t descr = 0;
    cusparseStatus = cusparseCreateMatDescr(&descr);

    CUDA_CHECK_STATUS(cusparseStatus);
#endif
    get_temp_buffer<double>(1024);
    m_tempRetBuffer = MyCudaAlloc<double>(4);
	
}
    
CUDAManager &CUDAManager::get_instance()
{
    static bool bInited = false;
    static CUDAManager ch;
    if(bInited==false)
    {
        ch.init();
        bInited = true;
    }
    return ch;
}
  
}

//#endif
