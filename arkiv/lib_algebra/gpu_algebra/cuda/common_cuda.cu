#define USE_CUDA
#include "../plugins/experimental/jitsg/jitsg_include.h"

GPU_KERNEL
void CUDA_vector_add_02_kernel (const int num,
        GPUMEM FPTYPE* x, 
        const FPTYPE beta, GPUMEM const FPTYPE *y)
{
    const int idx = GPU_X_INDEX;
	if (idx < num)
        x[idx] = beta*y[idx];
}

GPU_KERNEL
void CUDA_vector_add_2_kernel (const int num,
        const FPTYPE alpha, GPUMEM FPTYPE* x, 
        const FPTYPE beta, GPUMEM const FPTYPE *y)
{
    const int idx = GPU_X_INDEX;
	if (idx < num)
        x[idx] = alpha*x[idx] + beta*y[idx];
}

GPU_KERNEL
void CUDA_vector_add_03_kernel (const int num,
        GPUMEM FPTYPE* x, 
        const FPTYPE beta, GPUMEM const FPTYPE *y, 
        const FPTYPE gamma, GPUMEM const FPTYPE *z)
{
    const int idx = GPU_X_INDEX;
	if (idx < num)
        x[idx] = beta*y[idx] + gamma*z[idx];
}

GPU_KERNEL
void CUDA_vector_add_3_kernel (const int num, 
        const FPTYPE alpha, GPUMEM FPTYPE* x, 
        const FPTYPE beta, GPUMEM const FPTYPE *y, 
        const FPTYPE gamma, GPUMEM const FPTYPE *z)
{
    const int idx = GPU_X_INDEX;
	if (idx < num)
        x[idx] = alpha*x[idx] + beta*y[idx] + gamma*z[idx];
}

GPU_KERNEL
void CUDA_JacobiApply_Kernel (const GPUMEM FPTYPE *diagInv, GPUMEM FPTYPE *corr, const GPUMEM FPTYPE *defect, const int N)
{
	const int idx = GPU_X_INDEX;
	if(idx < N)
		corr[idx] = defect[idx] * diagInv[idx];
}


extern "C" bool 
CUDA_VecAdd2(const int len, FPTYPE alpha, FPTYPE *x, FPTYPE beta, const FPTYPE *y)
{
	if(alpha == 0.0)
    {
        EXECUTE_1D(len, CUDA_vector_add_02_kernel, (len, x, beta, y));
    }
    else
    {
        EXECUTE_1D(len, CUDA_vector_add_2_kernel, (len, alpha, x, beta, y));    
    }
	return true;
}

extern "C" bool 
CUDA_VecAdd3(const int len, FPTYPE alpha, FPTYPE *x, FPTYPE beta, const FPTYPE *y, FPTYPE gamma, const FPTYPE *z)
{   
	if(alpha == 0.0)
    {
        EXECUTE_1D(len, CUDA_vector_add_03_kernel, (len, x, beta, y, gamma, z));
    }
    else
    {
        EXECUTE_1D(len, CUDA_vector_add_3_kernel, (len, alpha, x, beta, y, gamma, z));
    }
	return true;	
}




GPU_KERNEL
void CUDA_VecAdd_2_Kernel (GPUMEM FPTYPE* dest, FPTYPE alpha1, const GPUMEM FPTYPE* v1, FPTYPE alpha2, const  GPUMEM FPTYPE* v2, const int N)
{
    const int idx = GPU_X_INDEX;
	if (idx < N)
        dest[idx] = alpha1*v1[idx] + alpha2*v2[idx];
}

GPU_KERNEL
void CUDA_VecAdd_3_Kernel (GPUMEM FPTYPE* dest, FPTYPE alpha1, const GPUMEM FPTYPE* v1, FPTYPE alpha2, const GPUMEM FPTYPE* v2, FPTYPE alpha3, const GPUMEM FPTYPE* v3, const int N)
{
    const int idx = GPU_X_INDEX;
	if (idx < N)
        dest[idx] = alpha1*v1[idx] + alpha2*v2[idx] + alpha3*v3[idx];
}

extern "C" bool
CUDA_VecAdd_2(FPTYPE *dest, FPTYPE alpha1, const FPTYPE *v1, FPTYPE alpha2, const FPTYPE *v2, const int N)
{
	EXECUTE_1D(N, CUDA_VecAdd_2_Kernel, (dest, alpha1, v1, alpha2, v2, N));
    return true;
}

extern "C" bool
CUDA_VecAdd_3(FPTYPE *dest, FPTYPE alpha1, FPTYPE *v1, FPTYPE alpha2, const FPTYPE *v2, FPTYPE alpha3, const FPTYPE *v3, const int N)
{
	EXECUTE_1D(N, CUDA_VecAdd_3_Kernel, (dest, alpha1, v1, alpha2, v2, alpha3, v3, N));
    return true;
}

extern "C" bool
CUDA_JacobiApply(const FPTYPE *diagInv, FPTYPE *corr, const FPTYPE *defect, const int N)
{
	EXECUTE_1D(N, CUDA_JacobiApply_Kernel, (diagInv, corr, defect, N));
    return true;
}


