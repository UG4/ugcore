/* 
 * File:   common_cuda.h
 * Author: mrupp
 *
 * Created on 18. Januar 2013, 17:18
 */

#ifndef COMMON_CUDA_H
#define	COMMON_CUDA_H

#define FPTYPE double

extern "C"
bool
CUDA_VecAdd2(const int len, FPTYPE alpha, FPTYPE *x, FPTYPE beta, const FPTYPE *y);

extern "C"
bool 
CUDA_VecAdd3(const int len, FPTYPE alpha,  FPTYPE *x, FPTYPE beta, const FPTYPE *y, FPTYPE gamma, const FPTYPE *z);

extern "C" bool
CUDA_VecAdd_2(FPTYPE *dest, FPTYPE alpha1, const FPTYPE *v1, FPTYPE alpha2, const FPTYPE *v2, const int N);

extern "C" bool
CUDA_VecAdd_3(FPTYPE *dest, FPTYPE alpha1, const FPTYPE *v1, FPTYPE alpha2, const FPTYPE *v2, FPTYPE alpha3, const FPTYPE *v3, const int N);

#endif	/* COMMON_CUDA_H */

