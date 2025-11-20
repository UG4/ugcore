/*
 * Copyright (c) 2009-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__CRS_ALGEBRA__VECTOR__
#define __H__UG__CRS_ALGEBRA__VECTOR__

#include "../cpu_algebra/vector.h"
#include "cuda/cuda_manager.h"
#include "cuda/common_cuda.h"

namespace ug{
///////////////////////////////////////////////////////////////////
//							GPUVector
///////////////////////////////////////////////////////////////////

/// \addtogroup crs_algebra
/// \{

//!
template <typename TValueType>
class GPUVector : public Vector<TValueType>
{
public:

	using value_type = TValueType ;
	using vector_type = GPUVector<TValueType> ;

	using super = Vector<TValueType>;
	using super::size;
	using super::resize;
	using super::reserve;


	//! constructor
	GPUVector() : Vector<TValueType>() {m_GPUState = ON_CPU; }

	//! constructor with length
	GPUVector(size_t _length) : Vector<TValueType>(_length) {m_GPUState = ON_CPU; }

	//! clones the vector (deep-copy) including values
	SmartPtr<vector_type> clone() const;

	//! clones the vector (deep-copy) excluding values
	SmartPtr<vector_type> clone_without_values() const;

	void resize(size_t newSize, bool bCopyValues=true)
	{
		UG_LOG(this << "GPUVector::resize(" << newSize << ")\n");
		assure_on_cpu();
		m_GPUState = ON_CPU;
		super::resize(newSize, bCopyValues);
	}
	void reserve(size_t newCapacity, bool bCopyValues=true)
	{
		UG_LOG(this << "GPUVector::reserve(" << newCapacity << ")\n");
		reserve(newCapacity, bCopyValues);
	}


	//! access element i of the vector
	inline value_type &operator [] (size_t i)
	{
		assure_on_cpu();
		m_GPUState = ON_CPU;
		return super::operator [] (i);
	}
	inline const value_type &operator [] (size_t i) const
	{
		assure_on_cpu();
		return super::operator [] (i);
	}


protected:
	//! virtual clone using covariant return type
	virtual vector_type* virtual_clone() const;

	//! virtual clone using covariant return type excluding values
	virtual vector_type* virtual_clone_without_values() const;

public:
	void assure_on_gpu()
	{
		if(on_gpu()) return;
		if(m_sizeOnGPU != size())
		{
			cudaFree(m_devValues);
			m_devValues = CudaCreateAndCopyToDevice(*this);
		}
		else
			CudaCpyToDevice(m_devValues, *this);
		m_GPUState = m_GPUState | ON_GPU;
	}

	void assure_on_gpu() const
	{
		const_cast<GPUVector<value_type> *>(this)->assure_on_gpu();
	}

	void assure_on_cpu()
	{
		if(on_cpu()) return;
		// do this before so CudaCpyToHost can access [0] as dest without
		// calling assure_on_cpu again.
		m_GPUState = m_GPUState | ON_CPU;
		CudaCpyToHost(*this, m_devValues);
	}

	void assure_on_cpu() const
	{
		const_cast<GPUVector<value_type> *>(this)->assure_on_cpu();
	}

	bool on_cpu()
	{
		return m_GPUState & ON_CPU;
	}

	bool on_gpu()
	{
		return m_GPUState & ON_GPU;
	}

private:
	enum GPU_STATE
	{
		ON_GPU = 1,
		ON_CPU = 2,
		ON_GPU_AND_CPU = 3
	};
	int m_GPUState;

public:
	double *get_dev_ptr()
	{
		assure_on_gpu();
		m_GPUState = ON_GPU; // not valid on CPU anymore
		return m_devValues;
	}
	const double *get_dev_ptr() const
	{
		assure_on_gpu();
		return m_devValues;
	}

public:
	inline void operator = (const GPUVector<value_type> &v)
	{
		CUDA_VecAdd2(size(), 0.0, get_dev_ptr(), 1.0, v.get_dev_ptr());
	}
	inline void operator += (const GPUVector<value_type> &v)
	{
		CUDA_VecAdd2(size(), 1.0, get_dev_ptr(), 1.0, v.get_dev_ptr());
	}
	inline void operator -= (const GPUVector<value_type> &v)
	{
		CUDA_VecAdd2(size(), 1.0, get_dev_ptr(), -1.0, v.get_dev_ptr());
	}

	inline void add(double alpha, const GPUVector<value_type> &v)
	{
		CUDA_VecAdd2(size(), 1.0, get_dev_ptr(), alpha, v.get_dev_ptr());
	}

	inline void operator *= (const number &a)
	{
		CUDA_VecAdd2(size(), 0.0, get_dev_ptr(), a, get_dev_ptr());
	}

	//! return sqrt(sum values[i]^2) (euclidian norm)
	inline double norm() const
	{
		double res=0;
		cublasDnrm2(CUDAHelper::get_cublasHandle(), size(), get_dev_ptr(), 1, &res);
		return res;
	}

	double dotprod(const GPUVector<value_type> &w) const
	{
		assert(size() == w.size());
		double res=0;
		cublasDdot(CUDAHelper::get_cublasHandle(), size(), get_dev_ptr(), 1, w.get_dev_ptr(), 1, &res);
		cudaThreadSynchronize();
		return res;
	}

private:
	double *m_devValues;
	size_t m_sizeOnGPU;
};

template<typename value_type>
GPUVector<value_type>* GPUVector<value_type>::virtual_clone() const
{
	return new GPUVector<value_type>(*this);
}

template<typename value_type>
SmartPtr<GPUVector<value_type> > GPUVector<value_type>::clone() const
{
	return SmartPtr<GPUVector<value_type> >(this->virtual_clone());
}

template<typename value_type>
GPUVector<value_type>* GPUVector<value_type>::virtual_clone_without_values() const
{
	return new GPUVector<value_type>(this->size());
}

template<typename value_type>
SmartPtr<GPUVector<value_type> > GPUVector<value_type>::clone_without_values() const
{
	return SmartPtr<GPUVector<value_type> >(this->virtual_clone_without_values());
}


/////////////////////////////////////////////////////////


// templated

// operations for vectors
//-----------------------------------------------------------------------------
// these functions execute vector operations by using the operations on the elements of the vector

// todo: change vector_t to TE_VEC<vector_t>


// VecScale: These function calculate dest = sum_i alpha_i v_i

//! calculates dest = alpha1*v1
template<typename T>
inline void VecScaleAssign(GPUVector<T> &dest, double alpha1, const GPUVector<T> &v1)
{
	UG_LOG("VecScaleAssign\n");
	for(size_t i=0; i<dest.size(); i++)
		VecScaleAssign(dest[i], alpha1, v1[i]);
}

//! sets dest = v1 entrywise
template<typename T>
inline void VecAssign(GPUVector<T> &dest, const GPUVector<T> &v1)
{
	UG_LOG("VecAssign\n");
	for(size_t i=0; i<dest.size(); i++)
		dest[i] = v1[i];
}

//! calculates dest = alpha1*v1 + alpha2*v2
template<typename T>
inline void VecScaleAdd(GPUVector<T> &dest, double alpha1, const GPUVector<T> &v1, double alpha2, const GPUVector<T> &v2)
{
	CUDA_VecAdd_2(dest.get_dev_ptr(), alpha1, v1.get_dev_ptr(), alpha2, v2.get_dev_ptr(), dest.size());
}

//! calculates dest = alpha1*v1 + alpha2*v2 + alpha3*v3
template<typename T>
inline void VecScaleAdd(GPUVector<T> &dest, double alpha1, const GPUVector<T> &v1, double alpha2, const GPUVector<T> &v2, double alpha3, const GPUVector<T> &v3)
{
	CUDA_VecAdd_3(dest.get_dev_ptr(), alpha1, v1.get_dev_ptr(), alpha2, v2.get_dev_ptr(), alpha3, v3.get_dev_ptr(), dest.size());
}


// VecProd

//! calculates s += scal<a, b>
template<typename T>
inline void VecProd(const GPUVector<T> &v1, const GPUVector<T> &v2, double &res)
{
//	UG_LOG("VecProd\n");
	assert(v1.size() == v2.size());
	cublasDdot(CUDAHelper::get_cublasHandle(), v1.size(), v1.get_dev_ptr(), 1, v2.get_dev_ptr(), 1, &res);
	cudaThreadSynchronize();
}

//! returns scal<a, b>
template<typename T>
inline double VecProd(const GPUVector<T> &v1, const GPUVector<T> &v2)
{
//	UG_LOG("VecProd\n");
	double res = 0;
	VecProd(v1, v2, res);
	return res;
}


//! calculates s += norm_2^2(a)
template<typename T>
inline void VecNormSquaredAdd(const GPUVector<T> &a, const GPUVector<T> &b, double &sum)
{
	UG_LOG("VecNormSA\n");
	for(int i=0; i<a.size(); i++)
		VecNormSquaredAdd(a[i], sum);
}

//! returns norm_2^2(a)
template<typename T>
inline double VecNormSquared(const GPUVector<T> &a, const GPUVector<T> &b)
{
	UG_LOG("VecNormS\n");
	double sum=0;
	VecNormSquaredAdd(a, sum);
	return sum;
}










// end group crs_algebra
/// \}

} // namespace ug

#endif
