/*
 *  arrayStorage.h
 *  flexamg
 *
 *  Created by Martin Rupp on 18.12.09.
 *  Copyright 2009 G-CSC, University of Frankfurt. All rights reserved.
 *
 *  variable (stl::vector-like) and fixed (T arr[SIZE];-like)
 *  arrays with same interface.
 */

#ifndef __H__UG__MARTIN_ALGEBRA__ARRAYSTORAGE__
#define __H__UG__MARTIN_ALGEBRA__ARRAYSTORAGE__

//#include "misc.h"

namespace ug{

////////////////////////////////////////////////////////////////////////////////
//! array with fixed storage in the struct
/*! This class can be used like variableArray, but has a fixed size.
	Use together with storage_traits:
	storage_traits<storage_type, double, 3, 3>::array_type is either
	fixed Array or variableArray.
 */
template<typename T, size_t n>
class fixedArray
{	
public:
	//! constructor. does nothing
	fixedArray(){}
	//! construct with size
	fixedArray(size_t n_) { UG_ASSERT(n_ == n, "fixed array set to " << n << ", tried " << n_ << "."); }
	
	//fixedArray(const fixedArray &other); // default copy constructor
	
public:
	//! access refernce of an element
	inline T &operator [] (size_t i) { return values[i]; }
	//! access const element
	inline const T &operator [] (size_t i) const { return values[i]; }
	
public:	
	//! compare with other fixedArray
	bool operator == (const fixedArray<T, n> &other) const
	{
		return memcmp(values, other.values, sizeof(T)*n)==0;
	}
	//! compare with other fixedArray
	bool operator != (const fixedArray<T, n> &other) const
	{
		return ! operator == (other);
	}
	
	//! swap with other array
	void swap(fixedArray<T,n> &other)
	{
		std::swap(values, other.values);
	}
	
public:
	//! get size
	inline size_t size() const 
	{ 
		return n;
	}
	
	//! change size (for fixed array trivial)
	inline void resize(size_t s, bool bZero=true)
	{
		UG_ASSERT(s <= n, "fixed Array is too small (max " << n << ", tried " << s << ".");		
		if(bZero)
			memset(values, 0, sizeof(T)*n);
	}
	
private:
	T values[n];	///< data
};

////////////////////////////////////////////////////////////////////////////////
//! array with variableStorage on the heap
/*!
 optimized for small memory consumption (smaller than stl::vector),
 needed because this array can be used millions of times inside sparse matrices
 \sa fixedArray
 */
template<typename T>
class variableArray
{
public:
	//! constructor
	variableArray()
	{
		values = 0;
		n = 0;
	}
	//! deconstrucotr
	~variableArray()
	{
		if(values) delete[] values;
	}
	
	//! construct with size
	variableArray(size_t n_)
	{
		values = 0;
		n = 0;
		resize(n_);
	}
	
	//! copy constructor
	variableArray(const variableArray<T> &other)
	{
		n = other.n;
		values = new T[n];
		memcpy(values, other.values, sizeof(T)*n);
	}
	
public:
	//! access refernce of an element
	inline T &operator [] (size_t i) { return values[i]; }
	//! access const element
	inline const T &operator [] (size_t i) const { return values[i]; }
	
public:
	//! compare with other variableArray
	bool operator == (const variableArray<T> &other) const
	{
		return memcmp(values, other.values, sizeof(T)*n)==0;
	}
	//! compare with other variableArray
	bool operator != (const variableArray<T> &other) const
	{
		return ! operator == (other);
	}
	
	//! assign other variableArray
	void operator = (const variableArray<T> &other)
	{
		if(n != other.n)
		{
			n = other.n;
			if(values) delete[] values;
			values = new T[n];			
		}
		memcpy(values, other.values, sizeof(T)*n);
	}
	
	void swap(variableArray &a)
	{
		std::swap(a.values, values);
		std::swap(a.n, n);
	}
	
public:
	inline size_t size() const 
	{ 
		return n;
	}

	inline void resize(size_t s, bool bZero=true)
	{
		if(s > n)
		{			
			T *m = new T[s];
			if(values)
			{
				memcpy(m, values, sizeof(T)*n);
				delete[] values;
				if(bZero)
					memset(m+n, 0, sizeof(T)*(s-n));
			}
			else if(bZero)
				memset(m, 0, sizeof(T)*s);			
			n = s;
			values = m;
		}		
	}	

private:
	T *values;
	size_t n;
};

////////////////////////////////////////////////////////////////////////////////
//! two-dimensional array with fixed Storage
/*!
 use accessor (r, c) or getAt(r, c) function
 \sa fixedArray
 */
template<typename T, size_t rows, size_t cols>
class fixedArray2
{
public:
	fixedArray2() { memset(values, 0, sizeof(double)*rows*cols); }
	fixedArray2(size_t rows_, size_t cols_) 
	{
		UG_ASSERT(rows_ == rows && cols == cols, "fixed Array! (rows " << rows << ", cols " << cols << " tried " << rows_ << ", " << cols_ << ".)");		
	}
	
	inline T &operator () (size_t r, size_t c) { ensure(r+1, c+1); return values[c + r*cols]; }
	inline const T &operator () (size_t r, size_t c) const { ensure(r+1, c+1); return values[c + r*cols]; }
	inline T &getAt (size_t r, size_t c) { ensure(r+1, c+1); return values[c + r*cols]; }
	inline const T &getAt (size_t r, size_t c) const { ensure(r+1, c+1); return values[c + r*cols]; }
	inline T &operator [] (size_t i) { UG_ASSERT(i<rows*cols && i >= 0, ""); return values[i]; }
	inline const T &operator [] (size_t i)  const { UG_ASSERT(i<rows*cols && i >= 0, ""); return values[i]; }
	size_t size() const { return rows*cols; }

	//! ensure size of rows*cols, but do not make smaller
	inline void ensure(size_t rows_, size_t cols_) const
	{
		UG_ASSERT(rows_ <= rows && cols <= cols, "fixed Array is too small (max rows " << rows << ", cols " << cols << " tried " << rows_ << ", " << cols_ << ".)");		
	}
	//! set size (bZero: if true, zero new mem)
	inline void resize(size_t rows_, size_t cols_, bool bZero=true)
	{
		UG_ASSERT(rows_ == rows && cols == cols, "fixed Array! (rows " << rows << ", cols " << cols << " tried " << rows_ << ", " << cols_ << ".)");		
		memset(values, 0, sizeof(double)*rows*cols);
	}
	
	size_t num_rows() const { return rows; }
	size_t num_cols() const { return cols; }
	
	void swap(fixedArray2<T, rows, cols> &other)
	{
		swap(values, other.values);
	}
	
private:
	T values[rows*cols];	
};

////////////////////////////////////////////////////////////////////////////////
//! two-dimensional array with variable Stroage
/*!
 use accessor (r, c) or getAt(r, c) function
 \sa variableArray
 */
template<typename T>
class variableArray2
{
public:
	variableArray2() : values() { cols = 0; }
	variableArray2(size_t rows_, size_t cols_) : values(rows_*cols_)
	{
		cols = cols_;
	}
	
	variableArray2(const variableArray2<T> &other)
	{
		values = other.values;
		cols = other.cols;
	}
public:
	void ensure(size_t r, size_t c)
	{
		if(c < num_cols())
		{
			if(r >= num_rows())
				values.resize(r*num_cols());
			return;
		}
		else if(r > num_rows() || c > num_cols())
		{
			size_t newc = num_cols(), newr = num_rows();
			if(c >= num_cols())
				newc = c;				
			if(r >= num_rows())
				newr = r;

			variableArray<T> values2(newr*newc);
			for(size_t i=0; i<newr*newc; i++) values2[i] = 0.0;
			// copy old data
			for(size_t ir=0; ir<num_rows(); ir++)
			{
				for(size_t ic=0; ic<num_cols(); ic++)
					values2[ic + ir*newc] = values[ic + ir*cols];				
			}
			values.swap(values2);
			cols = newc;
		}
	}
	
	inline void resize(size_t rows_, size_t cols_, bool bZero=true)
	{
		values.resize(rows_*cols_);
		cols = cols_;
		return;
		if(bZero)
		{
			for(size_t i=0; i<rows_*cols_; i++) values[i] = 0;
	//		memset(&values[0], 0, sizeof(T)*rows_*cols_);
		}
	}
	void swap(variableArray2<T> &other)
	{
		values.swap(other.values);
		std::swap(other.cols, cols);
	}
	
	inline T &operator () (size_t r, size_t c) { ensure(r+1, c+1); return values[c + r*cols]; }
	inline const T &operator () (size_t r, size_t c) const { return values[c + r*cols]; }
	inline T &operator [] (size_t i) { return values[i]; }
	inline const T &operator [] (size_t i) const { return values[i]; }	
	inline size_t num_rows() const {	return cols != 0 ? values.size()/cols : 0; }
	inline size_t num_cols() const { return cols; }
	inline size_t size() const { return values.size(); }
	
	
private:
	variableArray<T> values;
	size_t cols;
};


////////////////////////////////////////////////////////////////////////////////

//! helper class for selecting storage type
class fixedStorage
{
public:
	static const char *getType() { return "fix"; }
};

//! helper class for selecting storage type
class variableStorage 
{
public:
	static const char *getType() { return "var"; }
};

////////////////////////////////////////////////////////////////////////////////

//! storage traits: used to get corresponding one and two-dimensional arrays for a storage class
template<typename storage_type, typename value_type, size_t rows, size_t cols> class storage_traits;

//! fixed Storage traits
template<typename value_type, size_t rows, size_t cols>
class storage_traits<fixedStorage, value_type, rows, cols>
{
public:
	typedef fixedArray<value_type, rows> array_type;
	typedef fixedArray2<value_type, rows, cols> array2_type;	
};

//! variable storage traits
template <typename value_type, size_t rows, size_t cols>
class storage_traits<variableStorage, value_type, rows, cols>
{	
public:
	typedef variableArray<value_type> array_type;
	typedef variableArray2<value_type> array2_type;
};

} // namespace ug

#endif
