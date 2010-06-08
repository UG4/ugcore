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
#pragma once
//#include "misc.h"

namespace ug{

////////////////////////////////////////////////////////////////////////////////
//! array with fixed storage in the struct
/*! This class can be used like variableArray, but has a fixed size.
	Use together with storage_traits:
	storage_traits<storage_type, double, 3, 3>::array_type is either
	fixed Array or variableArray.
 */
template<typename T, int n>
class fixedArray
{	
public:
	//! constructor. does nothing
	fixedArray(){}
	//! construct with size
	fixedArray(int n_) { UG_ASSERT(n_ == n, "fixed array set to " << n << ", tried " << n_ << "."); }
	
	//fixedArray(const fixedArray &other); // default copy constructor
	
public:
	//! access refernce of an element
	inline T &operator [] (int i) { return values[i]; }
	//! access const element
	inline const T &operator [] (int i) const { return values[i]; }
	
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
	inline void setSize(int s, bool bZero=true)
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
	variableArray(int n_)
	{
		values = 0;
		n = 0;
		setSize(n_);
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
	inline T &operator [] (int i) { return values[i]; }
	//! access const element
	inline const T &operator [] (int i) const { return values[i]; }
	
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

	inline void setSize(int s, bool bZero=true)
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
	int n;
};

////////////////////////////////////////////////////////////////////////////////
//! two-dimensional array with fixed Storage
/*!
 use accessor (r, c) or getAt(r, c) function
 \sa fixedArray
 */
template<typename T, int rows, int cols>
class fixedArray2
{
public:
	fixedArray2() { memset(values, 0, sizeof(double)*rows*cols); }
	fixedArray2(int rows_, int cols_) 
	{
		UG_ASSERT(rows_ == rows && cols == cols, "fixed Array! (rows " << rows << ", cols " << cols << " tried " << rows_ << ", " << cols_ << ".)");		
	}
	
	inline T &operator () (int r, int c) { ensure(r+1, c+1); return values[c + r*cols]; }
	inline const T &operator () (int r, int c) const { ensure(r+1, c+1); return values[c + r*cols]; }
	inline T &getAt (int r, int c) { ensure(r+1, c+1); return values[c + r*cols]; }
	inline const T &getAt (int r, int c) const { ensure(r+1, c+1); return values[c + r*cols]; }
	inline T &operator [] (int i) { UG_ASSERT(i<rows*cols && i >= 0, ""); return values[i]; }
	inline const T &operator [] (int i)  const { UG_ASSERT(i<rows*cols && i >= 0, ""); return values[i]; }
	int size() const { return rows*cols; }

	//! ensure size of rows*cols, but do not make smaller
	inline void ensure(int rows_, int cols_) const
	{
		UG_ASSERT(rows_ <= rows && cols <= cols, "fixed Array is too small (max rows " << rows << ", cols " << cols << " tried " << rows_ << ", " << cols_ << ".)");		
	}
	//! set size (bZero: if true, zero new mem)
	inline void setSize(int rows_, int cols_, bool bZero=true)
	{
		UG_ASSERT(rows_ == rows && cols == cols, "fixed Array! (rows " << rows << ", cols " << cols << " tried " << rows_ << ", " << cols_ << ".)");		
		memset(values, 0, sizeof(double)*rows*cols);
	}
	
	int getRows() const { return rows; }
	int getCols() const { return cols; }
	
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
	variableArray2(int rows_, int cols_) : values(rows_*cols_)
	{
		cols = cols_;
	}
	
	variableArray2(const variableArray2<T> &other)
	{
		values = other.values;
		cols = other.cols;
	}
public:
	void ensure(int r, int c)
	{
		if(c < getCols())
		{
			if(r >= getRows())
				values.setSize(r*getCols());
			return;
		}
		else if(r > getRows() || c > getCols())
		{
			int newc = getCols(), newr = getRows();
			if(c >= getCols())
				newc = c;				
			if(r >= getRows())
				newr = r;

			variableArray<T> values2(newr*newc);
			for(int i=0; i<newr*newc; i++) values2[i] = 0.0;
			// copy old data
			for(int ir=0; ir<getRows(); ir++)
			{
				for(int ic=0; ic<getCols(); ic++)
					values2[ic + ir*newc] = values[ic + ir*cols];				
			}
			values.swap(values2);
			cols = newc;
		}
	}
	
	inline void setSize(int rows_, int cols_, bool bZero=true)
	{
		values.setSize(rows_*cols_);
		cols = cols_;
		return;
		if(bZero)
		{
			for(int i=0; i<rows_*cols_; i++) values[i] = 0;
	//		memset(&values[0], 0, sizeof(T)*rows_*cols_);
		}
	}
	void swap(variableArray2<T> &other)
	{
		values.swap(other.values);
		std::swap(other.cols, cols);
	}
	
	inline T &operator () (int r, int c) { ensure(r+1, c+1); return values[c + r*cols]; }
	inline const T &operator () (int r, int c) const { return values[c + r*cols]; }
	inline T &operator [] (int i) { return values[i]; }
	inline const T &operator [] (int i) const { return values[i]; }	
	inline int getRows() const {	return cols != 0 ? values.size()/cols : 0; }
	inline int getCols() const { return cols; }
	inline int size() const { return values.size(); }
	
	
private:
	variableArray<T> values;
	int cols;
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
template<typename storage_type, typename value_type, int rows, int cols> class storage_traits;

//! fixed Storage traits
template<typename value_type, int rows, int cols>
class storage_traits<fixedStorage, value_type, rows, cols>
{
public:
	typedef fixedArray<value_type, rows> array_type;
	typedef fixedArray2<value_type, rows, cols> array2_type;	
};

//! variable storage traits
template <typename value_type, int rows, int cols>
class storage_traits<variableStorage, value_type, rows, cols>
{	
public:
	typedef variableArray<value_type> array_type;
	typedef variableArray2<value_type> array2_type;
};

} // namespace ug
