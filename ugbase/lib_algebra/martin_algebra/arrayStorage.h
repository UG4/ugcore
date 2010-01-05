/*
 *  arrayStorage.h
 *  flexamg
 *
 *  Created by Martin Rupp on 18.12.09.
 *  Copyright 2009 . All rights reserved.
 *
 */
#pragma once


template<typename T, int n>
class fixedArray
{	
public:
	fixedArray(){}
	fixedArray(int n_) { ASSERT(n_ == n); }
	inline size_t size() const 
	{ 
		return n;
	}
	inline T &operator [] (int i) { return values[i]; }
	inline T operator [] (int i) const { return values[i]; }
	bool operator == (const fixedArray<T, n> &other) const
	{
		return memcmp(values, other.values, sizeof(T)*n)==0;
	}
	bool operator != (const fixedArray<T, n> &other) const
	{
		return ! operator == (other);
	}
	inline void setSize(int s, bool bZero=true)
	{
		ASSERT2(s <= n, "fixed Array is too small (max " << n << ", tried " << s << ".");		
		if(bZero)
			memset(values, 0, sizeof(T)*n);
	}
private:
	T values[n];	
};

template<typename T>
class variableArray
{
public:
	variableArray()
	{
		values = 0;
		n = 0;
	}
	variableArray(int n_)
	{
		values = 0;
		n = 0;
		setSize(n_);
	}
	
	variableArray(const variableArray<T> &other)
	{
		n = other.n;
		values = new T[n];
		memcpy(values, other.values, sizeof(T)*n);
	}
	
	bool operator == (const variableArray<T> &other) const
	{
		return memcmp(values, other.values, sizeof(T)*n)==0;
	}
	bool operator != (const variableArray<T> &other) const
	{
		return ! operator == (other);
	}
	
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
	
	~variableArray()
	{
		if(values) delete[] values;
	}
	inline size_t size() const 
	{ 
		return n;
	}
	inline T &operator [] (int i) { return values[i]; }
	inline T operator [] (int i) const { return values[i]; }
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
	
	void swap(variableArray &a)
	{
		std::swap(a.values, values);
		std::swap(a.n, n);
	}
private:
	T *values;
	int n;
};




template<typename T, int rows, int cols>
class fixedArray2
{
public:
	fixedArray2() { memset(values, 0, sizeof(double)*rows*cols); }
	fixedArray2(int rows_, int cols_) 
	{
		ASSERT2(rows_ == rows && cols == cols, "fixed Array! (rows " << rows << ", cols " << cols << " tried " << rows_ << ", " << cols_ << ".)");		
	}
	inline T &operator () (int r, int c) { ensure(r+1, c+1); return values[c + r*cols]; }
	inline T operator () (int r, int c) const { ensure(r+1, c+1); return values[c + r*cols]; }
	inline T &getAt (int r, int c) { ensure(r+1, c+1); return values[c + r*cols]; }
	inline T getAt (int r, int c) const { ensure(r+1, c+1); return values[c + r*cols]; }
	inline T &operator [] (int i) { ASSERT(i<rows*cols && i >= 0); return values[i]; }
	inline T operator [] (int i)  const { ASSERT(i<rows*cols && i >= 0); return values[i]; }	
	int size() const { return rows*cols; }

	inline void ensure(int rows_, int cols_) const
	{
		ASSERT2(rows_ <= rows && cols <= cols, "fixed Array is too small (max rows " << rows << ", cols " << cols << " tried " << rows_ << ", " << cols_ << ".)");		
	}
	inline void setSize(int rows_, int cols_, bool bZero=true)
	{
		ASSERT2(rows_ == rows && cols == cols, "fixed Array! (rows " << rows << ", cols " << cols << " tried " << rows_ << ", " << cols_ << ".)");		
		memset(values, 0, sizeof(double)*rows*cols);
	}
	
	int getRows() const { return rows; }
	int getCols() const { return cols; }
private:
	T values[rows*cols];	
};

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
	
	inline T &operator () (int r, int c) { ensure(r+1, c+1); return values[c + r*cols]; }
	inline T operator () (int r, int c) const { return values[c + r*cols]; }
	inline T &operator [] (int i) { return values[i]; }
	inline T operator [] (int i) const { return values[i]; }	
	int getRows() const {	return cols != 0 ? values.size()/cols : 0; }
	int getCols() const { return cols; }
	int size() const { return values.size(); }
	
	
private:
	variableArray<T> values;
	int cols;
};

class fixedStorage
{
public:
	static const char *getType() { return "fix"; }
};

class variableStorage 
{
public:
	static const char *getType() { return "var"; }
};

template<typename storage_type, typename value_type, int rows, int cols> class storage_traits;

template<typename value_type, int rows, int cols>
class storage_traits<fixedStorage, value_type, rows, cols>
{
public:
	typedef fixedArray<value_type, rows> array_type;
	typedef fixedArray2<value_type, rows, cols> array2_type;	
};

template <typename value_type, int rows, int cols>
class storage_traits<variableStorage, value_type, rows, cols>
{	
public:
	typedef variableArray<value_type> array_type;
	typedef variableArray2<value_type> array2_type;
};