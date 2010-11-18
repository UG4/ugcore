/*
 *  fixed_array.h
 *
 *  Created by Martin Rupp on 21.07.10.
 *  Copyright 2010 . All rights reserved.
 *
 */

#ifndef __H__UG__COMMON__FIXED_ARRAY_SPECIALIZATION_H__
#define __H__UG__COMMON__FIXED_ARRAY_SPECIALIZATION_H__

namespace ug{

template<typename T>
class FixedArray1<T, 1>
{
public:
	typedef T value_type;
	typedef size_t size_type;
	enum {Size = 1 };

public:
	FixedArray1() { values[0] = 0; }
	FixedArray1(const FixedArray1<T, 1> &other) { values[0] = other.values[0]; }

	// capacity
	inline size_type
	size() const { return Size; }

	inline bool
	resize(size_type newN) { assert(newN == Size); return true; }

	// Element access
	inline const T &
	operator[] (size_type i) const
	{
		return at(i);
	}

	inline T &
	operator[] (size_type i)
	{
		return at(i);
	}

	inline const T &
	at(size_type i) const { assert(i<Size); return values[i]; }

	inline T &
	at(size_type i) { assert(i<Size); return values[i]; }


	union
	{
		struct
		{
			value_type x;
		};

		value_type values[1];
	};
};



template<typename T>
class FixedArray1<T, 2>
{
public:
	typedef T value_type;
	typedef size_t size_type;
	enum {Size = 2 };

public:
	FixedArray1() { values[0] = 0.0; values[1] = 0.0; }
	FixedArray1(const FixedArray1<T, 2> &other) { values[0] = other.values[0]; values[1] = other.values[1]; }

	// capacity
	inline size_type
	size() const { return Size; }

	inline bool
	resize(size_type newN) { assert(newN == Size); return true; }

	// Element access
	inline const T &
	operator[] (size_type i) const
	{
		return at(i);
	}

	inline T &
	operator[] (size_type i)
	{
		return at(i);
	}

	inline const T &
	at(size_type i) const { assert(i<Size); return values[i]; }

	inline T &
	at(size_type i) { assert(i<Size); return values[i]; }

	union
	{
		struct
		{
			value_type x;
			value_type y;
		};

		value_type values[2];
	};
};


template<typename T>
class FixedArray1<T, 3>
{
public:
	typedef T value_type;
	typedef size_t size_type;
	enum {Size = 3 };

public:
	FixedArray1() { for(int i=0; i<Size; i++) values[i]=0.0; }
	FixedArray1(const FixedArray1<T, 3> &other) { values[0] = other.values[0]; values[1] = other.values[1]; values[2] = other.values[2]; }

	// capacity
	inline size_type
	size() const { return Size; }

	inline bool
	resize(size_type newN) { assert(newN == Size); return true; }

	// Element access
	inline const T &
	operator[] (size_type i) const
	{
		return at(i);
	}

	inline T &
	operator[] (size_type i)
	{
		return at(i);
	}

	inline const T &
	at(size_type i) const { assert(i<Size); return values[i]; }

	inline T &
	at(size_type i) { assert(i<Size); return values[i]; }


	union
	{
		struct
		{
			value_type x;
			value_type y;
			value_type z;
		};

		value_type values[3];
	};
};

template<typename T>
class FixedArray1<T, 4>
{
public:
	typedef T value_type;
	typedef size_t size_type;
	enum {Size = 4 };

public:
	FixedArray1() { for(int i=0; i<Size; i++) values[i]=0.0; }
	FixedArray1(const FixedArray1<T, 4> &other) { for(int i=0; i<Size; i++) values[i]=other.values[i]; }

	// capacity
	inline size_type
	size() const { return Size; }

	inline bool
	resize(size_type newN) { assert(newN == Size); return true; }

	// Element access
	inline const T &
	operator[] (size_type i) const
	{
		return at(i);
	}

	inline T &
	operator[] (size_type i)
	{
		return at(i);
	}

	inline const T &
	at(size_type i) const { assert(i<Size); return values[i]; }

	inline T &
	at(size_type i) { assert(i<Size); return values[i]; }

	union
	{
		struct
		{
			value_type x;
			value_type y;
			value_type z;
			value_type w;
		};

		value_type values[4];
	};
};
}
#endif // __H__UG__COMMON__FIXED_ARRAY_SPECIALIZATION_H__
