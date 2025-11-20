/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
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


#ifndef __H__UG__COMMON__FIXED_ARRAY_SPECIALIZATION_H__
#define __H__UG__COMMON__FIXED_ARRAY_SPECIALIZATION_H__

namespace ug{

template<typename T>
class FixedArray1<T, 1>
{
public:
	using value_type = T ;
	using size_type = size_t ;
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
	operator [] (size_type i) const
	{
		return at(i);
	}

	inline T &
	operator [] (size_type i)
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
	using value_type = T;
	using size_type =  size_t ;
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
	operator [] (size_type i) const
	{
		return at(i);
	}

	inline T &
	operator [] (size_type i)
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
	using value_type = T ;
	using size_type =size_t ;
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
	operator [] (size_type i) const
	{
		return at(i);
	}

	inline T &
	operator [] (size_type i)
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
	using value_type = T ;
	using size_type = size_t;
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
	operator [] (size_type i) const
	{
		return at(i);
	}

	inline T &
	operator [] (size_type i)
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
#endif