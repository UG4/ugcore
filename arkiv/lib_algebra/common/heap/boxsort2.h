/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
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


#ifndef __H__UG__LIB_DISC__AMG_SOLVER__BoxPriorityQueue2_H__
#define __H__UG__LIB_DISC__AMG_SOLVER__BoxPriorityQueue2_H__

#include <vector>

namespace ug{

/*#ifdef UGASSERT
#define MYASSERT(a, b) UGASSERT(a, b)
#else
#define MYASSERT(a, b) assert(a && ##b)
#endif

template<typename T>
inline size_t get_val(const T &t)
{
	return t.get_val();
}

inline size_t get_val(size_t i)
{
	return i;
}*/

//! maximal value for T::get_val(). Keep in mind that memory requirements are O(max get_val()).
const int BOXPRIORITYQUEUE2_MAXIMAL_VALUE = 500;

/**
 * \brief updateable priority queue class.
  * this priority queue works on an external array of elements T.
  * \note none of the elements of m_arr[0]..m_arr[size-1] are in the queue at the begining.
  * You can insert elements by using BoxPriorityQueue::insert_item(i);
 *
 * \param T type of elements in the maxheap. Has to support size_t T::get_val() const, and get_val() has to
 * be between 0 and and small number (for example, 500).
 *
 * Elements are put in "boxes" depending on their get_val()-value - this makes sorting extremely easy:
 * sorting and accessing is O(1), and inserting is maximal a) creation of enough boxes, and b) inserting
 * into a box, so it is O(maximal get_val()). This algorithm is like a "radical" radix exchange sort.
 * \note memory requirements are O(maximal get_val()) + O(size).
 * \note: in the current implementation, this is not a "stable" sort
 *
 */
template<typename T>
class BoxPriorityQueue2
{
private:
	BoxPriorityQueue2(const BoxPriorityQueue2<T> &other);

public:
	BoxPriorityQueue2()
	{

	}
	/** constructor
	 * \param	n		maximal number of elements
	 * \param	arr_	array with elements which are to compare. note that non of these are in the queue in the beginning
	 */
	BoxPriorityQueue2(size_t n, const T *arr_) : m_size(0), m_height(0)
	{
		create(n, arr_);
	}

	BoxPriorityQueue2(const std::vector<T> &v) : m_size(0), m_height(0)
	{
		create(v.size(), &v[0]);
	}
	//! deconstructor
	~BoxPriorityQueue2() = default;

	/** create
	 * \brief creates the queue on a external array arr_[0]..arr_[n-1]
	 */
	void create(size_t n, const T *arr_)
	{
		if(n != m_size)
		{
			m_values.resize(n);
			m_prev.resize(n);
			m_next.resize(n);
			m_size = n;
		}
		m_size = n;

		arr = arr_;

		reset();
	}

	void create(const std::vector<T> &v)
	{
		create(v.size(), &v[0]);
	}

	void reset()
	{
		m_height = 0;
		m_box.resize(0);
		for(size_t i=0; i<m_prev.size(); i++)
			m_prev[i] = -2;
	}

	//! insertItem
	//! inserts Item arr[i] (see constructor) into heap
	//! @param i index of item in arr
	void insert_item(size_t i)
	{
		UG_ASSERT(is_in(i), "item is already in");
		size_check(i);
		size_t val = get_val(arr[i]);
		UG_ASSERT(val < (size_t)BOXPRIORITYQUEUE2_MAXIMAL_VALUE, "T::get_val() has to be < " << BOXPRIORITYQUEUE_MAXIMAL_VALUE << " but is " << val);
		if(m_box.size() < val+1)
		{
			m_boxBegin.resize(val+1, -1);
			m_boxEnd.resize(val+1, -1);
		}

		if(m_boxBegin[val] == -1)
		{
			m_prev[i] = -1;
			m_next[i] = -1;
			m_boxBegin[val] = m_boxEnd[val] = i;
		}
		else
		{
			m_prev[i] = m_boxEnd[val];
			m_next[i] = -1;
			m_boxEnd[val] = i;
		}
		m_values[i] = val;
		m_height++;
	}

	// remove
	//! removes index i in arr
	void remove(size_t i)
	{
		size_check(i);
		UG_ASSERT(is_in(i), "item is not in");

		size_t val = m_values[i];
		if(m_prev[i] == -1)	m_boxBegin[val] = m_next[i];
		else				m_next[m_prev[i]] = m_next[i];
		if(m_next[i] == -1) m_boxEnd[val] = m_prev[i];
		else				m_prev[m_next[i]] = m_prev[i];

		m_prev[i] = -2;
		m_height--;
		while(m_boxBegin.size() && m_boxBegin.back() == -1)
		{
			m_boxBegin.pop_back();
			m_boxEnd.pop_back();
		}
	}

	// removeMax
	//! returns the index in arr of maximal element of the heap
	size_t remove_max()
	{
		UG_ASSERT(m_height > 0 && m_box.size() > 0, "queue empty! (height = " << m_height << ", boxsize = " << m_box.size() << ").");

		int i = m_boxBegin.back();
		remove(i);
		return i;
	}

	// get_max
	//! returns the index in m_arr of maximal element of the heap
	size_t get_max() const
	{
		UG_ASSERT(m_height > 0 || m_box.size() <= 0, "queue empty! (height = " << m_height << ", boxsize = " << m_box.size() << ").");
		return m_boxBegin.back();
	}

	//!
	//! @param i index in arr for which to update
	void update(size_t i)
	{
		if(is_in(i) && m_values[i] != get_val(arr[i]))
		{
			remove(i);
			insert_item(i);
		}
	}

	bool is_in(size_t i)
	{
		return m_prev[i] != -2;
	}

	//! returns size of external array
	size_t arr_size() const
	{
		return m_size;
	}

	//! returns nr of elements in the heap.
	size_t height() const
	{
		return m_height;
	}

	//!
	//! debug print output
	/*void print()
	{

		for(int i=0; i<m_size; i++)
		{
			cout << "Element " << i;
			if(m_posInBox[i] < 0)	cout << " not in box." << endl;
			else cout << " val = " << m_values[i] << " get_val = " <<  arr[i].get_val() << " m_posInBox= " << m_posInBox[i] << endl;
		}

		cout << m_box.size() << " boxs. content:" << endl;
		for(int i=0; i<m_box.size(); i++)
		{
			cout << "m_box[" << i << "]: ";
			for(int j=0; j<m_box[i].size(); j++)
			{
				cout << " " << m_box[i][j];
				if(m_values[m_box[i][j]] != i) cout << " <-err ";
				if( m_posInBox[m_box[i][j]] != j) cout << " <-err2 ";
			}
			cout << endl;
		}
		cout.flush();
	}*/

private:
	inline void size_check(size_t i)
	{
		UG_ASSERT(i < arr_size(), "accessing element " << i << ", but size of the array is only " << arr_size());
	}
	const T *arr;				//< pointer to array with elements of type T

	stdvector<stdvector<size_t> > m_box;

	stdvector<int> m_prev, m_next, m_boxBegin, m_boxEnd;
	stdvector<size_t> m_values; 							//< m_values[i] is the value used for sorting of element i
	size_t m_size;											//< maximal size of the PQ = size of array arr
	size_t m_height;
};

} // namespace ug

#endif
