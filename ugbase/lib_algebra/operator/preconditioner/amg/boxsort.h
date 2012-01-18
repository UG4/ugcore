/*
 *  BoxPriorityQueue.h
 *
 *  Created by Martin Rupp on 24.09.10.
 *  Copyright 2010 . All rights reserved.
 *
 */

#ifndef __H__UG__LIB_DISC__AMG_SOLVER__BoxPriorityQueue_H__
#define __H__UG__LIB_DISC__AMG_SOLVER__BoxPriorityQueue_H__

#include <vector>

namespace ug{

#ifdef UGASSERT
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
}

//! maximal value for T::get_val(). Keep in mind that memory requirements are O(max get_val()).
const int BOXPRIORITYQUEUE_MAXIMAL_VALUE = 500;

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
class BoxPriorityQueue
{
private:
	BoxPriorityQueue(const BoxPriorityQueue<T> &other);

public:
	BoxPriorityQueue() : arr(NULL), m_box(), m_posInBox(), m_values(), m_size(0), m_height(0)
	{

	}
	/** constructor
	 * \param	n		maximal number of elements
	 * \param	arr_	array with elements which are to compare. note that non of these are in the queue in the beginning
	 */
	BoxPriorityQueue(size_t n, const T *arr_) : m_size(0), m_height(0)
	{
		create(n, arr_);
	}

	BoxPriorityQueue(const std::vector<T> &v) : m_size(0), m_height(0)
	{
		create(v.size(), &v[0]);
	}
	//! deconstructor
	~BoxPriorityQueue()
	{
	}

	/** create
	 * \brief creates the queue on a external array arr_[0]..arr_[n-1]
	 */
	void create(size_t n, const T *arr_)
	{
		if(n != m_size)
		{
			m_values.resize(n);
			m_posInBox.resize(n);
			m_size = n;
		}
		m_height = 0;
		m_box.resize(0);
		arr = arr_;
		for(size_t i=0; i<n; i++) m_posInBox[i] = -1;
	}

	void create(const std::vector<T> &v)
	{
		create(v.size(), &v[0]);
	}

	void reset()
	{
		m_box.resize(0);
	}

	//! insertItem
	//! inserts Item arr[i] (see constructor) into heap
	//! @param i index of item in arr
	void insert_item(size_t i)
	{
		size_check(i);
		size_t val = get_val(arr[i]);
		UG_ASSERT(val < (size_t)BOXPRIORITYQUEUE_MAXIMAL_VALUE, "T::get_val() has to be < " << BOXPRIORITYQUEUE_MAXIMAL_VALUE << " but is " << val);
		if(m_box.size() < val+1)
		{
			//size_t s = m_box.size();
			m_box.resize(val+1);
			//			for(;s<val+1; s++)
			//				m_box[s].reserve(10);
		}


		m_box[val].push_back(i);
		m_values[i] = val;
		m_posInBox[i] = m_box[val].size()-1;
		m_height++;
	}

	// remove
	//! removes index i in arr
	void remove(size_t i)
	{
		size_check(i);

		if(m_posInBox[i] == -1) return;
		//UG_ASSERT(m_posInBox[i] != -1, "item " << i << " cannot be removed from queue, since it is not in the queue.");

		size_t val = m_values[i];

		// swap last of this box and i
		/*size_t last = m_box[val].back();
		if(last!=i)
			assert(m_posInBox[i] < (int)m_box[val].size()-1);
		m_box[val].pop_back();
		if(last != i)
		{

			m_box[val][m_posInBox[i]] = last;
			m_posInBox[last] = m_posInBox[i];
			}*/
		m_box[val].erase(m_box[val].begin()+m_posInBox[i]);
		for(size_t j=0; j<m_box[val].size(); j++)
		  {
		    m_posInBox[m_box[val][j]] = j;
		  }

		m_posInBox[i] = -1;

		while(m_box.size() && m_box.back().size() == 0)
			m_box.pop_back();
		m_height--;
	}

	// removeMax
	//! returns the index in arr of maximal element of the heap
	size_t remove_max()
	{
		UG_ASSERT(m_height > 0 && m_box.size() > 0, "queue empty! (height = " << m_height << ", boxsize = " << m_box.size() << ").");

		size_t val = m_box.size()-1;
		size_t i = m_box[val][0];
		remove(i);
		return i;
	}

	// get_max
	//! returns the index in m_arr of maximal element of the heap
	size_t get_max() const
	{
		UG_ASSERT(m_height > 0 || m_box.size() <= 0, "queue empty! (height = " << m_height << ", boxsize = " << m_box.size() << ").");

		size_t val = m_box.size()-1;
		return m_box[val][0];
	}

	//!
	//! @param i index in arr for which to update
	void update(size_t i)
	{
		if(m_posInBox[i] != -1 && m_values[i] != get_val(arr[i]))
		{
			remove(i);
			insert_item(i);
		}
	}

	bool is_in(size_t i)
	{
		return m_posInBox[i] != -1;
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

	stdvector<int> m_posInBox;	//< item is at box[m_values[i]]
	stdvector<size_t> m_values; 	//< m_values[i] is the value used for sorting of element i
	size_t m_size;				//< maximal size of the PQ = size of array arr
	size_t m_height;
};

} // namespace ug

#endif
