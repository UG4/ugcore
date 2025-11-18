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

#ifndef __H__UG__LIB_ALGEBRA__AMG_SOLVER__MAXHEAP_H__
#define __H__UG__LIB_ALGEBRA__AMG_SOLVER__MAXHEAP_H__

namespace ug{
// maxheap
//--------------
/**
 * \brief updateable priority queue class.
 * unlike most PQ implementations, we need a method to inform the PQ of updated elements
 * thats why we cannot use priority_queue from the STL.
 * maxheap works on an external array of elements T, note that none of the elements of m_arr[0]..m_arr[size-1]
 * are in the heap at the begining. You can insert elements by using maxheap::insert(i);
 * \note a radix heap would be faster for some applications.
 * \param T type of elements in the maxheap. Have to support bool operator <(const T &a, const T &b)
 */
template<typename T>
class maxheap
{
private:
	maxheap(const maxheap<T> *other);

public:
	maxheap()
	{
		m_arr = nullptr;
		m_height = 0;
		m_size = 0;
	}


	/** constructor
	 * \param	n		maximal number of elements
	 * \param	arr_	array with elements which are to compare. note that non of these are in the heap in the beginning
	 */
	maxheap(size_t n, const T *arr_)
	{
		m_size = 0;
		m_heap = nullptr;
		m_posinheap = nullptr;
		create(n, arr_);
	}

	maxheap(const std::vector<T> &v)
	{
		m_size = 0;
		m_heap = nullptr;
		m_posinheap = nullptr;
		create(v.size(), &v[0]);
	}
	//! deconstructor
	~maxheap()
	{
	}


	/** create
	 * \brief creates the heap on a external array arr_[0]..arr_[n-1]
	 */
	void create(size_t n, const T *arr_)
	{
		if(m_size != n)
		{
			m_heap.resize(n);
			m_posinheap.resize(n);
		}
		
		m_arr = arr_;
		m_height = 0;
		for(size_t i=0; i<n; i++) m_posinheap[i] = -1;
		for(size_t i=0; i<n; i++) m_heap[i] = -1;
		m_size = n;
		
	}
	void create(const std::vector<T> &v)
	{
		create(v.size(), &v[0]);
	}
	//! reset
	//! set m_height 0 (= remove all items from the heap)
	void reset()
	{
		m_height = 0;
	}

	//! insertItem
	//! inserts Item m_arr[i] (see constructor) into the heap
	//! \param i index of item in m_arr
	void insert_item(int i)
	{
		if(!is_in(i))
		{
			UG_ASSERT(m_height < m_size, "more elements added than there are in the external array. double adds?");
			m_posinheap[i] = m_height;
			m_heap[m_height] = i;
			m_height++;
		}
		upheap(i);
	}

	// remove
	//! removes index i in m_arr
	void remove(int i)
	{
		if(!is_in(i)) return;
		int j = m_heap[m_height-1];
		myswap(i, j);
		m_heap[m_height-1] = -1;
		m_height--;
		m_posinheap[i] = -1;

		downheap(j);
	}
	
	// remove_max
	//! returns the index in m_arr of maximal element of the m_heap and removes it from the m_heap
	int remove_max()
	{
		UG_ASSERT(m_height > 0, "m_heap already empty");
		int m = m_heap[0];
		remove(m);
		return m;
	}

	// get_max
	//! returns the index in m_arr of maximal element of the heap
	int get_max() const
	{
		return m_heap[0];
	}

	/**
	 * \brief update the heap position of array element i (m_arr[i])
	 * use this method if you have changed the value used in operator < of
	 * element i.
	 * \param i index in m_arr for which to update
	 */
	void update(int i)
	{
		if(!is_in(i)) return;
		if(m_arr[i] > m_arr[parent(i)])
			upheap(i);
		else
			downheap(i);
	}
	
	bool is_in(size_t i)
	{
		return m_posinheap[i] != -1;
	}


	//!
	//! debug print output
	void print() const
	{
		std::cout << "maxheap, size = " << m_size << ", height = " << m_height << std::endl;
		for(size_t i=0; i<m_height; i++)
		{
			std::cout << i << ": pos: " << m_heap[i] << " parent: " << parent(m_heap[i]) << " element: " << m_arr[m_heap[i]] <<
					(m_arr[m_heap[i]] > m_arr[parent(m_heap[i])] ? " ERR " : "") << std::endl;

		}
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

private:
	//!
	//! \param i index in m_arr for which to upheap
	//! \sa update
	void upheap(int i)
	{
		if(m_posinheap[i] == -1) return;
		while(m_arr[i] > m_arr[parent(i)])
			myswap(i, parent(i));
	}
	
	//!
	//! @param i index in m_arr for which to downheap
	void downheap(int i)
	{
		if(m_posinheap[i] == -1) return;
		while(true)
		{		
			const T &l = m_arr[leftchild(i)];
			const T &r = m_arr[rightchild(i)];
			const T &t = m_arr[i];
			if(l > t || r > t)
			{
				if(l > r)
					myswap(leftchild(i), i);
				else
					myswap(rightchild(i), i);
			}
			else
				break;
		}
	}
	
	// parent
	//! returns the index in m_arr of the parent of i
	//! all indices in m_arr, NOT in m_heap
	//! \param index index in m_arr of element
	//! \return returns index in m_arr of the parent
	int parent(int index) const
	{
		int p = m_posinheap[index];
		int parentpos = (p == 0 ? 0 : (p+1)/2 -1);
		return m_heap[parentpos];
	}
	
	
	// leftchild
	//! returns the index in m_arr of the left child of i
	//! all indices in m_arr, NOT in m_heap
	//! \param index index in m_arr of element
	//! \return returns index in m_arr of the left child
	int leftchild(int index) const
	{
		int p = m_posinheap[index];
		p = (p+1)*2 -1;
		if(p < (int)m_height)
			return m_heap[p];
		else return index;
	}
	
	// rightchild
	//! returns the index in m_arr of the right child of i
	//! all indices in m_arr, NOT in m_heap
	//! \param index index in m_arr of element
	//! \return returns index in m_arr of the right child
	int rightchild(int index) const
	{
		int p = m_posinheap[index];
		p = (p+1)*2 -1 +1;
		if(p < (int)m_height)
			return m_heap[p];
		else return index;
	}
	
	// myswap
	//! 
	//! swaps elements i and j of array m_posinheap and their counterparts in m_heap
	void myswap(int i, int j)
	{
		int posinheapi = m_posinheap[i];
		int posinheapj = m_posinheap[j];
		
		m_heap[posinheapi] = j;
		m_heap[posinheapj] = i;
		
		m_posinheap[i] = posinheapj;
		m_posinheap[j] = posinheapi;
	}
	
private:
	const T *m_arr;			//< pointer to array with elements of type T
	stdvector<int> m_heap;		//< m_heap of the elements m_heap[0] is the index of the largest element of m_arr[i] forall i=0..m_size-1 and m_posinheap[i] != -1
	stdvector<int> m_posinheap;	//< m_posinheap[i] is the position of element m_arr[i] in the m_heap. -1 if removed, otherwise m_posinheap[m_heap[i]] = i
	size_t m_height;		//< m_height of the m_heap, elements m_heap[0]..m_heap[m_height-1] are valid.
	size_t m_size;			//< maximal size of the m_heap = size of array m_arr
};
	
} // namespace ug

#endif // __H__UG__LIB_ALGEBRA__AMG_SOLVER__MAXHEAP_H__
