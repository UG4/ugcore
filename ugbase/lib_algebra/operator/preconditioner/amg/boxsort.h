/*
 *  BoxSort.h
 *
 *  Created by Martin Rupp on 24.09.10.
 *  Copyright 2010 . All rights reserved.
 *
 */

#ifndef __H__LIB_DISCRETIZATION__AMG_SOLVER__BOXSORT_H__
#define __H__LIB_DISCRETIZATION__AMG_SOLVER__BOXSORT_H__


namespace ug{

// note: in the current implementation, this is not a "stable" sort
/**
 * \brief updateable priority queue class.
  * maxheap works on an external array of elements T, note that none of the elements of m_arr[0]..m_arr[size-1]
 * are in the queue at the begining. You can insert elements by using BoxSort::insert(i);
 *
 * \param T type of elements in the maxheap. Have to support size_t T::get_val() const, and get_val() has to
 * be between 0 and and small number (for example, 100).
 */
template<typename T>
class BoxSort
{
private:
	BoxSort(const BoxSort<T> &other);

public:
	BoxSort() : arr(NULL), box(), posInBox(NULL), values(NULL), m_size(0), m_height(0)
	{

	}
	/** constructor
	 * \param	n		maximal number of elements
	 * \param	arr_	array with elements which are to compare. note that non of these are in the queue in the beginning
	 */
	BoxSort(size_t n, T *arr_) : arr(NULL), box(), posInBox(NULL), values(NULL), m_size(0), m_height(0)
	{
		create(n, arr_);
	}
	//! deconstructor
	~BoxSort()
	{
		delete [] values;
		delete [] posInBox;
	}

	/** create
	 * \brief creates the queue on a external array arr_[0]..arr_[n-1]
	 */
	void create(size_t n, T *arr_)
	{
		if(n != m_size)
		{
			if(posInBox) delete [] posInBox;
			if(values) delete [] values;
			posInBox = new int[n];
			values = new size_t[n];
			m_size = n;
		}
		m_height = 0;
		box.resize(0);
		arr = arr_;
		for(size_t i=0; i<n; i++) posInBox[i] = -1;
	}

	void reset()
	{
		box.resize(0);
	}

	//! insertItem
	//! inserts Item arr[i] (see constructor) into heap
	//! @param i index of item in arr
	void insert_item(size_t i)
	{
		size_check(i);
		size_t val = arr[i].get_val();
		UG_ASSERT(val < 5000, "boxsort is not meant for big values");
		if(box.size() < val+1)
		{
			//size_t s = box.size();
			box.resize(val+1);
			//			for(;s<val+1; s++)
			//				box[s].reserve(10);
		}


		box[val].push_back(i);
		values[i] = val;
		posInBox[i] = box[val].size()-1;
		m_height++;
	}

	// remove
	//! removes index i in arr
	void remove(size_t i)
	{
		size_check(i);

		UG_ASSERT(posInBox[i] != -1, "item " << i << " cannot be removed from queue, since it is not in the queue.");

		size_t val = values[i];

		// swap last of this box and i
		size_t last = box[val].back();
		if(last!=i)
			assert(posInBox[i] < (int)box[val].size()-1);
		box[val].pop_back();
		if(last != i)
		{

			box[val][posInBox[i]] = last;
			posInBox[last] = posInBox[i];
		}

		posInBox[i] = -1;

		while(box.size() && box.back().size() == 0)
			box.pop_back();
		m_height--;
	}

	// removeMax
	//! returns the index in arr of maximal element of the heap
	size_t remove_max()
	{
		UG_ASSERT(m_height > 0 && box.size() > 0, "queue empty! (height = " << m_height << ", boxsize = " << box.size() << ").");

		size_t val = box.size()-1;
		size_t i = box[val].back();
		box[val].pop_back();
		posInBox[i] = -1;
		while(box.size() && box.back().size() == 0)
			box.pop_back();
		m_height--;
		return i;
	}

	// get_max
	//! returns the index in m_arr of maximal element of the heap
	size_t get_max() const
	{
		UG_ASSERT(m_height > 0 || box.size() <= 0, "queue empty! (height = " << m_height << ", boxsize = " << box.size() << ").");

		size_t val = box.size()-1;
		return box[val].back();
	}

	//!
	//! @param i index in arr for which to update
	void update(size_t i)
	{
		if(posInBox[i] != -1 && values[i] != arr[i].get_val())
		{
			remove(i);
			insert_item(i);
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

	//!
	//! debug print output
	/*void print()
	{

		for(int i=0; i<m_size; i++)
		{
			cout << "Element " << i;
			if(posInBox[i] < 0)	cout << " not in box." << endl;
			else cout << " val = " << values[i] << " get_val = " <<  arr[i].get_val() << " posInBox= " << posInBox[i] << endl;
		}

		cout << box.size() << " boxs. content:" << endl;
		for(int i=0; i<box.size(); i++)
		{
			cout << "box[" << i << "]: ";
			for(int j=0; j<box[i].size(); j++)
			{
				cout << " " << box[i][j];
				if(values[box[i][j]] != i) cout << " <-err ";
				if( posInBox[box[i][j]] != j) cout << " <-err2 ";
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
	T *arr;				//< pointer to array with elements of type T

	stdvector<stdvector<size_t> > box;

	int *posInBox;	//< item is at box[values[i]]
	size_t *values;		//< values[i] is the value used for sorting of element i
	size_t m_size;		//< maximal size of the PQ = size of array arr
	size_t m_height;
};

} // namespace ug

#endif
