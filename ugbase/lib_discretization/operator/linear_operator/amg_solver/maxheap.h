/*
 *  maxheap.h
 *  flexamg
 *
 *  Created by Martin Rupp on 02.10.09.
 *  Copyright 2009 G-CSC. All rights reserved.
 *
 */
#pragma once

namespace ug{
// maxheap
//--------------
//! maxheap priority queue class. 
//! unlike most PQ implementations, we need a method to inform the PQ of updated elements
//! thats why we cannot use priority_queue from the STL.
//! surely a radix heap would be faster. 
//! problem is: we somehow want a STABLE sort, for that coarsening is in the direction of the numbering of the elements
//! but this is somehow artificial.
template<typename T>
class maxheap
{
public:
	//! constructor
	//! @param	n		maximal number of elements
	//! @param	arr_	array with elements which are to compare. note that non of these are in the heap in the beginning
	maxheap(int n, T *arr_)
	{
		arr = arr_;
		height = 0;
		heap = new int[n];
		posinheap = new int[n];
		for(int i=0; i<n; i++) posinheap[i] = -1;
		size = n;
	}
	//! deconstructor
	~maxheap()
	{
		delete [] heap;
		delete [] posinheap;
	}

	//! reset
	//! set height 0
	void create(int n, T *arr_)
	{
		delete[] heap;
		delete [] posinheap;
		
		arr = arr_;
		height = 0;
		heap = new int[n];
		posinheap = new int[n];
		for(int i=0; i<n; i++) posinheap[i] = -1;
		size = n;
		
	}
	
	void reset()
	{
		height = 0;
	}

	//! insertItem
	//! inserts Item arr[i] (see constructor) into heap
	//! @param i index of item in arr
	void insertItem(int i)
	{
		UG_ASSERT(height < size, "heap too small");
		posinheap[i] = height;
		heap[height] = i;
		height++;
		upheap(i);
	}

	// remove
	//! removes index i in arr
	void remove(int i)
	{
		if(posinheap[i] == -1) return;
		int j = heap[height-1];
		myswap(i, j);
		height--;
		downheap(j);
		posinheap[i] = -1;
	}
	
	// removeMax
	//! returns the index in arr of maximal element of the heap
	int removeMax()
	{
		UG_ASSERT(height > 0, "heap already empty");
		myswap(heap[0], heap[height-1]);
		height--;
		downheap(heap[0]);
		return heap[height];
	}

	//!
	//! @param i index in arr for which to update
	void update(int i)
	{
		if(posinheap[i] == -1) return;
		if(arr[i] > arr[parent(i)])
			upheap(i);
		else
			downheap(i);
	}
	
	//!
	//! @param i index in arr for which to upheap
	void upheap(int i)
	{
		if(posinheap[i] == -1) return;
		while(arr[i] > arr[parent(i)])
			myswap(i, parent(i));
	}
	
	//!
	//! @param i index in arr for which to downheap
	void downheap(int i)
	{
		if(posinheap[i] == -1) return;
		while(1)
		{		
			T &l = arr[leftchild(i)];
			T &r = arr[rightchild(i)];
			T &t = arr[i];
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
	
	//!
	//! debug print output
	void print()
	{
		for(int i=0; i<height; i++)
		{
			cout << i << ": pos: " << heap[i] << " parent: " << parent(heap[i]) << " rating: " << arr[heap[i]].rating << (arr[heap[i]] > arr[parent(heap[i])] ? " ERR " : "");// << endl;
			cout << (arr[heap[i]].rating > arr[parent(heap[i])].rating ? " r> " : "") << (&arr[heap[i]] > &arr[parent(heap[i])] ? " a> " : "") << endl;
		}
		
	}

private:		
	// parent
	//! returns the index in arr of the parent of i
	//! all indices in arr, NOT of the heap
	//! @param index index in arr of element
	//! @return returns index in arr of the parent
	int parent(int index)
	{
		int p = posinheap[index];
		int parentpos = (p == 0 ? 0 : (p+1)/2 -1);
		return heap[parentpos];
	}
	
	
	// leftchild
	//! returns the index in arr of the left child of i
	//! all indices in arr, NOT of the heap
	//! @param index index in arr of element
	//! @return returns index in arr of the left child
	int leftchild(int index)
	{
		int p = posinheap[index];
		p = (p+1)*2 -1;
		if(p < height)
			return heap[p];
		else return index;
	}
	
	// rightchild
	//! returns the index in arr of the right child of i
	//! all indices in arr, NOT of the heap
	//! @param index index in arr of element
	//! @return returns index in arr of the right child
	int rightchild(int index)
	{
		int p = posinheap[index];
		p = (p+1)*2 -1 +1;
		if(p < height)
			return heap[p];
		else return index;
	}
	
	// myswap
	//! 
	//! swaps elements i and j of array posinheap and their counterparts in heap
	void myswap(int i, int j)
	{
		int posinheapi = posinheap[i];
		int posinheapj = posinheap[j];
		
		heap[posinheapi] = j;
		heap[posinheapj] = i;
		
		posinheap[i] = posinheapj;
		posinheap[j] = posinheapi;
	}
	
			
	T *arr;				//< pointer to array with elements of type T
	int *heap;			//< heap of the elements heap[0] is the index of the largest element of arr[0]..arr[size]
	int *posinheap;		//< posinheap[i] is the position of element arr[i] in the heap. -1 if removed
	int height;			//< height of the heap
	int size;			//< maximal size of the heap = size of array arr
};
	
} // namespace ug
