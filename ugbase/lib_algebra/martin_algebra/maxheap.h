/*
 *  maxheap.h
 *  flexamg
 *
 *  Created by Martin Rupp on 02.10.09.
 *  Copyright 2009 G-CSC. All rights reserved.
 *
 */
#pragma once

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
	void reset()
	{
		for(int i=0; i<size; i++) posinheap[i] = -1;
		height = 0;
	}

	//! insertItem
	//! inserts Item arr[i] (see constructor) into heap
	//! @param i index of item in arr
	void insertItem(int i)
	{
		ASSERT2(height < size, "heap too small");
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
		ASSERT2(height > 0, "heap already empty");
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

/*
 
template<typename T>
class maxheap
	{
	public:
		maxheap(int n, T *arr_)
		{
			arr = arr_;
			height = 0;
			heap = new int[n];
			posinheap = new int[n];
			for(int i=0; i<n; i++) posinheap[i] = -1;
			size = n;
		}
		~maxheap()
		{
			delete [] heap;
			delete [] posinheap;
		}
		
		void reset()
		{
			for(int i=0; i<size; i++) posinheap[i] = -1;
			height = 0;
		}
		
		void insertItem(int i)
		{
			ASSERT2(height < size, "heap too small");
			posinheap[i] = height;
			heap[height] = i;
			height++;
			upheap(i);
		}
		
		int parent(int index)
		{
			int p = index;
			int parentpos = (p == 0 ? 0 : (p+1)/2 -1);
			return parentpos;
		}
		
		int leftchild(int index)
		{
			int p = index;
			p = (p+1)*2 -1;
			return p;
		}
		
		int rightchild(int index)
		{
			int p = index;
			p = (p+1)*2 -1 +1;
			return p;
		}
		
		void remove(int i)
		{
			if(posinheap[i] == -1) return;
			posinheap[i] = -1;
		}
		
		void assureNotDeleted(int i)
		{
			if(posinheap[heap[i]] == -1)
			{
				int L = leftchild(i); // leftchild
				if(L < height)
				{				
					assureNotDeleted(L);
					T &l = arr[heap[L]];
					
					int R = L+1;
					if(R < height)
					{
						assureNotDeleted(R);
						T &r = arr[heap[R]];
						
						if(l > r)
							myswap(L, i);
						else
							myswap(R, i);
					}
					else
						myswap(L, i);
				}
			}
		}
		
		int removeMax()
		{
			//ASSERT2(height > 0, "heap already empty");
			assureNotDeleted(0);
			posinheap[heap[0]] = -1;
			return heap[0];
		}
		
		void myswap(int i, int j)
		{		
			if(posinheap[heap[i]] != -1)
				posinheap[heap[i]] = j;
			if(posinheap[heap[j]] != -1)
				posinheap[heap[j]] = i;
			swap(heap[i], heap[j]);
		}
		
		void upheap(int i)
		{		
			if(posinheap[i] == -1) return;
			int posi = posinheap[i];
			int posp = parent(posi);
			while(arr[heap[posi]] > arr[heap[posp]])
			{
				myswap(posi, posp);
				posi = posp;
				posp = parent(posp);
			}
		}
		
		void print()
		{
			cout << "priority queue: " << endl;
			for(int i=0; i<height; i++)
			{
				cout << i << ": pos in element array: " << heap[i] << " parent element: " << arr[heap[parent(i)]] << " element: " << arr[heap[i]] 
				<< ( posinheap[heap[i]] == -1 ? "DELETED!" : "") << endl;
				
			}
		}
		
		T *arr;
		int *heap;
		int *posinheap;
		int height;	
		int size;
	};
*/


/*
 
 int main (int argc, char **argv) 
 {
 int bla[] = {2, 5, 7, 6, 1, 4, 3, 9, 8, 20};
 maxheap<int> pq(9, bla);
 for(int i=0; i<7; i++)
 {
 pq.insertItem(i);
 pq.print();
 cout << endl << endl;
 }
 pq.insertItem(7);
 bla[7] = -1;
 pq.remove(7);
 pq.insertItem(8);
 pq.print();
 pq.assureNotDeleted(0);
 pq.print();
 
 int i;
 do
 {
 
 i= pq.removeMax();
 
 cout << bla[i] << " ";
 bla[i] = -1;
 pq.assureNotDeleted(0);
 pq.print();
 } while(i != -1);
 return 3;
 }
 
 #if 0 // old
 */