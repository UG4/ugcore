/*
 *  maxheap.h
 *  flexamg
 *
 *  Created by Martin Rupp on 02.10.09.
 *  Copyright 2009 G-CSC. All rights reserved.
 *
 */
#pragma once
#if 0
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
#else
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
			int p = posinheap[index];
			int parentpos = (p == 0 ? 0 : (p+1)/2 -1);
			return heap[parentpos];
		}
		
		int leftchild(int index)
		{
			int p = posinheap[index];
			p = (p+1)*2 -1;
			if(p < height)
				return heap[p];
			else return index;
		}
		
		int rightchild(int index)
		{
			int p = posinheap[index];
			p = (p+1)*2 -1 +1;
			if(p < height)
				return heap[p];
			else return index;
		}
		
		void remove(int i)
		{
			if(posinheap[i] == -1) return;
			int j = heap[height-1];
			myswap(i, j);
			height--;
			downheap(j);
			posinheap[i] = -1;
		}
		
		int removeMax()
		{
			ASSERT2(height > 0, "heap already empty");
			myswap(heap[0], heap[height-1]);
			height--;
			downheap(heap[0]);
			return heap[height];
		}
		
		void myswap(int i, int j)
		{
			int posinheapi = posinheap[i];
			int posinheapj = posinheap[j];
			
			heap[posinheapi] = j;
			heap[posinheapj] = i;
			
			posinheap[i] = posinheapj;
			posinheap[j] = posinheapi;
		}
		
		void upheap(int i)
		{
			if(posinheap[i] == -1) return;
			while(arr[i] > arr[parent(i)])
				myswap(i, parent(i));
		}
		
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
		
		void print()
		{
			for(int i=0; i<height; i++)
			{
				cout << i << ": pos: " << heap[i] << " parent: " << parent(heap[i]) << " rating: " << arr[heap[i]].rating << (arr[heap[i]] > arr[parent(heap[i])] ? " ERR " : "");// << endl;
				cout << (arr[heap[i]].rating > arr[parent(heap[i])].rating ? " r> " : "") << (&arr[heap[i]] > &arr[parent(heap[i])] ? " a> " : "") << endl;
			}
			
		}
		
		T *arr;
		int *heap;
		int *posinheap;
		int height;	
		int size;
	};
#endif