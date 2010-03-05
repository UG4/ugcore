/*
 *  interface.h
 *  flexamg
 *
 *  Created by Martin Rupp on 24.02.10.
 *  Copyright 2010 G-CSC, University of Frankfurt. All rights reserved.
 *
 */
#pragma once

#include <iomanip>
typedef unsigned int uint;

template<int N>
class MultiIndex
{
public:
	typedef uint size_type;
	typedef uint single_index_type;
	
public:
	MultiIndex()
	{
		FORCE_CREATION { print(); }
		for(int i=0; i<N; i++) m_indices[i] = 0;
	}
	inline size_type num_index()
	{ return N; }
	
	inline single_index_type &operator[] (size_type i)
	{
		assert(i < N);
		return m_indices[i];
	}

	inline const single_index_type &operator[] (size_type i) const
	{
		assert(i < N);
		return m_indices[i];
	}
	
	void print()
	{
		cout << *this << endl;
	}
	
	friend ostream &operator << (ostream &out, MultiIndex<N> &index)
	{
		out << "[";
		for(int i=0; i<N; i++)
			out << index[i] << " ";
		out << "]";
		return out;
	}
	
private:
	single_index_type m_indices[N];
};



class IndexInfo
{
public:
	IndexInfo() : m_IndexInfoList()
	{
		m_isLeaf = false;
		m_num_comp = 0;
		m_num_index = 0;
	}

	void set_num_index(std::size_t n)
	{
		m_num_index = n;
		if(!is_leaf())
			m_IndexInfoList.resize(n);
	}
	
	void set_num_comp(std::size_t n)
	{		
		m_isLeaf = true;
		m_num_comp = n;
	}
	
	std::size_t num_index() const
	{
		return m_num_index;
	}
	
	std::size_t num_comp() const
	{
		return m_num_comp;
	}
	
	bool is_leaf() const
	{
		return m_isLeaf;
	}
	
	IndexInfo& get_index_info(uint i)
	{
		assert(i < m_IndexInfoList.size());
		assert(!is_leaf());
		return m_IndexInfoList[i];
	}
	const IndexInfo& get_index_info(uint i) const
	{
		assert(i < m_IndexInfoList.size());
		assert(!is_leaf());
		return m_IndexInfoList[i];
	}
	
	bool print_info(int offset = 0) const
	{
		using namespace std;
		
		if(!is_leaf())
		{
			cout  << setw(offset) << m_num_index << " Children: " << std::endl;
			for(uint i = 0; i < m_IndexInfoList.size(); ++i)
			{
				cout << setw(offset) << "#" << i <<": ";
				m_IndexInfoList[i].print_info(offset + 4);
			}
			return true;
		}
		else
		{
			cout  << setw(offset) << "is BigMatrixLeaf: num_index = " << m_num_index << ", num_comp = " << m_num_comp << "." << std::endl;
			return true;
		}
		return false;
	}
	
private:
	bool m_isLeaf;
	std::size_t m_num_comp; // for leaf
	std::size_t m_num_index;
	std::vector<IndexInfo> m_IndexInfoList;
};