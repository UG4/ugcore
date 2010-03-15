/*
 * multi_indices.h
 *
 *  Created on: 22.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__MULTI_INDEX__MULTI_INDICES__
#define __H__LIBDISCRETIZATION__MULTI_INDEX__MULTI_INDICES__

#include <vector>
#include <iostream>

namespace ug{

template<int N>
class MultiIndex
{
	public:
		typedef uint size_type;
		typedef uint single_index_type;

	public:
		inline size_type num_index()
		{ return N; }

		inline single_index_type& operator[] (size_type i)
		{
			assert(i < N);
			return m_indices[i];
		}

	private:
		single_index_type m_indices[N];
};

class IndexInfo
{
public:
	void set_leaf(bool b)
	{
		m_isLeaf = b;
		if(m_isLeaf)
		{
			m_num_index = 0;
			m_num_comp = 0;
			m_IndexInfoList.clear();
		}
		else
		{
			m_num_index = 0;
			m_IndexInfoList.resize(0);
		}
	}

	void set_num_index(std::size_t n)
	{
		m_num_index = n;
		if(!m_isLeaf)
			m_IndexInfoList.resize(n);
	}

	bool set_num_comp(std::size_t n)
	{
		if(!m_isLeaf) return false;
		m_num_comp = n;
		return true;
	}

	std::size_t num_index()
	{
		return m_num_index;
	}

	std::size_t num_comp()
	{
		return m_num_comp;
	}

	bool is_leaf()
	{
		return m_isLeaf;
	}

	IndexInfo& get_index_info(uint i)
	{
		assert(i < m_IndexInfoList.size());
		assert(!m_isLeaf);
		return m_IndexInfoList[i];
	}

	bool print_info(int offset = 4) const
	{
		using namespace std;

		if(!m_isLeaf)
		{
			std::cout << m_num_index << " Children: " << std::endl;
			for(uint i = 0; i < m_IndexInfoList.size(); ++i)
			{
				for(int j=0; j<offset; j++) std::cout << " ";
				std::cout << "#" << i <<": ";
				m_IndexInfoList[i].print_info(offset + 4);
			}
			return true;
		}
		else
		{
			std::cout << "is BigMatrixLeaf: num_index = " << m_num_index << ", num_comp = " << m_num_comp << "." << std::endl;
			return true;
		}
		return false;
	}

private:
	bool m_isLeaf;
	std::size_t m_num_index;
	std::size_t m_num_comp; // for leaf
	std::vector<IndexInfo> m_IndexInfoList;
};





}


#endif /* __H__LIBDISCRETIZATION__MULTI_INDEX__MULTI_INDICES__ */
