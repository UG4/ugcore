/*
 * operations_mat.h
 *
 *  Created on: 17.11.2010
 *      Author: mrupp
 */

#ifndef __H__UG__LIB_ALGEBRA__INDEX_VIEW_H__
#define __H__UG__LIB_ALGEBRA__INDEX_VIEW_H__

namespace ug
{


template<typename index_type>
class IndexView
{
public:
	IndexView(const std::vector<index_type> &indices)
	{
		local_to_global = indices;
		init();
	}

	IndexView(index_type *indices, size_t size)
	{
		local_to_global.resize(size);
		for(size_t i=0; i<size; i++)
			local_to_global[i] = indices[i];
		init();
	}

	typedef std::vector<index_type> iterator;

	iterator begin()
	{
		return local_to_global.begin();
	}

	iterator end()
	{
		return local_to_global.end();
	}

	bool is_in_view(index_type i)
	{
		// we need this function to be very fast, because it is used often in matrix-vector-multiplication
		// this uses memory, but is faster than searching the local_to_global array, even with binary search
		if(i >= is_in_locals.size()) return false;
		return is_in_locals[i];
	}

private:
	void init()
	{
		index_type maxindex = (*max_element(local_to_global.begin(), local_to_global.end()));
		is_in_locals.clear();
		is_in_locals.resize(maxindex+1, false);
		for(size_t i=0; i<local_to_global.size(); i++)
			is_in_locals[local_to_global[i]] = true;
	}
	std::vector<index_type> local_to_global;
	std::vector<bool> is_in_locals;
	//std::vector<int> global_to_lobal;
};
template<typename index_type>
class SliceIndexView
{
public:
	SliceIndexView(index_type from, index_type to) : m_from(from), m_to(to) { }

	struct iterator
	{
		iterator(index_type from) : i(from) { }
		index_type i;

		bool operator != (const iterator &other) const { return i != other.i; }
		bool operator == (const iterator &other) const { return i == other.i; }

		void operator ++() { ++i; }
		int operator *() const { return i; }
	};

	inline iterator begin()
	{
		return iterator(m_from);
	}

	inline iterator end()
	{
		return iterator(m_to);
	}

	bool is_in_view(index_type i)
	{
		return i >= m_from && i < m_to;
	}

private:
	index_type m_from;
	index_type m_to;
};

/*
template<typename vector_t, typename view_t>
class VectorView
{
	VectorView(vector_t &vector, view_t &view) : m_vector(vector), m_view(view)
	{	}

	typedef view_t::iterator index_iterator;
	index_iterator begin() { return m_view.begin(); }
	index_iterator end() { return m_view.end(); }

	typename vector_t::value_type operator[] (const index_iterator &it)
	{
		return vector_t::operator[] ( (*it) );
	}
}*/


/*
//! calculates dest = alpha1*v1
template<typename vector_t>
inline void VecScaleAssign(vector_t &dest, double alpha1, const vector_t &v1)
{
	for(dest::index_iterator it = dest.begin(); it != view.end(); ++it)
	{
	// so?
		size_t i = (*it);
		VecScaleAssign(dest[i], alpha1, v1[i]);
	// oder so
	  	VecScaleAssign(dest[it], alpha1, v1[it]);
	}
}

oder

template<typename vector_t, typename view_t>
inline void VecScaleAssign(vector_t &dest, double alpha1, const vector_t &v1, const view_t &view)
{
	for(view_t::index_iterator it = view.begin(); it != view.end(); ++it)
	{
		size_t i = (*it);
		VecScaleAssign(dest[i], alpha1, v1[i]);
	}
}
*/


template<typename index_type>
class BlockSliceIndexView
{
public:
	BlockSliceIndexView(index_type from, index_type to, index_type block) : m_from(from), m_to(to), m_block(block) { }

	struct iterator
	{
		iterator(index_type from, index_type block) : i(from), m_block(block) { }
		index_type i;
		index_type m_block;

		bool operator != (const iterator &other) const { return i != other.i; }
		bool operator == (const iterator &other) const { return i == other.i; }

		void operator ++() { ++i; }
		int index() const { return i; }
		int block() const { return block; }
	};

	inline iterator begin()
	{
		return iterator(m_from, m_block);
	}

	inline iterator end()
	{
		return iterator(m_to, m_block);
	}

	bool is_in_view(const iterator &i)
	{
		return is_in_view(i.i);
	}
	bool is_in_view(index_type i)
	{
		return i >= m_from && i < m_to;
	}

private:
	index_type m_from;
	index_type m_to;
};

/*
template<typename vector_t, typename view_t>
class BlockVectorView
{
	BlockVectorView(vector_t &vector, view_t &view) : m_vector(vector), m_view(view)
	{	}

	typedef view_t::iterator index_iterator;
	index_iterator begin() { return m_view.begin(); }
	index_iterator end() { return m_view.end(); }

	typename vector_t::value_type::value_type operator[] (const index_iterator &it)
	{
		return BlockRef(vector_t::operator[] ( it.index() ), it.block());
	}
}*/

//! calculates dest = alpha1*v1
/*template<typename vector_t>
inline void VecScaleAssign(vector_t &dest, double alpha1, const vector_t &v1)
{
	for(dest::index_iterator it = dest.begin(); it != view.end(); ++it)
	  	VecScaleAssign(dest[it], alpha1, v1[it]); // automatically on BlockRef(index, block) if used with BlockVectorView(vec, BlockSliceIndexView(10, 20, 2))
}
*/

} // namespace ug


#endif // __H__UG__LIB_ALGEBRA__INDEX_VIEW_H__
