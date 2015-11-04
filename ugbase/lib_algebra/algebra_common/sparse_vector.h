/*
 * sparse_vector.h
 *
 *  Created on: 13.05.2013
 *      Author: mrupp
 */

#ifndef SPARSE_VECTOR_H_
#define SPARSE_VECTOR_H_

#include <map>
namespace ug{

template<typename T>
class SparseVector
{
	size_t m_size;
	typedef std::map<size_t, T> container;
	container data;
public:
	typedef T value_type;
	class const_iterator : public container::const_iterator
	{
		using container::const_iterator::operator *;
	public:
		const_iterator(typename container::const_iterator it) : container::const_iterator(it) {}
		const T &value() const { return (operator *()).second; }
		size_t index() const { return (operator *()).first; }
	};

	class iterator : public container::iterator
	{
		using container::iterator::operator *;
	public:
		iterator(typename container::iterator it) : container::iterator(it) {}
		const T &value() const { return (operator *()).second; }
		T &value() { return (operator *())->second; }
		size_t index() const { return (operator *()).first; }
	};

	SparseVector(size_t s) : m_size(s)
	{
	}
	const_iterator begin() const
	{
		const_iterator c(data.begin());
		return c;
	}
	const_iterator end() const
	{
		return const_iterator(data.end());
	}

	const T &operator()(size_t c) const
	{
		assert(c < m_size);
		return data[c];
	}
	T &operator()(size_t c)
	{
		assert(c < m_size);
		return data[c];
	}
	bool has_connection(size_t c) const
	{
		return data.find(c) != data.end();
	}

	size_t size()
	{
		return m_size;
	}
	size_t num_connections() const
	{
		return data.size();
	}

	void print() const
	{
		for(const_iterator it=begin(); it != end(); ++it)
		{
			//if(BlockNorm(it.value()) == 0.0) continue;
			UG_LOG("(" << it.index() << " -> " << it.value() << ")");
		}

		UG_LOG("\n");
	}

};

}
#endif /* SPARSE_VECTOR_H_ */
