#include <vector>

#ifndef SMART_PTR_VECTOR_H
#define SMART_PTR_VECTOR_H

namespace ug{
template<typename T>
class SmartPtrVector : public std::vector<SmartPtr<T> >
{
	typedef std::vector<SmartPtr<T> > super;
private:
	SmartPtrVector(SmartPtrVector &);
public:
	SmartPtrVector() : super() {}
	SmartPtrVector(size_t size) : super(size)
	{
	}

	const T &operator () (size_t i) const
	{
		assert(i < super::size());
		return *super::operator[](i);
	}
	T &operator () (size_t i)
	{
		assert(i < super::size());
		return *super::operator[](i);
	}
};

/*
{
private:
	SmartPtrVector(SmartPtrVector &);
public:
	SmartPtrVector();
	SmartPtrVector(size_t size) : m_data(size)
	{
	}

	T &operator [] (size_t i)
	{
		assert(i < m_data.size());
		return *m_data[i];
	}

	const T &operator [] (size_t i) const
	{
		assert(i < m_data.size());
		return *m_data[i];
	}

	SmartPtr<T> get_smart_ptr(size_t i)
	{
		assert(i < m_data.size());
		return m_data[i];
	}


	void push_back(SmartPtr<T> t)
	{
		m_data.push_back(t);
	}

	size_t size() const
	{
		return m_data.size();
	}

	void clear()
	{
		m_data.clear();
	}

private:
	std::vector<SmartPtr<T> > m_data;
};*/




}


#endif
