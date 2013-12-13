#ifndef __H__LIB_ALGEBRA__STL_DEBUG__
#define __H__LIB_ALGEBRA__STL_DEBUG__

#include <vector>

#ifndef NDEBUG
#include "common/error.h"

namespace ug{

template<typename T, class Allocator = std::allocator<T> >
class stdvector : public std::vector<T, Allocator>
{
private:
	typedef std::vector<T, Allocator> super;
public:
	typedef typename super::size_type size_type;

private:

	inline void size_check(size_t i) const
	{
		UG_COND_THROW(i >= super::size(), "accessing element " << i << " but vector only has size " << super::size() << ".")
	}

public:
	explicit stdvector(const Allocator& a= Allocator()) : super(a) { }
	explicit stdvector(size_type n, const T& value = T(), const Allocator& a= Allocator()) : super(n, value, a) { }

	stdvector(const std::vector<T, Allocator> &x) : super(x) { }
	stdvector(const stdvector<T, Allocator> &x) : super(x) { }

	inline typename super::reference operator[] (size_t i)
	{
		size_check(i);
		return super::operator[](i);
	}
	inline typename super::const_reference operator[] (size_t i) const
	{
		size_check(i);
		return super::operator[](i);
	}
};

}
#else

#define stdvector std::vector

#endif


#endif // __H__LIB_ALGEBRA__LIB_ALGEBRA__
