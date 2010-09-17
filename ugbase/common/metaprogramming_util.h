// author: andreasvogel

#ifndef __H__COMMON__METAPROGRAMMING_UTIL__
#define __H__COMMON__METAPROGRAMMING_UTIL__

namespace ug {

template <int N>
struct Int2Type {
	enum{ value = N};
	typedef int value_type;
};



}

#endif /* __H__COMMON__METAPROGRAMMING_UTIL__ */
