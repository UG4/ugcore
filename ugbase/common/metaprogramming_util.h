// author: andreasvogel

#ifndef __H__COMMON__METAPROGRAMMING_UTIL__
#define __H__COMMON__METAPROGRAMMING_UTIL__

namespace ug {

template <int N>
struct Int2Type {
	enum{ value = N};
	typedef int value_type;
};

//////////////////////////////
// TypeList
//////////////////////////////

// empty type
struct EmptyType {};

// TypeList
template
<
  typename T1=EmptyType,
  typename T2=EmptyType,
  typename T3=EmptyType,
  typename T4=EmptyType,
  typename T5=EmptyType
> struct TypeList;

// implementation of TypeList
template
<
  typename T1,
  typename T2,
  typename T3,
  typename T4,
  typename T5
>
struct TypeList
{
  typedef T1 head;
  typedef TypeList< T2, T3, T4, T5 > tail;
  enum{length = tail::length+1};
};

// empty typelist specialization
template<>
struct TypeList< EmptyType, EmptyType, EmptyType, EmptyType >
{
  enum{length = 0};
};

//////////////////////////////
// TypeValueList
//////////////////////////////

// TypeList
template <typename TTypeList> struct TypeValueList
{
	typedef typename TTypeList::head head;
	typedef typename TTypeList::tail tail;

	head hd;
	TypeValueList<tail> tl;

	TypeValueList(const head& _hd,
				  const TypeValueList<tail>& typValList) :
		hd(_hd), tl(typValList)	{}

};

template <>
struct TypeValueList< TypeList<> > {};


}

#endif /* __H__COMMON__METAPROGRAMMING_UTIL__ */
