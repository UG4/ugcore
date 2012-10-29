// author: andreasvogel

#ifndef __H__COMMON__METAPROGRAMMING_UTIL__
#define __H__COMMON__METAPROGRAMMING_UTIL__

namespace ug {

template <int N>
struct Int2Type {
	enum{ value = N};
	typedef int value_type;
};

template <class T>
struct Pointer2Value{};

template <class T>
struct Pointer2Value<T*>{
	typedef T type;
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
  typename T5=EmptyType,
  typename T6=EmptyType,
  typename T7=EmptyType,
  typename T8=EmptyType,
  typename T9=EmptyType,
  typename T10=EmptyType,
  typename T11=EmptyType
> struct TypeList;

// implementation of TypeList
template
<
  typename T1,
  typename T2,
  typename T3,
  typename T4,
  typename T5,
  typename T6,
  typename T7,
  typename T8,
  typename T9,
  typename T10,
  typename T11
>
struct TypeList
{
  typedef T1 head;
  typedef TypeList< T2, T3, T4, T5, T6, T7, T8, T9, T10, T11 > tail;
  enum{length = tail::length+1};
};

// empty typelist specialization
template<>
struct TypeList< EmptyType, EmptyType, EmptyType, EmptyType,
				 EmptyType, EmptyType, EmptyType, EmptyType, EmptyType,
				 EmptyType, EmptyType>
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

	TypeValueList(head _hd,
				  TypeValueList<tail> typValList) :
		hd(_hd), tl(typValList)	{}

};

template <>
struct TypeValueList< TypeList<> > {};


//////////////////////////////
// Factorial
//////////////////////////////

/** returns the factorial of n
 * The struct value is n!
 */
template <size_t n>
struct Factorial
{
    enum{value = n*Factorial<n-1>::value};
};

template <>
struct Factorial<1>
{
    enum {value = 1};
};

//////////////////////////////
// Pow
//////////////////////////////

/** returns the power of n^d
 */
template <int n, size_t d>
struct Pow
{
    enum{value = n*Pow<n, d-1>::value};
};

template <int n>
struct Pow<n, 0>
{
    enum {value = 1};
};

//////////////////////////////
// BinomialCoefficient
//////////////////////////////

/** returns static value of binomial coefficient
 * The struct value is:
 *
 * 	  n!
 * ---------
 * k! (n-k)!
 */
template <size_t n, int k>
struct BinomialCoefficient
{
    enum { value = 	Factorial<n>::value/
    				(Factorial<k>::value*Factorial<n-k>::value) };
};

// end rekursion
template <size_t n>
struct BinomialCoefficient<n,0>
{
    enum { value = 1};
};
// end rekursion
template <size_t n>
struct BinomialCoefficient<n,-1>
{
    enum { value = 0};
};
// end rekursion
template <size_t n>
struct BinomialCoefficient<n,-2>
{
    enum { value = 0};
};
// end rekursion
template <size_t n>
struct BinomialCoefficient<n,-3>
{
    enum { value = 0};
};
// end rekursion
template <size_t n>
struct BinomialCoefficient<n,-4>
{
    enum { value = 0};
};

//////////////////////////////
// UniqueTypeID (sreiter)
//////////////////////////////
///	a singleton class that returns a new id for each type
class UniqueTypeIDProvider{
	public:
		static UniqueTypeIDProvider& inst(){
			static UniqueTypeIDProvider utid;
			return utid;
		}

		size_t new_id()	{return ++m_id;}

	private:
		UniqueTypeIDProvider() : m_id(0)	{}
		size_t m_id;
};

///	This method associated a unique unsigned integer value with each type.
template <class TType>
size_t GetUniqueTypeID()
{
	static size_t typeID = UniqueTypeIDProvider::inst().new_id();
	return typeID;
}

}

#endif /* __H__COMMON__METAPROGRAMMING_UTIL__ */
