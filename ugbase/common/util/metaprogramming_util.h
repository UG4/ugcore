/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#ifndef __H__COMMON__METAPROGRAMMING_UTIL__
#define __H__COMMON__METAPROGRAMMING_UTIL__

namespace ug {

/// \addtogroup ugbase_common_util
/// \{

template <int N>
struct Int2Type {
	enum{ value = N};

	using value_type = int;
};

template <typename T>
struct Pointer2Value{};

template <typename T>
struct Pointer2Value<T*>{
	using type = T;
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
  typename T11=EmptyType,
  typename T12=EmptyType
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
  typename T11,
  typename T12
>
struct TypeList
{
	using head = T1;
	using tail = TypeList< T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12 >;
	enum{length = tail::length+1};
};

// empty typelist specialization
template<>
struct TypeList< EmptyType, EmptyType, EmptyType, EmptyType,
				 EmptyType, EmptyType, EmptyType, EmptyType, EmptyType,
				 EmptyType, EmptyType, EmptyType>
{
  enum{length = 0};
};

//////////////////////////////
// TypeValueList
//////////////////////////////

// TypeList
template <typename TTypeList> struct TypeValueList
{
	using head = typename TTypeList::head;
	using tail = typename TTypeList::tail;

	head hd;
	TypeValueList<tail> tl;

	explicit TypeValueList() {
		std::cerr << "incomplete TVL\n";
	}

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
    static constexpr size_t value = n*Factorial<n-1>::value;
};

template <>
struct Factorial<1>
{
    static constexpr size_t value = 1;
};

//////////////////////////////
// Pow
//////////////////////////////

/** returns the power of n^d
 */
template <int n, size_t d>
struct Pow
{
     static constexpr int value = n*Pow<n, d-1>::value;
};

template <int n>
struct Pow<n, 0>
{
    static constexpr int value = 1;
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
    static constexpr size_t value = Factorial<n>::value/
    				(Factorial<k>::value*Factorial<n-k>::value);
};

// end rekursion
template <size_t n>
struct BinomialCoefficient<n,0>
{
    static constexpr size_t value = 1;
};
// end rekursion
template <size_t n>
struct BinomialCoefficient<n,-1>
{
    static constexpr size_t value = 0;
};
// end rekursion
template <size_t n>
struct BinomialCoefficient<n,-2>
{
    static constexpr size_t value = 0;
};
// end rekursion
template <size_t n>
struct BinomialCoefficient<n,-3>
{
    static constexpr size_t value = 0;
};
// end rekursion
template <size_t n>
struct BinomialCoefficient<n,-4>
{
    static constexpr size_t value = 0;
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
template <typename TType>
size_t GetUniqueTypeID()
{
	static size_t typeID = UniqueTypeIDProvider::inst().new_id();
	return typeID;
}

// end group ugbase_common_util
/// \}

}

#endif
