
#ifndef __H__UG_INTERFACE__UGBRIDGE__TYPELIST__
#define __H__UG_INTERFACE__UGBRIDGE__TYPELIST__

#include "common/metaprogramming_util.h"

namespace ug {

namespace interface{


template <typename TFunc>
struct func_traits {};

////////////////////////////////
// global function traits
////////////////////////////////

template <typename TRet>
struct func_traits <TRet (*) ()>
{
	typedef TRet result_type;
	typedef TypeList<> params_type;
	static TRet apply(TRet (*fp)(),  TypeValueList<params_type>& args)
	{
		return fp();
	};
};


template <typename TRet, typename P1>
struct func_traits <TRet (*) (P1)>
{
	typedef TRet result_type;
	typedef TypeList<P1> params_type;
	static TRet apply(TRet (*fp)(P1),  TypeValueList<params_type>& args)
	{
		return fp(args.head);
	};
};


template <typename TRet, typename T1, typename T2>
struct func_traits <TRet (*) (T1, T2)>
{
	typedef TRet result_type;
	typedef TypeList<T1, T2> params_type;
	static TRet apply(TRet (*fp)(T1, T2),  TypeValueList<params_type>& args)
	{
		return fp(args.head, args.tail.head);
	};
};

//todo: implement more ...

////////////////////////////////
// non-const method traits
////////////////////////////////

#define FUNC_TRAITS_GENERAL_NON_CONST_MEMBER \
	static const bool const_method = false;\
	typedef TClass class_type;\
	typedef TRet result_type


template <typename TClass, typename TRet>
struct func_traits <TRet (TClass::*) ()>
{
	FUNC_TRAITS_GENERAL_NON_CONST_MEMBER;
	typedef TypeList<> params_type;
	static TRet apply(TRet (TClass::*fp)(), TClass* obj, TypeValueList<params_type>& args)
	{
		return (obj->*fp)();
	};
};


template <typename TClass, typename TRet, typename P1>
struct func_traits <TRet (TClass::*) (P1)>
{
	FUNC_TRAITS_GENERAL_NON_CONST_MEMBER;
	typedef TypeList<P1> params_type;
	static TRet apply(TRet (TClass::*fp)(P1), TClass* obj, TypeValueList<params_type>& args)
	{
		return (obj->*fp)(args.head);
	};
};


template <typename TClass, typename TRet, typename T1, typename T2>
struct func_traits <TRet (TClass::*) (T1, T2)>
{
	FUNC_TRAITS_GENERAL_NON_CONST_MEMBER;
	typedef TypeList<T1, T2> params_type;
	static TRet apply(TRet (TClass::*fp)(T1, T2), TClass* obj, TypeValueList<params_type>& args)
	{
		return (obj->*fp)(args.head, args.tail.head);
	};
};


} // end namespace interface

} // end namespace ug

#endif /* __H__UG_INTERFACE__UGBRIDGE__TYPELIST__ */
