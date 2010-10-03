
#ifndef __H__UG_BRIDGE__FUNCTION_TRAITS__
#define __H__UG_BRIDGE__FUNCTION_TRAITS__

#include "common/metaprogramming_util.h"

namespace ug
{
namespace bridge
{

template <typename TFunc>
struct func_traits {};

////////////////////////////////
// global function traits
////////////////////////////////

template <typename TRet>
struct func_traits <TRet (*) ()>
{
	typedef TRet return_type;
	typedef TypeList<> params_type;
	static TRet apply(TRet (*fp)(),  TypeValueList<params_type>& args)
	{
		return fp();
	};
};


template <typename TRet, typename P1>
struct func_traits <TRet (*) (P1)>
{
	typedef TRet return_type;
	typedef TypeList<P1> params_type;
	static TRet apply(TRet (*fp)(P1),  TypeValueList<params_type>& args)
	{
		return fp(args.hd);
	};
};


template <typename TRet, typename T1, typename T2>
struct func_traits <TRet (*) (T1, T2)>
{
	typedef TRet return_type;
	typedef TypeList<T1, T2> params_type;
	static TRet apply(TRet (*fp)(T1, T2),  TypeValueList<params_type>& args)
	{
		return fp(args.hd, args.tl.hd);
	};
};

template <typename TRet, typename T1, typename T2, typename T3>
struct func_traits <TRet (*) (T1, T2, T3)>
{
	typedef TRet return_type;
	typedef TypeList<T1, T2, T3> params_type;
	static TRet apply(TRet (*fp)(T1, T2, T3),  TypeValueList<params_type>& args)
	{
		return fp(args.hd, args.tl.hd, args.tl.tl.hd);
	};
};

template <typename TRet, typename T1, typename T2, typename T3,
		  typename T4>
struct func_traits <TRet (*) (T1, T2, T3, T4)>
{
	typedef TRet return_type;
	typedef TypeList<T1, T2, T3, T4> params_type;
	static TRet apply(TRet (*fp)(T1, T2, T3, T4),  TypeValueList<params_type>& args)
	{
		return fp(args.hd, args.tl.hd, args.tl.tl.hd, args.tl.tl.tl.hd);
	};
};

template <typename TRet, typename T1, typename T2, typename T3,
		  typename T4, typename T5>
struct func_traits <TRet (*) (T1, T2, T3, T4, T5)>
{
	typedef TRet return_type;
	typedef TypeList<T1, T2, T3> params_type;
	static TRet apply(TRet (*fp)(T1, T2, T3, T4, T5),  TypeValueList<params_type>& args)
	{
		return fp(args.hd, args.tl.hd, args.tl.tl.hd, args.tl.tl.tl.hd,
				 args.tl.tl.tl.tl.hd);
	};
};

template <typename TRet, typename T1, typename T2, typename T3,
		  typename T4, typename T5, typename T6>
struct func_traits <TRet (*) (T1, T2, T3, T4, T5, T6)>
{
	typedef TRet return_type;
	typedef TypeList<T1, T2, T3> params_type;
	static TRet apply(TRet (*fp)(T1, T2, T3, T4, T5, T6),  TypeValueList<params_type>& args)
	{
		return fp(args.hd, args.tl.hd, args.tl.tl.hd, args.tl.tl.tl.hd,
				 args.tl.tl.tl.tl.hd, args.tl.tl.tl.tl.tl.hd);
	};
};

//todo: implement more ...

////////////////////////////////
// non-const method traits
////////////////////////////////

#define FUNC_TRAITS_GENERAL_NON_CONST_MEMBER \
	static const bool const_method = false;\
	typedef TClass class_type;\
	typedef TRet return_type


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
		return (obj->*fp)(args.hd);
	};
};


template <typename TClass, typename TRet, typename T1, typename T2>
struct func_traits <TRet (TClass::*) (T1, T2)>
{
	FUNC_TRAITS_GENERAL_NON_CONST_MEMBER;
	typedef TypeList<T1, T2> params_type;
	static TRet apply(TRet (TClass::*fp)(T1, T2), TClass* obj, TypeValueList<params_type>& args)
	{
		return (obj->*fp)(args.hd, args.tl.hd);
	};
};

template <typename TClass, typename TRet, typename T1, typename T2, typename T3>
struct func_traits <TRet (TClass::*) (T1, T2, T3)>
{
	FUNC_TRAITS_GENERAL_NON_CONST_MEMBER;
	typedef TypeList<T1, T2, T3> params_type;
	static TRet apply(TRet (TClass::*fp)(T1, T2, T3), TClass* obj, TypeValueList<params_type>& args)
	{
		return (obj->*fp)(args.hd, args.tl.hd, args.tl.tl.hd);
	};
};

////////////////////////////////
// const method traits
////////////////////////////////

#define FUNC_TRAITS_GENERAL_CONST_MEMBER \
	static const bool const_method = true;\
	typedef TClass class_type;\
	typedef TRet return_type


template <typename TClass, typename TRet>
struct func_traits <TRet (TClass::*) () const>
{
	FUNC_TRAITS_GENERAL_CONST_MEMBER;
	typedef TypeList<> params_type;
	static TRet apply(TRet (TClass::*fp)() const, const TClass* obj, TypeValueList<params_type>& args)
	{
		return (obj->*fp)();
	};
};


template <typename TClass, typename TRet, typename P1>
struct func_traits <TRet (TClass::*) (P1) const>
{
	FUNC_TRAITS_GENERAL_CONST_MEMBER;
	typedef TypeList<P1> params_type;
	static TRet apply(TRet (TClass::*fp)(P1) const, const TClass* obj, TypeValueList<params_type>& args)
	{
		return (obj->*fp)(args.hd);
	};
};


template <typename TClass, typename TRet, typename T1, typename T2>
struct func_traits <TRet (TClass::*) (T1, T2) const>
{
	FUNC_TRAITS_GENERAL_CONST_MEMBER;
	typedef TypeList<T1, T2> params_type;
	static TRet apply(TRet (TClass::*fp)(T1, T2) const, const TClass* obj, TypeValueList<params_type>& args)
	{
		return (obj->*fp)(args.hd, args.tl.hd);
	};
};

template <typename TClass, typename TRet, typename T1, typename T2, typename T3>
struct func_traits <TRet (TClass::*) (T1, T2, T3) const>
{
	FUNC_TRAITS_GENERAL_CONST_MEMBER;
	typedef TypeList<T1, T2, T3> params_type;
	static TRet apply(TRet (TClass::*fp)(T1, T2, T3) const, const TClass* obj, TypeValueList<params_type>& args)
	{
		return (obj->*fp)(args.hd, args.tl.hd, args.tl.tl.hd);
	};
};

} // end namespace bridge
} // end namespace ug

#endif /* __H__UG_BRIDGE__FUNCTION_TRAITS__ */
