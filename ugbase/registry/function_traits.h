/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
 * Authors: Sebastian Reiter, Andreas Vogel
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

#ifndef __H__UG_BRIDGE__FUNCTION_TRAITS__
#define __H__UG_BRIDGE__FUNCTION_TRAITS__

#include "common/util/metaprogramming_util.h"

/**
 * Maximum number of template arguments for a bridge traits class. These are
 * the maximal available function arguments.
 *
 * NOTE: IF YOU INCREASE THE NUMBER OF TEMPLATES BELOW, THIS NUMBER MUST BE
 * 		 ADJUSTED AS WELL.
 */
const int UG_REGISTRY_MAX_NUM_ARGS = 12;

namespace ug
{
namespace bridge
{

/// \addtogroup registry
/// \{

template <typename TFunc>
struct func_traits {};

////////////////////////////////////////////////////////////////////////////////
// global function traits
////////////////////////////////////////////////////////////////////////////////

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
	typedef TypeList<T1, T2, T3, T4, T5> params_type;
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
	typedef TypeList<T1, T2, T3, T4, T5, T6> params_type;
	static TRet apply(TRet (*fp)(T1, T2, T3, T4, T5, T6),  TypeValueList<params_type>& args)
	{
		return fp(args.hd, args.tl.hd, args.tl.tl.hd, args.tl.tl.tl.hd,
				 args.tl.tl.tl.tl.hd, args.tl.tl.tl.tl.tl.hd);
	};
};

template <typename TRet, typename T1, typename T2, typename T3,
		  typename T4, typename T5, typename T6, typename T7>
struct func_traits <TRet (*) (T1, T2, T3, T4, T5, T6, T7)>
{
	typedef TRet return_type;
	typedef TypeList<T1, T2, T3, T4, T5, T6, T7> params_type;
	static TRet apply(TRet (*fp)(T1, T2, T3, T4, T5, T6, T7),  TypeValueList<params_type>& args)
	{
		return fp(args.hd, args.tl.hd, args.tl.tl.hd, args.tl.tl.tl.hd,
				 args.tl.tl.tl.tl.hd, args.tl.tl.tl.tl.tl.hd, args.tl.tl.tl.tl.tl.tl.hd);
	};
};

template <typename TRet, typename T1, typename T2, typename T3,
		  typename T4, typename T5, typename T6, typename T7, typename T8>
struct func_traits <TRet (*) (T1, T2, T3, T4, T5, T6, T7, T8)>
{
	typedef TRet return_type;
	typedef TypeList<T1, T2, T3, T4, T5, T6, T7, T8> params_type;
	static TRet apply(TRet (*fp)(T1, T2, T3, T4, T5, T6, T7, T8),  TypeValueList<params_type>& args)
	{
		return fp(args.hd, args.tl.hd, args.tl.tl.hd, args.tl.tl.tl.hd,
				 args.tl.tl.tl.tl.hd, args.tl.tl.tl.tl.tl.hd, args.tl.tl.tl.tl.tl.tl.hd,
				 args.tl.tl.tl.tl.tl.tl.tl.hd);
	};
};

template <typename TRet, typename T1, typename T2, typename T3,
		  typename T4, typename T5, typename T6, typename T7, typename T8, typename T9>
struct func_traits <TRet (*) (T1, T2, T3, T4, T5, T6, T7, T8, T9)>
{
	typedef TRet return_type;
	typedef TypeList<T1, T2, T3, T4, T5, T6, T7, T8, T9> params_type;
	static TRet apply(TRet (*fp)(T1, T2, T3, T4, T5, T6, T7, T8, T9),  TypeValueList<params_type>& args)
	{
		return fp(args.hd, args.tl.hd, args.tl.tl.hd, args.tl.tl.tl.hd,
				 args.tl.tl.tl.tl.hd, args.tl.tl.tl.tl.tl.hd, args.tl.tl.tl.tl.tl.tl.hd,
				 args.tl.tl.tl.tl.tl.tl.tl.hd, args.tl.tl.tl.tl.tl.tl.tl.tl.hd);
	};
};

template <typename TRet, typename T1, typename T2, typename T3,
		  typename T4, typename T5, typename T6, typename T7, typename T8, typename T9,
		  typename T10>
struct func_traits <TRet (*) (T1, T2, T3, T4, T5, T6, T7, T8, T9, T10)>
{
	typedef TRet return_type;
	typedef TypeList<T1, T2, T3, T4, T5, T6, T7, T8, T9, T10> params_type;
	static TRet apply(TRet (*fp)(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10),  TypeValueList<params_type>& args)
	{
		return fp(args.hd, args.tl.hd, args.tl.tl.hd, args.tl.tl.tl.hd,
				 args.tl.tl.tl.tl.hd, args.tl.tl.tl.tl.tl.hd, args.tl.tl.tl.tl.tl.tl.hd,
				 args.tl.tl.tl.tl.tl.tl.tl.hd, args.tl.tl.tl.tl.tl.tl.tl.tl.hd,
				 args.tl.tl.tl.tl.tl.tl.tl.tl.tl.hd);
	};
};

template <typename TRet, typename T1, typename T2, typename T3,
		  typename T4, typename T5, typename T6, typename T7, typename T8, typename T9,
		  typename T10, typename T11>
struct func_traits <TRet (*) (T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11)>
{
	typedef TRet return_type;
	typedef TypeList<T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11> params_type;
	static TRet apply(TRet (*fp)(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11),  TypeValueList<params_type>& args)
	{
		return fp(args.hd, args.tl.hd, args.tl.tl.hd, args.tl.tl.tl.hd,
				 args.tl.tl.tl.tl.hd, args.tl.tl.tl.tl.tl.hd, args.tl.tl.tl.tl.tl.tl.hd,
				 args.tl.tl.tl.tl.tl.tl.tl.hd, args.tl.tl.tl.tl.tl.tl.tl.tl.hd,
				 args.tl.tl.tl.tl.tl.tl.tl.tl.tl.hd, args.tl.tl.tl.tl.tl.tl.tl.tl.tl.tl.hd);
	};
};

template <typename TRet, typename T1, typename T2, typename T3,
		  typename T4, typename T5, typename T6, typename T7, typename T8, typename T9,
		  typename T10, typename T11, typename T12>
struct func_traits <TRet (*) (T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12)>
{
	typedef TRet return_type;
	typedef TypeList<T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12> params_type;
	static TRet apply(TRet (*fp)(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12),  TypeValueList<params_type>& args)
	{
		return fp(args.hd, args.tl.hd, args.tl.tl.hd, args.tl.tl.tl.hd,
				 args.tl.tl.tl.tl.hd, args.tl.tl.tl.tl.tl.hd, args.tl.tl.tl.tl.tl.tl.hd,
				 args.tl.tl.tl.tl.tl.tl.tl.hd, args.tl.tl.tl.tl.tl.tl.tl.tl.hd,
				 args.tl.tl.tl.tl.tl.tl.tl.tl.tl.hd, args.tl.tl.tl.tl.tl.tl.tl.tl.tl.tl.hd,
         args.tl.tl.tl.tl.tl.tl.tl.tl.tl.tl.tl.hd);
	};
};

////////////////////////////////////////////////////////////////////////////////
// non-const method traits
////////////////////////////////////////////////////////////////////////////////

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

template <typename TClass, typename TRet, typename T1, typename T2, typename T3,
			typename T4>
struct func_traits <TRet (TClass::*) (T1, T2, T3, T4)>
{
	FUNC_TRAITS_GENERAL_NON_CONST_MEMBER;
	typedef TypeList<T1, T2, T3, T4> params_type;
	static TRet apply(TRet (TClass::*fp)(T1, T2, T3, T4), TClass* obj, TypeValueList<params_type>& args)
	{
		return (obj->*fp)(args.hd, args.tl.hd, args.tl.tl.hd, args.tl.tl.tl.hd);
	};
};

template <typename TClass, typename TRet, typename T1, typename T2, typename T3,
			typename T4, typename T5>
struct func_traits <TRet (TClass::*) (T1, T2, T3, T4, T5)>
{
	FUNC_TRAITS_GENERAL_NON_CONST_MEMBER;
	typedef TypeList<T1, T2, T3, T4, T5> params_type;
	static TRet apply(TRet (TClass::*fp)(T1, T2, T3, T4, T5), TClass* obj, TypeValueList<params_type>& args)
	{
		return (obj->*fp)(args.hd, args.tl.hd, args.tl.tl.hd, args.tl.tl.tl.hd, args.tl.tl.tl.tl.hd);
	};
};

template <typename TClass, typename TRet, typename T1, typename T2, typename T3,
			typename T4, typename T5, typename T6>
struct func_traits <TRet (TClass::*) (T1, T2, T3, T4, T5, T6)>
{
	FUNC_TRAITS_GENERAL_NON_CONST_MEMBER;
	typedef TypeList<T1, T2, T3, T4, T5, T6> params_type;
	static TRet apply(TRet (TClass::*fp)(T1, T2, T3, T4, T5, T6), TClass* obj, TypeValueList<params_type>& args)
	{
		return (obj->*fp)(args.hd, args.tl.hd, args.tl.tl.hd, args.tl.tl.tl.hd, args.tl.tl.tl.tl.hd,
							args.tl.tl.tl.tl.tl.hd);
	};
};

template <typename TClass, typename TRet, typename T1, typename T2, typename T3,
			typename T4, typename T5, typename T6, typename T7>
struct func_traits <TRet (TClass::*) (T1, T2, T3, T4, T5, T6, T7)>
{
	FUNC_TRAITS_GENERAL_NON_CONST_MEMBER;
	typedef TypeList<T1, T2, T3, T4, T5, T6, T7> params_type;
	static TRet apply(TRet (TClass::*fp)(T1, T2, T3, T4, T5, T6, T7), TClass* obj, TypeValueList<params_type>& args)
	{
		return (obj->*fp)(args.hd, args.tl.hd, args.tl.tl.hd, args.tl.tl.tl.hd, args.tl.tl.tl.tl.hd,
							args.tl.tl.tl.tl.tl.hd, args.tl.tl.tl.tl.tl.tl.hd);
	};
};

template <typename TClass, typename TRet, typename T1, typename T2, typename T3,
			typename T4, typename T5, typename T6, typename T7, typename T8>
struct func_traits <TRet (TClass::*) (T1, T2, T3, T4, T5, T6, T7, T8)>
{
	FUNC_TRAITS_GENERAL_NON_CONST_MEMBER;
	typedef TypeList<T1, T2, T3, T4, T5, T6, T7, T8> params_type;
	static TRet apply(TRet (TClass::*fp)(T1, T2, T3, T4, T5, T6, T7, T8), TClass* obj, TypeValueList<params_type>& args)
	{
		return (obj->*fp)(args.hd, args.tl.hd, args.tl.tl.hd, args.tl.tl.tl.hd, args.tl.tl.tl.tl.hd,
							args.tl.tl.tl.tl.tl.hd, args.tl.tl.tl.tl.tl.tl.hd, args.tl.tl.tl.tl.tl.tl.tl.hd);
	};
};

template <typename TClass, typename TRet, typename T1, typename T2, typename T3,
			typename T4, typename T5, typename T6, typename T7, typename T8, typename T9>
struct func_traits <TRet (TClass::*) (T1, T2, T3, T4, T5, T6, T7, T8, T9)>
{
	FUNC_TRAITS_GENERAL_NON_CONST_MEMBER;
	typedef TypeList<T1, T2, T3, T4, T5, T6, T7, T8, T9> params_type;
	static TRet apply(TRet (TClass::*fp)(T1, T2, T3, T4, T5, T6, T7, T8, T9), TClass* obj, TypeValueList<params_type>& args)
	{
		return (obj->*fp)(args.hd, args.tl.hd, args.tl.tl.hd, args.tl.tl.tl.hd, args.tl.tl.tl.tl.hd,
							args.tl.tl.tl.tl.tl.hd, args.tl.tl.tl.tl.tl.tl.hd, args.tl.tl.tl.tl.tl.tl.tl.hd,
							args.tl.tl.tl.tl.tl.tl.tl.tl.hd);
	};
};

////////////////////////////////////////////////////////////////////////////////
// const method traits
////////////////////////////////////////////////////////////////////////////////

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

template <typename TClass, typename TRet, typename T1, typename T2, typename T3,
			typename T4>
struct func_traits <TRet (TClass::*) (T1, T2, T3, T4) const>
{
	FUNC_TRAITS_GENERAL_CONST_MEMBER;
	typedef TypeList<T1, T2, T3, T4> params_type;
	static TRet apply(TRet (TClass::*fp)(T1, T2, T3, T4) const, const TClass* obj, TypeValueList<params_type>& args)
	{
		return (obj->*fp)(args.hd, args.tl.hd, args.tl.tl.hd, args.tl.tl.tl.hd);
	};
};

template <typename TClass, typename TRet, typename T1, typename T2, typename T3,
			typename T4, typename T5>
struct func_traits <TRet (TClass::*) (T1, T2, T3, T4, T5) const>
{
	FUNC_TRAITS_GENERAL_CONST_MEMBER;
	typedef TypeList<T1, T2, T3, T4, T5> params_type;
	static TRet apply(TRet (TClass::*fp)(T1, T2, T3, T4, T5) const, const TClass* obj, TypeValueList<params_type>& args)
	{
		return (obj->*fp)(args.hd, args.tl.hd, args.tl.tl.hd, args.tl.tl.tl.hd, args.tl.tl.tl.tl.hd);
	};
};

template <typename TClass, typename TRet, typename T1, typename T2, typename T3,
			typename T4, typename T5, typename T6>
struct func_traits <TRet (TClass::*) (T1, T2, T3, T4, T5, T6) const>
{
	FUNC_TRAITS_GENERAL_CONST_MEMBER;
	typedef TypeList<T1, T2, T3, T4, T5, T6> params_type;
	static TRet apply(TRet (TClass::*fp)(T1, T2, T3, T4, T5, T6) const, const TClass* obj, TypeValueList<params_type>& args)
	{
		return (obj->*fp)(args.hd, args.tl.hd, args.tl.tl.hd, args.tl.tl.tl.hd, args.tl.tl.tl.tl.hd,
				 	 	 args.tl.tl.tl.tl.tl.hd);
	};
};

template <typename TClass, typename TRet, typename T1, typename T2, typename T3,
			typename T4, typename T5, typename T6,  typename T7>
struct func_traits <TRet (TClass::*) (T1, T2, T3, T4, T5, T6, T7) const>
{
	FUNC_TRAITS_GENERAL_CONST_MEMBER;
	typedef TypeList<T1, T2, T3, T4, T5, T6, T7> params_type;
	static TRet apply(TRet (TClass::*fp)(T1, T2, T3, T4, T5, T6, T7) const, const TClass* obj, TypeValueList<params_type>& args)
	{
		return (obj->*fp)(args.hd, args.tl.hd, args.tl.tl.hd, args.tl.tl.tl.hd, args.tl.tl.tl.tl.hd,
				 	 	 args.tl.tl.tl.tl.tl.hd, args.tl.tl.tl.tl.tl.tl.hd);
	};
};

template <typename TClass, typename TRet, typename T1, typename T2, typename T3,
			typename T4, typename T5, typename T6,  typename T7, typename T8>
struct func_traits <TRet (TClass::*) (T1, T2, T3, T4, T5, T6, T7, T8) const>
{
	FUNC_TRAITS_GENERAL_CONST_MEMBER;
	typedef TypeList<T1, T2, T3, T4, T5, T6, T7, T8> params_type;
	static TRet apply(TRet (TClass::*fp)(T1, T2, T3, T4, T5, T6, T7, T8) const, const TClass* obj, TypeValueList<params_type>& args)
	{
		return (obj->*fp)(args.hd, args.tl.hd, args.tl.tl.hd, args.tl.tl.tl.hd, args.tl.tl.tl.tl.hd,
				 	 	 args.tl.tl.tl.tl.tl.hd, args.tl.tl.tl.tl.tl.tl.hd, args.tl.tl.tl.tl.tl.tl.tl.hd);
	};
};

template <typename TClass, typename TRet, typename T1, typename T2, typename T3,
			typename T4, typename T5, typename T6,  typename T7, typename T8, typename T9>
struct func_traits <TRet (TClass::*) (T1, T2, T3, T4, T5, T6, T7, T8, T9) const>
{
	FUNC_TRAITS_GENERAL_CONST_MEMBER;
	typedef TypeList<T1, T2, T3, T4, T5, T6, T7, T8, T9> params_type;
	static TRet apply(TRet (TClass::*fp)(T1, T2, T3, T4, T5, T6, T7, T8, T9) const, const TClass* obj, TypeValueList<params_type>& args)
	{
		return (obj->*fp)(args.hd, args.tl.hd, args.tl.tl.hd, args.tl.tl.tl.hd, args.tl.tl.tl.tl.hd,
				 	 	 args.tl.tl.tl.tl.tl.hd, args.tl.tl.tl.tl.tl.tl.hd, args.tl.tl.tl.tl.tl.tl.tl.hd,
				 	 	 args.tl.tl.tl.tl.tl.tl.tl.tl.hd);
	};
};

////////////////////////////////////////////////////////////////////////////////
// constructor function traits
////////////////////////////////////////////////////////////////////////////////

template <typename T, typename TTypelist>
struct constructor_traits;

template <typename T>
struct constructor_traits <T, TypeList<> >
{
	typedef TypeList<> params_type;
	static T* apply(TypeValueList<params_type>& args)
	{
		(void)args;
		return new T();
	};
};

template <typename T, typename T1>
struct constructor_traits <T, TypeList<T1> >
{
	typedef TypeList<T1> params_type;
	static T* apply(TypeValueList<params_type>& args)
	{
		return new T(args.hd);
	};
};

template <typename T, typename T1, typename T2>
struct constructor_traits <T, TypeList<T1, T2> >
{
	typedef TypeList<T1, T2> params_type;
	static T* apply(TypeValueList<params_type>& args)
	{
		return new T(args.hd, args.tl.hd);
	};
};

template <typename T, typename T1, typename T2, typename T3>
struct constructor_traits <T, TypeList<T1, T2, T3> >
{
	typedef TypeList<T1, T2, T3> params_type;
	static T* apply(TypeValueList<params_type>& args)
	{
		return new T(args.hd, args.tl.hd, args.tl.tl.hd);
	};
};

template <typename T, typename T1, typename T2, typename T3, typename T4>
struct constructor_traits <T, TypeList<T1, T2, T3, T4> >
{
	typedef TypeList<T1, T2, T3, T4> params_type;
	static T* apply(TypeValueList<params_type>& args)
	{
		return new T(args.hd, args.tl.hd, args.tl.tl.hd, args.tl.tl.tl.hd);
	};
};

template <typename T, typename T1, typename T2, typename T3, typename T4, typename T5>
struct constructor_traits <T, TypeList<T1, T2, T3, T4, T5> >
{
	typedef TypeList<T1, T2, T3, T4, T5> params_type;
	static T* apply(TypeValueList<params_type>& args)
	{
		return new T(args.hd, args.tl.hd, args.tl.tl.hd, args.tl.tl.tl.hd, args.tl.tl.tl.tl.hd);
	};
};

template <typename T, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
struct constructor_traits <T, TypeList<T1, T2, T3, T4, T5, T6> >
{
	typedef TypeList<T1, T2, T3, T4, T5, T6> params_type;
	static T* apply(TypeValueList<params_type>& args)
	{
		return new T(args.hd, args.tl.hd, args.tl.tl.hd, args.tl.tl.tl.hd, args.tl.tl.tl.tl.hd,
				args.tl.tl.tl.tl.tl.hd);
	};
};

template <typename T, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7>
struct constructor_traits <T, TypeList<T1, T2, T3, T4, T5, T6, T7> >
{
	typedef TypeList<T1, T2, T3, T4, T5, T6, T7> params_type;
	static T* apply(TypeValueList<params_type>& args)
	{
		return new T(args.hd, args.tl.hd, args.tl.tl.hd, args.tl.tl.tl.hd, args.tl.tl.tl.tl.hd,
				args.tl.tl.tl.tl.tl.hd, args.tl.tl.tl.tl.tl.tl.hd);
	};
};

template <typename T, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8>
struct constructor_traits <T, TypeList<T1, T2, T3, T4, T5, T6, T7, T8> >
{
	typedef TypeList<T1, T2, T3, T4, T5, T6, T7, T8> params_type;
	static T* apply(TypeValueList<params_type>& args)
	{
		return new T(args.hd, args.tl.hd, args.tl.tl.hd, args.tl.tl.tl.hd, args.tl.tl.tl.tl.hd,
				args.tl.tl.tl.tl.tl.hd, args.tl.tl.tl.tl.tl.tl.hd, args.tl.tl.tl.tl.tl.tl.tl.hd);
	};
};

// end group registry
/// \}

} // end namespace bridge
} // end namespace ug

#endif /* __H__UG_BRIDGE__FUNCTION_TRAITS__ */
