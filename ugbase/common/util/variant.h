/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
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

#ifndef __H__UG_variant__
#define __H__UG_variant__

#include <iostream>
#include <string>

#include "smart_pointer.h"
#include "common/types.h"
#include "common/ug_config.h"

#ifdef UG_FOR_LUA
#include "bindings/lua/lua_function_handle.h"
#include "bindings/lua/lua_table_handle.h"
#endif

namespace ug {

/// \addtogroup ugbase_common_types
/// \{

// todo: Change numbers in VRL
// Type values must not be changed! Bindings rely on the exact values.
// Append new types at the end and update bindings.
// If in doubt contact binding developers!
//	BUT NOW HAVE BEEN CHANGED (A.Vogel)
//	HERE ARE THE OLD VALUES, BELOW THE NEW
/*enum ParameterTypes
{
	PT_UNKNOWN = 0,
	PT_BOOL = 1,
	PT_INTEGER = 2,
	PT_NUMBER = 3,
	PT_CSTRING = 4,
	PT_STD_STRING = 5,
	PT_POINTER = 6,
	PT_CONST_POINTER = 7,
	PT_SMART_POINTER = 8,
	PT_CONST_SMART_POINTER = 9
};
*/


///	A variant can represent variables of different types.
/**	A variant will represent the value, which was assigned during construction
 * or in an assignment operation.
 * Use the different to_.. methods, to retrieve the encapsulated value. Take
 * care only to retrieve values, which are compatible with the encapsulated value.
 *
 * Number and boolean types are automatically converted to each other:
 * VT_BOOL, VT_INT, VT_FLOAT, VT_DOUBLE.
 *
 * A stdstring is automatically converted to a cstr, but not vice-versa.
 *
 * You may inspect the type represented by a variant by calling the type() method,
 * which will return a constant enumerated in Variant::Type.
 */
class UG_API Variant{
	public: // TODO?: hide enum
		enum Type{
			VT_INVALID = 0,
			VT_BOOL = 1,
			VT_INT = 2,
			VT_SIZE_T = 3,
			VT_FLOAT = 4,
			VT_DOUBLE = 5,
			VT_CSTRING = 6,
			VT_STDSTRING = 7,
			VT_POINTER = 8,
			VT_CONST_POINTER = 9,
			VT_SMART_POINTER = 10,
			VT_CONST_SMART_POINTER = 11
#ifdef UG_FOR_LUA
			,VT_LUA_FUNCTION_HANDLE = 12
			,VT_LUA_TABLE_HANDLE = 13
#endif
		};

	public:
		Variant();

		Variant(bool val);

		Variant(int val);

		Variant(size_t val);

		Variant(float val);

		Variant(double val);

		Variant(const char* val);

		Variant(const std::string& val);

		Variant(void* val);

		Variant(const void* val);

		Variant(const SmartPtr<void>& val);

		Variant(const ConstSmartPtr<void>& val);
#ifdef UG_FOR_LUA
		Variant(LuaFunctionHandle val);
		Variant(LuaTableHandle val);
#endif

		Variant(const Variant& v);

		~Variant();

		const Variant& operator =(const Variant& v);

		[[nodiscard]] inline Type type() const {return m_type;}

		template <typename T>
		inline static Type type() {return VT_INVALID;}

		[[nodiscard]] bool is_valid() const{return m_type;}
		[[nodiscard]] bool to_bool() const;
		[[nodiscard]] int to_int() const;
		[[nodiscard]] size_t to_size_t() const;
		[[nodiscard]] float to_float() const;
		[[nodiscard]] number to_number() const;
		[[nodiscard]] double to_double() const;
		[[nodiscard]] const char* to_c_string() const;
		[[nodiscard]] const std::string& to_std_string() const;
		[[nodiscard]] void* to_pointer() const;
		[[nodiscard]] const void* to_const_pointer() const;
		[[nodiscard]] SmartPtr<void> to_smart_pointer() const;
		[[nodiscard]] ConstSmartPtr<void> to_const_smart_pointer() const;
#ifdef UG_FOR_LUA
		[[nodiscard]] LuaFunctionHandle to_lua_function_handle() const;
		[[nodiscard]] LuaTableHandle to_lua_table_handle() const;
#endif

		template <typename T>
		inline T to() const;

	private:
	///	this method only performs assignment, no prior clean up
		void assign_variant(const Variant& v);

	///	returns the name of the current type
		[[nodiscard]] const char* type_name() const;

	private:
		union{
			bool					m_bool;
			int						m_int;
			size_t					m_size_t;
			float					m_float;
			double					m_double;
			const char*				m_cstring;
			std::string*			m_stdstring;
			void*					m_pointer;
			const void* 			m_constptr;
			SmartPtr<void>*			m_smartptr;
			ConstSmartPtr<void>* 	m_constsmartptr;
#ifdef UG_FOR_LUA
			LuaFunctionHandle		m_luafcthandle;
			LuaTableHandle		m_luatblhandle;
#endif
		};

		Type	m_type;
};

template <> inline bool 				Variant::to<bool>() const 					{return to_bool();}
template <> inline int 					Variant::to<int>() const 					{return to_int();}
template <> inline size_t 				Variant::to<size_t>() const 				{return to_size_t();}
template <> inline float 				Variant::to<float>() const 					{return to_float();}
template <> inline double 				Variant::to<double>() const 				{return to_double();}
template <> inline const char* 			Variant::to<const char*>() const 			{return to_c_string();}
template <> inline std::string 			Variant::to<std::string>() const 			{return to_std_string();}
template <> inline const std::string& 	Variant::to<const std::string&>() const 	{return to_std_string();}
template <> inline void* 				Variant::to<void*>() const 					{return to_pointer();}
template <> inline const void* 		   	Variant::to<const void*>() const 			{return to_const_pointer();}
template <> inline SmartPtr<void> 	   	Variant::to<SmartPtr<void> >() const 		{return to_smart_pointer();}
template <> inline ConstSmartPtr<void> 	Variant::to<ConstSmartPtr<void> >() const 	{return to_const_smart_pointer();}
#ifdef UG_FOR_LUA
template <> inline LuaFunctionHandle  Variant::to<LuaFunctionHandle>() const 		{return to_lua_function_handle();}
template <> inline LuaTableHandle 	Variant::to<LuaTableHandle>() const 		{return to_lua_table_handle();}
#endif

template <> inline Variant::Type Variant::type<bool>() 					{return VT_BOOL;}
template <> inline Variant::Type Variant::type<int>() 					{return VT_INT;}
template <> inline Variant::Type Variant::type<size_t>() 				{return VT_SIZE_T;}
template <> inline Variant::Type Variant::type<float>() 				{return VT_FLOAT;}
template <> inline Variant::Type Variant::type<double>() 				{return VT_DOUBLE;}
template <> inline Variant::Type Variant::type<const char*>() 			{return VT_CSTRING;}
template <> inline Variant::Type Variant::type<std::string>() 			{return VT_STDSTRING;}
template <> inline Variant::Type Variant::type<const std::string&>() 	{return VT_STDSTRING;}
template <> inline Variant::Type Variant::type<void*>() 				{return VT_POINTER;}
template <> inline Variant::Type Variant::type<const void*>() 			{return VT_CONST_POINTER;}
template <> inline Variant::Type Variant::type<SmartPtr<void> >() 		{return VT_SMART_POINTER;}
template <> inline Variant::Type Variant::type<ConstSmartPtr<void> >() 	{return VT_CONST_SMART_POINTER;}
#ifdef UG_FOR_LUA
template <> inline Variant::Type Variant::type<LuaFunctionHandle>() 	{return VT_LUA_FUNCTION_HANDLE;}
template <> inline Variant::Type Variant::type<LuaTableHandle>() 	{return VT_LUA_TABLE_HANDLE;}
#endif


// end group ugbase_common_types
/// \}

}//	end of namespace

UG_API std::ostream& operator << (std::ostream& outStream, const ug::Variant& v);

#endif
