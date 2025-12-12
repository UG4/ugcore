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

#include "variant.h"

#include <iostream>

#include "common/error.h"


#define untested() ( std::cerr <<  "@@#\n@@@:"<< __FILE__ << ":"<< __LINE__ \
          <<":" << __func__ << "\n" )

using namespace std;

namespace ug {

Variant::Variant() :
	m_type(VT_INVALID)
{}

Variant::Variant(bool val) :
	m_bool(val),
	m_type(VT_BOOL)
{}

Variant::Variant(int val) :
	m_int(val),
	m_type(VT_INT)
{}

Variant::Variant(size_t val) :
	m_size_t(val),
	m_type(VT_SIZE_T)
{}

Variant::Variant(float val) :
	m_float(val),
	m_type(VT_FLOAT)
{}

Variant::Variant(double val) :
	m_double(val),
	m_type(VT_DOUBLE)
{}

Variant::Variant(const char* val) :
	m_cstring(val),
	m_type(VT_CSTRING)
{}

Variant::Variant(const std::string& val) :
	m_stdstring(new string(val)),
	m_type(VT_STDSTRING)
{}

Variant::Variant(void* val) :
	m_pointer(val),
	m_type(VT_POINTER)
{}

Variant::Variant(const void* val) :
	m_constptr(val),
	m_type(VT_CONST_POINTER)
{}

Variant::Variant(const SmartPtr<void>& val) :
	m_smartptr(new SmartPtr<void>(val)),
	m_type(VT_SMART_POINTER)
{}

Variant::Variant(const ConstSmartPtr<void>& val) :
	m_constsmartptr(new ConstSmartPtr<void>(val)),
	m_type(VT_CONST_SMART_POINTER)
{}

#ifdef UG_FOR_LUA
Variant::Variant(LuaFunctionHandle val) :
	m_luafcthandle(val),
	m_type(VT_LUA_FUNCTION_HANDLE)
{}
Variant::Variant(LuaTableHandle val) :
	m_luatblhandle(val),
	m_type(VT_LUA_TABLE_HANDLE)
{}
#endif

Variant::Variant(const Variant& v) : m_type(VT_INVALID)
{
	assign_variant(v);
}

Variant::~Variant()
{
	if(m_type == VT_STDSTRING)
		delete m_stdstring;

	if(m_type == VT_SMART_POINTER)
		delete m_smartptr;
	if(m_type == VT_CONST_SMART_POINTER)
		delete m_constsmartptr;
}

const Variant& Variant::operator = (const Variant& v)
{
//	if the variant encapsulates an std::string or a smartptr,
//	we first have to delete the	old instance
	if(m_type == VT_STDSTRING)
		delete m_stdstring;
	if(m_type == VT_SMART_POINTER)
		delete m_smartptr;
	if(m_type == VT_CONST_SMART_POINTER)
		delete m_constsmartptr;

//	now assign the new value
	assign_variant(v);

//	done. return reference to this
	return *this;
}

void Variant::assign_variant(const Variant& v)
{
#ifdef UG_FOR_LUA
	// need std::variant behaviour in order to be able to
	// work with objects.
	if(m_type == VT_LUA_TABLE_HANDLE && v.m_type != m_type){ untested();
		m_luatblhandle.~LuaTableHandle();
	}else if(m_type == VT_LUA_TABLE_HANDLE && v.m_type == m_type){ untested();
		// assign.
		m_luatblhandle = v.m_luatblhandle;
	}else if(v.m_type == VT_LUA_TABLE_HANDLE && v.m_type != m_type){
		// create in place, reuse memory.
		new (&m_luatblhandle) LuaTableHandle(v.m_luatblhandle);
	}else{
	}
#endif
	switch(v.m_type){
		case VT_BOOL:
			m_bool = v.m_bool;
			break;
		case VT_INT:
			m_int = v.m_int;
			break;
		case VT_SIZE_T:
			m_size_t = v.m_size_t;
			break;
		case VT_FLOAT:
			m_float = v.m_float;
			break;
		case VT_DOUBLE:
			m_double = v.m_double;
			break;
		case VT_CSTRING:
			m_cstring = v.m_cstring;
			break;
		case VT_STDSTRING:
			m_stdstring = new string;
			*m_stdstring = *v.m_stdstring;
			break;
		case VT_POINTER:
			m_pointer = v.m_pointer;
			break;
		case VT_CONST_POINTER:
			m_constptr = v.m_constptr;
			break;
		case VT_SMART_POINTER:
			m_smartptr = new SmartPtr<void>(*v.m_smartptr);
			break;
		case VT_CONST_SMART_POINTER:
			m_constsmartptr = new ConstSmartPtr<void>(*v.m_constsmartptr);
			break;
#ifdef UG_FOR_LUA
		case VT_LUA_FUNCTION_HANDLE:
			m_luafcthandle = v.m_luafcthandle;
			break;
		case VT_LUA_TABLE_HANDLE:
			break;
#endif
		default:
			break;
	}

	m_type = v.m_type;
}

bool Variant::to_bool() const
{
	switch(m_type){
		case VT_BOOL:	return m_bool;
		case VT_INT:	return m_int != 0;
		case VT_SIZE_T:	return m_size_t != 0;
		case VT_FLOAT:	return m_float != 0;
		case VT_DOUBLE:	return m_double != 0;
		default: break;
	}

	UG_THROW("Variant: can't convert " << type_name() << " to bool.");
//	this should never be reached
	return false;
}

int Variant::to_int() const
{
	switch(m_type){
		case VT_BOOL:	return (int)m_bool;
		case VT_INT:	return m_int;
		case VT_SIZE_T:	return static_cast<int>(m_size_t);
		case VT_FLOAT:	return static_cast<int>(m_float);
		case VT_DOUBLE:	return static_cast<int>(m_double);
		default: break;
	}

	UG_THROW("Variant: can't convert " << type_name() << " to int.");
//	this should never be reached
	return 0;
}

size_t Variant::to_size_t() const
{
	// note: we only allow save casts, i.e. not int->size_t
	switch(m_type){
		case VT_BOOL:	return (size_t)m_bool;
		case VT_SIZE_T:	return m_size_t;
		default: break;
	}

	UG_THROW("Variant: can't convert " << type_name() << " to int.");
//	this should never be reached
	return 0;
}

float Variant::to_float() const
{
	switch(m_type){
		case VT_BOOL:	return (float)m_bool;
		case VT_INT:	return static_cast<float>(m_int);
		case VT_SIZE_T:	return static_cast<float>(m_size_t);
		case VT_FLOAT:	return m_float;
		case VT_DOUBLE:	return static_cast<float>(m_double);
		default: break;
	}

	UG_THROW("Variant: can't convert " << type_name() << " to float.");
//	this should never be reached
	return 0;
}

number Variant::to_number() const
{
	switch(m_type){
		case VT_BOOL:	return m_bool;
		case VT_INT:	return m_int;
		case VT_SIZE_T:	return m_size_t;
		case VT_FLOAT:	return m_float;
		case VT_DOUBLE:	return m_double;
		default: break;
	}

	UG_THROW("Variant: can't convert " << type_name() << " to number.");
//	this should never be reached
	return 0;
}

double Variant::to_double() const
{
	switch(m_type){
		case VT_BOOL:	return m_bool;
		case VT_INT:	return m_int;
		case VT_SIZE_T:	return m_size_t;
		case VT_FLOAT:	return m_float;
		case VT_DOUBLE:	return m_double;
		default: break;
	}

	UG_THROW("Variant: can't convert " << type_name() << " to double.");
//	this should never be reached
	return 0;
}

const char* Variant::to_c_string() const
{
	if(m_type == VT_CSTRING)
		return m_cstring;

	if(m_type == VT_STDSTRING)
		return m_stdstring->c_str();

	UG_THROW("Variant: can't convert " << type_name() << " to cstring.");
//	this should never be reached
	return 0;
}

const std::string& Variant::to_std_string() const
{
	if(m_type == VT_STDSTRING)
		return *m_stdstring;

	UG_THROW("Variant: can't convert " << type_name() << " to stdstring.");
//	this should never be reached
	return *m_stdstring;
}

void* Variant::to_pointer() const
{
	if(m_type == VT_POINTER)
		return m_pointer;

	UG_THROW("Variant: can't convert " << type_name() << " to pointer.");
}

const void* Variant::to_const_pointer() const
{
	if(m_type == VT_CONST_POINTER)
		return m_constptr;

	if(m_type == VT_POINTER)
		return const_cast<const void*>(m_pointer);

	UG_THROW("Variant: can't convert " << type_name() << " to const pointer.");
}

SmartPtr<void> Variant::to_smart_pointer() const
{
	if(m_type == VT_SMART_POINTER)
		return *m_smartptr;

	UG_THROW("Variant: can't convert " << type_name() << " to smart pointer.");
}

ConstSmartPtr<void> Variant::to_const_smart_pointer() const
{
	if(m_type == VT_CONST_SMART_POINTER)
		return *m_constsmartptr;

	if(m_type == VT_SMART_POINTER)
		return *m_smartptr;

	UG_THROW("Variant: can't convert " << type_name() << " to const smart pointer.");
}

#ifdef UG_FOR_LUA
LuaFunctionHandle Variant::to_lua_function_handle() const
{
	return m_luafcthandle;
}
LuaTableHandle Variant::to_lua_table_handle() const
{
	return m_luatblhandle;
}
#endif

const char* Variant::type_name() const
{
	switch(m_type){
		case VT_INVALID:	return "invalid_type";
		case VT_BOOL:		return "bool";
		case VT_SIZE_T:		return "size_t";
		case VT_INT:		return "int";
		case VT_FLOAT:		return "float";
		case VT_DOUBLE:		return "double";
		case VT_CSTRING:	return "cstring";
		case VT_STDSTRING:	return "stdstring";
		case VT_POINTER:	return "pointer";
		case VT_CONST_POINTER:	return "constpointer";
		case VT_SMART_POINTER:	return "smartpointer";
		case VT_CONST_SMART_POINTER:	return "constsmartpointer";

		default: return "unknown_type";
	}
}

}//	end of namespace

std::ostream& operator << (std::ostream& outStream, const ug::Variant& v)
{
	using namespace ug;
	// TODO?: hide enum
	switch(v.type()){
		case Variant::VT_BOOL:
			if(v.to_bool())
				outStream << "true";
			else
				outStream << "false";
			break;
		case Variant::VT_INT:
			outStream << v.to_int();
			break;
		case Variant::VT_SIZE_T:
			outStream << v.to_size_t();
			break;
		case Variant::VT_FLOAT:
			outStream << v.to_float();
			break;
		case Variant::VT_DOUBLE:
			outStream << v.to_double();
			break;
		case Variant::VT_CSTRING:
			outStream << v.to_c_string();
			break;
		case Variant::VT_STDSTRING:
			outStream << v.to_std_string();
			break;
		default:
			outStream << "?";
	}
	return outStream;
}
