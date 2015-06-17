// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 04.02.2012 (d.m.y())
//
// this code was originally created for a private project of mine,
// may however be freely used in ug.

#include "common/error.h"
#include "variant.h"

using namespace std;

namespace ug{

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
#endif

Variant::Variant(const Variant& v)
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

const Variant& Variant::operator=(const Variant& v)
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
		case VT_SIZE_T:	return (int)m_size_t;
		case VT_FLOAT:	return (int)m_float;
		case VT_DOUBLE:	return (int)m_double;
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
		case VT_INT:	return (float)m_int;
		case VT_SIZE_T:	return (float)m_size_t;
		case VT_FLOAT:	return m_float;
		case VT_DOUBLE:	return (float)m_double;
		default: break;
	}

	UG_THROW("Variant: can't convert " << type_name() << " to float.");
//	this should never be reached
	return 0;
}

number Variant::to_number() const
{
	switch(m_type){
		case VT_BOOL:	return (number)m_bool;
		case VT_INT:	return (number)m_int;
		case VT_SIZE_T:	return (number)m_size_t;
		case VT_FLOAT:	return (number)m_float;
		case VT_DOUBLE:	return (number)m_double;
		default: break;
	}

	UG_THROW("Variant: can't convert " << type_name() << " to number.");
//	this should never be reached
	return 0;
}

double Variant::to_double() const
{
	switch(m_type){
		case VT_BOOL:	return (double)m_bool;
		case VT_INT:	return (double)m_int;
		case VT_SIZE_T:	return (double)m_size_t;
		case VT_FLOAT:	return (double)m_float;
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

std::ostream& operator<< (std::ostream& outStream, const ug::Variant& v)
{
	using namespace ug;
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
