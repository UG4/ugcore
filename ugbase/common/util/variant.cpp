// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 04.02.2012 (d.m.y)
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


Variant::Variant(const Variant& v)
{
	assign_variant(v);
}

Variant::~Variant()
{
	if(m_type == VT_STDSTRING)
		delete m_stdstring;
}

const Variant& Variant::operator=(const Variant& v)
{
//	if the variant encapsulates an std::string, we first have to delete the
//	old instance
	if(m_type == VT_STDSTRING)
		delete m_stdstring;

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
		case VT_FLOAT:	return (int)m_float;
		case VT_DOUBLE:	return (int)m_double;
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
		case VT_FLOAT:	return m_float;
		case VT_DOUBLE:	return (float)m_double;
		default: break;
	}

	UG_THROW("Variant: can't convert " << type_name() << " to float.");
//	this should never be reached
	return 0;
}

double Variant::to_double() const
{
	switch(m_type){
		case VT_BOOL:	return (double)m_bool;
		case VT_INT:	return (double)m_int;
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
//	this should never be reached
	return 0;
}

const char* Variant::type_name() const
{
	switch(m_type){
		case VT_INVALID:	return "invalid_type";
		case VT_BOOL:		return "bool";
		case VT_INT:		return "int";
		case VT_FLOAT:		return "float";
		case VT_DOUBLE:		return "double";
		case VT_CSTRING:	return "cstring";
		case VT_STDSTRING:	return "stdstring";
		case VT_POINTER:	return "pointer";

		default: return "unknown_type";
	}
}

}//	end of namespace
