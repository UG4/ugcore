// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 04.02.2012 (d.m.y())
//
// this code was originally created for a private project of mine,
// may however be freely used in ug.

#ifndef __H__UG_variant__
#define __H__UG_variant__

#include <string>
#include "smart_pointer.h"
#include "common/types.h"
#include "common/ug_config.h"

namespace ug{

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
	public:
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

		Variant(const Variant& v);

		~Variant();

		const Variant& operator=(const Variant& v);

		inline Type type() const	{return m_type;}

		template <typename T>
		inline static Type type() {return VT_INVALID;}

		bool to_bool() const;
		int to_int() const;
		size_t to_size_t() const;
		float to_float() const;
		number to_number() const;
		double to_double() const;
		const char* to_c_string() const;
		const std::string& to_std_string() const;
		void* to_pointer() const;
		const void* to_const_pointer() const;
		SmartPtr<void> to_smart_pointer() const;
		ConstSmartPtr<void> to_const_smart_pointer() const;

		template <typename T>
		inline T to() const;

	private:
	///	this method only performs assignment, no prior clean up
		void assign_variant(const Variant& v);

	///	returns the name of the current type
		const char* type_name() const;

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


// end group ugbase_common_types
/// \}

}//	end of namespace

UG_API std::ostream& operator<< (std::ostream& outStream, const ug::Variant& v);

#endif
