// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 04.02.2012 (d.m.y)
//
// this code was originally created for a private project of mine,
// may however be freely used in ug.

#ifndef __H__UG_variant__
#define __H__UG_variant__

#include <string>

namespace ug{

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
class Variant{
	public:
		enum Type{
			VT_INVALID,
			VT_BOOL,
			VT_INT,
			VT_FLOAT,
			VT_DOUBLE,
			VT_CSTRING,
			VT_STDSTRING,
			VT_POINTER
		};

	public:
		Variant();
		Variant(bool val);
		Variant(int val);
		Variant(float val);
		Variant(double val);
		Variant(const char* val);
		Variant(const std::string& val);
		Variant(void* val);

		Variant(const Variant& v);

		~Variant();

		const Variant& operator=(const Variant& v);

		inline Type type() const	{return m_type;}

		bool to_bool() const;
		int to_int() const;
		float to_float() const;
		double to_double() const;
		const char* to_c_string() const;
		const std::string& to_std_string() const;
		void* to_pointer() const;

	private:
	///	this method only performs assignment, no prior clean up
		void assign_variant(const Variant& v);

	///	returns the name of the current type
		const char* type_name() const;

	private:
		union{
			bool			m_bool;
			int				m_int;
			float			m_float;
			double			m_double;
			const char*		m_cstring;
			std::string*	m_stdstring;
			void*			m_pointer;
		};

		Type	m_type;
};

}//	end of namespace

#endif
