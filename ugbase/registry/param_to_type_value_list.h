//	authors: Sebastian Reiter, Andreas Vogel

#ifndef __H__UG_BRIDGE__PARAM_TO_TYPE_VALUE_LIST__
#define __H__UG_BRIDGE__PARAM_TO_TYPE_VALUE_LIST__

#include "parameter_stack.h"
#include "class_name_provider.h"
#include <string>
#include <vector>
#include <cstring>

namespace ug
{
namespace bridge
{

//////////////////////////////
//////////////////////////////
// ParameterStack DataStack
//////////////////////////////
//////////////////////////////

///	Generic PLStack structure. Raises a compile error if used directly by ParameterStackToTypeValueList
/**	The methods in this class should give an idea, on how the methods in
 * a concrete PLStack class should look like. The methods have to be public
 * there. For the generic PLStack, all members are private to enforce a
 * compile error, if no specialization for a given TData type exists.
 */
template <typename TData>
struct PLStack
{
	private:
		static void write(ParameterStack& ps, TData data, int index);
		static void push(ParameterStack& ps);
		static TData read(const ParameterStack& ps, int index);
};

//////////////////////////////
// build-in data types
//////////////////////////////
template <>
struct PLStack<bool>
{
	static void push(ParameterStack& ps)
	{
		ps.push_bool();
	}
	static void write(ParameterStack& ps, bool data, int index)
	{
		ps.set_bool(index, data);
	}
	static bool read(const ParameterStack& ps, int index)
	{
		return ps.to_bool(index);
	}
};

template <>
struct PLStack<int>
{
	static void push(ParameterStack& ps)
	{
		ps.push_integer();
	}
	static void write(ParameterStack& ps, int data, int index)
	{
		ps.set_integer(index, data);
	}
	static int read(const ParameterStack& ps, int index)
	{
		return ps.to_integer(index);
	}
};

template <>
struct PLStack<size_t>
{
	static void push(ParameterStack& ps)
	{
		ps.push_integer();
	}
	static void write(ParameterStack& ps, size_t data, int index)
	{
		ps.set_integer(index, (int)data);
	}
	static size_t read(const ParameterStack& ps, int index)
	{
		return (size_t)ps.to_integer(index);
	}
};

/* on SOME systems, size_t is not equal to unsigned int...
template <>
struct PLStack<unsigned int>
{
	static void push(ParameterStack& ps)
	{
		ps.push_integer();
	}
	static void write(ParameterStack& ps, unsigned int data, int index)
	{
		ps.set_integer(index, (int)data);
	}
	static unsigned int read(const ParameterStack& ps, int index)
	{
		return (unsigned int)ps.to_integer(index);
	}
};
*/
template <>
struct PLStack<float>
{
	static void push(ParameterStack& ps)
	{
		ps.push_number();
	}
	static void write(ParameterStack& ps, float data, int index)
	{
		ps.set_number(index, data);
	}
	static float read(const ParameterStack& ps, int index)
	{
		return ps.to_number(index);
	}
};

template <>
struct PLStack<double>
{
	static void push(ParameterStack& ps)
	{
		ps.push_number();
	}
	static void write(ParameterStack& ps, double data, int index)
	{
		ps.set_number(index, data);
	}
	static double read(const ParameterStack& ps, int index)
	{
		return ps.to_number(index);
	}
};

template <>
struct PLStack<const char*>
{
	static void push(ParameterStack& ps)
	{
		ps.push_cstring();
	}
	static void write(ParameterStack& ps, const char* data, int index)
	{
		ps.set_cstring(index, data);
	}
	static const char* read(const ParameterStack& ps, int index)
	{
		return ps.to_cstring(index);
	}
};

template <>
struct PLStack<std::string>
{
	static void push(ParameterStack& ps)
	{
		ps.push_std_string();
	}
	static void write(ParameterStack& ps, const std::string& data, int index)
	{
		ps.set_std_string(index, data.c_str());
	}
	static std::string read(const ParameterStack& ps, int index)
	{
		return ps.to_std_string(index);
	}
};

template <>
struct PLStack<const std::string&>
{
	static void push(ParameterStack& ps)
	{
		ps.push_std_string();
	}
	static void write(ParameterStack& ps, const std::string& data, int index)
	{
		ps.set_std_string(index, data.c_str());
	}
	static const std::string& read(const ParameterStack& ps, int index)
	{
		return ps.to_std_string(index);
	}
};


//////////////////////////////
// classes
//////////////////////////////

template <typename T>
struct PLStack<SmartPtr<std::vector<T> > >
{
	static void push(ParameterStack& ps)
	{
		ps.push_smart_pointer_std_vector<T>();
	}
	static void write(ParameterStack& ps, const SmartPtr<std::vector<T> >& data, int index)
	{
		ps.set_smart_pointer(index, data);
	}
	static SmartPtr<std::vector<T> > read(const ParameterStack& ps, int index)
	{
		return ps.to_smart_pointer<std::vector<T> >(index);
	}
};





template <typename TClass>
struct PLStack<TClass*>
{
	static void push(ParameterStack& ps)
	{
		ps.push_pointer<TClass>();
	}
	static void write(ParameterStack& ps, TClass* data, int index)
	{
		ps.set_pointer(index, data);
	}
	static TClass* read(const ParameterStack& ps, int index)
	{
		return ps.to_pointer<TClass>(index);
	}
};

template <typename TClass>
struct PLStack<const TClass*>
{
	static void push(ParameterStack& ps)
	{
		ps.push_const_pointer<TClass>();
	}
	static void write(ParameterStack& ps, const TClass* data, int index)
	{
		ps.set_const_pointer(index, data);
	}
	static const TClass* read(const ParameterStack& ps, int index)
	{
		return ps.to_const_pointer<TClass>(index);
	}
};

template <typename TClass>
struct PLStack<TClass&>
{
	static void push(ParameterStack& ps)
	{
		ps.push_pointer<TClass>();
	}
	static void write(ParameterStack& ps, TClass& data, int index)
	{
		ps.set_pointer(index, &data);
	}
	static TClass& read(const ParameterStack& ps, int index)
	{
		return *ps.to_pointer<TClass>(index);
	}
};

template <typename TClass>
struct PLStack<const TClass&>
{
	static void push(ParameterStack& ps)
	{
		ps.push_const_pointer<TClass>();
	}
	static void write(ParameterStack& ps, const TClass& data, int index)
	{
		ps.set_const_pointer(index, &data);
	}
	static const TClass& read(const ParameterStack& ps, int index)
	{
		return *ps.to_const_pointer<TClass>(index);
	}
};

template <typename TClass>
struct PLStack<SmartPtr<TClass> >
{
	static void push(ParameterStack& ps)
	{
		ps.push_smart_pointer<TClass>();
	}
	static void write(ParameterStack& ps, const SmartPtr<TClass>& data, int index)
	{
		ps.set_smart_pointer(index, data);
	}
	static SmartPtr<TClass> read(const ParameterStack& ps, int index)
	{
		return ps.to_smart_pointer<TClass>(index);
	}
};

template <typename TClass>
struct PLStack<ConstSmartPtr<TClass> >
{
	static void push(ParameterStack& ps)
	{
		ps.push_const_smart_pointer<TClass>();
	}
	static void write(ParameterStack& ps, const ConstSmartPtr<TClass>& data, int index)
	{
		ps.set_const_smart_pointer(index, data);
	}
	static ConstSmartPtr<TClass> read(const ParameterStack& ps, int index)
	{
		return ps.to_const_smart_pointer<TClass>(index);
	}
};

//////////////////////////////
// ParameterStackToTypeValueList
//////////////////////////////

// ParameterStackToTypeValueList
template <typename TTypeList, int index = 0>
struct ParameterStackToTypeValueList;

// specialization for empty TypeValueList
template <int index>
struct ParameterStackToTypeValueList<TypeList<>, index> :
	public TypeValueList<TypeList<> >
{
	ParameterStackToTypeValueList(const ParameterStack& in) {}
};

// implementation for non-empty TypeValueList
template <typename TTypeList, int index>
struct ParameterStackToTypeValueList :
	public TypeValueList<TTypeList>
{
	typedef	 typename TTypeList::head	head;
	typedef	 typename TTypeList::tail	tail;

	ParameterStackToTypeValueList(const ParameterStack& in) :
		TypeValueList<TTypeList>(	PLStack<head>::read(in, index),
									ParameterStackToTypeValueList<tail, index+1>(in))
			{}
};


//////////////////////////////
// WriteTypeValueToParameterStackTop
//////////////////////////////

// WriteTypeValueToParameterStackTop
/*
template <typename TType>
void PushTypeValueToParameterStack(TType& val, ParameterStack& out)
{
	PLStack<TType>::push(out);
	PLStack<TType>::write(out, val, -1);
};
*/

//////////////////////////////
// CreateParameterStack
//////////////////////////////

// implementation for non-empty
template <typename TTypeList>
struct CreateParameterStack
{
	typedef	typename TTypeList::head	head;
	typedef	typename TTypeList::tail	tail;

	static void create(ParameterStack& in)
	{
	    // add parameter so list
		PLStack<head>::push(in);

		// recursive call for remaining types and names
		CreateParameterStack<tail>::create(in);
	}
};

// specialization for empty typelist
template <>
struct CreateParameterStack<TypeList<> >
{
	static void create(ParameterStack& in)
	{
		// do nothing and end recursion
	}
};


} // end namespace bridge

} // end namespace ug

#endif /* __H__UG_BRIDGE__PARAM_TO_TYPE_VALUE_LIST__ */
