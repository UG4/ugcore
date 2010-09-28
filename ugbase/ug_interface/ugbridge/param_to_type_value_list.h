
#ifndef __H__UG_INTERFACE__UGBRIDGE__PARAM_TO_TYPE_VALUE_LIST__
#define __H__UG_INTERFACE__UGBRIDGE__PARAM_TO_TYPE_VALUE_LIST__

#include "parameter_stack.h"
#include <string>

namespace ug {

namespace interface{


template <typename TClass>
struct ClassNameProvider
{
	static void set_name(const char* name) {m_name = name;}
	static const char* name() {return m_name;}

	private:
		static const char* m_name;
};

template <typename TClass>
const char* ClassNameProvider<TClass>::m_name = "";


//////////////////////////////
//////////////////////////////
// ParameterStack DataStack
//////////////////////////////
//////////////////////////////

template <typename TData>
struct PLStack
{
	private:
		static void push(ParameterStack& ps);
		static void write(ParameterStack& ps, TData data, int index);
		static TData read(const ParameterStack& ps, int index);
};

//////////////////////////////
// build-in data types
//////////////////////////////

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
struct PLStack<const char*>
{
	static void push(ParameterStack& ps)
	{
		ps.push_string();
	}
	static void write(ParameterStack& ps, const char* data, int index)
	{
		ps.set_string(index, data);
	}
	static const char* read(const ParameterStack& ps, int index)
	{
		return ps.to_string(index);
	}
};

template <>
struct PLStack<std::string>
{
	static void push(ParameterStack& ps)
	{
		ps.push_string();
	}
	static void write(ParameterStack& ps, const std::string& data, int index)
	{
		ps.set_string(index, data.c_str());
	}
	static std::string read(const ParameterStack& ps, int index)
	{
		return std::string(ps.to_string(index));
	}
};

template <>
struct PLStack<const std::string&>
{
	static void push(ParameterStack& ps)
	{
		ps.push_string();
	}
	static void write(ParameterStack& ps, const std::string& data, int index)
	{
		ps.set_string(index, data.c_str());
	}
	static std::string read(const ParameterStack& ps, int index)
	{
		return std::string(ps.to_string(index));
	}
};


//////////////////////////////
// classes
//////////////////////////////

template <typename TClass>
struct PLStack<TClass*>
{
	static void push(ParameterStack& ps)
	{
		ps.push_pointer(ClassNameProvider<TClass>::name());
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
struct PLStack<TClass&>
{
	static void push(ParameterStack& ps)
	{
		ps.push_reference(ClassNameProvider<TClass>::name());
	}
	static void write(ParameterStack& ps, TClass& data, int index)
	{
		ps.set_reference(index, data);
	}
	static TClass& read(const ParameterStack& ps, int index)
	{
		return ps.to_reference<TClass>(index);
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
	typedef	typename TTypeList::head	head;
	typedef	typename TTypeList::tail	tail;

	ParameterStackToTypeValueList(const ParameterStack& in) :
		TypeValueList<TTypeList>(	PLStack<head>::read(in, index),
									ParameterStackToTypeValueList<tail, index+1>(in))
			{}
};


//////////////////////////////
// WriteTypeValueToParameterStackTop
//////////////////////////////

// WriteTypeValueToParameterStackTop
template <typename TType>
void PushTypeValueToParameterStack(TType& val, ParameterStack& out)
{
	PLStack<TType>::push(out);
	PLStack<TType>::write(out, val, -1);
};


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


} // end namespace interface

} // end namespace ug

#endif /* __H__UG_INTERFACE__UGBRIDGE__PARAM_TO_TYPE_VALUE_LIST__ */
