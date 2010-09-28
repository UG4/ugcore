
#ifndef __H__UG_INTERFACE__UGBRIDGE__PARAM_TO_TYPE_VALUE_LIST__
#define __H__UG_INTERFACE__UGBRIDGE__PARAM_TO_TYPE_VALUE_LIST__

#include "parameter_stack.h"

namespace ug {

namespace interface{

//////////////////////////////
// ParameterStack DataStack
//////////////////////////////

template <typename TData>
struct PLStack
{
	private:
		static void push(ParameterStack& ps);
		static void write(ParameterStack& ps, TData data, int index);
		static TData read(const ParameterStack& ps, int index);
};

template <>
struct PLStack<void>
{
	static void add(ParameterStack& pl, const char* name)
	{
		// do nothing
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
void WriteTypeValueToParameterStackTop(TType& val, ParameterStack& out)
{
	PLStack<TType>::write(out, 0);
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
	static void create(ParameterStack& in, const char* paramValNames, const char* delimiter)
	{
		// do nothing and end recursion
	}
};


} // end namespace interface

} // end namespace ug

#endif /* __H__UG_INTERFACE__UGBRIDGE__PARAM_TO_TYPE_VALUE_LIST__ */
