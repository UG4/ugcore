
#ifndef __H__UG_INTERFACE__UGBRIDGE__PARAM_TO_TYPE_VALUE_LIST__
#define __H__UG_INTERFACE__UGBRIDGE__PARAM_TO_TYPE_VALUE_LIST__

#include "../interface_base.h"


namespace ug {

namespace interface{

//////////////////////////////
// ParameterList DataStack
//////////////////////////////

template <typename TData>
struct PLStack
{
	private:
		static void add(ParameterList& pl, const char* name);
		static void write(const ParameterList& pl, TData data, int index);
		static TData read(const ParameterList& pl, int index);
};

template <>
struct PLStack<void>
{
	static void add(ParameterList& pl, const char* name)
	{
		// do nothing
	}
};


template <>
struct PLStack<double>
{
	static void add(ParameterList& pl, const char* name)
	{
		pl.add_double(name);
	}
	static void write(const ParameterList& pl, double data, int index)
	{
		pl.set_double(index, data);
	}
	static double read(const ParameterList& pl, int index)
	{
		return pl.to_double(index);
	}
};


//////////////////////////////
// ParamToTypeValueList
//////////////////////////////

// ParamToTypeValueList
template <typename TTypeList, int index = 0>
struct ParamToTypeValueList;

// specialization for empty TypeValueList
template <int index>
struct ParamToTypeValueList<TypeList<>, index> :
	public TypeValueList<TypeList<> >
{
	ParamToTypeValueList(const ParameterList& in) {}
};

// implementation for non-empty TypeValueList
template <typename TTypeList, int index>
struct ParamToTypeValueList :
	public TypeValueList<TTypeList>
{
	typedef	typename TTypeList::head	head;
	typedef	typename TTypeList::tail	tail;

	ParamToTypeValueList(const ParameterList& in) :
		TypeValueList<TTypeList>(	PLStack<head>::read(in, index),
									ParamToTypeValueList<tail, index+1>(in))
			{}
};

//////////////////////////////
// ParamToTypeValueList
//////////////////////////////

// implementation for non-empty
template <typename TTypeList>
struct CreateParameterList
{
	typedef	typename TTypeList::head	head;
	typedef	typename TTypeList::tail	tail;

	static void create(ParameterList& in, const char* paramValNames, const char* delimiter)
	{
		std::string names(paramValNames);
		std::string delim(delimiter);

		// find first param name
		string::size_type begin = names.find_first_not_of(delim, 0);
	    string::size_type end   = names.find_first_of(delim, begin);

	    // add parameter so list
		PLStack<head>::add(in, names.substr(begin, end - begin).c_str());

		// recursive call for remaining types and names
		CreateParameterList<tail>::create(in, names.substr(end, string::npos).c_str(), delimiter);
	}
};

// specialization for empty typelist
template <>
struct CreateParameterList<TypeList<> >
{
	static void create(ParameterList& in, const char* paramValNames, const char* delimiter)
	{
		// do nothing and end recursion
	}
};


} // end namespace interface

} // end namespace ug

#endif /* __H__UG_INTERFACE__UGBRIDGE__PARAM_TO_TYPE_VALUE_LIST__ */
