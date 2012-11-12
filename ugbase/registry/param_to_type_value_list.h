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

////////////////////////////////////////////////////////////////////////////////
// ParameterStackToTypeValueList
////////////////////////////////////////////////////////////////////////////////

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
		TypeValueList<TTypeList>(	in.to<head>(index),
									ParameterStackToTypeValueList<tail, index+1>(in))
			{}
};

////////////////////////////////////////////////////////////////////////////////
// CreateParameterInfo
////////////////////////////////////////////////////////////////////////////////

// implementation for non-empty
template <typename TTypeList>
struct CreateParameterInfo
{
	typedef	typename TTypeList::head	head;
	typedef	typename TTypeList::tail	tail;

	static void create(ParameterInfo& in)
	{
	    // add parameter so list
		in.push_type<head>();

		// recursive call for remaining types and names
		CreateParameterInfo<tail>::create(in);
	}
};

// specialization for empty typelist
template <>
struct CreateParameterInfo<TypeList<> >
{
	static void create(ParameterInfo& in)
	{
		// do nothing and end recursion
	}
};

///	Creation of return value
template <typename TRet>
struct CreateParameterInfoOut {
	static void create(ParameterInfo& stack){
		CreateParameterInfo<TypeList<TRet> >::create(stack);
}};

///	Creation of void return value (template specialization)
template <>
struct CreateParameterInfoOut<void>{
	static void create(ParameterInfo& stack){}
};


} // end namespace bridge

} // end namespace ug

#endif /* __H__UG_BRIDGE__PARAM_TO_TYPE_VALUE_LIST__ */
