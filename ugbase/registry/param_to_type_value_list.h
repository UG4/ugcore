/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
 * Authors: Sebastian Reiter, Andreas Vogel
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

#ifndef __H__UG_BRIDGE__PARAM_TO_TYPE_VALUE_LIST__
#define __H__UG_BRIDGE__PARAM_TO_TYPE_VALUE_LIST__

#include "parameter_stack.h"
#include "class_name_provider.h"
#include <string>
#include <vector>
#include <cstring>
#include "function_traits.h"

namespace ug
{
namespace bridge
{

/// \addtogroup registry
/// \{

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
	using head = typename TTypeList::head;
	using tail = typename TTypeList::tail;

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
	using head = typename TTypeList::head;
	using tail = typename TTypeList::tail;

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

template <>
struct CreateParameterInfoOut<CustomReturn> {
	static void create(ParameterInfo& stack){}		
};

///	Creation of void return value (template specialization)
template <>
struct CreateParameterInfoOut<void>{
	static void create(ParameterInfo& stack){}
};

// end group registry
/// \}

} // end namespace bridge

} // end namespace ug

#endif