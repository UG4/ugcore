/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel
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

#include <iostream>
#include <sstream>
#include <string>

// include bridge
#include "bridge/bridge.h"
#include "bridge/util.h"

#include "lua_user_data.h"

using namespace std;

namespace ug
{
	

///	returns true if callback exists
bool CheckLuaCallbackName(const char* name)
{
//	get lua state
	lua_State* L = ug::script::GetDefaultLuaState();

//	obtain a reference
	lua_getglobal(L, name);

//	check if reference is valid
	if(lua_isnil(L, -1)) return false;
	else return true;
}


LuaUserNumberNumberFunction::LuaUserNumberNumberFunction()
{
	m_L = ug::script::GetDefaultLuaState();
	m_callbackRef = LUA_NOREF;
}

void LuaUserNumberNumberFunction::set_lua_callback(const char* luaCallback)
{
	m_callbackName = luaCallback;
//	store the callback function in the registry and obtain a reference.
	lua_getglobal(m_L, m_callbackName.c_str());

//	make sure that the reference is valid
	if(lua_isnil(m_L, -1)){
		UG_THROW("ERROR in LuaUserNumberNumberFunction::set_lua_callback(...):"
				"Specified callback does not exist: " << m_callbackName);
	}

	m_callbackRef = luaL_ref(m_L, LUA_REGISTRYINDEX);
}

number LuaUserNumberNumberFunction::operator() (int numArgs, ...) const
{
//	push the callback function on the stack
	lua_rawgeti(m_L, LUA_REGISTRYINDEX, m_callbackRef);

	va_list ap;
	va_start(ap, numArgs);

	for(int i = 0; i < numArgs; ++i)
	{
		number val = va_arg(ap, number);
		lua_pushnumber(m_L, val);
//		UG_LOG("Push value i=" << i << ": " << val<<"\n");
	}

	va_end(ap);


	if(lua_pcall(m_L, numArgs, 1, 0) != 0)
	{
		UG_THROW("ERROR in 'LuaUserNumberNumberFunction::operator(...)': Error while "
				 "running callback '" << m_callbackName << "',"
				 " lua message: "<< lua_tostring(m_L, -1) << "\n");
	}

	number c = ReturnValueToNumber(m_L, -1);
	lua_pop(m_L, 1);

	return c;
}




namespace bridge{
namespace LuaUserData{


template <typename TData, int dim>
void RegisterLuaUserDataType(Registry& reg, string type, string grp)
{
	string suffix = GetDimensionSuffix<dim>();
	string tag = GetDimensionTag<dim>();

//	LuaUser"Type"
	{
		typedef ug::LuaUserData<TData, dim> T;
		typedef CplUserData<TData, dim> TBase;
		string name = string("LuaUser").append(type).append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.template add_constructor<void (*)(const char*)>("Callback")
			.template add_constructor<void (*)(LuaFunctionHandle)>("handle")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, string("LuaUser").append(type), tag);
	}

//	LuaCondUser"Type"
	{
		typedef ug::LuaUserData<TData, dim, bool> T;
		typedef CplUserData<TData, dim, bool> TBase;
		string name = string("LuaCondUser").append(type).append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.template add_constructor<void (*)(const char*)>("Callback")
			.template add_constructor<void (*)(LuaFunctionHandle)>("handle")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, string("LuaCondUser").append(type), tag);
	}
}

/**
 * Class exporting the functionality. All functionality that is to
 * be used in scripts or visualization must be registered here.
 */
struct Functionality
{

/**
 * Function called for the registration of Dimension dependent parts.
 * All Functions and Classes depending on the Dimension
 * are to be placed here when registering. The method is called for all
 * available Dimension types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <int dim>
static void Dimension(Registry& reg, string grp)
{
	string suffix = GetDimensionSuffix<dim>();
	string tag = GetDimensionTag<dim>();

	RegisterLuaUserDataType<number, dim>(reg, "Number", grp);
	RegisterLuaUserDataType<MathVector<dim>, dim>(reg, "Vector", grp);
	RegisterLuaUserDataType<MathMatrix<dim,dim>, dim>(reg, "Matrix", grp);

//	LuaUserFunctionNumber
	{
		typedef LuaUserFunction<number, dim, number> T;
		typedef DependentUserData<number, dim> TBase;
		string name = string("LuaUserFunctionNumber").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.template add_constructor<void (*)(const char*, int)>("LuaCallbackName#NumberOfArguments")
			.template add_constructor<void (*)(const char*, int, bool)>("LuaCallbackName#NumberOfArguments#PosTimeFlag")
			.template add_constructor<void (*)(LuaFunctionHandle, int)>("LuaCallbackName#NumberOfArguments")
			.template add_constructor<void (*)(LuaFunctionHandle, int, bool)>("LuaCallbackName#NumberOfArguments#PosTimeFlag")			
			.add_method("set_deriv", static_cast<void (T::*)(size_t, const char*)>(&T::set_deriv))
			.add_method("set_deriv", static_cast<void (T::*)(size_t, LuaFunctionHandle)>(&T::set_deriv))
			.add_method("set_input", static_cast<void (T::*)(size_t, SmartPtr<CplUserData<number, dim> >)>(&T::set_input))
			.add_method("set_input", static_cast<void (T::*)(size_t, number)>(&T::set_input))
			.add_method("set_input_and_deriv", &T::set_input_and_deriv)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "LuaUserFunctionNumber", tag);
	}

//	LuaUserFunctionMatrixNumber
	{
		typedef LuaUserFunction<MathMatrix<dim,dim>, dim, number> T;
		typedef DependentUserData<MathMatrix<dim,dim>, dim> TBase;
		string name = string("LuaUserFunctionMatrixNumber").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.template add_constructor<void (*)(const char*, int)>("LuaCallbackName#NumberOfArguments")
			.template add_constructor<void (*)(const char*, int, bool)>("LuaCallbackName#NumberOfArguments#PosTimeFlag")
			.template add_constructor<void (*)(LuaFunctionHandle, int)>("LuaCallbackName#NumberOfArguments")
			.template add_constructor<void (*)(LuaFunctionHandle, int, bool)>("LuaCallbackName#NumberOfArguments#PosTimeFlag")			
			.add_method("set_deriv", static_cast<void (T::*)(size_t, const char*)>(&T::set_deriv))
			.add_method("set_input", static_cast<void (T::*)(size_t, SmartPtr<CplUserData<number, dim> >)>(&T::set_input))
			.add_method("set_input", static_cast<void (T::*)(size_t, number)>(&T::set_input))
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "LuaUserFunctionMatrixNumber", tag);
	}

//	LuaUserFunctionVectorNumber
	{
		typedef LuaUserFunction<MathVector<dim>, dim, number > T;
		typedef DependentUserData<MathVector<dim>, dim> TBase;
		string name = string("LuaUserFunctionVectorNumber").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.template add_constructor<void (*)(const char*, int)>("LuaCallbackName#NumberOfArguments")
			.template add_constructor<void (*)(const char*, int, bool)>("LuaCallbackName#NumberOfArguments#PosTimeFlag")
			.template add_constructor<void (*)(LuaFunctionHandle, int)>("LuaCallbackName#NumberOfArguments")
			.template add_constructor<void (*)(LuaFunctionHandle, int, bool)>("LuaCallbackName#NumberOfArguments#PosTimeFlag")			
			.add_method("set_deriv", static_cast<void (T::*)(size_t, const char*)>(&T::set_deriv))
			.add_method("set_input", static_cast<void (T::*)(size_t, SmartPtr<CplUserData<number, dim> >)>(&T::set_input))
			.add_method("set_input", static_cast<void (T::*)(size_t, number)>(&T::set_input))
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "LuaUserFunctionVectorNumber", tag);
	}
/*
//	LuaUserFunctionNumberVector
	{
		typedef LuaUserFunction<number, dim, MathVector<dim> > T;
		typedef DependentUserData<number, dim> TBase;
		string name = string("LuaUserFunctionNumberVector").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.template add_constructor<void (*)(const char*, int)>("LuaCallbackName, NumberOfArguments")
			.template add_constructor<void (*)(const char*, int, bool)>("LuaCallbackName, NumberOfArguments, PosTimeFlag")
			.add_method("set_deriv", &T::set_deriv)
			.add_method("set_input", static_cast<void (T::*)(size_t, SmartPtr<CplUserData<number, dim> >)>(&T::set_input))
			.add_method("set_input", static_cast<void (T::*)(size_t, number)>(&T::set_input))
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "LuaUserFunctionNumberVector", tag);
	}
*/
}

/**
 * Function called for the registration of Domain and Algebra independent parts.
 * All Functions and Classes not depending on Domain and Algebra
 * are to be placed here when registering.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
static void Common(Registry& reg, string grp)
{

//	LuaUserNumberNumberFunction
	{
		typedef LuaUserNumberNumberFunction T;
		reg.add_class_<T>("LuaUserNumberNumberFunction", grp)
			.add_constructor()
			.add_method("set_lua_callback", &T::set_lua_callback)
			.set_construct_as_smart_pointer(true);
	}

//	LuaFunctionNumber
	{
		typedef LuaFunction<number, number> T;
		typedef IFunction<number, number> TBase;
		reg.add_class_<T, TBase>("LuaFunctionNumber", grp)
			.add_constructor()
			.add_method("set_lua_callback", &T::set_lua_callback)
			.set_construct_as_smart_pointer(true);
	}

} // Common
}; // end Functionality
}// end LuaUserData

void RegisterLuaUserData(Registry& reg, string grp)
{
	typedef LuaUserData::Functionality Functionality;

	try{
		RegisterCommon<Functionality>(reg,grp);
		RegisterDimensionDependent<Functionality>(reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

}// namespace bridge
}// namespace ug
