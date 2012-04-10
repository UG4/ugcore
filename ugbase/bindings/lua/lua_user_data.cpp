// author: Andreas Vogel

#include <iostream>
#include <sstream>

#include "registry/registry.h"
#include "bridge/bridge.h"

#include "lua_user_data.h"

using namespace std;

namespace ug
{



LuaUserNumberNumberFunction::LuaUserNumberNumberFunction()
{
	m_L = ug::script::GetDefaultLuaState();
	m_callbackRef = LUA_NOREF;
}

void LuaUserNumberNumberFunction::set_lua_callback(const char* luaCallback)
{
	m_callbackName = luaCallback;
//	store the callback function in the registry and obtain a reference.
	lua_getglobal(m_L, m_callbackName);

//	make sure that the reference is valid
	if(lua_isnil(m_L, -1)){
		UG_THROW_FATAL("ERROR in LuaUserNumberNumberFunction::set_lua_callback(...):"
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
		std::stringstream ss;
		ss << "ERROR in 'LuaUserNumberNumberFunction::operator(...)': Error while "
				"running callback '" << m_callbackName << "',"
				" lua message: "<< lua_tostring(m_L, -1) << "\n";
		throw(UGError(true, ss.str().c_str()));
	}

	number c = ReturnValueToNumber(m_L, -1);
	lua_pop(m_L, 1);

	return c;
}




namespace bridge
{

///////////////////////////////////////////////////////////////////////////////
// Registration
///////////////////////////////////////////////////////////////////////////////

template <typename TData, int dim>
bool RegisterLuaUserDataType(Registry& reg, string type, const char* parentGroup)
{
	string grp = std::string(parentGroup);

	string dimSuffix = GetDomainSuffix<dim>();
	string dimTag = GetDomainTag<dim>();

//	LuaUser"Type"
	{
		typedef LuaUserData<TData, dim> T;
		typedef IPData<TData, dim> TBase;
		typedef boost::function<void (TData& res, const MathVector<dim>& x,number time)> TBase2;
		string name = string("LuaUser").append(type).append(dimSuffix);
		reg.add_class_<T, TBase, TBase2>(name, grp)
			.template add_constructor<void (*)(const char*)>("Callback")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, string("LuaUser").append(type), dimTag);
	}

//	LuaCondUser"Type"
	{
		typedef LuaUserData<TData, dim, bool> T;
		typedef IPData<TData, dim, bool> TBase;
		typedef boost::function<bool (TData& res, const MathVector<dim>& x,number time)> TBase2;
		string name = string("LuaCondUser").append(type).append(dimSuffix);
		reg.add_class_<T, TBase, TBase2>(name, grp)
			.template add_constructor<void (*)(const char*)>("Callback")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, string("LuaCondUser").append(type), dimTag);
	}

//	LuaBoundary"Type"
	{
		typedef LuaBoundaryData<TData, dim> T;
		typedef boost::function<bool (TData& res, const MathVector<dim>& x,number time)> TBase;
		string name = string("LuaBoundary").append(type).append(dimSuffix);
		reg.add_class_<T, TBase>(name, grp)
			.template add_constructor<void (*)(const char*)>("Callback")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, string("LuaBoundary").append(type), dimTag);
	}

	return true;
}

template <int dim>
void RegisterLuaUserData(ug::bridge::Registry& reg, const char* parentGroup)
{
	string grp = std::string(parentGroup);

	string dimSuffix = GetDomainSuffix<dim>();
	string dimTag = GetDomainTag<dim>();

	RegisterLuaUserDataType<number, dim>(reg, "Number", parentGroup);
	RegisterLuaUserDataType<MathVector<dim>, dim>(reg, "Vector", parentGroup);
	RegisterLuaUserDataType<MathMatrix<dim,dim>, dim>(reg, "Matrix", parentGroup);

//	LuaUserFunctionNumber
	{
		typedef LuaUserFunction<number, dim, number> T;
		typedef DataLinkerEqualData<number, dim, number> TBase;
		string name = string("LuaUserFunctionNumber").append(dimSuffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.add_method("set_lua_value_callback", &T::set_lua_value_callback)
			.add_method("set_lua_deriv_callback", &T::set_lua_deriv_callback)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "LuaUserFunctionNumber", dimTag);
	}

//	LuaUserFunctionMatrixNumber
	{
		typedef LuaUserFunction<MathMatrix<dim,dim>, dim, number> T;
		typedef DataLinkerEqualData<MathMatrix<dim,dim>, dim, number> TBase;
		string name = string("LuaUserFunctionMatrixNumber").append(dimSuffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.add_method("set_lua_value_callback", &T::set_lua_value_callback)
			.add_method("set_lua_deriv_callback", &T::set_lua_deriv_callback)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "LuaUserFunctionMatrixNumber", dimTag);
	}

}

void RegisterLuaUserData(Registry& reg, const char* parentGroup)
{

//	LuaUserNumberNumberFunction
	{
		typedef LuaUserNumberNumberFunction T;
		reg.add_class_<T>("LuaUserNumberNumberFunction", parentGroup)
			.add_constructor()
			.add_method("set_lua_callback", &T::set_lua_callback)
			.set_construct_as_smart_pointer(true);
	}

//	LuaFunctionNumber
	{
		typedef LuaFunction<number, number> T;
		typedef IFunction<number, number> TBase;
		reg.add_class_<T, TBase>("LuaFunctionNumber", parentGroup)
			.add_constructor()
			.add_method("set_lua_callback", &T::set_lua_callback)
			.set_construct_as_smart_pointer(true);
	}

#ifdef UG_DIM_1
	RegisterLuaUserData<1>(reg, parentGroup);
#endif
#ifdef UG_DIM_2
	RegisterLuaUserData<2>(reg, parentGroup);
#endif
#ifdef UG_DIM_3
	RegisterLuaUserData<3>(reg, parentGroup);
#endif
}

} // end namespace
} // end namepace
