/*
 * Copyright (c) 2013-2014:  G-CSC, Goethe University Frankfurt
 * Author: Martin Rupp
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

#ifndef LUA_TRAITS_H
#define	LUA_TRAITS_H

#include "common/util/number_util.h"
#include "externals/lua/lua.h"

namespace ug{
///////////////////////////////////////////////////////////////////////////////
// Lua Traits
///////////////////////////////////////////////////////////////////////////////


/// Helper to access a return value on the stack.
/**	If the value can't be converted to a number, an error is thrown*/
inline number ReturnValueToNumber(lua_State* L, int index){
	if(!lua_isnumber(L, index))
		UG_THROW("ReturnValueToNumber: Data passed from Lua: "
						"Can't convert return value to number!");
	number a = lua_tonumber(L, index);
	UG_COND_THROW(IsFiniteAndNotTooBig(a) == false, a << " not finite");
	return a;

}

/// Helper to access a return value on the stack.
/**	If the value can't be converted to a boolean, an error is thrown*/
inline number ReturnValueToBool(lua_State* L, int index){
	if(!lua_isboolean(L, index))
		UG_THROW("ReturnValueToBool: Data passed from Lua: "
						"Can't convert return value to boolean!");

	return lua_toboolean(L, index);
}

/// Helper to access a return value on the stack.
/**	If the value can't be converted to a integer, an error is thrown*/
inline int ReturnValueToInteger(lua_State* L, int index){
	if(!lua_isnumber(L, index))
		UG_THROW("ReturnValueToBool: Data passed from Lua: "
						"Can't convert return value to integer!");

	return lua_tointeger(L, index);
}


/// Lua Traits to push/pop on lua stack
template <typename TData>
struct lua_traits;

template <>
struct lua_traits<void>
{
	static constexpr int size = 0;

	static void push(lua_State*	L, const bool&){}

	static void read(lua_State* L, bool&, int index = -1){}

	static void do_return(const bool&) {}

	static bool check(lua_State* L, int index = -1){return true;}

	static std::string signature()
	{
		return "void";
	}

	static std::string name()
	{
		return "void";
	}
};

template <>
struct lua_traits<bool>
{
	static constexpr int size = 1;

	static void push(lua_State*	L, const bool& b)
	{
		lua_pushboolean(L, b);
	}

	static void read(lua_State* L, bool& b, int index = -1)
	{
		b =  ReturnValueToBool(L, index);
	}

	static bool check(lua_State* L, int index = -1)
	{
		return lua_isboolean(L, index);
	}

	static bool do_return(const bool& b) {return b;}

	static std::string signature()
	{
		return "bool";
	}

	static std::string name()
	{
		return "Bool";
	}
};

template <>
struct lua_traits<int>
{
	static constexpr int size = 1;

	static void push(lua_State*	L, const int& c)
	{
		lua_pushinteger(L, c);
	}

	static void read(lua_State* L, int& c, int index = -1)
	{
		c = ReturnValueToInteger(L, index);
	}

	static bool check(lua_State* L, int index = -1)
	{
		return lua_isnumber(L, index);
	}

	static std::string signature()
	{
		return "integer";
	}

	static std::string name()
	{
		return "Integer";
	}
};

template <>
struct lua_traits<number>
{
	static constexpr int size = 1;

	static void push(lua_State*	L, const number& c)
	{
		lua_pushnumber(L, c);
	}

	static void read(lua_State* L, number& c, int index = -1)
	{
		c = ReturnValueToNumber(L, index);
	}
	
	static void read(number &c, double ret[1], bool *dummy)
	{
		c = ret[1];
	}
	static void read(number &c, double ret[1], void *dummy)
	{
		c = ret[0];
	}	

	static bool check(lua_State* L, int index = -1)
	{
		return lua_isnumber(L, index);
	}

	static std::string signature()
	{
		return "number";
	}

	static std::string name()
	{
		return "Number";
	}
};

template <std::size_t dim>
struct lua_traits< MathVector<dim> >
{
	static constexpr int size = dim;

	static void push(lua_State*	L, const MathVector<dim>& x)
	{
		for(size_t i = 0; i < dim; ++i)
			lua_pushnumber(L, x[i]);
	}

	static void read(lua_State* L, MathVector<dim>& x, int index = -1)
	{
		for(size_t i = 0; i < dim; ++i){
				x[dim-1-i] = ReturnValueToNumber(L, index--);
		}
	}
	
	static void read(MathVector<dim>& x, const double ret[dim], bool *dummy)
	{
		for(size_t i=0; i<dim; i++)
			x[i] = ret[i+1];
	}
	
	static void read(MathVector<dim>& x, const double ret[dim], void *dummy)
	{
		for(size_t i=0; i<dim; i++)
			x[i] = ret[i];
	}

	static bool check(lua_State* L, int index = -1)
	{
		for(size_t i = 0; i < dim; ++i){
			if(!lua_isnumber(L, index--)) return false;
		}
		return true;
	}

	static std::string signature()
	{
		static char cmp[] = {'x', 'y', 'z'};
		std::stringstream ss;
		for(size_t i = 0; i < dim; ++i){
			if(i != 0) ss << ", ";
			ss << "v" << cmp[i];
		}
		return ss.str();
	}

	static std::string name()
	{
		std::stringstream ss;
		ss << "Vector";
		return ss.str();
	}
};


template <std::size_t dim>
struct lua_traits< MathMatrix<dim, dim> >
{
	static constexpr int size = dim*dim;

	static void push(lua_State*	L, const MathMatrix<dim, dim>& D)
	{
		for(size_t i = 0; i < dim; ++i){
			for(size_t j = 0; j < dim; ++j){
				lua_pushnumber(L, D[i][j]);
			}
		}

	}

	static void read(lua_State* L, MathMatrix<dim, dim>& D, int index = -1)
	{
		for(size_t i = 0; i < dim; ++i){
			for(size_t j = 0; j < dim; ++j){
				D[dim-1-i][dim-1-j] = ReturnValueToNumber(L, index--);
			}
		}
	}

	static void read(MathMatrix<dim, dim>& D, const double ret[size], bool *dummy)
	{
		for(size_t i = 0; i < dim; ++i){
			for(size_t j = 0; j < dim; ++j){
				D[i][j] = ret[i*dim+j+1];
			}
		}
	}
	static void read(MathMatrix<dim, dim>& D, const double ret[size], void *dummy)
	{
		for(size_t i = 0; i < dim; ++i){
			for(size_t j = 0; j < dim; ++j){
				D[i][j] = ret[i*dim+j];
			}
		}
	}

	static bool check(lua_State* L, int index = -1)
	{
		for(size_t i = 0; i < dim; ++i){
			for(size_t j = 0; j < dim; ++j){
				if(!lua_isnumber(L, index--)) return false;
			}
		}
		return true;
	}

	static std::string signature()
	{
		static char cmp[] = {'x', 'y', 'z'};
		std::stringstream ss;
		for(size_t i = 0; i < dim; ++i){
			for(size_t j = 0; j < dim; ++j){
				if(i != 0 || j != 0) ss << ", ";
				ss << "D" << cmp[i] << cmp[j];
			}
		}
		return ss.str();
	}

	static std::string name()
	{
		std::stringstream ss;
		ss << "Matrix";
		return ss.str();
	}
};
	
/* for debugging callbacks
 * 
		if(compare(out, out2)==false)
		{
			UG_LOG("\neval_value------------\n" << m_luaC.name << "\n"  << out << "\n--\n" << out2 << "\n------------\n");
			for(int i=0; i<lua_traits<TData>::size+1; i++)
				UG_LOG(ret[i] << "\n");
		}

inline bool compare(double &a, double &b)
{
	return a == b;
}

inline bool compare(MathVector<3, double> &x, MathVector<3, double> &y)
{
	for(int i=0; i<3; i++)
		if(x[i] != y[i]) return false;
	return true;
}

inline bool compare(MathVector<2, double> &x, MathVector<2, double> &y)
{
	for(int i=0; i<2; i++)
		if(x[i] != y[i]) return false;
	return true;
}

inline bool compare(MathMatrix<3, 3, double> &A, MathMatrix<3, 3, double> &B)
{
	for(int j=0; j<3; j++)
		for(int i=0; i<3; i++)
			if(A[i][j] != B[i][j])
				return false;
	return true;
	
}
inline bool compare(MathMatrix<2, 2, double> &A, MathMatrix<2, 2, double> &B)
{
	for(int j=0; j<2; j++)
		for(int i=0; i<2; i++)
			if(A[i][j] != B[i][j])
				return false;
	return true;
	
}*/

}
#endif	/* LUA_TRAITS_H */


