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

#ifndef __H__UG_BRIDGE__GLOBAL_FUNCTION__
#define __H__UG_BRIDGE__GLOBAL_FUNCTION__

#include <string>
#include <vector>
#include <sstream>
#include "parameter_stack.h"
#include "function_traits.h"
#include "param_to_type_value_list.h"
#include "common/ug_config.h"
#include <iostream>

#ifdef PROFILE_BRIDGE
#ifndef UG_PROFILER
	#error "You need to define UG_PROFILER to use PROFILE_BRIDGE"
#endif
#include "common/profiler/runtime_profile_info.h"

#endif

namespace ug
{
namespace bridge
{

/// \addtogroup registry
/// \{

///	Exception throw, if method name has not been given
struct UG_REGISTRY_ERROR_FunctionOrMethodNameMissing {};

/// Base class for function/method export
class UG_API ExportedFunctionBase
{
	public:
		ExportedFunctionBase(const std::string& funcName, const std::string& funcOptions,
		                     const std::string& retValInfos, const std::string& paramInfos,
		                     const std::string& tooltip, const std::string& help);

	///	name of function
		const std::string& name() const {return m_name;}

	///	name of function
		const std::string& options() const {return m_methodOptions;}

	/// name of return value
		const std::string& return_name() const {return return_info(0);}

	///	type info of return type
		const std::string& return_info(size_t i) const {return m_vRetValInfo.at(i);}

	/// type info of return value
		const std::vector<std::string>& return_info_vec() const {return m_vRetValInfo;}

	/// number of parameters.
		size_t num_parameter() const {return m_vvParamInfo.size();}

	///	number of info strings for one parameter
		size_t num_infos(size_t i) const {return m_vvParamInfo.at(i).size();}

	/// class name of parameter i
		const char* parameter_class_name(size_t i) const	{return params_in().class_name((int)i);}
		
	/// name of parameter i
		const std::string& parameter_name(size_t i) const {return parameter_info(i, 0);}

	///	type info of all parameters
		const std::string& parameter_info(size_t i, size_t j) const	{return m_vvParamInfo.at(i).at(j);}

	/// type info of i th parameters
		const std::vector<std::string>& parameter_info_vec(size_t i) const {return m_vvParamInfo.at(i);}

	///	whole string of all type infos for of all parameters
		const std::string& parameter_info_string() const {return m_paramInfos;}

	/// gives some information to the exported functions
		const std::string& tooltip() const {return m_tooltip;}

	/// help informations
		const std::string& help() const {return m_help;}

	/// parameter list for input values
		const ParameterInfo& params_in() const	{return m_paramsIn;}

	/// parameter list for input values
		const ParameterInfo& params_out() const {return m_paramsOut;}

	// todo: we export non-const here, since we can not make ExportedClass<TClass> a friend
	/// non-const export of param list
		ParameterInfo& params_in() {return m_paramsIn;}

	/// returns true if all parameters of the function are correctly declared
		bool check_consistency(std::string classname = "") const;

	protected:
		template <typename TFunc>
		void create_parameter_stack()
		{
		////////////////////////////////////////////////
		//	Create parameter stack for PARAMETERS
		////////////////////////////////////////////////
			using params_type = typename func_traits<TFunc>::params_type;
			CreateParameterInfo<params_type>::create(m_paramsIn);

		//	arbitrary choosen minimum number of infos exported
		//	(If values non given we set them to an empty string)
			const size_t MinNumInfos = 3; // for "name | style | options"

		//	Fill missing Parameter
			m_vvParamInfo.resize(m_paramsIn.size());
		//	resize missing infos for each parameter
			for(int i = 0; i < (int)m_vvParamInfo.size(); ++i)
				for(size_t j = m_vvParamInfo.at(i).size(); j < MinNumInfos; ++j)
					m_vvParamInfo.at(i).emplace_back("");

		////////////////////////////////////////////////
		//	Create parameter stack for RETURN VALUES
		////////////////////////////////////////////////
			using return_type = typename func_traits<TFunc>::return_type;
			CreateParameterInfoOut<return_type>::create(m_paramsOut);

		//	resize missing infos for return value
			for(size_t j = m_vRetValInfo.size(); j < MinNumInfos; ++j)
				m_vRetValInfo.emplace_back("");
		}

	// 	help function to tokenize the parameter string
		void tokenize(const std::string& str, std::vector<std::string>& tokens,
		              const char delimiter);

	protected:
		std::string m_name;
		std::string m_methodOptions;

		std::string m_retValInfos; // string with Infos about return type
		std::vector<std::string> m_vRetValInfo; // tokenized Infos

	// 	string with Infos about parameter
		std::string m_paramInfos;

	// 	tokenized strings for each Parameter and each Info (name | #style | options | ...)
		std::vector<std::vector<std::string> > m_vvParamInfo;

		std::string m_tooltip;
		std::string m_help;

		ParameterInfo m_paramsIn;
		ParameterInfo m_paramsOut;
};


///This class describes a wrapper for a c++ - function, that is exported by ug
class UG_API ExportedFunction : public ExportedFunctionBase
{
	public:
	//	all c++ functions are wrapped by a proxy function of the following type
		using ProxyFunc = void(*)(void* func, const ParameterStack& in, ParameterStack& out);

		template <typename TFunc>
		ExportedFunction(	TFunc f, ProxyFunc pf,
							const std::string& name, const std::string& funcOptions,
							const std::string& group,
							const std::string& retValInfos, const std::string& paramInfos,
							const std::string& tooltip, const std::string& help)
			: ExportedFunctionBase(name, funcOptions, retValInfos,
			                       paramInfos, tooltip, help),
			  m_group(group), m_func((void*)f), m_proxy_func(pf)
#ifdef PROFILE_BRIDGE
			  ,m_dpi(ExportedFunctionBase::name().c_str(), true, "registry", false)
#endif
		{
			create_parameter_stack<TFunc>();
		}

	/// executes the function
		void execute(const ParameterStack& paramsIn, ParameterStack& paramsOut) const
		{
#ifdef PROFILE_BRIDGE
			m_dpi.beginNode();
#endif
			m_proxy_func(m_func, paramsIn, paramsOut);

#ifdef PROFILE_BRIDGE
			m_dpi.endNode();
#endif

		}

	///	return groups
		const std::string& group() const {return m_group;}

	protected:
	/// save groups
		std::string m_group;

	/// pointer to to function
		void* m_func;

	/// proxy function
		ProxyFunc m_proxy_func;

#ifdef PROFILE_BRIDGE
		mutable RuntimeProfileInfo m_dpi;
#endif
};

////////////////////////////////////////////////////////////////////////
//	ExportedFunctionGroup (sreiter)
///	Groups of Functions - useful to realize overloaded functions
class UG_API ExportedFunctionGroup
{
	public:
	///	constructor
		ExportedFunctionGroup(const std::string& name) : m_name(name){}

	///	destructor
		~ExportedFunctionGroup()
		{
			for(size_t i = 0; i < m_overloads.size(); ++i)
				delete m_overloads[i].m_func;
		}

	///	name of function group
		const std::string& name() const {return m_name;}

	///	adds an overload. Returns false if the overload already existed.
		template <typename TFunc>
		ExportedFunction*
		add_overload (	TFunc f, ExportedFunction::ProxyFunc pf,
		                const std::string& funcOptions, const std::string& group,
		                const std::string& retValInfos, const std::string& paramInfos,
		                const std::string& tooltip, const std::string& help)
		{
			size_t typeID = GetUniqueTypeID<TFunc>();

		//	make sure that the overload didn't exist
			if(get_overload_by_type_id(typeID))return nullptr;

		//	create a new overload
			auto* func = new ExportedFunction(f, pf, m_name,
												funcOptions, group, retValInfos,
												paramInfos, tooltip, help);

			m_overloads.push_back(Overload(func, typeID));
			return func;
		}

		size_t num_overloads() const {return m_overloads.size();}

		ExportedFunction* get_overload(size_t index)
			{return m_overloads.at(index).m_func;}

		const ExportedFunction* get_overload(size_t index) const
			{return m_overloads.at(index).m_func;}

		template <typename TType>
		ExportedFunction* get_overload_by_type()
		{
			size_t typeID = GetUniqueTypeID<TType>();
			return get_overload_by_type_id(typeID);
		}

		template <typename TType>
		const ExportedFunction* get_overload_by_type() const
		{
			size_t typeID = GetUniqueTypeID<TType>();
			return get_overload_by_type_id(typeID);
		}

		ExportedFunction* get_overload_by_type_id(size_t typeID)
		{
			for(size_t i = 0; i < m_overloads.size(); ++i){
				if(m_overloads[i].m_typeID == typeID)
					return m_overloads[i].m_func;
			}
			return nullptr;
		}

		const ExportedFunction* get_overload_by_type_id(size_t typeID) const
		{
			for(size_t i = 0; i < m_overloads.size(); ++i){
				if(m_overloads[i].m_typeID == typeID)
					return m_overloads[i].m_func;
			}
			return nullptr;
		}

		size_t get_overload_type_id(size_t index) const
			{return m_overloads.at(index).m_typeID;}

	private:
		struct Overload{
			Overload()	= default;
			Overload(ExportedFunction* func, size_t typeID) :
				m_func(func), m_typeID(typeID) {}
			ExportedFunction* 	m_func;
			size_t				m_typeID;
		};

		std::string m_name;
		std::vector<Overload>	m_overloads;
};

/**
 * Function that forwards for function pointer for a signature
 */
template <typename TFunc, typename TRet = typename func_traits<TFunc>::return_type>
struct FunctionProxy
{
	static void apply(void* func, const ParameterStack& in, ParameterStack& out)
	{
		using params_type = typename func_traits<TFunc>::params_type;
		TFunc fp = (TFunc) func;

	//  convert parameter stack
		ParameterStackToTypeValueList<params_type> args(in);

	//  apply
		TRet res = func_traits<TFunc>::apply(fp, args);

	//  write result
		out.push<TRet>(res);
	}
};

// specialization if no return value given
template <typename TFunc>
struct FunctionProxy<TFunc, void>
{
	static void apply(void* func, const ParameterStack& in, ParameterStack& out)
	{
		using params_type = typename func_traits<TFunc>::params_type;
		TFunc fp = (TFunc) func;

	//  convert parameter stack
		ParameterStackToTypeValueList<params_type> args(in);

	//  apply
		func_traits<TFunc>::apply(fp, args);
	}
};

// end group registry
/// \}

} // end namespace bridge
} // end namespace ug


#endif