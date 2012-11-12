//	authors: Sebastian Reiter, Andreas Vogel

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
#include "common/profiler/dynamic_profiling.h"

#endif

namespace ug
{
namespace bridge
{

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
			typedef typename func_traits<TFunc>::params_type params_type;
			CreateParameterInfo<params_type>::create(m_paramsIn);

		//	arbitrary choosen minimum number of infos exported
		//	(If values non given we set them to an empty string)
			const size_t MinNumInfos = 3; // for "name | style | options"

		//	Fill missing Parameter
			m_vvParamInfo.resize(m_paramsIn.size());

		//	resize missing infos for each parameter
			for(int i = 0; i < (int)m_vvParamInfo.size(); ++i)
				for(size_t j = m_vvParamInfo.at(i).size(); j < MinNumInfos; ++j)
					m_vvParamInfo.at(i).push_back(std::string(""));

		////////////////////////////////////////////////
		//	Create parameter stack for RETURN VALUES
		////////////////////////////////////////////////
			typedef typename func_traits<TFunc>::return_type return_type;
			CreateParameterInfoOut<return_type>::create(m_paramsOut);

		//	resize missing infos for return value
			for(size_t j = m_vRetValInfo.size(); j < MinNumInfos; ++j)
				m_vRetValInfo.push_back(std::string(""));
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
		typedef void (*ProxyFunc)(void* func, const ParameterStack& in, ParameterStack& out);

		template <typename TFunc>
		ExportedFunction(	TFunc f, ProxyFunc pf,
							const std::string& name, const std::string& funcOptions,
							const std::string& group,
							const std::string& retValInfos, const std::string& paramInfos,
							const std::string& tooltip, const std::string& help)
			: ExportedFunctionBase(name, funcOptions, retValInfos,
			                       paramInfos, tooltip, help),
			  m_group(group), m_func((void*)f), m_proxy_func(pf)
		{
#ifdef PROFILE_BRIDGE
			m_dpi.init(ExportedFunctionBase::name().c_str(), true, "registry", false);
#endif
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
			m_dpi.endCurNode();
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
		mutable DynamicProfileInformation m_dpi;
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
		template <class TFunc>
		bool add_overload(TFunc f, ExportedFunction::ProxyFunc pf,
		                  const std::string& funcOptions, const std::string& group,
		                  const std::string& retValInfos, const std::string& paramInfos,
		                  const std::string& tooltip, const std::string& help)
		{
			size_t typeID = GetUniqueTypeID<TFunc>();

		//	make sure that the overload didn't exist
			if(get_overload_by_type_id(typeID))return false;

		//	create a new overload
			ExportedFunction* func = new ExportedFunction(f, pf, m_name,
												funcOptions, group, retValInfos,
												paramInfos, tooltip, help);

			m_overloads.push_back(Overload(func, typeID));
			return true;
		}

		size_t num_overloads() const {return m_overloads.size();}

		ExportedFunction* get_overload(size_t index)
			{return m_overloads.at(index).m_func;}

		const ExportedFunction* get_overload(size_t index) const
			{return m_overloads.at(index).m_func;}

		template <class TType>
		ExportedFunction* get_overload_by_type()
		{
			size_t typeID = GetUniqueTypeID<TType>();
			return get_overload_by_type_id(typeID);
		}

		template <class TType>
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
			return NULL;
		}

		const ExportedFunction* get_overload_by_type_id(size_t typeID) const
		{
			for(size_t i = 0; i < m_overloads.size(); ++i){
				if(m_overloads[i].m_typeID == typeID)
					return m_overloads[i].m_func;
			}
			return NULL;
		}

		size_t get_overload_type_id(size_t index) const
			{return m_overloads.at(index).m_typeID;}

	private:
		struct Overload{
			Overload()	{}
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
		typedef typename func_traits<TFunc>::params_type params_type;
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
		typedef typename func_traits<TFunc>::params_type params_type;
		TFunc fp = (TFunc) func;

	//  convert parameter stack
		ParameterStackToTypeValueList<params_type> args(in);

	//  apply
		func_traits<TFunc>::apply(fp, args);
	}
};


} // end namespace bridge
} // end namespace ug


#endif /* __H__UG_BRIDGE__GLOBAL_FUNCTION__ */
