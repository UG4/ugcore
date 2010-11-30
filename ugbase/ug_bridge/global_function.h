
#ifndef __H__UG_BRIDGE__GLOBAL_FUNCTION__
#define __H__UG_BRIDGE__GLOBAL_FUNCTION__

#include <string>
#include <vector>
#include <sstream>
#include "parameter_stack.h"
#include "function_traits.h"
#include "param_to_type_value_list.h"
#include "ug_bridge/class_helper.h"
#include <iostream>

namespace ug
{
namespace bridge
{

template <typename TRet>
struct CreateParameterOutStack
{
	static void create(ParameterStack& stack)
	{
		CreateParameterStack<TypeList<TRet> >::create(stack);
	}
};


template <>
struct CreateParameterOutStack<void>
{
	static void create(ParameterStack& stack)
	{
	}
};

struct UG_REGISTRY_ERROR_FunctionOrMethodNameMissing {};

/** Base class for function/method export
 */
class ExportedFunctionBase
{
	public:
		ExportedFunctionBase(	const char* funcInfos, const char* retValInfos, const char* paramInfos,
								const char* tooltip, const char* help)
		: m_funcInfos(funcInfos), m_retValInfos(retValInfos), m_paramInfos(paramInfos),
		  m_tooltip(tooltip), m_help(help)
		{
		//	get name and visualization options of function
			std::vector<std::string> vFuncInfoTmp;
			tokenize(m_funcInfos, vFuncInfoTmp, '|');

		//	set name
			if(vFuncInfoTmp.size() >= 1) m_name = trim(vFuncInfoTmp[0]);
			else throw(UG_REGISTRY_ERROR_FunctionOrMethodNameMissing());

		//	set options if given
			if(vFuncInfoTmp.size() >= 2) m_methodOptions = trim(vFuncInfoTmp[1]);
			else m_methodOptions = "";

		//	other fields are neglected

		//	Tokenize string for return value (separated by '|')
			tokenize(m_retValInfos, m_vRetValInfo, '|');

		//	Tokenize string for parameters into infos per one parameter (separated by '#')
			std::vector<std::string> vParamInfoTmp;
			tokenize(m_paramInfos, vParamInfoTmp, '#');
			m_vvParamInfo.resize(vParamInfoTmp.size());

		//	Tokenite each info-string of one parameter into single infos (separated by '|')
			for(size_t i = 0; i < vParamInfoTmp.size(); ++i)
			{
				tokenize(vParamInfoTmp[i], m_vvParamInfo[i], '|');
			}
		};

	///	name of function
		const std::string& name() const 							{return m_name;}

	///	name of function
		const std::string& options() const 							{return m_methodOptions;}

	/// name of return value
		const std::string& return_name() const 						{return return_info(0);}

	///	type info of return type
		const std::string& return_info(size_t i) const				{return m_vRetValInfo.at(i);}

	/// type info of return value
		const std::vector<std::string>& return_info_vec() const {return m_vRetValInfo;}

	/// number of parameters.
		size_t num_parameter() const 								{return m_vvParamInfo.size();}

	///	number of info strings for one parameter
		size_t num_infos(size_t i) const 							{return m_vvParamInfo.at(i).size();}

	/// name of parameter i
		const std::string& parameter_name(size_t i) const 			{return parameter_info(i, 0);}

	///	type info of all parameters
		const std::string& parameter_info(size_t i, size_t j) const	{return m_vvParamInfo.at(i).at(j);}

	/// type info of i th parameters
		const std::vector<std::string>& parameter_info_vec(size_t i) const {return m_vvParamInfo.at(i);}

	///	whole string of all type infos for of all parameters
		const std::string& parameter_info_string() const			{return m_paramInfos;}

	/// gives some information to the exported functions
		const std::string& tooltip() const 							{return m_tooltip;}

	/// help informations
		const std::string& help() const 							{return m_help;}

	/// parameter list for input values
		const ParameterStack& params_in() const						{return m_paramsIn;}

	/// parameter list for input values
		const ParameterStack& params_out() const					{return m_paramsOut;}

		// todo: we export non-const here, since we can not make ExportedClass_<TClass> a friend
	/// non-const export of param list
		ParameterStack& params_in() 								{return m_paramsIn;}

		// returns false if parameters of the function are undeclared (foreward-declared) classes
		bool check_consistency(const char *classname=NULL) const
		{
			int function_found = 0;
			for(int j=0; j<params_in().size(); j++)
			{
				if(params_in().is_parameter_undeclared(j))
				{

					if(function_found == 0)
					{
						function_found++;
						UG_LOG("Function ");
						PrintFunctionInfo(*this, false, classname);
						UG_LOG(": parameter " << j);
					}
					else
					{	UG_LOG(", " << j);	}
				}
			}
			for(int j=0; j<params_out().size(); j++)
			{
				if(params_out().is_parameter_undeclared(j))
				{

					if(function_found == 0)
					{
						function_found++;
						UG_LOG("Function ");
						PrintFunctionInfo(*this, false, classname);
						UG_LOG(": return value " << j);
					}
					else
					{	UG_LOG(", return value " << j);	}
				}
			}

			if(function_found)
			{
				UG_LOG(": undeclared class.\n");
				return true;
			}
			else return false;
		}

	protected:
		template <typename TFunc>
		void create_parameter_stack()
		{
			typedef typename func_traits<TFunc>::params_type params_type;
			CreateParameterStack<params_type>::create(m_paramsIn);

		//	arbitrary choosen minimum number of infos exported
		//	(If values non given we set them to an empty string)
			const size_t MinNumInfos = 3; // for "name | style | options"

		//	Fill missing Parameter
			m_vvParamInfo.resize(m_paramsIn.size());

		//	resize missing infos for each parameter
			for(int i = 0; i < (int)m_vvParamInfo.size(); ++i)
				for(size_t j = m_vvParamInfo.at(i).size(); j < MinNumInfos; ++j)
					m_vvParamInfo.at(i).push_back(std::string(""));

			typedef typename func_traits<TFunc>::return_type return_type;
			CreateParameterOutStack<return_type>::create(m_paramsOut);

		//	resize missing infos for return value
			for(size_t j = m_vRetValInfo.size(); j < MinNumInfos; ++j)
				m_vRetValInfo.push_back(std::string(""));
		}

		// help function to tokenize the parameter string
		void tokenize(const std::string& str, std::vector<std::string>& tokens, const char delimiter)
		{
			tokens.clear();
			std::stringstream tokenstream;
			tokenstream << str;
			std::string token;

			while ( std::getline (tokenstream, token, delimiter ) )
			{
				tokens.push_back(trim(token));
			}
		}

		std::string trim(const std::string& str)
		{
			const size_t start = str.find_first_not_of(" \t");
			const size_t end = str.find_last_not_of(" \t");
			if(start == std::string::npos || end == std::string::npos) return "";
			return str.substr(start, end - start + 1);
		}

	protected:
		std::string m_funcInfos;
		std::string m_name;
		std::string m_methodOptions;

		std::string m_retValInfos; // string with Infos about return type
		std::vector<std::string> m_vRetValInfo; // tokenized Infos

		// string with Infos about parameter
		std::string m_paramInfos;

		// tokenized strings for each Parameter and each Info (name |�style | options | ...)
		std::vector<std::vector<std::string> > m_vvParamInfo;

		std::string m_tooltip;
		std::string m_help;

		ParameterStack m_paramsIn;
		ParameterStack m_paramsOut;
};

/** function exported from ug
 * This class describes a wrapper for a c++ - function, that is exported by ug
 */
class ExportedFunction : public ExportedFunctionBase
{
	// all c++ functions are wrapped by a proxy function of the following type
	typedef void (*ProxyFunc)(void* func, const ParameterStack& in, ParameterStack& out);

	public:
		template <typename TFunc>
		ExportedFunction(	TFunc f, ProxyFunc pf,
							const char* name, const char* group,
							const char* retValInfos, const char* paramInfos,
							const char* tooltip, const char* help)
			: ExportedFunctionBase( name , retValInfos, paramInfos, tooltip, help),
			  m_group(group), m_func((void*)f), m_proxy_func(pf)
		{
			create_parameter_stack<TFunc>();
		}

	/// executes the function
		void execute(const ParameterStack& paramsIn, ParameterStack& paramsOut) const
		{
			m_proxy_func(m_func, paramsIn, paramsOut);
		}

	///	return groups
		const std::string& group() const {return m_group;}

	protected:
		// save groups
		std::string m_group;

		// pointer to to function
		void* m_func;

		// proxy function
		ProxyFunc m_proxy_func;
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
		typedef typename func_traits<TFunc>::return_type return_type;
		return_type res = func_traits<TFunc>::apply(fp, args);

	//  write result
		//PushTypeValueToParameterStack(res, out);
		PLStack<return_type>::push(out);
		PLStack<return_type>::write(out, res, -1);
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
