
#ifndef __H__UG_BRIDGE__GLOBAL_FUNCTION__
#define __H__UG_BRIDGE__GLOBAL_FUNCTION__

#include <string>
#include <vector>
#include "parameter_stack.h"
#include "function_traits.h"
#include "param_to_type_value_list.h"

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


/** Base class for function/method export
 */
class ExportedFunctionBase
{
	public:
		ExportedFunctionBase(	const char* name, const char* retValName, const char* paramValNames,
								const char* tooltip, const char* help)
		: m_name(name), m_retValName(retValName), m_paramValNames(paramValNames),
		  m_tooltip(tooltip), m_help(help)
		{
			tokenize(m_paramValNames, m_vParamValNames, ",");
		};

	///	name of function
		const std::string& name() const 							{return m_name;}

	/// name of return value
		const std::string& return_name() const 						{return m_retValName;}

	/// string of all return parameters
		const std::string& parameter_names_string() const 			{return m_paramValNames;}

	/// number of parameters.
		size_t num_parameter() const 								{return m_vParamValNames.size();}

	/// name of parameter i
		const std::string& parameter_name(size_t i) const 			{return m_vParamValNames.at(i);}

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

	protected:
		template <typename TFunc>
		void create_parameter_stack()
		{
			typedef typename func_traits<TFunc>::params_type params_type;
			CreateParameterStack<params_type>::create(m_paramsIn);
			for(int i = (int)m_vParamValNames.size(); i < m_paramsIn.size(); ++i)
				m_vParamValNames.push_back(std::string(""));

			typedef typename func_traits<TFunc>::return_type return_type;
			CreateParameterOutStack<return_type>::create(m_paramsOut);
		}
		
		// help function to tokenize the parameter string
		void tokenize(const std::string& str, std::vector<std::string>& tokens, const std::string& delimiters)
		{
			using namespace std;
			tokens.clear();
		    // Skip delimiters at beginning.
		    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
		    // Find first "non-delimiter".
		    string::size_type pos     = str.find_first_of(delimiters, lastPos);

		    while (string::npos != pos || string::npos != lastPos)
		    {
		        // Found a token, add it to the vector.
		        tokens.push_back(str.substr(lastPos, pos - lastPos));
		        // Skip delimiters.  Note the "not_of"
			    lastPos = str.find_first_not_of(delimiters, pos);
			    // Find next "non-delimiter"
			    pos = str.find_first_of(delimiters, lastPos);
		    }
		}

	protected:
		std::string m_name;
		std::string m_retValName;
		std::string m_paramValNames;
		std::vector<std::string> m_vParamValNames;
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
							const char* name, const char* retValName, const char* paramValNames,
							const char* tooltip, const char* help)
			: ExportedFunctionBase( name , retValName, paramValNames, tooltip, help),
			  m_func((void*)f), m_proxy_func(pf)
		{
			create_parameter_stack<TFunc>();
		}

	/// executes the function
		void execute(const ParameterStack& paramsIn, ParameterStack& paramsOut) const
		{
			m_proxy_func(m_func, paramsIn, paramsOut);
		}

	protected:
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
