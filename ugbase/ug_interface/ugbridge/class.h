
#ifndef __H__UG_INTERFACE__UGBRIDGE__CLASS__
#define __H__UG_INTERFACE__UGBRIDGE__CLASS__

#include "function_traits.h"
#include "global_function.h"

namespace ug {

namespace interface{

// dummy for unknown class name (used at initialization)
const char* UnknownClassName;


template <typename TClass>
struct ClassNameProvider
{
	//  name of class of type 'TClass'
		static const char* name () 					{return m_name;}
	//  set name of class of type 'TClass'
		static void set_name (const char* name)		{m_name = name;}
	//  return if class is const
		static bool is_const ()						{return false;}

	private:
		static const char *m_name;
};

// Initialization of name to unknown type
template <typename TClass>
const char* ClassNameProvider<TClass>::m_name = UnknownClassName;

// Specialization for const classes
template <typename TClass>
struct ClassNameProvider<const TClass>
{
	//  name of class of type 'TClass'
		static const char* name () 					{return ClassNameProvider<TClass>::name();}
	//  set name of class of type 'TClass'
		static void set_name (const char* name)		{ClassNameProvider<TClass>::set_name(name);}
	//  return if class is const
		static bool is_const ()						{return true;}
};

template <typename TClass, typename TMethod>
struct ProxyMethod
{
	static void apply(void* method, TClass* obj, const ParameterList& in, const ParameterList& out)
	{
	//  cast to method pointer
		TMethod* mp = (TMethod) method;

	//  get parameter
		typedef typename func_traits<TMethod>::params_type params_type;
		ParamToTypeValueList<params_type> args(in);

	//  apply method
		func_traits<TMethod>::apply(mp, obj, args);
	}
};



/** function exported from ug
 * This class describes a wrapper for a c++ - function, that is exported by ug
 */
template <typename TClass>
class ExportedMethod : public ExportedFunctionBase
{
	// make Registry a friend
	friend class InterfaceRegistry;

	// all c++ functions are wrapped by a proxy function of the following type
	typedef void (*ProxyFunc)(void* func, TClass* obj, const ParameterList& in, const ParameterList& out);

	public:
		ExportedMethod(	void* f, ProxyFunc pf,
						const char* name, const char* retValName, const char* paramValNames,
						const char* tooltip, const char* help)
		: ExportedFunctionBase(f, name , retValName, paramValNames, tooltip, help),
		  m_proxy_func(pf)
		{}

	/// executes the function
		void execute(TClass* obj) const
		{
			m_proxy_func(m_func, obj, m_paramsIn, m_paramsOut);
		}

	private:
		ProxyFunc m_proxy_func;
};

class IExportedInstance {};

/** Memory for Instances
 *
 */
template <typename TClass>
class ExportedInstance : public IExportedInstance
{
	public:
		ExportedInstance() {}

	protected:
		TClass m_instance;
};


/** Base class for exported Classes
 *
 */
class IExportedClass
{
	public:
		virtual IExportedInstance* create() = 0;

		virtual ~IExportedClass() {};
};

template <typename TClass>
class ExportedClass_ : public IExportedClass
{
	public:
		ExportedClass_ ();
		ExportedClass_ (const char *name)
		{
		//  remember class name for this type
			ClassNameProvider<TClass>::set_name(name);
		}

		// constructor registration
		template <typename TFunc>
		ExportedClass_<TClass>& constructor ();

		// Method registration
		template <typename TMethod>
		ExportedClass_<TClass>& method (	const char* methodName, TMethod func,
											const char* retValName = "", const char* paramValNames = "",
											const char* tooltip = "", const char* help = "")
		{
			//  create new exported function
				m_vMethod.push_back(new ExportedMethod<TClass>(	(void*) func, &ProxyMethod<TClass, TMethod>::apply,
																methodName, retValName, paramValNames,
																tooltip, help));

			//  create parameter in list
				ParameterList& in = m_vMethod.back()->params_in();
				typedef typename func_traits<TMethod>::params_type params_type;
				CreateParameterList<params_type>::create(in, paramValNames, ",");

			//  create parameter out list
				ParameterList& out = m_vMethod.back()->params_out();
				typedef typename func_traits<TMethod>::result_type result_type;
				CreateParameterList<TypeList<result_type> >::create(out, retValName, ",");

				return *this;
		}

		////////////////////////
		// memory management
		////////////////////////
		virtual IExportedInstance* create()
		{
			m_vInstance.push_back(new ExportedInstance<TClass>());
			return m_vInstance.back();
		}

		virtual ~ExportedClass_()
		{
		//  delete instances
			for(size_t i = 0; i < m_vInstance.size(); ++i)
			{
				delete m_vInstance[i];
			}
		}

	private:
		std::vector<ExportedMethod<TClass>*> m_vMethod;

		std::vector<IExportedInstance*> m_vInstance;
};

} // end namespace interface

} // end namespace ug


#endif /* __H__UG_INTERFACE__UGBRIDGE__CLASS__ */
