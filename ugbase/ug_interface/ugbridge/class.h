
#ifndef __H__UG_INTERFACE__UGBRIDGE__CLASS__
#define __H__UG_INTERFACE__UGBRIDGE__CLASS__

#include "function_traits.h"
#include "global_function.h"

namespace ug {

namespace interface{

// predeclaration
template <typename TClass>
class ExportedClass_;

/** function exported from ug
 * This class describes a wrapper for a c++ - function, that is exported by ug
 */
class ExportedMethod : public ExportedFunctionBase
{
	// make Registry a friend
	friend class InterfaceRegistry;

	// all c++ functions are wrapped by a proxy function of the following type
	typedef void (*ProxyFunc)(void* func, IExportedInstance& obj, const ParameterList& in, const ParameterList& out);

	public:
		ExportedMethod(	void* f, ProxyFunc pf,
						const char* name, const char* retValName, const char* paramValNames,
						const char* tooltip, const char* help)
		: ExportedFunctionBase(f, name , retValName, paramValNames, tooltip, help),
		  m_proxy_func(pf)
		{}

	/// executes the function
		void execute(IExportedInstance& obj) const
		{
			m_proxy_func(m_func, obj, m_paramsIn, m_paramsOut);
		}

	private:
		ProxyFunc m_proxy_func;
};

template <typename TClass, typename TMethod>
struct ProxyMethod
{
	static void apply(void* method, IExportedInstance& obj, const ParameterList& in, const ParameterList& out)
	{
	//  cast to method pointer
		TMethod* mptr = (TMethod) method;

	//  cast object to type
		TClass* objPtr = (TClass*) (obj.get_raw_instance());

	//  get parameter
		typedef typename func_traits<TMethod>::params_type params_type;
		ParamToTypeValueList<params_type> args(in);

	//  apply method
		func_traits<TMethod>::apply(mptr, objPtr, args);
	}
};


/** Base class for exported Classes
 *
 */
class IExportedClass
{
	public:
	//  name of class
		virtual const std::string& name() const = 0;

	//  number of method of the class
		virtual size_t num_methods() const = 0;

	//  get exported method
		virtual const ExportedMethod& get_method(size_t i) const = 0;

	//  create an instance
		virtual void* create() = 0;

	//  virtual destructor
		virtual ~IExportedClass() {};
};

template <typename TClass>
class ExportedClass_ : public IExportedClass
{
	private:
	//  disallow creation
		ExportedClass_ () {};
		ExportedClass_ (const ExportedClass_& other);

	//  singleton provider
		static ExportedClass_<TClass>& inst()
		{
			static ExportedClass_<TClass> inst;
			return inst;
		}

	public:
	//  get already created instance
		static ExportedClass_<TClass>& get_inst() {return inst();}

	//  get instance and set name (can only be called with equal name)
		static ExportedClass_<TClass>& get_inst(const char* name)
		{
			// todo: Error handling
			if(std::string(name) == "")
				UG_ASSERT(0, "You must specify a name for a class.");
			if(std::string(name) != m_name)
				UG_ASSERT(0, "Registering a name twice.");

			m_name = std::string(name);
			return inst();
		}

	//  name of class
		virtual const std::string& name() const {return m_name;}

	//  number of registered methods
		virtual size_t num_methods() const {return m_vMethod.size();}

	//  get exported method
		virtual const ExportedMethod& get_method(size_t i) const {return &m_vMethod.at(i);}

	//  Method registration
		template <typename TMethod>
		ExportedClass_<TClass>& method (	const char* methodName, TMethod func,
											const char* retValName = "", const char* paramValNames = "",
											const char* tooltip = "", const char* help = "")
		{
			//  create new exported function
				m_vMethod.push_back(new ExportedMethod(	(void*) func, &ProxyMethod<TClass, TMethod>::apply,
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
		virtual void* create()
		{
			new TClass();
		}
		virtual ~ExportedClass_()
		{
		
		//  delete methods
			for(size_t i = 0; i < m_vMethod.size(); ++i)
			{
				delete m_vMethod[i];
			}
		}

	private:
		static std::string m_name;

		std::vector<ExportedMethod*> m_vMethod;
};

// Initialization of name to unknown type
template <typename TClass>
std::string ExportedClass_<TClass>::m_name = std::string("");

} // end namespace interface

} // end namespace ug


#endif /* __H__UG_INTERFACE__UGBRIDGE__CLASS__ */
