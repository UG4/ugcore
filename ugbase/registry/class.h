
#ifndef __H__UG_BRIDGE__CLASS__
#define __H__UG_BRIDGE__CLASS__

#include <cstdlib>
#include <cstring>
#include <string>
#include <boost/type_traits.hpp>

#include "parameter_stack.h"
#include "function_traits.h"
#include "global_function.h"
#include "common/common.h"
#include "error.h"
#include "registry_util.h"

#ifdef PROFILE_BRIDGE
#ifndef UG_PROFILER
	#error "You need to define UG_PROFILER to use PROFILE_BRIDGE"
#endif
#include "common/profiler/runtime_profile_info.h"
#else
#endif


namespace ug
{
namespace bridge
{

/// \addtogroup registry
/// \{

class MethodPtrWrapper
{
	public:
		template <typename TMethod>
		MethodPtrWrapper(TMethod m)
		{
			size = sizeof(TMethod);
			data = malloc(size);
			memcpy(data, &m, size);
		}
		
		MethodPtrWrapper(const MethodPtrWrapper& mpw)
		{
			size = mpw.size;
			data = malloc(size);
			memcpy(data, mpw.get_raw_ptr(), size);
		}
		
		~MethodPtrWrapper()	{free(data);}
		
		void* get_raw_ptr() const {return data;}
		
	protected:
		void*	data;
		int		size;
};

///	Performs a reinterpret cast on the given pointer, then calls delete on it
template <class TClass> void CastAndDelete(const void* ptr)
{
	delete reinterpret_cast<const TClass*>(ptr);
}


/** function exported from ug
 * This class describes a wrapper for a c++ - function, that is exported by ug
 */
class ExportedMethod : public ExportedFunctionBase
{
	public:
	// all c++ functions are wrapped by a proxy function of the following type
		typedef void (*ProxyFunc)(const MethodPtrWrapper& func, void* obj,
								  const ParameterStack& in, ParameterStack& out);

	public:
		ExportedMethod(	const MethodPtrWrapper& m, ProxyFunc pf,
						const std::string& name, const std::string& className,
						const std::string& methodOptions,
						const std::string& retValInfos, const std::string& paramInfos,
						const std::string& tooltip, const std::string& help)
		: ExportedFunctionBase(name, methodOptions, retValInfos, paramInfos, tooltip, help),
		  m_ptrWrapper(m), m_proxy_func(pf), m_className(className)
#ifdef PROFILE_BRIDGE
		,m_dpi((m_className + ":" + name).c_str(), true, "registry", false)
#endif
		{
		}

	/// executes the function
		void execute(void* obj, const ParameterStack& paramsIn, ParameterStack& paramsOut) const
		{
#ifdef PROFILE_BRIDGE
			m_dpi.beginNode();
#endif
			m_proxy_func(m_ptrWrapper, obj, paramsIn, paramsOut);

#ifdef PROFILE_BRIDGE
			m_dpi.endNode();
#endif
		}

	/// \todo: replace this method with a better integrated way.
		template <typename TFunc>
		void create_parameter_stack()
		{
			ExportedFunctionBase::create_parameter_stack<TFunc>();
		}
		
	///	returns the class name this method belongs to
		const std::string& class_name() const {return m_className;}

	private:
	/// pointer to function (stored in a wrapper)
		MethodPtrWrapper m_ptrWrapper;
	
	/// proxy function to call method
		ProxyFunc m_proxy_func;

	/// name of class this method belongs to
		std::string m_className;

#ifdef PROFILE_BRIDGE
		mutable RuntimeProfileInfo m_dpi;
#endif
};

////////////////////////////////////////////////////////////////////////////////
//	ExportedMethodGroup (sreiter)
////////////////////////////////////////////////////////////////////////////////

///	Groups of methods - useful to realize overloaded methods
class ExportedMethodGroup
{
	typedef ExportedMethod::ProxyFunc ProxyFunc;

	public:
		ExportedMethodGroup(const std::string& name) : m_name(name)
		{}

		~ExportedMethodGroup()
		{
			for(size_t i = 0; i < m_overloads.size(); ++i)
				delete m_overloads[i].m_func;
		}

	///	name of function group
		const std::string& name() const {return m_name;}

	///	adds an overload. Returns false if the overload already existed.
		template <class TFunc>
		bool add_overload(	const TFunc& m, ProxyFunc pf,
		                  	const std::string& className,
							const std::string& methodOptions, const std::string& retValInfos,
							const std::string& paramInfos, const std::string& tooltip,
							const std::string& help)
		{
			size_t typeID = GetUniqueTypeID<TFunc>();

		//	make sure that the overload didn't exist
			if(get_overload_by_type_id(typeID))
				return false;

		//	create a new overload
			ExportedMethod* func = new ExportedMethod(MethodPtrWrapper(m),
												pf, m_name, className,
												methodOptions, retValInfos,
												paramInfos, tooltip, help);

			m_overloads.push_back(Overload(func, typeID));


		//  create parameter in list
			func->create_parameter_stack<TFunc>();

			return true;
		}

		size_t num_overloads() const {return m_overloads.size();}

		ExportedMethod* get_overload(size_t index) {return m_overloads.at(index).m_func;}

		const ExportedMethod* get_overload(size_t index) const {return m_overloads.at(index).m_func;}

		template <class TType>
		ExportedMethod* get_overload_by_type()
		{
			size_t typeID = GetUniqueTypeID<TType>();
			return get_overload_by_type_id(typeID);
		}

		template <class TType>
		const ExportedMethod* get_overload_by_type() const
		{
			size_t typeID = GetUniqueTypeID<TType>();
			return get_overload_by_type_id(typeID);
		}

		ExportedMethod* get_overload_by_type_id(size_t typeID)
		{
			for(size_t i = 0; i < m_overloads.size(); ++i){
				if(m_overloads[i].m_typeID == typeID)
					return m_overloads[i].m_func;
			}
			return NULL;
		}

		const ExportedMethod* get_overload_by_type_id(size_t typeID) const
		{
			for(size_t i = 0; i < m_overloads.size(); ++i){
				if(m_overloads[i].m_typeID == typeID)
					return m_overloads[i].m_func;
			}
			return NULL;
		}

		size_t get_overload_type_id(size_t index) const {return m_overloads.at(index).m_typeID;}

	private:
		struct Overload{
			Overload()	{}
			Overload(ExportedMethod* func, size_t typeID)
				: m_func(func), m_typeID(typeID)
			{}
			ExportedMethod* 	m_func;
			size_t				m_typeID;
		};

		std::string m_name;
		std::vector<Overload>	m_overloads;
};

template <typename TClass, typename TMethod,
		  typename TRet = typename func_traits<TMethod>::return_type>
struct MethodProxy
{
	static void apply(const MethodPtrWrapper& method, void* obj,
	                  const ParameterStack& in, ParameterStack& out)
	{
	//  cast to method pointer
		TMethod mptr = *(TMethod*) method.get_raw_ptr();

	//  cast object to type
		TClass* objPtr = (TClass*) (obj);

	//  get parameter
		typedef typename func_traits<TMethod>::params_type params_type;
		ParameterStackToTypeValueList<params_type> args(in);

	//  apply method
		TRet res = func_traits<TMethod>::apply(mptr, objPtr, args);

	//  write return value
		out.push<TRet>(res);
	}
};

template <typename TClass, typename TMethod>
struct MethodProxy<TClass, TMethod, void>
{
	static void apply(const MethodPtrWrapper& method, void* obj,
	                  const ParameterStack& in, ParameterStack& out)
	{
	//  cast to method pointer
		TMethod mptr = *(TMethod*) method.get_raw_ptr();

	//  cast object to type
		TClass* objPtr = (TClass*) (obj);

	//  get parameter
		typedef typename func_traits<TMethod>::params_type params_type;
		ParameterStackToTypeValueList<params_type> args(in);

	//  apply method
		func_traits<TMethod>::apply(mptr, objPtr, args);
	}
};

////////////////////////////////////////////////////////////////////////////////
// Exported Constructor
////////////////////////////////////////////////////////////////////////////////

/// describing information for constructor
class UG_API ExportedConstructor
{
	public:
	// all c++ functions are wrapped by a proxy function of the following type
		typedef void* (*ProxyFunc)(const ParameterStack& in);

	public:
		ExportedConstructor(ProxyFunc pf,
		                    const std::string& className, const std::string& options,
		                    const std::string& paramInfos,
		                    const std::string& tooltip, const std::string& help);

	/// executes the function
		void* create(const ParameterStack& paramsIn) const
		{
#ifdef PROFILE_BRIDGE
			m_dpi.beginNode();
#endif

			void *pRet = m_proxy_func(paramsIn);

#ifdef PROFILE_BRIDGE
			m_dpi.endNode();
#endif
			return pRet;
		}

	///	options
		const std::string& options() const {return m_options;}

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

	/// non-const export of param list
		ParameterInfo& params_in() {return m_paramsIn;}

	/// returns true if all parameters of the function are correctly declared
		bool check_consistency(std::string classname) const;

		template <typename TFunc>
		void create_parameter_stack()
		{
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
		}

	protected:
	// 	help function to tokenize the parameter string
		void tokenize(const std::string& str, std::vector<std::string>& tokens,
		              const char delimiter);

	protected:
	private:
	/// proxy function to call method
		ProxyFunc m_proxy_func;

	///	name of class constructed
		std::string m_className;

	///	options
		std::string m_options;

	/// string with Infos about parameter
		std::string m_paramInfos;

	///	tokenized strings for each Parameter and each Info (name | style | options | ...)
		std::vector<std::vector<std::string> > m_vvParamInfo;

		std::string m_tooltip;
		std::string m_help;

		ParameterInfo m_paramsIn;

#ifdef PROFILE_BRIDGE
		mutable RuntimeProfileInfo m_dpi;
#endif
};

template <typename TClass, typename TMethod>
struct ConstructorProxy
{
	static void* create(const ParameterStack& in)
	{
	//  get parameter
		typedef typename func_traits<TMethod>::params_type params_type;
		ParameterStackToTypeValueList<params_type> args(in);

	//  apply method
		TClass* newInst = constructor_traits<TClass, params_type>::apply(args);

	//  return new pointer
		return (void*) newInst;
	}
};

template <typename TClass>
void DestructorProxy(void* obj)
{
	TClass* pObj = (TClass*)obj;
	delete pObj;
}

////////////////////////////////////////////////////////////////////////////////
// Interface Exported Class
////////////////////////////////////////////////////////////////////////////////

/// Base class for exported Classes
class IExportedClass
{
	public:
		typedef void (*DeleteFunction)(const void*);

	public:
	///  name of class
		virtual const std::string& name() const = 0;

	///	get groups
		virtual const std::string& group() const = 0;

	///	name node of class
		virtual const ClassNameNode& class_name_node() const = 0;

	/// get tooltip
		virtual const std::string& tooltip() const = 0;

	///	name-list of class hierarchy
		virtual const std::vector<const char*>* class_names() const = 0;

	///  number of method of the class
		virtual size_t num_methods() const = 0;

	///	number of registered const-methods
		virtual size_t num_const_methods() const = 0;
		
	///  get exported method
		virtual const ExportedMethod& get_method(size_t i) const = 0;

	/// get exported const-method
		virtual const ExportedMethod& get_const_method(size_t i) const = 0;
		
	///	returns the number of overloads of a method
		virtual size_t num_overloads(size_t funcInd) const = 0;

	///	returns the number of overloads of a const method
		virtual size_t num_const_overloads(size_t funcInd) const = 0;

	///	returns the i-th overload of a method
		virtual const ExportedMethod& get_overload(size_t funcInd, size_t oInd) const = 0;

	///	returns the i-th overload of a const method
		virtual const ExportedMethod& get_const_overload(size_t funcInd, size_t oInd) const = 0;

	///	returns the i-th method group (all overloads of the i-th function)
		virtual const ExportedMethodGroup& get_method_group(size_t ind) const = 0;

	///	returns the i-th method group (all overloads of the i-th function)
		virtual const ExportedMethodGroup& get_const_method_group(size_t ind) const = 0;

		virtual const ExportedMethodGroup* get_exported_method_group(const std::string& name) const = 0;

		virtual const ExportedMethodGroup* get_const_exported_method_group(const std::string& name) const = 0;

	/**  can we create instances of this class
	 *	(i.e. the class does not contain pure virtual functions)*/
		virtual bool is_instantiable() const = 0;

	///	number of registered constructors
		virtual size_t num_constructors() const = 0;

	///	get exported constructor
		virtual const ExportedConstructor& get_constructor(size_t i) const = 0;

	///	true if the class shall be wrapped in a SmartPtr on construction
		virtual bool construct_as_smart_pointer() const = 0;

	///	destructur for object
		virtual void destroy(void* obj) const = 0;

	///	returns a function which will call delete on the object
		virtual DeleteFunction get_delete_function() const = 0;

	///	returns false is consistency-check failed
		virtual bool check_consistency() const;

	///  virtual destructor
		virtual ~IExportedClass() {};
};

///	A base implementation with non-template methods.
/**	Speeds up compilation times.*/
class ExportedClassBaseImpl : public IExportedClass
{
	private:
	//  disallow
		ExportedClassBaseImpl();
		ExportedClassBaseImpl(const ExportedClassBaseImpl& other);

	public:
		ExportedClassBaseImpl(const std::string& tooltip);

	/// destructor
		virtual ~ExportedClassBaseImpl();

	/// tooltip
		virtual const std::string& tooltip() const;

	/// number of registered methods (overloads are not counted)
		virtual size_t num_methods() const;

	///	number of registered const-methods (overloads are not counted)
		virtual size_t num_const_methods() const;

	/// returns the first overload of an exported function
		virtual const ExportedMethod& get_method(size_t i) const;

	/// returns the first overload of an exported const function
		virtual const ExportedMethod& get_const_method(size_t i) const;

	///	returns the number of overloads of a method
		virtual size_t num_overloads(size_t funcInd) const;

	///	returns the number of overloads of a const method
		virtual size_t num_const_overloads(size_t funcInd) const;

	///	returns the i-th overload of a method
		virtual const ExportedMethod& get_overload(size_t funcInd, size_t oInd) const;

	///	returns the i-th overload of a const method
		virtual const ExportedMethod& get_const_overload(size_t funcInd, size_t oInd) const;

	///	returns the i-th method group (all overloads of the i-th function)
		virtual const ExportedMethodGroup& get_method_group(size_t ind) const;

	///	returns the i-th method group (all overloads of the i-th function)
		virtual const ExportedMethodGroup& get_const_method_group(size_t ind) const;

		virtual const ExportedMethodGroup* get_exported_method_group(const std::string& name) const;

		virtual const ExportedMethodGroup* get_const_exported_method_group(const std::string& name) const;

	///	number of registered constructors
		virtual size_t num_constructors() const;

	///	get exported constructor
		virtual const ExportedConstructor& get_constructor(size_t i) const;

	///	returns whether the class shall be wrapped in a SmartPtr on construction
		virtual bool construct_as_smart_pointer() const;

	///	sets whether the class shall be wrapped in a SmartPtr
	/**	Returns a reference to this, so that it can be used in a chained call.*/
		virtual void set_construct_as_smart_pointer(bool enable);

	/// is instantiable
		virtual bool is_instantiable() const;

	///	destructur for object
		virtual void destroy(void* obj) const;

	protected:
	///	returns if a constructor overload is registered
		bool constructor_type_id_registered(size_t typeID);

	/// returns true if methodname is already used by a method in this class
		bool constmethodname_registered(const std::string& name);

	/// returns true if methodname is already used by a method in this class
		bool methodname_registered(const std::string& name);

		ExportedMethodGroup* get_exported_method_group(const std::string& name);

		ExportedMethodGroup* get_const_exported_method_group(const std::string& name);

	protected:
		struct ConstructorOverload{
			ConstructorOverload()	{}
			ConstructorOverload(ExportedConstructor* func, size_t typeID)
				: m_constructor(func), m_typeID(typeID)
			{}
			ExportedConstructor* 	m_constructor;
			size_t					m_typeID;
		};

		std::vector<ConstructorOverload> m_vConstructor;

		typedef void (*DestructorFunc)(void*);
		DestructorFunc m_destructor;

		std::vector<ExportedMethodGroup*> m_vMethod;
		std::vector<ExportedMethodGroup*> m_vConstMethod;
		std::string m_tooltip;

		bool	m_constructAsSmartPtr;
};


///	This template class represents real c++ classes in the registry.
template <typename TClass>
class ExportedClass : public ExportedClassBaseImpl
{
	private:
	//  disallow
		ExportedClass () {};
		ExportedClass (const ExportedClass& other);

	public:
	//  contructor
		ExportedClass(const std::string& name, const std::string& group,
		               const std::string& tooltip)
				: ExportedClassBaseImpl(tooltip)
		{
			ClassNameProvider<TClass>::set_name(name, group, true);
			ClassNameProvider<const TClass>::set_name(name, group, true);
		}

	/// destructor
		virtual ~ExportedClass(){}

	/// name of class
		virtual const std::string& name() const {return ClassNameProvider<TClass>::name();}

	///	name node of class
		virtual const ClassNameNode& class_name_node() const {return ClassNameProvider<TClass>::class_name_node();}

	///	get groups
		virtual const std::string& group() const {return ClassNameProvider<TClass>::group();}

	//\todo: remove this method, use class name nodes instead
	///	class-hierarchy
		virtual const std::vector<const char*>* class_names() const	{return &ClassNameProvider<TClass>::names();}

	/// Method registration
	/**
	 * @param methodName the name of the method to appear in the registry
	 * @param func pointer to the member function (e.g. &MyClass::my_method)
	 * @param retValInfos string documenting what the function returns (optional)
	 * @param paramInfos string documenting the parameters of the function
	 * seperate parameters with an # e.g. "x#y#z" (don't specify the type of the values)  (optional)
	 * @param toolTip small documentation for the function (optional)
	 * @param help help string for the function
	 * @sa \ref pageUG4Registry
	 * @sa \ref secSTHowToSpecifyParameterInformation
	 */
		template <typename TMethod>
		ExportedClass<TClass>& add_method (std::string methodName, TMethod func,
		                                    std::string retValInfos = "", std::string paramInfos = "",
		                                    std::string tooltip = "", std::string help = "")
		{
		//	At this point the method name contains parameters (name|param1=...).
		//todo: they should be removed and specified with an extra parameter.

			std::string strippedMethodName = methodName;
			std::string methodOptions;
			std::string::size_type pos = strippedMethodName.find("|");
			if(pos != std::string::npos){
				methodOptions = strippedMethodName.substr(pos + 1, strippedMethodName.length() - pos);
				strippedMethodName = strippedMethodName.substr(0, pos);
			}

		//	trim whitespaces
			strippedMethodName = TrimString(strippedMethodName);
			methodOptions = TrimString(methodOptions);

		// 	check that name is not empty
			if(strippedMethodName.empty())
			{
				UG_THROW_REGISTRY_ERROR(strippedMethodName,
				        "Trying to register empty method name.");
			}
			
			// check that name does not contain illegal characters
			if (!IsValidRegistryIdentifier(strippedMethodName)) {
				UG_THROW_REGISTRY_ERROR(strippedMethodName,
				"Trying to register method '"<< strippedMethodName << "' that"
				<< " contains illegal characters. "
				<< GetRegistryIdentifierMessage());
			}

		//	if the method is already in use, we have to add an overload.
		//	If not, we have to create a new method group
			ExportedMethodGroup* methodGrp = NULL;
			if(func_traits<TMethod>::const_method)
				methodGrp = get_const_exported_method_group(strippedMethodName);
			else
				methodGrp = get_exported_method_group(strippedMethodName);

			if(!methodGrp){
				methodGrp = new ExportedMethodGroup(strippedMethodName);
				if(func_traits<TMethod>::const_method)
					m_vConstMethod.push_back(methodGrp);
				else
					m_vMethod.push_back(methodGrp);
			}

			bool success = methodGrp->add_overload(func, &MethodProxy<TClass, TMethod>::apply,
												   ClassNameProvider<TClass>::name(),
												   methodOptions, retValInfos, paramInfos,
												   tooltip, help);

			if(!success)
			{
				UG_THROW_REGISTRY_ERROR(name(),
				"Trying to register method name '" << strippedMethodName
				<< "' to class '" << name() << "', but another method with this name "
				<< " and the same function signature is already registered for this class.");
			}

			return *this;
		}

	/// Make default constructor accessible
		ExportedClass<TClass>& add_constructor ()
		{
		//	add also in new style
			add_constructor<void (*)()>();

		//  remember constructor proxy
			m_destructor = &DestructorProxy<TClass>;

			return *this;
		}

	/// constructor registration
	/**
	 * @param paramInfos string documenting the parameters of the function
	 * seperate parameters with an # e.g. "x#y#z" (don't specify the type of the values)  (optional)
	 * @param toolTip small documentation for the function (optional)
	 * @param help help string for the function
	 * @param options style option of the constructor itself (for visualisation)
	 * @sa \ref pageUG4Registry
	 */
		template <typename TFunc>
		ExportedClass<TClass>& add_constructor(std::string paramInfos = "",
		                                       std::string tooltip = "", std::string help = "",
		                                       std::string options = "")
		{
		//	return-type must be void
			if(!(boost::is_void< typename func_traits<TFunc>::return_type >::value))
			{
				UG_THROW_REGISTRY_ERROR(name(),
				"Trying to register constructor of class "
				 <<name()<<"with non-void return value in signature function.");
			}

		//	type id of constructor
			size_t typeID = GetUniqueTypeID<TFunc>();

		//	make sure that the overload didn't exist
			if(constructor_type_id_registered(typeID))
			{
				UG_THROW_REGISTRY_ERROR(name(),
				"Trying to register constructor of class "
				 <<name()<<" with same signature twice.");
			}

		//	create new exported constructor
			ExportedConstructor* expConstr
				= new ExportedConstructor(	&ConstructorProxy<TClass, TFunc>::create,
				                          	ClassNameProvider<TClass>::name(),
											options, paramInfos, tooltip, help);

		//	create parameter stack
			expConstr->create_parameter_stack<TFunc>();

		//	rememeber it
			m_vConstructor.push_back(ConstructorOverload(expConstr, typeID));

		//	done
			return *this;
		}

	///	return pointer to the delete method
		virtual DeleteFunction get_delete_function() const
		{
			return CastAndDelete<TClass>;
		}
};

// end group registry
/// \}

} // end namespace bridge
} // end namespace ug


#endif /* __H__UG_BRIDGE__CLASS__ */
