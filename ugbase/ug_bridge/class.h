
#ifndef __H__UG_BRIDGE__CLASS__
#define __H__UG_BRIDGE__CLASS__

#include <cstdlib>
#include <cstring>
#include "parameter_stack.h"
#include "function_traits.h"
#include "global_function.h"
#include "common/common.h"

namespace ug
{
namespace bridge
{

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

/** function exported from ug
 * This class describes a wrapper for a c++ - function, that is exported by ug
 */
class ExportedMethod : public ExportedFunctionBase
{
	// all c++ functions are wrapped by a proxy function of the following type
	typedef void (*ProxyFunc)(const MethodPtrWrapper& func, void* obj, const ParameterStack& in, ParameterStack& out);

	public:
		ExportedMethod(	const MethodPtrWrapper& m, ProxyFunc pf,
						const char* name, const char* retValName, const char* paramValNames,
						const char* tooltip, const char* help,
						const char* retValInfoType, const char* paramValInfoType)
		: ExportedFunctionBase(name , retValName, paramValNames, tooltip, help, retValInfoType, paramValInfoType),
		  m_ptrWrapper(m), m_proxy_func(pf)
		{}

	/// executes the function
		void execute(void* obj, const ParameterStack& paramsIn, ParameterStack& paramsOut) const
		{
			m_proxy_func(m_ptrWrapper, obj, paramsIn, paramsOut);
		}

	/// \todo: replace this method with a better integrated way.
		template <typename TFunc>
		void create_parameter_stack()
		{
			ExportedFunctionBase::create_parameter_stack<TFunc>();
		}
		
	private:
		// pointer to function (stored in a wrapper)
		MethodPtrWrapper m_ptrWrapper;
	
		// proxy function to call method
		ProxyFunc m_proxy_func;
};

template <typename TClass, typename TMethod, typename TRet = typename func_traits<TMethod>::return_type>
struct MethodProxy
{
	static void apply(const MethodPtrWrapper& method, void* obj, const ParameterStack& in, ParameterStack& out)
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
		//PushTypeValueToParameterStack(res, out);
		PLStack<TRet>::push(out);
		PLStack<TRet>::write(out, res, -1);
	}
};

template <typename TClass, typename TMethod>
struct MethodProxy<TClass, TMethod, void>
{
	static void apply(const MethodPtrWrapper& method, void* obj, const ParameterStack& in, ParameterStack& out)
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

template <typename TClass>
TClass* ConstructorProxy() {return new TClass();}

/** Base class for exported Classes
 *
 */
class IExportedClass
{
	public:
	///  name of class
		virtual const char* name() const = 0;

	///	get groups
		virtual const std::string& group() const = 0;

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
		
	/**  can we create instances of this class
	 *	(i.e. the class does not contain pure virtual functions)*/
		virtual bool is_instantiable() const = 0;

	/**  create an instance
	 *  returns NULL id we cannot create instances of this type*/
		virtual void* create() const = 0;

	///  virtual destructor
		virtual ~IExportedClass() {};
};


template <typename TClass>
class ExportedClass_ : public IExportedClass
{
	private:
	//  disallow
		ExportedClass_ () {};
		ExportedClass_ (const ExportedClass_& other);

	public:
	//  contructor
		ExportedClass_(const char* name, const char* group = "") : m_constructor(NULL)
		{
			ClassNameProvider<TClass>::set_name(name, group, true);
			ClassNameProvider<const TClass>::set_name(name, group, true);
		}

	/// name of class
		virtual const char* name() const {return ClassNameProvider<TClass>::name();}

	///	get groups
		virtual const std::string& group() const {return ClassNameProvider<TClass>::group();}

	///	class-hierarchy
		virtual const std::vector<const char*>* class_names() const	{return &ClassNameProvider<TClass>::names();}

	/// number of registered methods
		virtual size_t num_methods() const {return m_vMethod.size();}

	///	number of registered const-methods
		virtual size_t num_const_methods() const {return m_vConstMethod.size();}
		
	/// get exported method
		virtual const ExportedMethod& get_method(size_t i) const {return *m_vMethod.at(i);}

	/// get exported const-method
		virtual const ExportedMethod& get_const_method(size_t i) const {return *m_vConstMethod.at(i);}

	/// Method registration
		template <typename TMethod>
		ExportedClass_<TClass>& add_method (	const char* methodName, TMethod func,
												const char* retValName = "", const char* paramValNames = "",
												const char* tooltip = "", const char* help = "",
												const char* retValInfoType = "", const char* paramValInfoType = "")
		{
		//	check that funcName is not already used
			bool bUsed = false;
			if(func_traits<TMethod>::const_method)
				bUsed = constmethodname_registered(methodName);
			else
				bUsed = methodname_registered(methodName);
			if(bUsed)
			{
				std::cout << "### Registry ERROR: Trying to register method name '" << methodName
						<< "' to class '" << name() << "', but another method with this name is already"
						<< " registered for this class."
						<< "\n### Please change register process. Aborting ..." << std::endl;
				exit(1);
			}


		//  create new exported function
			ExportedMethod* nMethod = NULL;
			nMethod = new ExportedMethod(	MethodPtrWrapper(func), &MethodProxy<TClass, TMethod>::apply,
													methodName, retValName, paramValNames,
													tooltip, help,
													retValInfoType, paramValInfoType);
			
			try{
		//  create parameter in list
				nMethod->create_parameter_stack<TMethod>();
			}
			catch(ug::bridge::UG_ERROR_ClassUnknownToRegistry ex)
			{
				UG_LOG("###  Registering method '" << methodName << "' for class '");
				UG_LOG( name() << "' requires argument of user-defined type,\n");
				UG_LOG("###  that has not yet been registered to this Registry.\n");
				exit(1);
			}

			if(func_traits<TMethod>::const_method)
				m_vConstMethod.push_back(nMethod);
			else
				m_vMethod.push_back(nMethod);

			return *this;
		}

	/// Make constructor accessible
	//  We use a pointer to ConstructorProxy since abstract base classes can not be created but registered.
	//  Each class that is instantiable must register its constructor
		ExportedClass_<TClass>& add_constructor ()
		{
		//  remember constructor proxy
			m_constructor = &ConstructorProxy<TClass>;

			return *this;
		}

	/// is instantiable
		virtual bool is_instantiable() const {return m_constructor != NULL;}

	/// create new instance of class
		virtual void* create() const
		{
			if(m_constructor != NULL)
				return (*m_constructor)();
			else
				return NULL;
		}

	/// destructor
		virtual ~ExportedClass_()
		{
		//  delete methods
			for(size_t i = 0; i < m_vMethod.size(); ++i)
				delete m_vMethod[i];

			for(size_t i = 0; i < m_vConstMethod.size(); ++i)
				delete m_vConstMethod[i];
		}

	protected:
		// returns true if methodname is already used by a method in this class
		bool constmethodname_registered(const char* name)
		{
			for(size_t i = 0; i < m_vConstMethod.size(); ++i)
			{
			//  compare strings
				if(strcmp(name, (m_vConstMethod[i]->name()).c_str()) == 0)
					return true;
			}
			return false;
		}
		// returns true if methodname is already used by a method in this class
		bool methodname_registered(const char* name)
		{
			for(size_t i = 0; i < m_vMethod.size(); ++i)
			{
			//  compare strings
				if(strcmp(name, (m_vMethod[i]->name()).c_str()) == 0)
					return true;
			}
			return false;
		}


	private:
		typedef TClass* (*ConstructorFunc)();
		ConstructorFunc m_constructor;

		std::vector<ExportedMethod*> m_vMethod;
		std::vector<ExportedMethod*> m_vConstMethod;
};

} // end namespace bridge
} // end namespace ug


#endif /* __H__UG_BRIDGE__CLASS__ */
