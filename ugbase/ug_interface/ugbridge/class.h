
#ifndef __H__UG_INTERFACE__UGBRIDGE__CLASS__
#define __H__UG_INTERFACE__UGBRIDGE__CLASS__

#include "parameter_stack.h"
#include "function_traits.h"
#include "global_function.h"
#include "common/common.h"

namespace ug {

namespace interface{

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
	typedef void (*ProxyFunc)(MethodPtrWrapper& func, void* obj, const ParameterStack& in, ParameterStack& out);

	public:
		ExportedMethod(	const MethodPtrWrapper& m, ProxyFunc pf,
						const char* name, const char* retValName, const char* paramValNames,
						const char* tooltip, const char* help)
		: ExportedFunctionBase(NULL, name , retValName, paramValNames, tooltip, help),
		  m_ptrWrapper(m), m_proxy_func(pf)
		{}

	/// executes the function
		void execute(void* obj, const ParameterStack& paramsIn, ParameterStack& paramsOut)
		{
			m_proxy_func(m_ptrWrapper, obj, paramsIn, paramsOut);
		}

	private:
		MethodPtrWrapper m_ptrWrapper;
	
		ProxyFunc m_proxy_func;
};

template <typename TClass, typename TMethod>
struct ProxyMethod
{
	static void apply(MethodPtrWrapper& method, void* obj, const ParameterStack& in, ParameterStack& out)
	{
	//  cast to method pointer
		TMethod mptr = *(TMethod*) method.get_raw_ptr();

	//  cast object to type
		TClass* objPtr = (TClass*) (obj);

	//  get parameter
		typedef typename func_traits<TMethod>::params_type params_type;
		typedef typename func_traits<TMethod>::result_type result_type;
		ParameterStackToTypeValueList<params_type> args(in);

	//  apply method
		result_type res = func_traits<TMethod>::apply(mptr, objPtr, args);

	//  write return value
		PushTypeValueToParameterStack(res, out);
	}
};


/** Base class for exported Classes
 *
 */
class IExportedClass
{
	public:
	//  name of class
		virtual const char* name() const = 0;

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
struct ClassNameProvider
{
	void set_name(const char* name) {m_name = name;}
	const char* name() const {return m_name;}

	private:
		const char* m_name;
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
		ExportedClass_(const char* name) : m_name(name) {}

	//  name of class
		virtual const char* name() const {return m_name;}

	//  number of registered methods
		virtual size_t num_methods() const {return m_vMethod.size();}

	//  get exported method
		virtual const ExportedMethod& get_method(size_t i) const {return *m_vMethod.at(i);}

	//  Method registration
		template <typename TMethod>
		ExportedClass_<TClass>& add_method (	const char* methodName, TMethod func,
												const char* retValName = "", const char* paramValNames = "",
												const char* tooltip = "", const char* help = "")
		{
			//  create new exported function
				m_vMethod.push_back(new ExportedMethod(	MethodPtrWrapper(func), &ProxyMethod<TClass, TMethod>::apply,
														methodName, retValName, paramValNames,
														tooltip, help));

			//  create parameter in list
				ParameterStack& in = m_vMethod.back()->params_in();
				typedef typename func_traits<TMethod>::params_type params_type;
				CreateParameterStack<params_type>::create(in);

			//  create parameter out list
				ParameterStack& out = m_vMethod.back()->params_out();
				typedef typename func_traits<TMethod>::result_type result_type;
				CreateParameterStack<TypeList<result_type> >::create(out);

				return *this;
		}

		////////////////////////
		// memory management
		////////////////////////
		virtual void* create()
		{
			return new TClass();
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
		const char* m_name;

		std::vector<ExportedMethod*> m_vMethod;
};

} // end namespace interface

} // end namespace ug


#endif /* __H__UG_INTERFACE__UGBRIDGE__CLASS__ */
