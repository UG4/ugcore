
#ifndef __H__UG_BRIDGE__CLASS__
#define __H__UG_BRIDGE__CLASS__

#include <cstdlib>
#include <cstring>
#include <string>
#include "parameter_stack.h"
#include "function_traits.h"
#include "global_function.h"
#include "common/common.h"

namespace ug
{
namespace bridge
{


struct UG_REGISTRY_ERROR_RegistrationFailed
{
		UG_REGISTRY_ERROR_RegistrationFailed(const char* name_)
			: name(name_)
		{}
		std::string name;
};

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
template <class TClass> void CastAndDelete(void* ptr)
{
	delete reinterpret_cast<TClass*>(ptr);
}

/** function exported from ug
 * This class describes a wrapper for a c++ - function, that is exported by ug
 */
class ExportedMethod : public ExportedFunctionBase
{
	public:
	// all c++ functions are wrapped by a proxy function of the following type
		typedef void (*ProxyFunc)(const MethodPtrWrapper& func, void* obj, const ParameterStack& in, ParameterStack& out);

	public:
		ExportedMethod(	const MethodPtrWrapper& m, ProxyFunc pf,
						const char* name, const char* className,
						const char* methodOptions,
						const char* retValInfos, const char* paramInfos,
						const char* tooltip, const char* help)
		: ExportedFunctionBase(name, methodOptions, retValInfos, paramInfos, tooltip, help),
		  m_ptrWrapper(m), m_proxy_func(pf), m_className(className)
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
		
	///	returns the class name this method belongs to
		const char* class_name() const {return m_className;}

	private:
		// pointer to function (stored in a wrapper)
		MethodPtrWrapper m_ptrWrapper;
	
		// proxy function to call method
		ProxyFunc m_proxy_func;

		// name of class this method belongs to
		const char* m_className;
};

////////////////////////////////////////////////////////////////////////
//	ExportedMethodGroup (sreiter)
///	Groups of methods - useful to realize overloaded methods
class ExportedMethodGroup
{
	typedef ExportedMethod::ProxyFunc ProxyFunc;

	public:
		ExportedMethodGroup(const char* name) : m_name(name)
		{}

		~ExportedMethodGroup()
		{
			for(size_t i = 0; i < m_overloads.size(); ++i)
				delete m_overloads[i].m_func;
		}

	///	name of function group
		const std::string& name() const 							{return m_name;}

	///	adds an overload. Returns false if the overload already existed.
		template <class TFunc>
		bool add_overload(	const TFunc& m, ProxyFunc pf,
		                  	const char* className,
							const char* methodOptions, const char* retValInfos,
							const char* paramInfos, const char* tooltip,
							const char* help)
		{
			size_t typeID = GetUniqueTypeID<TFunc>();

		//	make sure that the overload didn't exist
			if(get_overload_by_type_id(typeID))
				return false;

		//	create a new overload
			ExportedMethod* func = new ExportedMethod(MethodPtrWrapper(m),
												pf, m_name.c_str(), className,
												methodOptions, retValInfos,
												paramInfos, tooltip, help);

			m_overloads.push_back(Overload(func, typeID));


		//  create parameter in list
			func->create_parameter_stack<TFunc>();

			return true;
		}

		size_t num_overloads() const
			{return m_overloads.size();}

		ExportedMethod* get_overload(size_t index)
			{return m_overloads.at(index).m_func;}

		const ExportedMethod* get_overload(size_t index) const
			{return m_overloads.at(index).m_func;}

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

		size_t get_overload_type_id(size_t index) const
			{return m_overloads.at(index).m_typeID;}

	private:
		struct Overload{
			Overload()	{}
			Overload(ExportedMethod* func, size_t typeID) : m_func(func), m_typeID(typeID) {}
			ExportedMethod* 	m_func;
			size_t				m_typeID;
		};

		std::string m_name;
		std::vector<Overload>	m_overloads;
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

template <typename TClass>
void DestructorProxy(void* obj)
{
	TClass* pObj = (TClass*)obj;
	delete pObj;
}

/** Base class for exported Classes
 *
 */
class IExportedClass
{
	public:
		typedef void (*DeleteFunction)(void*);

	public:
	///  name of class
		virtual const char* name() const = 0;

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

	/**  can we create instances of this class
	 *	(i.e. the class does not contain pure virtual functions)*/
		virtual bool is_instantiable() const = 0;

	/**  create an instance
	 *  returns NULL id we cannot create instances of this type*/
		virtual void* create() const = 0;

	///	destructur for object
		virtual void destroy(void* obj) const = 0;

	///	returns a function which will call delete on the object
		virtual DeleteFunction get_delete_function() const = 0;

	///	returns false is consistency-check failed
		virtual bool check_consistency() const
		{
		//	get class name vector of all parents
			const std::vector<const char*>* vClassNames = class_names();

		//	check if class name vector correct
			if(vClassNames==NULL)
			{
				UG_LOG("ERROR in 'IExportedClass::check_consistency':"
						" Class name vector of parent classes missing for "
						"class '"<<this->name()<<"'.\n");
				return false;
			}

		//	loop all base classes
			for(size_t i = 0; i < (*vClassNames).size(); ++i)
			{
			//	get name of base class
				const char* baseName = (*vClassNames)[i];

			//	check the name
				if(baseName == NULL || strlen(baseName) == 0 || baseName[0] == '[')
				{
					if(i>0){
					UG_LOG("ERROR in 'IExportedClass::check_consistency':"
							" base class "<<i<<" of class '"<<this->name()<<
							"' has not been named.\n");
						return false;
					}
					else{
					UG_LOG("ERROR in 'IExportedClass::check_consistency':"
							" Class '"<<this->name()<<"' has not been named.\n");
						return false;
					}
				}
			}

		//	everything ok
			return true;
		}

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
		ExportedClass_(const char* name, const char* group = "", const char *tooltip = "")
				: m_constructor(NULL), m_destructor(NULL), m_tooltip(tooltip)
		{
			ClassNameProvider<TClass>::set_name(name, group, true);
			ClassNameProvider<const TClass>::set_name(name, group, true);
		}

	/// name of class
		virtual const char* name() const {return ClassNameProvider<TClass>::name();}

	///	name node of class
		virtual const ClassNameNode& class_name_node() const {return ClassNameProvider<TClass>::class_name_node();}

	/// tooltip
		virtual const std::string& tooltip() const{return m_tooltip;}

	///	get groups
		virtual const std::string& group() const {return ClassNameProvider<TClass>::group();}

	///	class-hierarchy
		virtual const std::vector<const char*>* class_names() const	{return &ClassNameProvider<TClass>::names();}

	/// number of registered methods (overloads are not counted)
		virtual size_t num_methods() const {return m_vMethod.size();}

	///	number of registered const-methods (overloads are not counted)
		virtual size_t num_const_methods() const {return m_vConstMethod.size();}
		
	/// returns the first overload of an exported function
		virtual const ExportedMethod& get_method(size_t i) const {return *m_vMethod.at(i)->get_overload(0);}

	/// returns the first overload of an exported const function
		virtual const ExportedMethod& get_const_method(size_t i) const {return *m_vConstMethod.at(i)->get_overload(0);}

	///	returns the number of overloads of a method
		virtual size_t num_overloads(size_t funcInd) const			{return m_vMethod.at(funcInd)->num_overloads();}

	///	returns the number of overloads of a const method
		virtual size_t num_const_overloads(size_t funcInd) const	{return m_vConstMethod.at(funcInd)->num_overloads();}

	///	returns the i-th overload of a method
		virtual const ExportedMethod& get_overload(size_t funcInd, size_t oInd) const	{return *m_vMethod.at(funcInd)->get_overload(oInd);}

	///	returns the i-th overload of a const method
		virtual const ExportedMethod& get_const_overload(size_t funcInd, size_t oInd) const	{return *m_vConstMethod.at(funcInd)->get_overload(oInd);}

	///	returns the i-th method group (all overloads of the i-th function)
		virtual const ExportedMethodGroup& get_method_group(size_t ind) const		{return *m_vMethod.at(ind);}

	///	returns the i-th method group (all overloads of the i-th function)
		virtual const ExportedMethodGroup& get_const_method_group(size_t ind) const	{return *m_vConstMethod.at(ind);}


	/// Method registration
		template <typename TMethod>
		ExportedClass_<TClass>& add_method (	const char* methodName, TMethod func,
												const char* retValInfos = "", const char* paramInfos = "",
												const char* tooltip = "", const char* help = "")
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
			{
				const size_t start = strippedMethodName.find_first_not_of(" \t");
				const size_t end = strippedMethodName.find_last_not_of(" \t");
				if(start != std::string::npos && end != std::string::npos)
					strippedMethodName = strippedMethodName.substr(start, end - start + 1);
			}
			{
				const size_t start = methodOptions.find_first_not_of(" \t");
				const size_t end = methodOptions.find_last_not_of(" \t");
				if(start != std::string::npos && end != std::string::npos)
					methodOptions = methodOptions.substr(start, end - start + 1);
			}


		// 	check that name is not empty
			if(strippedMethodName.empty())
			{
				UG_LOG("### Registry ERROR: Trying to register empty method name."
						<< "\n### Please change register process. Aborting ..." << std::endl);
				throw(UG_REGISTRY_ERROR_RegistrationFailed(strippedMethodName.c_str()));
			}

		//	if the method is already in use, we have to add an overload.
		//	If not, we have to create a new method group
			ExportedMethodGroup* methodGrp = NULL;
			if(func_traits<TMethod>::const_method)
				methodGrp = get_const_exported_method_group(strippedMethodName.c_str());
			else
				methodGrp = get_exported_method_group(strippedMethodName.c_str());

			if(!methodGrp){
				methodGrp = new ExportedMethodGroup(strippedMethodName.c_str());
				if(func_traits<TMethod>::const_method)
					m_vConstMethod.push_back(methodGrp);
				else
					m_vMethod.push_back(methodGrp);
			}

			bool success = methodGrp->add_overload(func, &MethodProxy<TClass, TMethod>::apply,
												   ClassNameProvider<TClass>::name(),
													methodOptions.c_str(), retValInfos, paramInfos,
													tooltip, help);

			if(!success)
			{
				UG_LOG("### Registry ERROR: Trying to register method name '" << strippedMethodName
						<< "' to class '" << name() << "', but another method with this name "
						<< " and the same function signature is already registered for this class."
						<< "\n### Please change register process. Aborting ..." << std::endl);
				throw(UG_REGISTRY_ERROR_RegistrationFailed(name()));
			}

			return *this;
		}

	/// Make constructor accessible
	//  We use a pointer to ConstructorProxy since abstract base classes can not be created but registered.
	//  Each class that is instantiable must register its constructor
		ExportedClass_<TClass>& add_constructor ()
		{
		//  remember constructor proxy
			m_constructor = &ConstructorProxy<TClass>;

			//  remember constructor proxy
			m_destructor = &DestructorProxy<TClass>;

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

	///	destructur for object
		virtual void destroy(void* obj) const
		{
			if(m_destructor != NULL)
				(*m_destructor)(obj);
		}

	///	return pointer to the delete method
		virtual DeleteFunction get_delete_function() const
		{
			return CastAndDelete<TClass>;
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

		ExportedMethodGroup* get_exported_method_group(const char* name)
		{
			for(size_t i = 0; i < m_vMethod.size(); ++i)
			{
			//  compare strings
				if(strcmp(name, (m_vMethod[i]->name()).c_str()) == 0)
					return m_vMethod[i];
			}
			return NULL;
		}

		ExportedMethodGroup* get_const_exported_method_group(const char* name)
		{
			for(size_t i = 0; i < m_vConstMethod.size(); ++i)
			{
			//  compare strings
				if(strcmp(name, (m_vConstMethod[i]->name()).c_str()) == 0)
					return m_vConstMethod[i];
			}
			return NULL;
		}

	private:
		typedef TClass* (*ConstructorFunc)();
		ConstructorFunc m_constructor;

		typedef void (*DestructorFunc)(void*);
		DestructorFunc m_destructor;

		std::vector<ExportedMethodGroup*> m_vMethod;
		std::vector<ExportedMethodGroup*> m_vConstMethod;
		std::string m_tooltip;
};

} // end namespace bridge
} // end namespace ug


#endif /* __H__UG_BRIDGE__CLASS__ */
