
#ifndef __H__UG_INTERFACE__UGBRIDGE__REGISTRY__
#define __H__UG_INTERFACE__UGBRIDGE__REGISTRY__

#include <vector>

#include "global_function.h"
#include "class.h"
#include "param_to_type_value_list.h"
#include "parameter_stack.h"

namespace ug {

namespace interface{

// InterfaceRegistry
/** registers functions and classes that are exported to scripts and visualizations
 *
 */
class InterfaceRegistry {
	public:
		InterfaceRegistry()	{}
	//////////////////////
	// global functions
	//////////////////////

	/**	References the template function proxy_function<TFunc> and stores
	 * it with the FuntionWrapper.
	 */
		template<class TFunc>
		InterfaceRegistry& add_function(const char* funcName, TFunc func,
										const char* retValName = "", const char* paramValNames = "",
										const char* tooltip = "", const char* help = "")
		{
		//  create new exported function
			m_vFunction.push_back(new ExportedFunction(	(void*) func, &FunctionProxy<TFunc>::apply,
														funcName, retValName, paramValNames,
														tooltip, help));
	
		//  create parameter in list
			ParameterStack& in = m_vFunction.back()->params_in();
			typedef typename func_traits<TFunc>::params_type params_type;
			CreateParameterStack<params_type>::create(in);

			return *this;
		}

	/// number of function registered at the Registry
		size_t num_functions()							{return m_vFunction.size();}

	/// returns an exported function
		ExportedFunction& get_function(size_t ind){return *m_vFunction.at(ind);}

	///////////////////
	// classes
	///////////////////

	/** Register a class at this registry
	 * This function registers any class
	 */
		template <typename TClass>
		ExportedClass_<TClass>& add_class_(const char *className)
		{
		//	todo: check whether a class with the specified name exists aready.
			ExportedClass_<TClass>* newClass = new ExportedClass_<TClass>(className);

			m_vClass.push_back(newClass);

			return *newClass;
		}

	/** Register a class at this registry
	 * This function registers any class together with its base class
	 */
		template <typename TClass, typename TBaseClass>
		ExportedClass_<TClass>& add_class_(const char *className)
		{
			ExportedClass_<TClass>* newClass = new ExportedClass_<TClass>(className);
			m_vClass.push_back(newClass);

			// set base class names
			ClassNameProvider<TClass>::template set_name<TBaseClass>(className);

			return *newClass;
		}

	/// number of classes registered at the Registry
		size_t num_classes()							{return m_vClass.size();}

	/// returns an exported function
		const IExportedClass& get_class(size_t ind)	{return *m_vClass.at(ind);}


	/// destructor
		~InterfaceRegistry()
		{
		//  delete registered functions
			for(size_t i = 0; i < m_vFunction.size(); ++i)
			{
				if(m_vFunction[i] != NULL)
					delete m_vFunction[i];
			}
		//  delete registered classes
			for(size_t i = 0; i < m_vClass.size(); ++i)
			{
				if(m_vClass[i] != NULL)
					delete m_vClass[i];
			}
		}

	private:
		InterfaceRegistry(const InterfaceRegistry& reg)	{}
		
		std::vector<ExportedFunction*>	m_vFunction;

		std::vector<IExportedClass*> m_vClass;
};

} // end namespace interface

} // end namespace ug


#endif /* __H__UG_INTERFACE__UGBRIDGE__REGISTRY__ */
