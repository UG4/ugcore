
#ifndef __H__UG_BRIDGE__REGISTRY__
#define __H__UG_BRIDGE__REGISTRY__

#include <vector>

#include "global_function.h"
#include "class.h"
#include "param_to_type_value_list.h"
#include "parameter_stack.h"

namespace ug
{
namespace bridge
{

// Registry
/** registers functions and classes that are exported to scripts and visualizations
 *
 */
class Registry {
	public:
		Registry()	{}
	//////////////////////
	// global functions
	//////////////////////

	/**	References the template function proxy_function<TFunc> and stores
	 * it with the FuntionWrapper.
	 */
		template<class TFunc>
		Registry& add_function(const char* funcName, TFunc func, const char* group = "",
								const char* retValName = "", const char* paramValNames = "",
								const char* tooltip = "", const char* help = "")
		{
		//  create new exported function
			m_vFunction.push_back(new ExportedFunction(	func, &FunctionProxy<TFunc>::apply,
														funcName, group,
														retValName, paramValNames,
														tooltip, help));
	
			return *this;
		}

	/// number of function registered at the Registry
		size_t num_functions() const						{return m_vFunction.size();}

	/// returns an exported function
		ExportedFunction& get_function(size_t ind) 			{return *m_vFunction.at(ind);}

	///////////////////
	// classes
	///////////////////

	/** Register a class at this registry
	 * This function registers any class
	 */
		template <typename TClass>
		ExportedClass_<TClass>& add_class_(const char* className, const char* group = "")
		{
		//	todo: check whether a class with the specified name exists aready.
			ExportedClass_<TClass>* newClass = new ExportedClass_<TClass>(className, group);

			m_vClass.push_back(newClass);

			return *newClass;
		}

	/** Register a class at this registry
	 * This function registers any class together with its base class
	 */
		template <typename TClass, typename TBaseClass>
		ExportedClass_<TClass>& add_class_(const char* className, const char* group = "")
		{
			ExportedClass_<TClass>* newClass = new ExportedClass_<TClass>(className, group);

			// set base class names
			try
			{
				ClassNameProvider<TClass>::template set_name<TBaseClass>(className, group);
			}
			catch(ug::bridge::UG_ERROR_ClassUnknownToRegistry* ex)
			{
				std::cout << "Registering class '" << className << "', that derives"
						 << " from class, that has " 
						 << "not yet been registered to this Registry.\n";
				exit(1);
			}

			m_vClass.push_back(newClass);
			return *newClass;
		}

	/// number of classes registered at the Registry
		size_t num_classes() const						{return m_vClass.size();}

	/// returns an exported function
		const IExportedClass& get_class(size_t ind)	const {return *m_vClass.at(ind);}


	/// destructor
		~Registry()
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
		Registry(const Registry& reg)	{}
		
		std::vector<ExportedFunction*>	m_vFunction;

		std::vector<IExportedClass*> m_vClass;
};

} // end namespace registry

} // end namespace ug


#endif /* __H__UG_BRIDGE__REGISTRY__ */
