
#ifndef __H__UG_INTERFACE__UGBRIDGE__REGISTRY__
#define __H__UG_INTERFACE__UGBRIDGE__REGISTRY__

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
			m_vFunction.push_back(new ExportedFunction(	(void*) func, &ProxyFunction<TFunc>,
														funcName, retValName, paramValNames,
														tooltip, help));

		//  create parameter in list
			ParameterStack& in = m_vFunction.back()->params_in();
			typedef typename func_traits<TFunc>::params_type params_type;
			CreateParameterStack<params_type>::create(in, paramValNames, ",");

		//  create parameter out list
			ParameterStack& out = m_vFunction.back()->params_out();
			typedef typename func_traits<TFunc>::result_type result_type;
			CreateParameterStack<TypeList<result_type> >::create(out, retValName, ",");

			return *this;
		}

	/// number of function registered at the Registry
		size_t num_functions()							{return m_vFunction.size();}

	/// returns an exported function
		const ExportedFunction& get_function(size_t ind){return *m_vFunction.at(ind);}

	///////////////////
	// classes
	///////////////////

	/** Register a class at this registry
	 * This function registers any class
	 */
		template <typename TClass>
		ExportedClass_<TClass>& add_class_(const char *className)
		{
			m_vClass.push_back(ExportedClass_<TClass>::get_inst(className));

			return *m_vClass.back();
		}

	/// number of classes registered at the Registry
		size_t num_classes()							{return m_vClass.size();}

	/// returns an exported function
		const IExportedClass& get_class(size_t ind)	{return *m_vClass.at(ind);}


	/// destructor
		~InterfaceRegistry()
		{
			for(size_t i = 0; i < m_vFunction.size(); ++i)
			{
				if(m_vFunction[i] != NULL)
					delete m_vFunction[i];
			}
		}

	private:
		std::vector<ExportedFunction*>	m_vFunction;

		std::vector<IExportedClass*> m_vClass;
};

} // end namespace interface

} // end namespace ug


#endif /* __H__UG_INTERFACE__UGBRIDGE__REGISTRY__ */
