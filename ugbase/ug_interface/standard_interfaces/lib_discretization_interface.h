/*
 * lib_discretization_interface.h
 *
 *  Created on: 22.09.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG_INTERFACE__STANDARD_INTERFACES__LIB_DISCRETIZATION_INTERFACE__
#define __H__UG_INTERFACE__STANDARD_INTERFACES__LIB_DISCRETIZATION_INTERFACE__

#include "../ug_interface.h"
#include "../registry.h"
#include "lib_discretization/lib_discretization.h"
#include "lib_grid_interface.h"

namespace ug
{
namespace interface
{

template <int dim>
class DomainObject : public ObjectBase<DomainObject<dim>, IObject>
{
	typedef MultiGrid grid_type;
	typedef MGSubsetHandler subset_handler_type;
	typedef Domain<dim, grid_type, subset_handler_type> domain_type;

	public:
		static const char* static_type_name()
		{
			static char name[9];
			static bool init = false;
			if(!init)
			{
				sprintf(name, "Domain%id", dim);
				init = true;
			}

			return name;
		}

		DomainObject() : m_gridObj(&m_domain.get_grid()), m_shObj(&m_domain.get_subset_handler())
		{
			typedef DomainObject THIS;

			{//	get_grid method
				MethodDesc& md = add_method("get_grid", &THIS::get_grid);
				md.params_out().add_object("grid", "MultiGrid");
			}

			{//	get_subset_handler method
				MethodDesc& md = add_method("get_subset_handler", &THIS::get_subset_handler);
				md.params_out().add_object("subsethandler", "MGSubsetHandler");
			}

		}

		~DomainObject(){}

	public:
		domain_type& get_domain() {return m_domain;}

	protected:
		void get_grid(ParameterList& in, ParameterList& out)
		{
			out.set_object(0, &m_gridObj);
		}

		void get_subset_handler(ParameterList& in, ParameterList& out)
		{
			out.set_object(0, &m_shObj);
		}

	private:
		Domain<dim, grid_type, subset_handler_type> m_domain;
		MultiGridObject m_gridObj;
		MGSubsetHandlerObject m_shObj;
};

} // end namespace interface
} // end namespace ug

#endif /* __H__UG_INTERFACE__STANDARD_INTERFACES__LIB_DISCRETIZATION_INTERFACE__ */
