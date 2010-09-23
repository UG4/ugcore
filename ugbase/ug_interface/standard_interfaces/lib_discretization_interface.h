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
#include "common/static_assert.h"
#include "lib_discretization/lib_discretization.h"
#include "lib_grid_interface.h"

namespace ug
{
namespace interface
{

template <int dim>
class MGDomainObject : public
	ObjectBase_ClassWrapper<MGDomainObject<dim>,
							Domain<dim, MultiGrid, MGSubsetHandler>,
							IObject>
{
	typedef MultiGrid grid_type;
	typedef MGSubsetHandler subset_handler_type;
	typedef Domain<dim, grid_type, subset_handler_type> domain_type;

	public:
		static const char* static_type_name()
		{
		   UG_STATIC_ASSERT(dim < 10, DIMENSION_TOO_HIGH);
		   static char name[12];
		   static int dummy = sprintf(name, "MGDomain%id", dim);
		   dummy = 0; // avoid warnings
		   return name;
		}

		MGDomainObject()
		{
			typedef MGDomainObject THIS;

			{//	get_grid method
				MethodDesc& md = add_method("get_grid", &THIS::get_grid);
				md.params_out().add_object("multi_grid", "MultiGrid");
			}

			{//	get_subset_handler method
				MethodDesc& md = add_method("get_subset_handler", &THIS::get_subset_handler);
				md.params_out().add_object("mg_subset_handler", "MGSubsetHandler");
			}

		}

	public:
		domain_type& get_domain() {return *this->get_inst();}

	protected:
		void get_grid(ParameterList& in, ParameterList& out)
		{
			domain_type& dom = get_domain();
			out.set_object(0, MultiGridObject::create_with_existing_instance(&dom.get_grid()));
		}

		void get_subset_handler(ParameterList& in, ParameterList& out)
		{
			domain_type& dom = get_domain();
			out.set_object(0, MGSubsetHandlerObject::create_with_existing_instance(&dom.get_subset_handler()));
		}
};

} // end namespace interface
} // end namespace ug

#endif /* __H__UG_INTERFACE__STANDARD_INTERFACES__LIB_DISCRETIZATION_INTERFACE__ */
