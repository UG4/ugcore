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

////////////////////////////////////////////////////////
// P1 Conform Function Pattern
////////////////////////////////////////////////////////

class P1ConformFunctionPatternObject : public
		ObjectBase_ClassWrapper<P1ConformFunctionPatternObject,
								P1ConformFunctionPattern,
								IObject>
{
	public:
		static const char* static_type_name()	{return "P1ConformFunctionPattern";}

		P1ConformFunctionPatternObject()
		{
			typedef P1ConformFunctionPatternObject THIS;

			{//	add_discrete_function
				MethodDesc& md = add_method("add_discrete_function", &THIS::add_discrete_function);
				md.params_in().add_string("name");
				md.params_in().add_string("type");
				md.params_in().add_int("dim");
				md.params_out().add_int("success");
			}

			{//	lock
				MethodDesc& md = add_method("lock", &THIS::lock);
				md.params_out().add_int("success");
			}
		}

	public:
		P1ConformFunctionPattern& get_function_pattern() {return *this->get_inst();}

	protected:
		void add_discrete_function(ParameterList& in, ParameterList& out)
		{
			string func_name = in.to_string(0);
			string func_type = in.to_string(1);
			int dim = in.to_int(2);

			P1ConformFunctionPattern& pattern = get_function_pattern();

			if(func_type != "P1")
			{
				UG_LOG("Only P1 Functions supported.\n");
				out.set_int(0, 1);
				return;
			}

			bool res = pattern.add_discrete_function(func_name, LSFS_LAGRANGEP1, dim);

			out.set_int(0, res);
		}

		void lock(ParameterList& in, ParameterList& out)
		{
			P1ConformFunctionPattern& pattern = get_function_pattern();

			out.set_int(0, pattern.lock());
		}
};

////////////////////////////////////////////////////////
// Approximation Space
////////////////////////////////////////////////////////

template <int dim, typename TAlgebra>
class ApproximationSpaceObject : public
	ObjectBase_ClassWrapper<ApproximationSpaceObject<dim, TAlgebra>,
							ApproximationSpace<Domain<dim, MultiGrid, MGSubsetHandler>, P1ConformDoFDistribution, TAlgebra>,
							IObject>
{
	typedef MultiGrid grid_type;
	typedef MGSubsetHandler subset_handler_type;
	typedef Domain<dim, grid_type, subset_handler_type> domain_type;
	typedef ApproximationSpace<domain_type, P1ConformDoFDistribution, TAlgebra> approximation_space_type;

	public:
		static const char* static_type_name()
		{
		   UG_STATIC_ASSERT(dim < 10, DIMENSION_TOO_HIGH);
		   static char name[22];
		   static int dummy = sprintf(name, "ApproximationSpace%id", dim);
		   dummy = 0; // avoid warnings
		   return name;
		}

		ApproximationSpaceObject()
		{
			typedef ApproximationSpaceObject THIS;

			{//	set_domain method
				MethodDesc& md = add_method("set_domain", &THIS::set_domain);
				char name[12]; sprintf(name, "MGDomain%id", dim);
				md.params_in().add_object(name, name);
			}

			{//	set_function_pattern method
				MethodDesc& md = add_method("set_function_pattern", &THIS::set_function_pattern);
				md.params_in().add_object("P1ConformFunctionPattern", "P1ConformFunctionPattern");
			}
		}

	public:
		approximation_space_type& get_approximation_space() {return *this->get_inst();}

	protected:
		void set_domain(ParameterList& in, ParameterList& out)
		{
			// unpack Object
			MGDomainObject<dim>* domainObj = static_cast<MGDomainObject<dim>*>(in.to_object(0));

			// get c++ types
			domain_type& domain = domainObj->get_domain();
			approximation_space_type& approxSpace = get_approximation_space();

			approxSpace.assign_domain(domain);
		}

		void set_function_pattern(ParameterList& in, ParameterList& out)
		{
			// unpack Object
			P1ConformFunctionPatternObject* funcPatternObj = static_cast<P1ConformFunctionPatternObject*>(in.to_object(0));

			// get c++ types
			P1ConformFunctionPattern& funcPattern = funcPatternObj->get_function_pattern();
			approximation_space_type& approxSpace = get_approximation_space();

			approxSpace.assign_function_pattern(funcPattern);
		}
};


} // end namespace interface
} // end namespace ug

#endif /* __H__UG_INTERFACE__STANDARD_INTERFACES__LIB_DISCRETIZATION_INTERFACE__ */
