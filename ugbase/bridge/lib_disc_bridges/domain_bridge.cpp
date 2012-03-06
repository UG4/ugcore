/*
 * domain_bridge.cpp
 *
 *  Created on: 22.09.2010
 *      Author: andreasvogel
 */

// extern headers
#include <iostream>
#include <sstream>
#include <string>

// include bridge
#include "../bridge.h"
#include "registry/registry.h"

// lib_algebra includes
#include "lib_algebra/cpu_algebra_types.h"
#include "lib_algebra/operator/operator_util.h"
#include "lib_algebra/operator/operator_interface.h"
#include "lib_algebra/operator/operator_inverse_interface.h"

// lib_disc includes
#include "lib_disc/domain.h"
#include "lib_disc/spatial_disc/domain_disc.h"
#include "lib_disc/function_spaces/grid_function.h"
#include "lib_disc/function_spaces/approximation_space.h"
#include "lib_disc/function_spaces/grid_function_util.h"
#include "lib_disc/function_spaces/interpolate.h"
#include "lib_disc/function_spaces/integrate.h"
#include "lib_disc/function_spaces/integrateDraft.h"
#include "lib_disc/function_spaces/error_indicator.h"
#include "lib_disc/dof_manager/cuthill_mckee.h"
#include "lib_disc/dof_manager/lexorder.h"

using namespace std;

namespace ug{
namespace bridge{

template <typename TDomain, typename TAlgebra>
void Register__Algebra_Domain(Registry& reg, string parentGroup)
{
//	typedef
	static const int dim = TDomain::dim;
	typedef typename TAlgebra::vector_type vector_type;
	typedef typename TAlgebra::matrix_type matrix_type;
	typedef ApproximationSpace<TDomain> approximation_space_type;
	typedef GridFunction<TDomain, SurfaceDoFDistribution, TAlgebra> TFct;

//	suffix and tag
	string dimAlgSuffix = GetDomainSuffix<TDomain>();
	dimAlgSuffix.append(GetAlgebraSuffix<TAlgebra>());

	string dimAlgTag = GetDomainTag<TDomain>();
	dimAlgTag.append(GetAlgebraTag<TAlgebra>());

//	group string
	string approxGrp = parentGroup; approxGrp.append("/ApproximationSpace");
	string domDiscGrp = parentGroup; domDiscGrp.append("/SpatialDisc");

//	DomainDiscretization
	{
		typedef IDomainDiscretization<TAlgebra> TBase;
		typedef DomainDiscretization<TDomain, TAlgebra> T;
		string name = string("DomainDiscretization").append(dimAlgSuffix);
		reg.add_class_<T, TBase>(name, domDiscGrp)
			.template add_constructor<void (*)(SmartPtr<ApproximationSpace<TDomain> >)>("ApproximationSpace")
			.add_method("add", static_cast<void (T::*)(IDomainConstraint<TDomain, TAlgebra>&)>(&T::add), "", "Post Process")
			.add_method("add", static_cast<void (T::*)(IDomainElemDisc<TDomain>&)>(&T::add), "", "Element Discretization")
			.add_method("add", static_cast<void (T::*)(IDiscretizationItem<TDomain, TAlgebra>&)>(&T::add), "", "DiscItem");
		reg.add_class_to_group(name, "DomainDiscretization", dimAlgTag);
	}

//	IDiscretizationItem
	{
		typedef IDiscretizationItem<TDomain, TAlgebra> T;
		string name = string("IDiscretizationItem").append(dimAlgSuffix);
		reg.add_class_<T>(name, domDiscGrp);
		reg.add_class_to_group(name, "IDiscretizationItem", dimAlgTag);
	}

//	GridFunction
	{
		string name = string("GridFunction").append(dimAlgSuffix);
		reg.add_class_<TFct, vector_type>(name, approxGrp)
			.template add_constructor<void (*)(SmartPtr<approximation_space_type>)>("ApproximationSpace")
			.add_method("assign", static_cast<void (TFct::*)(const vector_type&)>(&TFct::assign),
						"Success", "Vector")
			.add_method("clone", &TFct::clone);
		reg.add_class_to_group(name, "GridFunction", dimAlgTag);
	}

	string grp = parentGroup; grp.append("");

//	InterpolateFunction
	{
		typedef bool (*fct_type)(
				const boost::function<void (number& res,const MathVector<dim>& x, number time)>&,
				TFct&, const char*, number);
		reg.add_function("InterpolateFunction",
						 static_cast<fct_type>(&InterpolateFunction<TFct>),
						 grp);

		typedef bool (*fct_type_subset)(
				const boost::function<void (number& res,const MathVector<dim>& x, number time)>&,
				TFct&, const char*, number, const char*);
		reg.add_function("InterpolateFunction",
						 static_cast<fct_type_subset>(&InterpolateFunction<TFct>),
						 grp);
	}

//	L2Error
	{
		typedef number (*fct_type)(
				const boost::function<void (number& res,const MathVector<dim>& x, number time)>&,
				TFct&, const char*, number);
		reg.add_function("L2Error",
						 static_cast<fct_type>(&L2Error<TFct>),
						 grp);

		typedef number (*fct_type_subset)(
				const boost::function<void (number& res,const MathVector<dim>& x, number time)>&,
				TFct&, const char*, number, const char*);
		reg.add_function("L2Error",
						 static_cast<fct_type_subset>(&L2Error<TFct>),
						 grp);
	}

//	L2ErrorDraft
	{
		typedef number (*fct_type)(
				const boost::function<void (number& res,const MathVector<dim>& x, number time)>&,
				TFct&, const char*, number, int, const char*);
		reg.add_function("L2ErrorDraft",
						 static_cast<fct_type>(&L2ErrorDraft<TFct>),
						 grp);
	}

//	L2Norm
	{
		typedef number (*fct_type)(TFct&, const char*, int, const char*);


		reg.add_function("L2Norm",
						 static_cast<fct_type>(&L2Norm<TFct>),
						 grp);
	}

//	MarkForRefinement_GradientIndicator
	{
		string grp("ug4/Refinement/");
		reg.add_function("MarkForRefinement_GradientIndicator",
						 &MarkForRefinement_GradientIndicator<TDomain, SurfaceDoFDistribution, TAlgebra>, grp);
	}
}


template <typename TAlgebra>
static void Register__Algebra(Registry& reg, string parentGroup)
{
//	get group string
	string grp = parentGroup; grp.append("/Discretization");

	try{
#ifdef UG_DIM_1
		Register__Algebra_Domain<Domain1d, TAlgebra>(reg, grp);
#endif
#ifdef UG_DIM_2
		Register__Algebra_Domain<Domain2d, TAlgebra>(reg, grp);
#endif
#ifdef UG_DIM_3
		Register__Algebra_Domain<Domain3d, TAlgebra>(reg, grp);
#endif
	}
	catch(UG_REGISTRY_ERROR_RegistrationFailed ex)
	{
		UG_LOG("### ERROR in Register__Algebra: "
				"Registration failed (using name " << ex.name << ").\n");
		UG_THROW_FATAL("Registration failed.");
	}
}

template <typename TDomain>
void RegisterLibDiscDomain__Domain(Registry& reg, string parentGroup)
{
//	typedef
//	static const int dim = TDomain::dim;
	typedef ApproximationSpace<TDomain> approximation_space_type;

//	group string
	string approxGrp = parentGroup; approxGrp.append("/ApproximationSpace");

//	suffix and tag
	string dimSuffix = GetDomainSuffix<TDomain>();
	string dimTag = GetDomainTag<TDomain>();

//	Order Cuthill-McKee
	{
		reg.add_function("OrderCuthillMcKee", static_cast<void (*)(approximation_space_type&, bool)>(&OrderCuthillMcKee), approxGrp);
	}

//	Order lexicographically
	{
		reg.add_function("OrderLex", static_cast<void (*)(approximation_space_type&, const char*)>(&OrderLex<TDomain>), approxGrp);
	}
}

bool RegisterLibDisc_Domain(Registry& reg, string parentGroup)
{
	reg.add_class_<GridLevel>("GridLevel")
		.add_constructor();

	try{
#ifdef UG_CPU_1
	Register__Algebra<CPUAlgebra>(reg, parentGroup);
#endif
#ifdef UG_CPU_2
	Register__Algebra<CPUBlockAlgebra<2> >(reg, parentGroup);
#endif
#ifdef UG_CPU_3
	Register__Algebra<CPUBlockAlgebra<3> >(reg, parentGroup);
#endif
#ifdef UG_CPU_4
	Register__Algebra<CPUBlockAlgebra<4> >(reg, parentGroup);
#endif
#ifdef UG_CPU_VAR
	Register__Algebra<CPUVariableBlockAlgebra >(reg, parentGroup);
#endif

#ifdef UG_DIM_1
	RegisterLibDiscDomain__Domain<Domain1d>(reg, parentGroup);
#endif
#ifdef UG_DIM_2
	RegisterLibDiscDomain__Domain<Domain2d>(reg, parentGroup);
#endif
#ifdef UG_DIM_3
	RegisterLibDiscDomain__Domain<Domain3d>(reg, parentGroup);
#endif
	}
	catch(UG_REGISTRY_ERROR_RegistrationFailed ex)
	{
		UG_LOG("### ERROR in RegisterLibDiscDomain__Domain: "
				"Registration failed (using name " << ex.name << ").\n");
		UG_THROW_FATAL("Registration failed.");
	}

	return true;
}

}//	end of namespace ug
}//	end of namespace interface
