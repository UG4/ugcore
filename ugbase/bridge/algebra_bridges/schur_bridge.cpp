/*
 * Copyright (c) 2014-2015:  G-CSC, Goethe University Frankfurt
 * Author: Martin Rupp
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

// extern headers
#include <iostream>
#include <sstream>
#include <string>

// include bridge
#include "bridge/bridge.h"
#include "bridge/util.h"
#include "bridge/util_algebra_dependent.h"

// preconditioner
#include "lib_algebra/lib_algebra.h"
#ifdef UG_PARALLEL
#include "lib_algebra/operator/preconditioner/schur/schur.h"
#include "lib_algebra/operator/preconditioner/schur/schur_precond.h"
#include "lib_algebra/operator/preconditioner/schur/schur_complement_inverse.h"
#endif

#include "../util_overloaded.h"
using namespace std;

namespace ug{

namespace bridge{
namespace Schur{

/**
 * \defgroup precond_bridge Preconditioner Bridge
 * \ingroup algebra_bridge
 * \{
 */

/**
 * Class exporting the functionality. All functionality that is to
 * be used in scripts or visualization must be registered here.
 */



struct Functionality
{


/**
 * Function called for the registration of Algebra dependent parts.
 * All Functions and Classes depending on Algebra
 * are to be placed here when registering. The method is called for all
 * available Algebra types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <typename TAlgebra>
static void Algebra(Registry& reg, string grp)
{
	string suffix = GetAlgebraSuffix<TAlgebra>();
	string tag = GetAlgebraTag<TAlgebra>();

//	typedefs for this algebra

#ifdef UG_PARALLEL
typedef typename TAlgebra::vector_type vector_type;

	{
		string name = string("ISchurComplementInverse").append(suffix);
		reg.add_class_< ISchurComplementInverse<TAlgebra> >(name, grp) ;
		reg.add_class_to_group(name, "ISchurComplementInverse", tag);
	}

	// 	Schur complement preconditioner
	{
		typedef SchurPrecond<TAlgebra> T;
		typedef IPreconditioner<TAlgebra> TBase;
		string name = string("SchurComplement").append(suffix);
		reg.add_class_<T,TBase>(name, grp, "Schur complement preconditioner")
			.add_constructor()
			.add_method("set_dirichlet_solver", &T::set_dirichlet_solver, "","Dirichlet solver")
			.add_method("set_skeleton_solver", &T::set_skeleton_solver, "","Skeleton solver")
			.add_method("num_global_skeleton", &T::num_global_skeleton)
			.add_method("num_local_skeleton", &T::num_local_skeleton)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "SchurComplement", tag);
	}



	{
		typedef SchurInverseWithOperator<TAlgebra> T;
		typedef ISchurComplementInverse<TAlgebra> TBase;
		string name = string("SchurInverseWithOperator").append(suffix);
		reg.add_class_<T,TBase>(name, grp, "SchurInverseWithOperator")
			.ADD_CONSTRUCTOR( (SmartPtr<ILinearOperatorInverse<vector_type> > ) )("linOpInverse")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "SchurInverseWithOperator", tag);
	}

	{
		typedef SchurInverseWithFullMatrix<TAlgebra> T;
		typedef ISchurComplementInverse<TAlgebra> TBase;
		string name = string("SchurInverseWithFullMatrix").append(suffix);
		reg.add_class_<T,TBase>(name, grp, "SchurInverseWithFullMatrix")
			.ADD_CONSTRUCTOR( (SmartPtr<ILinearOperatorInverse<vector_type> > ) )("linOpInverse")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "SchurInverseWithFullMatrix", tag);
	}


	{
		typedef SchurInverseWithAGammaGamma<TAlgebra> T;
		typedef ISchurComplementInverse<TAlgebra> TBase;
		string name = string("SchurInverseWithAGammaGamma").append(suffix);
		reg.add_class_<T,TBase>(name, grp, "SchurInverseWithAGammaGamma")
			.ADD_CONSTRUCTOR( (SmartPtr<IPreconditionedLinearOperatorInverse<vector_type> > ) )("precLinOpInv")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "SchurInverseWithAGammaGamma", tag);
	}


	{
		typedef SchurInverseWithAutoFullMatrix<TAlgebra> T;
		typedef ISchurComplementInverse<TAlgebra> TBase;
		string name = string("SchurInverseWithAutoFullMatrix").append(suffix);
		reg.add_class_<T,TBase>(name, grp, "SchurInverseWithAutoFullMatrix")
			.ADD_CONSTRUCTOR( (SmartPtr<ILinearOperatorInverse<vector_type> > ) )("linOpInverse")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "SchurInverseWithAutoFullMatrix", tag);
	}



#endif
}


}; // end Functionality

// end group precond_bridge
/// \}

}// end Preconditioner

/// \addtogroup precond_bridge
void RegisterBridge_Schur(Registry& reg, string grp)
{
	grp.append("/Algebra/Preconditioner");
	typedef Schur::Functionality Functionality;

	try{
		RegisterAlgebraDependent<Functionality>(reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

} // namespace bridge
} // namespace ug
