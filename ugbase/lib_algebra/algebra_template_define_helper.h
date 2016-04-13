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

#ifndef ALGEBRA_TEMPLATE_DEFINE_HELPER_H_
#define ALGEBRA_TEMPLATE_DEFINE_HELPER_H_

#include "cpu_algebra_types.h"

#ifdef UG_CPU_1
	#define UG_ALGEBRA_CPP_TEMPLATE_DEFINE_1(TheTemplateClassType) \
		template class TheTemplateClassType<ug::CPUAlgebra>;
#else
	#define UG_ALGEBRA_CPP_TEMPLATE_DEFINE_1(TheTemplateClassType)
#endif

#ifdef UG_CPU_2
	#define UG_ALGEBRA_CPP_TEMPLATE_DEFINE_2(TheTemplateClassType) \
		template class TheTemplateClassType<ug::CPUBlockAlgebra<2> >;
#else
	#define UG_ALGEBRA_CPP_TEMPLATE_DEFINE_2(TheTemplateClassType)
#endif

#ifdef UG_CPU_3
	#define UG_ALGEBRA_CPP_TEMPLATE_DEFINE_3(TheTemplateClassType) \
		template class TheTemplateClassType<ug::CPUBlockAlgebra<3> >;
#else
	#define UG_ALGEBRA_CPP_TEMPLATE_DEFINE_3(TheTemplateClassType)
#endif

#ifdef UG_CPU_4
	#define UG_ALGEBRA_CPP_TEMPLATE_DEFINE_4(TheTemplateClassType) \
		template class TheTemplateClassType<ug::CPUBlockAlgebra<4> >;
#else
	#define UG_ALGEBRA_CPP_TEMPLATE_DEFINE_4(TheTemplateClassType)
#endif

#ifdef UG_CPU_5
	#define UG_ALGEBRA_CPP_TEMPLATE_DEFINE_5(TheTemplateClassType) \
		template class TheTemplateClassType<ug::CPUBlockAlgebra<5> >;
#else
	#define UG_ALGEBRA_CPP_TEMPLATE_DEFINE_5(TheTemplateClassType)
#endif

#ifdef UG_CPU_6
	#define UG_ALGEBRA_CPP_TEMPLATE_DEFINE_6(TheTemplateClassType) \
		template class TheTemplateClassType<ug::CPUBlockAlgebra<6> >;
#else
	#define UG_ALGEBRA_CPP_TEMPLATE_DEFINE_6(TheTemplateClassType)
#endif

#ifdef UG_CPU_VAR
	#define UG_ALGEBRA_CPP_TEMPLATE_DEFINE_VAR(TheTemplateClassType) \
		template class TheTemplateClassType<ug::CPUVariableBlockAlgebra >;
#else
	#define UG_ALGEBRA_CPP_TEMPLATE_DEFINE_VAR(TheTemplateClassType)
#endif

#ifdef UG_GPU
	#define UG_ALGEBRA_CPP_TEMPLATE_DEFINE_GPU(TheTemplateClassType) \
		template class TheTemplateClassType<GPUAlgebra>;
#else
	#define UG_ALGEBRA_CPP_TEMPLATE_DEFINE_GPU(TheTemplateClassType)
#endif

/**
 * use this with templates of the form
 * template<typename TALgebra> class MyClass;
 *
 * UG_ALGEBRA_CPP_TEMPLATE_DEFINE_ALL(MyClass)
 * will generate the code
 * template class TheTemplateClassType<ug::CPUAlgebra>;
 * template class TheTemplateClassType<ug::CPUBlockAlgebra<2> >;
 * template class TheTemplateClassType<ug::CPUBlockAlgebra<3> >;
 * ...
 */
#define UG_ALGEBRA_CPP_TEMPLATE_DEFINE_ALL(TheTemplateClassType) \
		UG_ALGEBRA_CPP_TEMPLATE_DEFINE_1(TheTemplateClassType)\
		UG_ALGEBRA_CPP_TEMPLATE_DEFINE_2(TheTemplateClassType) \
		UG_ALGEBRA_CPP_TEMPLATE_DEFINE_3(TheTemplateClassType) \
		UG_ALGEBRA_CPP_TEMPLATE_DEFINE_4(TheTemplateClassType) \
		UG_ALGEBRA_CPP_TEMPLATE_DEFINE_5(TheTemplateClassType) \
		UG_ALGEBRA_CPP_TEMPLATE_DEFINE_6(TheTemplateClassType) \
		UG_ALGEBRA_CPP_TEMPLATE_DEFINE_VAR(TheTemplateClassType) \
		UG_ALGEBRA_CPP_TEMPLATE_DEFINE_GPU(TheTemplateClassType)


#endif /* ALGEBRA_TEMPLATE_DEFINE_HELPER_H_ */
