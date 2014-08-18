/*
 * algebra_template_define_helper.h
 *
 *  Created on: 18.08.2014
 *      Author: mrupp
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
		UG_ALGEBRA_CPP_TEMPLATE_DEFINE_VAR(TheTemplateClassType) \
		UG_ALGEBRA_CPP_TEMPLATE_DEFINE_GPU(TheTemplateClassType)


#endif /* ALGEBRA_TEMPLATE_DEFINE_HELPER_H_ */
