/*
 * algebra_extensions.cpp
 *
 *  Created on: 04.03.2011
 *      Author: kosta
 */

/**	This file is only temporary. Methods which are declared here
 * should be reworked and moved into the main ug-section.
 * Please take care to not introduce unwanted dependencies.
 */
#include "ug.h"
#include "algebra_extensions.h"
#include "ug_bridge/ug_bridge.h"
#include "lib_algebra/algebra_chooser.h"

#include "ug_script/user_data/user_data.h"

#include <iostream>
#include <sstream>

namespace ug
{
extern enum_AlgebraType g_AlgebraType;
namespace bridge
{

template <typename TVector>
void KostaUpdate(TVector& xOut, const TVector& xOld, const TVector& v_m, const number dt,
			const LuaUserNumberNumberFunction& alpha, const LuaUserNumberNumberFunction& beta)
{
	UG_ASSERT(xOut.size() == xOld.size(), "Vector size does not match with first vector.");

	UG_ASSERT(xOut.size() == v_m.size(), "Vector size does not match with second vector.");


	for(size_t i = 0; i < xOut.size(); ++i)
	{
		const number V_m = v_m[i];
		const number x_old = xOld[i];

		xOut[i] = (x_old + dt*alpha(1, V_m))/(1 + dt*(alpha(1, V_m) + beta(1, V_m)));

	}
}

template <typename TVector>
void HhFlux(TVector& hhFlux, const TVector& nGate, const TVector& mGate, const TVector& hGate,
			const TVector& vm, const TVector& injection)
{
	UG_ASSERT(hhFlux.size() == nGate.size(), "Vector size does not match with first vector.");

	UG_ASSERT(hhFlux.size() == mGate.size(), "Vector size does not match with first vector.");
	UG_ASSERT(hhFlux.size() == hGate.size(), "Vector size does not match with first vector.");
	UG_ASSERT(hhFlux.size() == vm.size(), "Vector size does not match with first vector.");
	UG_ASSERT(hhFlux.size() == injection.size(), "Vector size does not match with first vector.");


	for(size_t i = 0; i < hhFlux.size(); i++)
	{
		const number V_m = vm[i];
		const number n = nGate[i];
		const number m = mGate[i];
		const number h = hGate[i];
		const number inj = injection[i];

		hhFlux[i] = (120*m*m*m*h*(V_m - 50) + 36*n*n*n*n*(V_m + 77) + 0.3*(V_m + 54.4) - inj);
	//	UG_LOG("hhFlux[i]="<<hhFlux[i]);
	}
}


bool RegisterAlgebraExtensions(Registry& reg, const char* parentGroup)
{
//	get group string
	std::stringstream groupString; groupString << parentGroup << "/AlgebraExtensions";
	std::string grp = groupString.str();
	reg.add_function("KostaUpdate", &KostaUpdate<CPUAlgebra::vector_type>, grp.c_str());
	reg.add_function("HhFlux", &HhFlux<CPUAlgebra::vector_type>, grp.c_str());
	return true;
}



}
}
