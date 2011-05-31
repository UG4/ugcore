/*
 * lib_disc_bridge_discs.cpp
 *
 *  Created on: 22.09.2010
 *      Author: andreasvogel
 */

// extern headers
#include <iostream>
#include <sstream>

// include bridge
#include "../../../ug_bridge.h"
#include "../../../registry.h"

// other files
#include "navier_stokes_bridge.h"
#include "density_driven_flow_bridge.h"
#include "convection_diffusion_bridge.h"
#include "constant_equation_bridge.h"
#include "thermohaline_flow_bridge.h"

#include "lib_discretization/spatial_discretization/disc_helper/conv_shape_interface.h"
#include "lib_discretization/spatial_discretization/disc_helper/conv_shape.h"

namespace ug
{

namespace bridge
{

/////////////////////////////////////////////////////////////////////////////
// Convection Shapes
/////////////////////////////////////////////////////////////////////////////

template <int dim>
void RegisterConvectionShapes(Registry& reg, const char* parentGroup)
{
	std::stringstream grpSS; grpSS << parentGroup << "/" << dim << "d";
	std::string grp = grpSS.str();

//	IConvectionShapes
	{
		typedef IConvectionShapes<dim> T;
		std::stringstream ss; ss << "IConvectionShapes" << dim << "d";
		reg.add_class_<T>(ss.str().c_str(), grp.c_str());
	}

//	ConvectionShapesNoUpwind
	{
		typedef ConvectionShapesNoUpwind<dim> T;
		typedef IConvectionShapes<dim> TBase;
		std::stringstream ss; ss << "NoUpwind" << dim << "d";
		reg.add_class_<T, TBase>(ss.str().c_str(), grp.c_str())
			.add_constructor();
	}

//	ConvectionShapesFullUpwind
	{
		typedef ConvectionShapesFullUpwind<dim> T;
		typedef IConvectionShapes<dim> TBase;
		std::stringstream ss; ss << "FullUpwind" << dim << "d";
		reg.add_class_<T, TBase>(ss.str().c_str(), grp.c_str())
			.add_constructor();
	}

//	ConvectionShapesWeightedUpwind
	{
		typedef ConvectionShapesWeightedUpwind<dim> T;
		typedef IConvectionShapes<dim> TBase;
		std::stringstream ss; ss << "WeightedUpwind" << dim << "d";
		reg.add_class_<T, TBase>(ss.str().c_str(), grp.c_str())
			.add_method("set_weight", &T::set_weight)
			.add_constructor();
	}

//	ConvectionShapesPartialUpwind
	{
		typedef ConvectionShapesPartialUpwind<dim> T;
		typedef IConvectionShapes<dim> TBase;
		std::stringstream ss; ss << "PartialUpwind" << dim << "d";
		reg.add_class_<T, TBase>(ss.str().c_str(), grp.c_str())
			.add_constructor();
	}
}

bool RegisterDynamicLibDiscInterfaceDiscs(Registry& reg, int algebra_type, const char* parentGroup)
{
	bool bReturn = true;

	bReturn &= RegisterDynamicNavierStokesDisc(reg, algebra_type, parentGroup);
	bReturn &= RegisterDynamicDensityDrivenFlowDisc(reg, algebra_type, parentGroup);
	bReturn &= RegisterDynamicConvectionDiffusionDisc(reg, algebra_type, parentGroup);
	bReturn &= RegisterDynamicConstantEquationDisc(reg, algebra_type, parentGroup);
	bReturn &= RegisterDynamicThermohalineFlowDisc(reg, algebra_type, parentGroup);

#ifdef UG_DIM_1
	RegisterConvectionShapes<1>(reg, parentGroup);
#endif
#ifdef UG_DIM_2
	RegisterConvectionShapes<2>(reg, parentGroup);
#endif
#ifdef UG_DIM_3
	RegisterConvectionShapes<3>(reg, parentGroup);
#endif

	return bReturn;
}

}//	end of namespace ug
}//	end of namespace interface
