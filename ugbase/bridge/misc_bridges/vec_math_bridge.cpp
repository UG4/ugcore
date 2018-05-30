/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
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

#include <string>
#include <sstream>

#include "bridge/bridge.h"
#include "bridge/util.h"
#include "common/math/ugmath.h"
#include "bridge/suffix_tag.h"

using namespace std;

namespace ug{

/// \defgroup vecmath_bridge VecMath Bridge
/// \ingroup misc_bridge
/// \{

////////////////////////////////////////////////////////////////////////////////
//	Some methods which assist in creating vectors
static SmartPtr<vector1> MakeVec(number x)
{
	return SmartPtr<vector1>(new vector1(x));
}

static SmartPtr<vector2> MakeVec(number x, number y)
{
	return SmartPtr<vector2>(new vector2(x, y));
}

static SmartPtr<vector3> MakeVec(number x, number y, number z)
{
	return SmartPtr<vector3>(new vector3(x, y, z));
}

static SmartPtr<vector4> MakeVec(number x, number y, number z, number w)
{
	return SmartPtr<vector4>(new vector4(x, y, z, w));
}

////////////////////////////////////////////////////////////////////////////////
namespace bridge{
////////////////////////////////////////////////////////////////////////////////
/*
template <int dim>
static void RegisterBridge_VecMath(Registry& reg, string grp)
{
	typedef MathVector<dim, number> vec_type;
	
	try
	{
		string suffix = GetDimensionSuffix<dim>();
		string tag = GetDimensionTag<dim>();

	//	register the class
		string vecName = "Vec";
		vecName.append(suffix);

		reg.add_class_<vec_type>(vecName, grp)
			.add_method("coord",
					static_cast<typename vec_type::value_type (vec_type::*)(size_t) const>(&vec_type::coord));
		reg.add_class_to_group(vecName, "Vec", tag);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}
*/
////////////////////////////////////////////////////////////////////////////////
static void RegisterVecMathBridge_DimIndep(Registry& reg, string grp)
{
	try
	{
		{
			typedef MathVector<1, number> vec_type;
			reg.add_class_<vec_type>("Vec1d", grp)
				.add_constructor()
				.add_constructor<void (*)(number)>()
				.add_method("set_coord", &vec_type::set_coord, "", "index # value",
				            "sets the value of the coordinate with the given index")
				.add_method("coord",
						static_cast<vec_type::value_type (vec_type::*)(size_t) const>(&vec_type::coord));
			//reg.add_class_to_group("Vec1d", "Vec", GetDimensionTag<1>());
		}
		{
			typedef MathVector<2, number> vec_type;
			reg.add_class_<vec_type>("Vec2d", grp)
				.add_constructor()
				.add_constructor<void (*)(number, number)>()
				.add_method("set_coord", &vec_type::set_coord, "", "index # value",
				            "sets the value of the coordinate with the given index")
				.add_method("coord",
						static_cast<vec_type::value_type (vec_type::*)(size_t) const>(&vec_type::coord));
			//reg.add_class_to_group("Vec2d", "Vec", GetDimensionTag<2>());
		}
		{
			typedef MathVector<3, number> vec_type;
			reg.add_class_<vec_type>("Vec3d", grp)
				.add_constructor()
				.add_constructor<void (*)(number, number, number)>()
				.add_method("set_coord", &vec_type::set_coord, "", "index # value",
				            "sets the value of the coordinate with the given index")
				.add_method("coord",
						static_cast<vec_type::value_type (vec_type::*)(size_t) const>(&vec_type::coord));
			//reg.add_class_to_group("Vec3d", "Vec", GetDimensionTag<3>());
		}
		{
			typedef MathVector<4, number> vec_type;
			reg.add_class_<vec_type>("Vec4d", grp)
				.add_constructor()
				.add_constructor<void (*)(number, number, number, number)>()
				.add_method("set_coord", &vec_type::set_coord, "", "index # value",
				            "sets the value of the coordinate with the given index")
				.add_method("coord",
						static_cast<vec_type::value_type (vec_type::*)(size_t) const>(&vec_type::coord));
			//reg.add_class_to_group("Vec4d", "Vec", GetDimensionTag<4>());
		}

	//	register make-methods
		reg.add_function("MakeVec", static_cast<SmartPtr<vector1> (*)(number)>(
							&MakeVec), grp)
			.add_function("MakeVec", static_cast<SmartPtr<vector2> (*)(number, number)>(
							&MakeVec), grp)
			.add_function("MakeVec", static_cast<SmartPtr<vector3> (*)(number, number, number)>(
							&MakeVec), grp)
			.add_function("MakeVec", static_cast<SmartPtr<vector4> (*)(number, number, number, number)>(
							&MakeVec), grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}


////////////////////////////////////////////////////////////////////////////////
void RegisterBridge_VecMath(Registry& reg, string parentGroup)
{
//	get group string
	string grp = parentGroup; grp.append("/Util/VecMath");

//	RegisterBridge_VecMath<1>(reg, grp);
//	RegisterBridge_VecMath<2>(reg, grp);
//	RegisterBridge_VecMath<3>(reg, grp);
//	RegisterBridge_VecMath<4>(reg, grp);
	RegisterVecMathBridge_DimIndep(reg, grp);
}

// end group vecmath_bridge
/// \}

}// end of namespace
}// end of namespace
