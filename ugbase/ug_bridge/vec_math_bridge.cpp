// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 14.07.2011 (m,d,y)
 
#include "ug_bridge.h"
#include "common/math/ugmath.h"

namespace ug{

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
template <int dim>
static bool RegisterVecMathBridge(Registry& reg, const char* parentGroup)
{
	typedef MathVector<dim, number> vec_type;

	try
	{
	//	the dimension suffix
		std::stringstream ss;	ss << dim << "d";
		std::string dimSuffix = ss.str();

	//	the dimension tag
		std::string dimTag = "dim=";
		dimTag.append(dimSuffix);

	//	get group string
		std::string grp = parentGroup;
		grp.append("/VecMath");

	//	register the class
		std::string vecName = "Vector";
		vecName.append(dimSuffix);

		reg.add_class_<vec_type>(vecName.c_str(), grp.c_str())
			.add_method("coord",
					static_cast<const number& (vec_type::*)(size_t) const>(&vec_type::coord));
		reg.add_class_to_group(vecName.c_str(), "Vector_d", dimTag.c_str());
	}
	catch(UG_REGISTRY_ERROR_RegistrationFailed ex)
	{
		UG_LOG("### ERROR in RegisterVecMathBridge: "
				"Registration failed (using name " << ex.name << ").\n");
		return false;
	}

	return true;
}

////////////////////////////////////////////////////////////////////////////////
static bool RegisterVecMathBridge_DimIndep(Registry& reg, const char* parentGroup)
{
	try
	{
	//	get group string
		std::stringstream groupString; groupString << parentGroup << "/VecMath";
		std::string strGrp = groupString.str();
		const char* grp = strGrp.c_str();

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
	catch(UG_REGISTRY_ERROR_RegistrationFailed ex)
	{
		UG_LOG("### ERROR in RegisterVecMathBridge_DimIndep: "
				"Registration failed (using name " << ex.name << ").\n");
		return false;
	}

	return true;
}


////////////////////////////////////////////////////////////////////////////////
bool RegisterVecMathBridge(Registry& reg, const char* parentGroup)
{
	bool bSuccess = true;

	bSuccess &= RegisterVecMathBridge<1>(reg, parentGroup);
	bSuccess &= RegisterVecMathBridge<2>(reg, parentGroup);
	bSuccess &= RegisterVecMathBridge<3>(reg, parentGroup);
	bSuccess &= RegisterVecMathBridge<4>(reg, parentGroup);
	bSuccess &= RegisterVecMathBridge_DimIndep(reg, parentGroup);

	return bSuccess;
}

}// end of namespace
}// end of namespace
