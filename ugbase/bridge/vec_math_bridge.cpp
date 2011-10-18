// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 14.07.2011 (m,d,y)
 
#include <string>
#include <sstream>

#include "bridge.h"
#include "common/math/ugmath.h"

using namespace std;

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
static bool RegisterVecMathBridge(Registry& reg, string grp)
{
	typedef MathVector<dim, number> vec_type;

	try
	{
	//	the dimension suffix
		stringstream ss;	ss << dim << "d";
		string dimSuffix = ss.str();

	//	the dimension tag
		string dimTag = "dim=";
		dimTag.append(dimSuffix);

	//	register the class
		string vecName = "Vector";
		vecName.append(dimSuffix);

		reg.add_class_<vec_type>(vecName, grp)
			.add_method("coord",
					static_cast<const number& (vec_type::*)(size_t) const>(&vec_type::coord));
		reg.add_class_to_group(vecName, "Vector_d", dimTag);
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
static bool RegisterVecMathBridge_DimIndep(Registry& reg, string grp)
{
	try
	{
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
bool RegisterVecMathBridge(Registry& reg, string parentGroup)
{
//	get group string
	string grp = parentGroup; grp.append("/Util/VecMath");

	bool bSuccess = true;

	bSuccess &= RegisterVecMathBridge<1>(reg, grp);
	bSuccess &= RegisterVecMathBridge<2>(reg, grp);
	bSuccess &= RegisterVecMathBridge<3>(reg, grp);
	bSuccess &= RegisterVecMathBridge<4>(reg, grp);
	bSuccess &= RegisterVecMathBridge_DimIndep(reg, grp);

	return bSuccess;
}

}// end of namespace
}// end of namespace
