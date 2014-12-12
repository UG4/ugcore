// author: Andreas Vogel

#include <iostream>
#include <sstream>
#include <boost/function.hpp>
#include <cmath>

// include bridge
#include "bridge/bridge.h"
#include "bridge/util.h"

#include "common/common.h"
#include "lib_disc/spatial_disc/user_data/const_user_data.h"
#include "lib_disc/spatial_disc/user_data/std_glob_pos_data.h"
#include "lib_disc/spatial_disc/user_data/linker/linker.h"
#include "lib_disc/spatial_disc/user_data/linker/scale_add_linker.h"
#include "lib_disc/spatial_disc/user_data/linker/inverse_linker.h"
#include "lib_disc/spatial_disc/user_data/linker/darcy_velocity_linker.h"
#include "lib_disc/spatial_disc/user_data/user_function.h"

using namespace std;

namespace ug{

template <class TData, int dim>
void PrintUserDataValue(const UserData<TData, dim>& ud, const MathVector<dim>& globIP,
						number time, int si)
{
	TData data;
	ud(data, globIP, time, si);
	UG_LOG(data);
}

template <class TData, int dim>
void PrintCondUserDataValue(const UserData<TData, dim, bool>& ud, const MathVector<dim>& globIP,
							number time, int si)
{
	TData data;
	bool ret = ud(data, globIP, time, si);
	UG_LOG(ret << ", " << data);
}


class RotatingCone2d
	: public StdGlobPosData<RotatingCone2d, number, 2>
{
	public:
		RotatingCone2d(double _eps, double _cx, double _cy, double _ax, double _ay,
		             double _nu, double _delta)
			: eps(_eps), cx(_cx), cy(_cy), ax(_ax), ay(_ay), nu(_nu), delta(_delta)
		{};

		inline void evaluate(number& val, const MathVector<2>& pos, number time, int si) const
		{
			const number t = time;
			const number x = pos[0];
			const number y = pos[1];

			const number xRot = cos(nu*t) * (x-cx) - sin(nu*t) * (y-cy);
			const number yRot = sin(nu*t) * (x-cx) + cos(nu*t) * (y-cy);

			const number expo = -((xRot - ax)*(xRot - ax) + (yRot - ay)*(yRot - ay)) / (delta + 4*eps*t);
			const number scale = delta/(delta+4*eps*t);

			val = scale * exp(expo);
		}

	protected:
		double eps, cx, cy, ax, ay, nu, delta;
};

class RotatingVelocity2d
	: public StdGlobPosData<RotatingVelocity2d, MathVector<2>, 2>
{
	public:
		RotatingVelocity2d(double _cx, double _cy, double _nu)
			: cx(_cx), cy(_cy), nu(_nu)
		{};

		inline void evaluate(MathVector<2>& val, const MathVector<2>& pos, number time, int si) const
		{
			const number x = pos[0];
			const number y = pos[1];

			val[0] = nu*(y - cx);
			val[1] = nu*(cy - x);
		}

	protected:
		double cx, cy, nu;
};


namespace bridge{

/**
 * \defgroup userdata_bridge User Data Bridge
 * \ingroup disc_bridge
 * \{
 */
template <typename TData, int dim>
void RegisterUserDataType(Registry& reg, string grp)
{
	string dimSuffix = GetDimensionSuffix<dim>();
	string dimTag = GetDimensionTag<dim>();

	string type = user_data_traits<TData>::name();

//	User"Type"
//	NOTE: For better readability this class is named User"Type"
//	      in vrl and lua. E.g. UserNumber, UserVector, ...
	{
		typedef UserData<TData, dim> T;
		typedef UserDataInfo TBase1;
		string name = string("User").append(type).append(dimSuffix);
		reg.add_class_<T,TBase1>(name, grp)
			.add_method("get_dim", &T::get_dim)
			.add_method("type", &T::type);
		reg.add_class_to_group(name, string("User").append(type), dimTag);
		reg.add_function("PrintUserDataValue", &PrintUserDataValue<TData, dim>,
						 grp, "", "userData#position#time#subsetIndex",
						 "Prints the value of the given user data at the given global position at the given time on the given subset.");
	}

//	CondUser"Type"
//	NOTE: For better readability this class is named CondUser"Type"
//	 	  in vrl and lua. E.g. CondUserNumber, CondUserVector, ...
	{
		typedef UserData<TData, dim, bool> T;
		typedef UserDataInfo TBase1;
		string name = string("CondUser").append(type).append(dimSuffix);
		reg.add_class_<T,TBase1>(name, grp);
		reg.add_class_to_group(name, string("CondUser").append(type), dimTag);
		reg.add_function("PrintUserDataValue", &PrintCondUserDataValue<TData, dim>,
						 grp, "", "userData#position#time#subsetIndex",
						 "Prints the value of the given user data at the given global position at the given time on the given subset.");
	}

//	CplUser"Type"
	{
		typedef CplUserData<TData, dim> T;
		typedef UserData<TData,dim> TBase1;
		string name = string("CplUser").append(type).append(dimSuffix);
		reg.add_class_<T,TBase1>(name, grp)
			.add_method("get_dim", &T::get_dim)
			.add_method("type", &T::type);
		reg.add_class_to_group(name, string("CplUser").append(type), dimTag);
	}

//	CondCplUser"Type"
	{
		typedef CplUserData<TData, dim, bool> T;
		typedef UserData<TData,dim,bool> TBase1;
		string name = string("CondCplUser").append(type).append(dimSuffix);
		reg.add_class_<T,TBase1>(name, grp);
		reg.add_class_to_group(name, string("CondCplUser").append(type), dimTag);
	}

//	DependentUserData"Type"
	{
		typedef DependentUserData<TData, dim> T;
		typedef CplUserData<TData, dim> TBase;
		string name = string("DependentUserData").append(type).append(dimSuffix);
		reg.add_class_<T, TBase >(name, grp);
		reg.add_class_to_group(name, string("DependentUserData").append(type), dimTag);
	}

//	ScaleAddLinker"Type"
	{
		typedef ScaleAddLinker<TData, dim, number> T;
		typedef DependentUserData<TData, dim> TBase;
		string name = string("ScaleAddLinker").append(type).append(dimSuffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_method("add", static_cast<void (T::*)(SmartPtr<CplUserData<number,dim> > , SmartPtr<CplUserData<TData,dim> >)>(&T::add))
			.add_method("add", static_cast<void (T::*)(number , SmartPtr<CplUserData<TData,dim> >)>(&T::add))
			.add_method("add", static_cast<void (T::*)(SmartPtr<CplUserData<number,dim> > , number)>(&T::add))
			.add_method("add", static_cast<void (T::*)(number,number)>(&T::add))
			.add_constructor()
			.template add_constructor<void (*)(const ScaleAddLinker<TData, dim, number>&)>()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, string("ScaleAddLinker").append(type), dimTag);
	}

}

// end group userdata_bridge
/// \}

namespace UserDataBridge{

/// \addtogroup userdata_bridge
/// \{

/**
 * Class exporting the functionality. All functionality that is to
 * be used in scripts or visualization must be registered here.
 */
struct Functionality
{

/**
 * Function called for the registration of Dimension dependent parts.
 * All Functions and Classes depending on the Dimension
 * are to be placed here when registering. The method is called for all
 * available Dimension types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <int dim>
static void Dimension(Registry& reg, string grp)
{
	string dimSuffix = GetDimensionSuffix<dim>();
	string dimTag = GetDimensionTag<dim>();

	RegisterUserDataType<number, dim>(reg, grp);
	RegisterUserDataType<MathVector<dim>, dim>(reg, grp);
	RegisterUserDataType<MathMatrix<dim,dim>, dim>(reg, grp);
	RegisterUserDataType<MathTensor<4,dim>, dim>(reg, grp);

//	ConstUserNumber
	{
		typedef ConstUserNumber<dim> T;
		typedef CplUserData<number, dim> TBase;
		string name = string("ConstUserNumber").append(dimSuffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.template add_constructor<void (*)(number)>("Value")
			.add_method("set", &T::set, "", "Value")
			.add_method("print", &T::print)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "ConstUserNumber", dimTag);
	}

//	ConstUserVector
	{
		typedef ConstUserVector<dim> T;
		typedef CplUserData<MathVector<dim>, dim> TBase;
		string name = string("ConstUserVector").append(dimSuffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.template add_constructor<void (*)(number)>("Values")
			.add_method("set_all_entries", &T::set_all_entries)
			.add_method("set_entry", &T::set_entry)
			.add_method("print", &T::print)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "ConstUserVector", dimTag);
	}

//	ConstUserMatrix
	{
		typedef ConstUserMatrix<dim> T;
		typedef CplUserData<MathMatrix<dim, dim>, dim> TBase;
		string name = string("ConstUserMatrix").append(dimSuffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.template add_constructor<void (*)(number)>("Diagonal Value")
			.add_method("set_diag_tensor", &T::set_diag_tensor)
			.add_method("set_all_entries", &T::set_all_entries)
			.add_method("set_entry", &T::set_entry)
			.add_method("print", &T::print)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "ConstUserMatrix", dimTag);
	}

//	ScaleAddLinkerMatrixVector
	{
		typedef MathVector<dim> TData;
		typedef MathMatrix<dim,dim> TDataScale;
		typedef ScaleAddLinker<TData, dim, TDataScale> T;
		typedef DependentUserData<TData, dim> TBase;
		string name = string("ScaleAddLinkerVectorMatrix").append(dimSuffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_method("add", static_cast<void (T::*)(SmartPtr<CplUserData<TDataScale,dim> > , SmartPtr<CplUserData<TData,dim> >)>(&T::add))
			.add_method("add", static_cast<void (T::*)(number , SmartPtr<CplUserData<TData,dim> >)>(&T::add))
			.add_method("add", static_cast<void (T::*)(SmartPtr<CplUserData<TDataScale,dim> > , number)>(&T::add))
			.add_method("add", static_cast<void (T::*)(number,number)>(&T::add))
			.add_constructor()
			.template add_constructor<void (*)(const ScaleAddLinker<TData, dim, TDataScale>&)>()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, string("ScaleAddLinkerVectorMatrix"), dimTag);
	}

//	ScaleAddLinkerVectorVector
	{
		typedef MathVector<dim> TData;
		typedef MathVector<dim> TDataScale;
		typedef ScaleAddLinker<TData, dim, TDataScale, number> T;
		typedef DependentUserData<number, dim> TBase;
		string name = string("ScaleAddLinkerVectorVector").append(dimSuffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_method("add", static_cast<void (T::*)(SmartPtr<CplUserData<TDataScale,dim> > , SmartPtr<CplUserData<TData,dim> >)>(&T::add))
			.add_method("add", static_cast<void (T::*)(number , SmartPtr<CplUserData<TData,dim> >)>(&T::add))
			.add_method("add", static_cast<void (T::*)(SmartPtr<CplUserData<TDataScale,dim> > , number)>(&T::add))
			.add_method("add", static_cast<void (T::*)(number,number)>(&T::add))
			.add_constructor()
			.template add_constructor<void (*)(const ScaleAddLinker<TData, dim, TDataScale,number>&)>()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, string("ScaleAddLinkerVectorVector"), dimTag);
	}

//	DarcyVelocityLinker
	{
		typedef DarcyVelocityLinker<dim> T;
		typedef DependentUserData<MathVector<dim>, dim> TBase;
		string name = string("DarcyVelocityLinker").append(dimSuffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_method("set_density", &T::set_density)
			.add_method("set_gravity", &T::set_gravity)
			.add_method("set_permeability", static_cast<void (T::*)(number)>(&T::set_permeability))
			.add_method("set_permeability", static_cast<void (T::*)(SmartPtr<CplUserData<MathMatrix<dim,dim>,dim> >)>(&T::set_permeability))
			.add_method("set_pressure_gradient", &T::set_pressure_gradient)
			.add_method("set_viscosity", static_cast<void (T::*)(number)>(&T::set_viscosity))
			.add_method("set_viscosity", static_cast<void (T::*)(SmartPtr<CplUserData<number,dim> >)>(&T::set_viscosity))
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "DarcyVelocityLinker", dimTag);
	}

//	InverseLinker"Type"
	{
		typedef InverseLinker<dim> T;
		typedef DependentUserData<number,dim> TBase;

			string name = string("InverseLinker").append(dimSuffix);
			reg.add_class_<T,TBase>(name, grp)
			.add_method("divide", static_cast<void (T::*)(SmartPtr<CplUserData<number,dim> > , SmartPtr<CplUserData<number,dim> >)>(&T::divide))
			.add_method("divide", static_cast<void (T::*)(number , SmartPtr<CplUserData<number,dim> >)>(&T::divide))
			.add_method("divide", static_cast<void (T::*)(SmartPtr<CplUserData<number,dim> > , number)>(&T::divide))
			.add_method("divide", static_cast<void (T::*)(number,number)>(&T::divide))
			.add_constructor()
			.template add_constructor<void (*)(const InverseLinker<dim>&)>()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, string("InverseLinker"), dimTag);

	}

}

/**
 * Function called for the registration of Domain and Algebra independent parts.
 * All Functions and Classes not depending on Domain and Algebra
 * are to be placed here when registering.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
static void Common(Registry& reg, string grp)
{
	reg.add_class_<IFunction<number> >("IFunctionNumber", grp);

//	UserDataInfo
	{
		reg.add_class_<UserDataInfo>("UserDataInfo", grp);
	}

	{
		typedef CplUserData<number, 2> TBase1;
		reg.add_class_<RotatingCone2d, TBase1>("RotatingCone2d", grp)
				.add_constructor<void (*)(double,double,double,double,double,double,double)>()
				.set_construct_as_smart_pointer(true);
		typedef CplUserData<MathVector<2>, 2> TBase2;
		reg.add_class_<RotatingVelocity2d, TBase2>("RotatingVelocity2d", grp)
				.add_constructor<void (*)(double,double,double)>()
				.set_construct_as_smart_pointer(true);
	}
}

}; // end Functionality

// end group userdata_bridge
/// \}

}// end UserData

/// \addtogroup userdata_bridge
void RegisterBridge_UserData(Registry& reg, string grp)
{
//	get group string
	grp.append("/Discretization/SpatialDisc/UserData");
	typedef UserDataBridge::Functionality Functionality;

	try{
		RegisterCommon<Functionality>(reg,grp);
		RegisterDimensionDependent<Functionality>(reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

} // end namespace
} // end namepace
