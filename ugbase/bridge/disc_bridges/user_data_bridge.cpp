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
#include "lib_disc/spatial_disc/user_data/common_user_data/common_user_data.h"

#include "lib_disc/spatial_disc/user_data/common_user_data/rotating_cone.h"

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
			.template add_constructor<void (*)(const std::vector<number>&)>("Values")
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
			.add_method("set_gravity", &T::set_gravity)
			.add_method("set_permeability", static_cast<void (T::*)(number)>(&T::set_permeability))
			.add_method("set_permeability", static_cast<void (T::*)(SmartPtr<CplUserData<MathMatrix<dim,dim>,dim> >)>(&T::set_permeability))
			.add_method("set_pressure_gradient", &T::set_pressure_gradient)
			.add_method("set_viscosity", static_cast<void (T::*)(number)>(&T::set_viscosity))
			.add_method("set_viscosity", static_cast<void (T::*)(SmartPtr<CplUserData<number,dim> >)>(&T::set_viscosity))
			.add_method("set_density", static_cast<void (T::*)(number)>(&T::set_density))
			.add_method("set_density", static_cast<void (T::*)(SmartPtr<CplUserData<number,dim> >)>(&T::set_density))
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


	{
		typedef ug::LognormalRandomField<MathMatrix<dim, dim>, dim> T;
		typedef CplUserData<MathMatrix<dim, dim>, dim> TBase;
		string name = string("LognormalRandomField").append(dimSuffix);;
		reg.add_class_<T, TBase>(name, grp)
			//.add_constructor()
			.template add_constructor<void (*)()>("LognormalRandomField")
			.template add_constructor<void (*)(size_t N, double mean_f, double sigma_f, double sigma)>("LognormalRandomField", "N#mean_f#sigma_f#sigma")
			.add_method("set_config",  &T::set_config, "", "N#mean_f#sigma_f#sigma")
			.add_method("set_no_exp",  &T::set_no_exp, "", "", "use this for display of log of the field")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, string("LognormalRandomField"), dimTag);
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

#ifdef UG_DIM_2
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
#endif
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
