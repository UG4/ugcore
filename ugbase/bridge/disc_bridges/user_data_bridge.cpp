/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel
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

#include <iostream>
#include <sstream>
#include <boost/function.hpp>
#include <cmath>

// include bridge
#include "bridge/bridge.h"
#include "bridge/util.h"
#include "bridge/util_domain_dependent.h"

// common
#include "common/common.h"

// lib disc
#include "lib_disc/spatial_disc/user_data/const_user_data.h"
#include "lib_disc/spatial_disc/user_data/std_glob_pos_data.h"
#include "lib_disc/spatial_disc/user_data/common_user_data/common_user_data.h"
#include "lib_disc/spatial_disc/user_data/common_user_data/raster_user_data.h"


#include "lib_disc/spatial_disc/user_data/linker/linker.h"
#include "lib_disc/spatial_disc/user_data/linker/scale_add_linker.h"
#include "lib_disc/spatial_disc/user_data/linker/inverse_linker.h"
#include "lib_disc/spatial_disc/user_data/linker/darcy_velocity_linker.h"
#include "lib_disc/spatial_disc/user_data/linker/bingham_viscosity_linker.h"
#include "lib_disc/spatial_disc/user_data/linker/projection_linker.h"
#include "lib_disc/spatial_disc/user_data/linker/adapter.h"
#include "lib_disc/spatial_disc/user_data/linker/interval_linker.h"
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


namespace bridge{

/**
 * \defgroup userdata_bridge User Data Bridge
 * \ingroup disc_bridge
 * \{
 */

template <typename TData, int dim, typename TTraits=user_data_traits<TData> >
void RegisterUserDataTypeA(Registry& reg, string grp)
{
	string dimSuffix = GetDimensionSuffix<dim>();
	string dimTag = GetDimensionTag<dim>();

	string type = TTraits::name();

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


}

template <typename TData, int dim, typename TTraits=user_data_traits<TData> >
void RegisterUserDataTypeB(Registry& reg, string grp)
{
	string dimSuffix = GetDimensionSuffix<dim>();
	string dimTag = GetDimensionTag<dim>();

	string type = TTraits::name();

	//	ScaleAddLinker"Type"
	{
			typedef ScaleAddLinker<TData, dim, number> T;
			typedef DependentUserData<TData, dim> TBase;
			string name = string("ScaleAddLinker").append(type).append(dimSuffix);
			reg.add_class_<T, TBase>(name, grp)
				.add_method("add", static_cast<void (T::*)(SmartPtr<CplUserData<number,dim> > , SmartPtr<CplUserData<TData,dim> >)>(&T::add))
				.add_method("add", static_cast<void (T::*)(number, SmartPtr<CplUserData<TData,dim> >)>(&T::add))
				.add_method("add", static_cast<void (T::*)(SmartPtr<CplUserData<number,dim> > , number)>(&T::add))
				.add_method("add", static_cast<void (T::*)(number,number)>(&T::add))
				.add_constructor()
				.template add_constructor<void (*)(const ScaleAddLinker<TData, dim, number>&)>()
				.set_construct_as_smart_pointer(true);
			reg.add_class_to_group(name, string("ScaleAddLinker").append(type), dimTag);
	}

}

template <typename TData, int dim, typename TTraits=user_data_traits<TData> >
void RegisterUserDataType(Registry& reg, string grp)
{
	RegisterUserDataTypeA<TData, dim, TTraits>(reg, grp);
	RegisterUserDataTypeB<TData, dim, TTraits>(reg, grp);
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
	/*
	typedef std::tuple<number, number> MathPair;
		typedef std::tuple<number, number, number> MathTriple;
	template <>
		struct user_data_traits<MathPair>{static std::string name() 	{return "Pair";}};*/
template <int dim>
static void Dimension(Registry& reg, string grp)
{
	string dimSuffix = GetDimensionSuffix<dim>();
	string dimTag = GetDimensionTag<dim>();

	RegisterUserDataType<number, dim>(reg, grp);
	RegisterUserDataType<MathVector<dim>, dim>(reg, grp);
	RegisterUserDataType<MathMatrix<dim,dim>, dim>(reg, grp);
	RegisterUserDataType<MathTensor<4,dim>, dim>(reg, grp);




	// RegisterUserDataType<MathPair, dim>(reg, grp);

	if (dim!=2)
	{
		// Register pair (corresponds to vector for dim==2)
		struct pair_traits {static std::string name() 	{return "Pair";}};
		RegisterUserDataTypeA<MathVector<2>, dim, pair_traits>(reg, grp); // Pair
	}
/*
	if (dim!=3)
	{
		// Register triple (corresponds to vector for dim==3)
		struct triple_traits {static std::string name() {return "Triple";}};
		RegisterUserDataType<MathVector<3>, dim, triple_traits, false>(reg, grp);
	}

*/

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
			.add_method("get", &T::get)
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

	// UserVectorEntryAdapter
	{
		string name = string("UserVectorEntryAdapter").append(dimSuffix);
		typedef UserVectorEntryAdapter<dim> T;
		typedef CplUserData<number, dim> TCplData;
		typedef typename T::encapsulated_type TEncaps;

		reg.template add_class_<T, TCplData>(name, grp)
		   .template add_constructor<void (*)() >("")
			.add_method("set_vector", &T::set_vector)
			.set_construct_as_smart_pointer(true);

		reg.add_class_to_group(name, "UserVectorEntryAdapter", dimTag);

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
#ifdef UG_FOR_LUA
			.add_method("set_permeability", static_cast<void (T::*)(const char* fctName)>(&T::set_permeability))
			.add_method("set_permeability", static_cast<void (T::*)(LuaFunctionHandle fct)>(&T::set_permeability))
#endif
			.add_method("set_pressure_gradient", &T::set_pressure_gradient)
			.add_method("set_viscosity", static_cast<void (T::*)(number)>(&T::set_viscosity))
			.add_method("set_viscosity", static_cast<void (T::*)(SmartPtr<CplUserData<number,dim> >)>(&T::set_viscosity))
			.add_method("set_density", static_cast<void (T::*)(number)>(&T::set_density))
			.add_method("set_density", static_cast<void (T::*)(SmartPtr<CplUserData<number,dim> >)>(&T::set_density))
			.add_method("set_derivative_mask",  &T::set_derivative_mask)
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "DarcyVelocityLinker", dimTag);
	}

//	BinghamViscosityLinker
	{
		typedef BinghamViscosityLinker<dim> T;
		typedef DependentUserData<number, dim> TBase;
		string name = string("BinghamViscosityLinker").append(dimSuffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_method("set_velocity_gradient", &T::set_velocity_gradient)
			.add_method("set_yield_stress", static_cast<void (T::*)(number)>(&T::set_yield_stress))
			.add_method("set_yield_stress", static_cast<void (T::*)(SmartPtr<CplUserData<number,dim> >)>(&T::set_yield_stress))
			.add_method("set_viscosity", static_cast<void (T::*)(number)>(&T::set_viscosity))
			.add_method("set_viscosity", static_cast<void (T::*)(SmartPtr<CplUserData<number,dim> >)>(&T::set_viscosity))
			.add_method("set_density", static_cast<void (T::*)(number)>(&T::set_density))
			.add_method("set_density", static_cast<void (T::*)(SmartPtr<CplUserData<number,dim> >)>(&T::set_density))
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "BinghamViscosityLinker", dimTag);
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

//	ProjectionLinker
	{
		typedef ProjectionLinker<dim> T;
		typedef DependentUserData<MathVector<dim>, dim> TBase;
		string name = string("ProjectionLinker").append(dimSuffix);
		reg.add_class_<T, TBase>(name, grp)
		   .template add_constructor<void (*) (SmartPtr<CplUserData<MathVector<dim>, dim> >)>()
		   .set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "ProjectionLinker", dimTag);
	}

//	LognormalRandomField
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

//	Inverse-distance-weighting interpolation
	{
		typedef IDWUserData<dim, number> T;
		typedef CplUserData<number, dim> TBase;
		string name = string("IDWUserData").append(dimSuffix);
		reg.add_class_<T, TBase>(name, grp)
		   .add_method("load_data_from", static_cast<void (T::*)(const char*)>(&T::load_data_from), "loads data from a file", "file name")
		   .add_method("set_order", static_cast<void (T::*)(number)>(&T::set_order), "sets order of the IDW-interpolation", "order")
		   .add_method("set_radius", static_cast<void (T::*)(number)>(&T::set_radius), "sets radius of the neighbourhood for the IDW-interpolation", "radius")
		   .add_constructor()
		   .template add_constructor<void(*)(number,number)> ("order#radius")
		   .set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "IDWUserData", dimTag);
	}

	// Composite user data (number)
	{
			string name = string("CompositeUserNumber").append(dimSuffix);
			typedef CompositeUserData<number, dim, void> T;
			reg.add_class_<T,typename T::base_type>(name, grp)
				.template add_constructor<void (*)()>()
				.template add_constructor<void (*)(bool)>("continuous")
				.add_method("add", static_cast<void (T::*)(int, typename T::ref_type)>(&T::add), "assign a user data object to a subset index", "si#userdata")
				.add_method("add", static_cast<void (T::*)(ConstSmartPtr<ISubsetHandler>, const char *, typename T::ref_type)>(&T::add), "assign a user data object to subsets by names", "names#userdata")
				.add_method("has", &T::has)
				.add_method("get", &T::get)
				.add_method("is_coupled", &T::is_coupled)
				.add_method("get_coupled", &T::get_coupled)
				.set_construct_as_smart_pointer(true);

			reg.add_class_to_group(name, "CompositeUserNumber", dimTag);
	}

	// Composite user data (vector)
	{
			string name = string("CompositeUserVector").append(dimSuffix);
			typedef CompositeUserData<MathVector<dim>, dim, void> T;
			reg.add_class_<T,typename T::base_type>(name, grp)
				.template add_constructor<void (*)()>()
				.template add_constructor<void (*)(bool)>("continuous")
				.add_method("add", static_cast<void (T::*)(int, typename T::ref_type)>(&T::add), "assign a user data object to a subset index", "si#userdata")
				.add_method("add", static_cast<void (T::*)(ConstSmartPtr<ISubsetHandler>, const char *, typename T::ref_type)>(&T::add), "assign a user data object to subsets by names", "names#userdata")
				.add_method("has", &T::has)
				.add_method("get", &T::get)
				.add_method("is_coupled", &T::is_coupled)
				.add_method("get_coupled", &T::get_coupled)
				.set_construct_as_smart_pointer(true);

			reg.add_class_to_group(name, "CompositeUserVector", dimTag);
	}

//	Interval filter linker
	{
		typedef IntervalNumberLinker<dim> T;
		typedef DependentUserData<number, dim> TBase;
		string name = string("IntervalNumberLinker").append(dimSuffix);
		reg.add_class_<T, TBase>(name, grp)
		   .add_method("set_default", static_cast<void (T::*)(number)>(&T::set_default), "sets the value out of the interval", "value")
		   .template add_constructor<void (*) (SmartPtr<CplUserData<number, dim> >, MathVector<dim>&, MathVector<dim>&)>("data#left#right")
		   .template add_constructor<void (*) (SmartPtr<CplUserData<number, dim> >, std::vector<number>, std::vector<number>)>("data#left#right")
		   .set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "IntervalNumberLinker", dimTag);
	}

//	GlobAttachmentElementUserData
	{
		typedef GlobAttachmentElementUserData<dim, number> T;
		typedef CplUserData<number, dim> TBase;
		string name = string("GlobAttachmentElementNumberData").append(dimSuffix);
		reg.add_class_<T, TBase>(name, grp)
			.template add_constructor<void (*)(SmartPtr<Grid>, const char *) >("AttachmentName")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "GlobAttachmentElementNumberData", dimTag);
	}
	
}

/**
 * Function called for the registration of Domain dependent parts.
 * All Functions and Classes depending on the Domain
 * are to be placed here when registering. The method is called for all
 * available Domain types, based on the current build options.
 *
 * @param reg	registry
 * @param grp	group for sorting of functionality
 */
template <typename TDomain>
static void Domain(Registry& reg, string grp)
{
	static const int dim = TDomain::dim;
	
	string suffix = GetDomainSuffix<TDomain>();
	string tag = GetDomainTag<TDomain>();
	
//	User data of a subset indicator (1 in the subset, 0 everywhere else)
	{
		string name = string("SubsetIndicatorUserData").append(suffix);
		typedef SubsetIndicatorUserData<TDomain> T;
		typedef UserData<number, dim> TBase;
		
		reg.add_class_<T, TBase> (name, grp)
			.template add_constructor<void (*)(ConstSmartPtr<TDomain>, const char*)>("Domain#Subsets")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "SubsetIndicatorUserData", tag);
	}

//	User data of a value indicator (1 in the subset, 0 everywhere else)
	{
		string name = string("ValueIndicatorUserData").append(suffix);
		typedef ValueIndicatorUserData<TDomain> T;
		typedef UserData<number, dim> TBase;
		
		reg.add_class_<T, TBase> (name, grp)
			.template add_constructor<void (*)(SmartPtr<TBase>, number, bool)>("Domain#threshold#greater")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "ValueIndicatorUserData", tag);
	}

// EdgeOrientation (vector)
	{
		string name = string("EdgeOrientation").append(suffix);
		typedef EdgeOrientation<TDomain> T;
		typedef CplUserData<MathVector<dim>, dim> TBase;

		reg.add_class_<T,TBase>(name, grp)
		   .template add_constructor<void (*)(SmartPtr<TDomain> ) >("")
		   .set_construct_as_smart_pointer(true);

		reg.add_class_to_group(name, "EdgeOrientation", tag);

	}
	
//	User data for evaluation of full-dimensional vector fields on hypersurfaces
	{
		string name = string("OutNormCmp").append(suffix);
		typedef OutNormCmp<TDomain> T;
		typedef UserData<MathVector<dim>, dim> TBase;
		
		reg.add_class_<T, TBase> (name, grp)
			.template add_constructor<void (*)(SmartPtr<TDomain>, SmartPtr<TBase>, const char*)>("Domain#Data#Subsets")
			.template add_constructor<void (*)(SmartPtr<TDomain>, SmartPtr<TBase>)>("Domain#Data")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "OutNormCmp", tag);
	}
   

	{ // Raster
		string name = string("HomoRasterUserData").append(suffix);
		typedef HomoRasterUserData<TDomain> T;
		typedef UserData<number, dim> TBase;
		
		reg.add_class_<T, TBase> (name, grp)
			.template add_constructor<void (*)(typename T::input_type)>("Raster")
			.add_method("set_order", &T::set_order)
			.add_method("set_domain", &T::set_domain)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "HomoRasterUserData", tag);
	
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
		reg.add_class_<UserDataInfo>("UserDataInfo", grp)
			.add_method("set_obj_name", &UserDataInfo::set_obj_name)
			.add_method("obj_name", &UserDataInfo::obj_name);
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


template <class TValue, int dim>
static void RegisterRasterUserData(Registry& reg, string name, string grp)
{

	string suffix = GetDimensionSuffix<dim>();
	string tag = GetDimensionTag<dim>();

	{ // Raster
		typedef RasterUserData<dim> T;
		typedef typename T::base_type TBase;
		string fullName = name + suffix;

		reg.add_class_<T,TBase>(fullName, grp)
			.template add_constructor<void (*)(typename T::input_type)>("RasterNumberData")
			.add_method("set_order", &T::set_order)
			.add_method("set_scale", &T::set_scale)
			.set_construct_as_smart_pointer(true);

		reg.add_class_to_group(fullName, name, tag);
	}

}




/// \addtogroup userdata_bridge
void RegisterBridge_UserData(Registry& reg, string grp)
{
//	get group string
	grp.append("/Discretization/SpatialDisc/UserData");
	typedef UserDataBridge::Functionality Functionality;

	try{
		RegisterCommon<Functionality>(reg,grp);
		RegisterDimensionDependent<Functionality>(reg,grp);
		RegisterDomainDependent<Functionality>(reg,grp);

#ifdef UG_DIM_2 // only for 2D/3D
		RegisterRasterUserData<number, 2>(reg, "RasterNumberData", grp);
#endif

#ifdef UG_DIM_3
		RegisterRasterUserData<number, 3>(reg, "RasterNumberData", grp);
#endif
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

} // end namespace
} // end namepace
