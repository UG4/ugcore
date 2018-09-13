/*
 * Copyright (c) 2014-2017:  G-CSC, Goethe University Frankfurt
 * Author: Arne Naegel
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

#ifndef __H__UG__LIB_DISC__FUNCTION_SPACE__METRIC_SPACES_H_
#define __H__UG__LIB_DISC__FUNCTION_SPACE__METRIC_SPACES_H_

// C++ includes
#include <vector>
#include <cmath>

// UG includes
#include "common/common.h"
#include "common/util/smart_pointer.h"

#include "lib_algebra/lib_algebra.h"
#include "lib_algebra/operator/debug_writer.h"

#include "lib_disc/function_spaces/integrate.h"
#include "lib_disc/function_spaces/grid_function.h"
#include "lib_disc/function_spaces/grid_function_user_data.h"
#include "lib_disc/function_spaces/integrate.h"

#include "lib_disc/spatial_disc/user_data/linker/darcy_velocity_linker.h"  // required for Darcy velocity linker

namespace ug {
/// Abstract base class for (algebraic) vectors
template<typename TVector>
class IBanachSpace
{

public:
	virtual ~IBanachSpace() {}

	/// euclidean norm (default)
	virtual double norm(TVector &x)
	{ return x.norm(); }

	virtual double distance(TVector& x, TVector& y)
	{
		SmartPtr<TVector> delta = x.clone();
		*delta -= y;
		return norm(*delta);
	}
};



/// Abstract base class for grid functions
template<typename TGridFunction>
class IGridFunctionSpace  :
		public IBanachSpace<typename TGridFunction::vector_type>
{

public:
	typedef typename TGridFunction::vector_type vector_type;
	typedef TGridFunction grid_function_type;

	/// DTOR
	virtual ~IGridFunctionSpace() {}

	/// norm (for grid functions)
	virtual double norm(TGridFunction& x)
	{ return x.norm(); }

	/// distance (for grid functions)
	virtual double distance(TGridFunction& x, TGridFunction& y)
	{
		SmartPtr<TGridFunction> delta = x.clone();
		*delta -= y;
		return norm(*delta);
	}

	/// norm (for vectors)
	virtual double norm(vector_type &x)
	{
		TGridFunction* gfX=dynamic_cast< TGridFunction*>(&x);
		UG_ASSERT(gfX!=NULL, "Huhh: GridFunction required!");
		return norm(*gfX);
	}

	/// distance (for vectors)
	virtual double distance(vector_type &x, vector_type &y)
	{ return distance(static_cast<TGridFunction &>(x), static_cast<TGridFunction &>(y)); }

	virtual double scaling() const
	{ return 1.0; }

	virtual std::string config_string() const
	 { return std::string("IGridFunctionSpace"); }
};



/** Auxiliary class for providing weights*/
template <typename W>
class IObjectWithWeights
{
public:
	typedef W weight_type;

	IObjectWithWeights()
	: m_spWeight(SPNULL) {}

	IObjectWithWeights(ConstSmartPtr<weight_type> spW)
	: m_spWeight(spW) {}

	/// for weighted norms
	void set_weight(ConstSmartPtr<weight_type> spWeight)
	{ m_spWeight = spWeight; }

	ConstSmartPtr<weight_type> get_weight()
	{ return m_spWeight; }

protected:
	ConstSmartPtr<weight_type> m_spWeight;
};


/** Auxiliary class for time dependence - SHOULD be replaced by product space*/
/*
class ITimeData {

public:
	virtual ~ITimeData() {};

	/// characteristic time
	virtual void update_time_data(number dt) = 0;

	virtual bool is_time_dependent() const
	{ return false; }

};

*/

//! Estimate the error (based on the difference between two grid functions)
template <typename TGridFunction>
class IComponentSpace :
		public IGridFunctionSpace<TGridFunction>
{
protected:
	std::string m_fctNames;
	const char* m_ssNames;
	int m_quadorder;

public:
	typedef IGridFunctionSpace<TGridFunction> base_type;
	static const int dim=TGridFunction::dim;

	IComponentSpace(const char *fctNames)
	: m_fctNames(fctNames), m_ssNames(NULL), m_quadorder(3){}

	IComponentSpace(const char *fctNames, int order)
	: m_fctNames(fctNames), m_ssNames(NULL), m_quadorder(order) {}

	/*IComponentSpace(const char *fctNames, int order, number scale)
	: m_fctNames(fctNames), m_ssNames(NULL), m_quadorder(order), m_scale(scale) {}
*/
	IComponentSpace(const char *fctNames, const char* ssNames, int order)
	: m_fctNames(fctNames), m_ssNames(ssNames), m_quadorder(order) {}

	virtual ~IComponentSpace() {};

	// per convention, norm must return sqrt of norm2
	virtual double norm(TGridFunction& uFine)
	{ return sqrt(norm2(uFine)); }

	// per convention, distance must return sqrt of distance2
	virtual double distance(TGridFunction& uFine, TGridFunction& uCoarse)
	{ return sqrt(distance2(uFine, uCoarse)); }

	virtual double norm2(TGridFunction& uFine)  = 0;
	virtual double distance2(TGridFunction& uFine, TGridFunction& uCoarse)  = 0;


public:
	/// print config string
	virtual std::string config_string() const
	{
		std::stringstream ss;

		if (this->m_ssNames)
		ss << this->m_fctNames << ", " << this->m_ssNames << ", " << this->m_quadorder
			<< ", type=" << std::endl;

		else
			ss << this->m_fctNames << ",  (no name), " << this->m_quadorder
						<< ", type=" << std::endl;
		return ss.str();
	}

};



//! Wrapper class for time dependence.
template <typename TGridFunction>
class TimeDependentSpace :
		public IGridFunctionSpace<TGridFunction>
{
public:
	typedef IGridFunctionSpace<TGridFunction> base_type;
	typedef IComponentSpace<TGridFunction> comp_space_type;

	/// time dependent CTOR
	TimeDependentSpace(SmartPtr<comp_space_type> spSpace, number tScale)
	: m_spSpatialSpace(spSpace), m_tScale(tScale) {};

	/// DTOR
	virtual ~TimeDependentSpace() {};

	using base_type::norm;
	using base_type::distance;
	using base_type::scaling;

	/// scaling (OVERRIDE)
	double scaling() const
	{ return (m_spSpatialSpace->scaling()*m_tScale); }

	/// characteristic time
	void update_time_data(number tScale)
	{  m_tScale = tScale; }

	/// print config string
	std::string config_string() const
	{
		std::stringstream ss;
		ss << "TimeDependentSpace for " <<  std::endl;
		ss <<  m_spSpatialSpace->config_string() << std::endl;
		return ss.str();
	}

protected:
	SmartPtr<comp_space_type> m_spSpatialSpace;
	number m_tScale;
};


/** Evaluates difference between two grid functions in L2 norm */
template <typename TGridFunction>
class L2ComponentSpace :
		public IComponentSpace<TGridFunction>,
		public IObjectWithWeights<typename L2DistIntegrand<TGridFunction>::weight_type >
{
public:
	typedef IComponentSpace<TGridFunction> base_type;
	typedef typename L2DistIntegrand<TGridFunction>::weight_type weight_type;
	typedef IObjectWithWeights<weight_type> weighted_obj_type;

	/// CTOR
	L2ComponentSpace(const char *fctNames)
	: base_type(fctNames), weighted_obj_type(make_sp(new ConstUserNumber<TGridFunction::dim>(1.0))) {};

	L2ComponentSpace(const char *fctNames, int order)
	: base_type(fctNames, order), weighted_obj_type(make_sp(new ConstUserNumber<TGridFunction::dim>(1.0))) {};

	L2ComponentSpace(const char *fctNames,  int order, double weight, const char* ssNames=0)
	: base_type(fctNames, ssNames, order), weighted_obj_type(make_sp(new ConstUserNumber<TGridFunction::dim>(weight))) {};

	L2ComponentSpace(const char *fctNames, int order, ConstSmartPtr<weight_type> spWeight, const char* ssNames=0)
	: base_type(fctNames, ssNames, order), weighted_obj_type(spWeight) {};

	/// DTOR
	~L2ComponentSpace() {};

	using IComponentSpace<TGridFunction>::norm;
	using IComponentSpace<TGridFunction>::distance;

	/// for weighted norms
	using weighted_obj_type::set_weight;
	using weighted_obj_type::get_weight;
	using weighted_obj_type::m_spWeight;

	/// \copydoc IComponentSpace<TGridFunction>::norm
	double norm2(TGridFunction& uFine)
	{ return L2Norm2(uFine, base_type::m_fctNames.c_str(), base_type::m_quadorder, base_type::m_ssNames, weighted_obj_type::m_spWeight); }

	/// \copydoc IComponentSpace<TGridFunction>::distance
	double distance2(TGridFunction& uFine, TGridFunction& uCoarse)
	{ return L2Distance2(uFine, base_type::m_fctNames.c_str(), uCoarse, base_type::m_fctNames.c_str(),
		base_type::m_quadorder, base_type::m_ssNames, weighted_obj_type::m_spWeight);}



};

/** Evaluates difference between two grid functions in L2 norm */
template <typename TGridFunction>
class L2QuotientSpace :
		public IComponentSpace<TGridFunction>,
		public IObjectWithWeights<typename L2DistIntegrand<TGridFunction>::weight_type >
{
public:
	typedef IComponentSpace<TGridFunction> base_type;
	typedef typename L2DistIntegrand<TGridFunction>::weight_type weight_type;
	typedef IObjectWithWeights<weight_type> weighted_obj_type;

	/// CTOR
	L2QuotientSpace(const char *fctNames)
	: base_type(fctNames), weighted_obj_type(make_sp(new ConstUserNumber<TGridFunction::dim>(1.0))) {};

	L2QuotientSpace(const char *fctNames, int order)
	: base_type(fctNames, order), weighted_obj_type(make_sp(new ConstUserNumber<TGridFunction::dim>(1.0))) {};

	L2QuotientSpace(const char *fctNames,  int order, double weight, const char* ssNames=0)
	: base_type(fctNames, ssNames, order), weighted_obj_type(make_sp(new ConstUserNumber<TGridFunction::dim>(weight))) {};

	L2QuotientSpace(const char *fctNames, int order, ConstSmartPtr<weight_type> spWeight, const char* ssNames=0)
	: base_type(fctNames, ssNames, order), weighted_obj_type(spWeight) {};

	/// DTOR
	~L2QuotientSpace() {};

	using IComponentSpace<TGridFunction>::norm;
	using IComponentSpace<TGridFunction>::distance;

	/// for weighted norms
	using weighted_obj_type::set_weight;
	using weighted_obj_type::get_weight;
	using weighted_obj_type::m_spWeight;

	/// \copydoc IComponentSpace<TGridFunction>::norm
	double norm2(TGridFunction& u)
	{
		typedef ConstUserNumber<TGridFunction::dim> MyConstUserData;

		SmartPtr<MyConstUserData> spConst= make_sp(new MyConstUserData(1.0));
		number Meas = Integral(spConst, u, base_type::m_ssNames, 0.0, 1, "best");
		number uAvg = Integral(u,  base_type::m_fctNames.c_str(), base_type::m_ssNames, base_type::m_quadorder);

		std::cerr << "Average:=" << uAvg <<"/" << Meas << " = " << uAvg/Meas << std::endl;
		SmartPtr<MyConstUserData> spAvg = make_sp(new MyConstUserData(uAvg/Meas));
		double qnorm = L2Error(spAvg, u, base_type::m_fctNames.c_str(), 0.0, base_type::m_quadorder, base_type::m_ssNames);

		return qnorm*qnorm;
	}

	/// \copydoc IComponentSpace<TGridFunction>::distance
	double distance2(TGridFunction& uFine, TGridFunction& uCoarse)
	{
		typedef ConstUserNumber<TGridFunction::dim> MyConstUserData;

		SmartPtr<MyConstUserData> spConst= make_sp(new MyConstUserData(1.0));
		number Meas = Integral(spConst, uFine, base_type::m_ssNames, 0.0, 1, "best");
		number avgFine = Integral(uFine, base_type::m_fctNames.c_str(), base_type::m_ssNames, base_type::m_quadorder);
		number avgCoarse = Integral(uCoarse, base_type::m_fctNames.c_str(), base_type::m_ssNames, base_type::m_quadorder);

		std::cerr << "Average:=(" << avgFine << "-" << avgCoarse<<")/" << Meas << " = " << (avgFine-avgCoarse)/Meas << std::endl;
		return L2Distance2(uFine, base_type::m_fctNames.c_str(),
						  uCoarse, base_type::m_fctNames.c_str(),
						  base_type::m_quadorder, base_type::m_ssNames, weighted_obj_type::m_spWeight,
						  (avgFine-avgCoarse)/Meas);
	}



};


/** Evaluates distance between two grid functions in H1 semi-norm */
template <typename TGridFunction>
class H1SemiComponentSpace
: public IComponentSpace<TGridFunction>,
  public IObjectWithWeights<typename H1SemiDistIntegrand<TGridFunction>::weight_type >
{
public:
	typedef IComponentSpace<TGridFunction> base_type;
	typedef typename H1SemiDistIntegrand<TGridFunction>::weight_type weight_type;
	typedef IObjectWithWeights<weight_type> weighted_obj_type;


	H1SemiComponentSpace(const char *fctNames)
	: base_type(fctNames), weighted_obj_type(make_sp(new ConstUserMatrix<TGridFunction::dim>(1.0))) {};

	H1SemiComponentSpace(const char *fctNames, int order)
	: base_type(fctNames, order), weighted_obj_type(make_sp(new ConstUserMatrix<TGridFunction::dim>(1.0))) {};

	H1SemiComponentSpace(const char *fctNames, int order, number weight, const char* ssNames=0)
	: base_type(fctNames, ssNames, order), weighted_obj_type(make_sp(new ConstUserMatrix<TGridFunction::dim>(weight)))  {};

	H1SemiComponentSpace(const char *fctNames, int order, ConstSmartPtr<weight_type> spWeight, const char* ssNames=0)
	: base_type(fctNames, ssNames, order), weighted_obj_type(spWeight) {};

	H1SemiComponentSpace(const char *fctNames, int order, const char* ssNames, ConstSmartPtr<weight_type> spWeight)
	: base_type(fctNames, ssNames, order), weighted_obj_type(spWeight) {};

	/// DTOR
	~H1SemiComponentSpace() {};


	/// \copydoc IComponentSpace<TGridFunction>::norm
	/// \copydoc IComponentSpace<TGridFunction>::distance
	using IComponentSpace<TGridFunction>::norm;
	using IComponentSpace<TGridFunction>::distance;

	/// \copydoc IComponentSpace<TGridFunction>::norm2
	double norm2(TGridFunction& uFine)
	{ return H1SemiNorm2<TGridFunction>(uFine, base_type::m_fctNames.c_str(), base_type::m_quadorder, NULL, weighted_obj_type::m_spWeight); }

	/// \copydoc IComponentSpace<TGridFunction>::distance2
	double distance2(TGridFunction& uFine, TGridFunction& uCoarse)
	{ return H1SemiDistance2<TGridFunction>(uFine, base_type::m_fctNames.c_str(), uCoarse, base_type::m_fctNames.c_str(), base_type::m_quadorder, m_spWeight); }

	/// for weighted norms
	using weighted_obj_type::set_weight;
	using weighted_obj_type::get_weight;
	using weighted_obj_type::m_spWeight;

};



/*
//! Evaluates distance w.r.t kinetic energy
template <typename TGridFunction>
class KineticEnergyComponentSpace
: public IComponentSpace<TGridFunction>,
  public IObjectWithWeights<typename L2Integrand<TGridFunction>::weight_type >
{
public:
	typedef IComponentSpace<TGridFunction> base_type;
	typedef typename L2Integrand<TGridFunction>::weight_type weight_type;
	typedef IObjectWithWeights<weight_type> weighted_obj_type;
	typedef DarcyVelocityLinker<TGridFunction::dim> velocity_type;


	KineticEnergyComponentSpace(SmartPtr<velocity_type> spVelocity, const char *fctNames)
	: base_type(fctNames), weighted_obj_type(make_sp(new ConstUserNumber<TGridFunction::dim>(1.0))), m_spVelocity(spVelocity) {};

	KineticEnergyComponentSpace(SmartPtr<velocity_type> spVelocity, const char *fctNames, int order)
	: base_type(fctNames, order), weighted_obj_type(make_sp(new ConstUserNumber<TGridFunction::dim>(1.0))), m_spVelocity(spVelocity)  {};

	KineticEnergyComponentSpace(SmartPtr<velocity_type> spVelocity, const char *fctNames, int order, number weight, const char* ssNames=0)
	: base_type(fctNames, ssNames, order), weighted_obj_type(make_sp(new ConstUserNumber<TGridFunction::dim>(weight))), m_spVelocity(spVelocity)   {};

	KineticEnergyComponentSpace(SmartPtr<velocity_type> spVelocity, const char *fctNames, int order, ConstSmartPtr<weight_type> spWeight, const char* ssNames=0)
	: base_type(fctNames, ssNames, order), weighted_obj_type(spWeight), m_spVelocity(spVelocity)  {};

	KineticEnergyComponentSpace(SmartPtr<velocity_type> spVelocity, const char *fctNames, int order, const char* ssNames, ConstSmartPtr<weight_type> spWeight)
	: base_type(fctNames, ssNames, order), weighted_obj_type(spWeight), m_spVelocity(spVelocity) {};

	/// DTOR
	~KineticEnergyComponentSpace() {};


	using IComponentSpace<TGridFunction>::m_quadorder;
	using IComponentSpace<TGridFunction>::norm;
	using IComponentSpace<TGridFunction>::distance;

	/// \copydoc IComponentSpace<TGridFunction>::norm2
	double norm2(TGridFunction& u)
	{
		const char* subsets = NULL;
		UserDataIntegrandSq<MathVector<TGridFunction::dim>, TGridFunction> integrand2(m_spVelocity, &u, 0.0);
		return IntegrateSubsets(integrand2, u, subsets, m_quadorder);
	}

	/// \copydoc IComponentSpace<TGridFunction>::distance2
	double distance2(TGridFunction& uFine, TGridFunction& uCoarse)
	{
		// UG_THROW("Not implemented!");
		const char* subsets = NULL;
		UserDataDistIntegrandSq<MathVector<TGridFunction::dim>, TGridFunction> integrand2(m_spVelocity, uFine, 0, uCoarse, 0);
		return IntegrateSubsets(integrand2, uFine, subsets, m_quadorder);
	}

	/// for weighted norms
	using weighted_obj_type::set_weight;
	using weighted_obj_type::get_weight;
	using weighted_obj_type::m_spWeight;

	void set_velocity(SmartPtr<velocity_type> spVelocity)
	{ m_spVelocity = spVelocity;}

protected:
	SmartPtr<velocity_type> m_spVelocity;

};
*/

//!
/** Evaluates distance between two grid functions in H1 semi-norm */
template <typename TGridFunction>
class H1EnergyComponentSpace
: public IComponentSpace<TGridFunction>,
  public IObjectWithWeights<typename H1EnergyDistIntegrand<TGridFunction>::weight_type >
{
public:
	typedef IComponentSpace<TGridFunction> base_type;
	typedef typename H1SemiDistIntegrand<TGridFunction>::weight_type weight_type;
	typedef IObjectWithWeights<weight_type> weighted_obj_type;
	typedef DarcyVelocityLinker<TGridFunction::dim> velocity_type;

	H1EnergyComponentSpace(const char *fctNames)
	: base_type(fctNames), weighted_obj_type(make_sp(new ConstUserMatrix<TGridFunction::dim>(1.0))), m_spVelocity(SPNULL) {};

	H1EnergyComponentSpace(const char *fctNames, int order)
	: base_type(fctNames, order), weighted_obj_type(make_sp(new ConstUserMatrix<TGridFunction::dim>(1.0))), m_spVelocity(SPNULL) {};

	H1EnergyComponentSpace(const char *fctNames, int order, number weight, const char* ssNames=0)
	: base_type(fctNames, ssNames, order), weighted_obj_type(make_sp(new ConstUserMatrix<TGridFunction::dim>(weight))), m_spVelocity(SPNULL)  {};

	H1EnergyComponentSpace(const char *fctNames, int order, ConstSmartPtr<weight_type> spWeight, const char* ssNames=0)
	: base_type(fctNames, ssNames, order), weighted_obj_type(spWeight), m_spVelocity(SPNULL) {};

	H1EnergyComponentSpace(const char *fctNames, int order, const char* ssNames, ConstSmartPtr<weight_type> spWeight)
	: base_type(fctNames, ssNames, order), weighted_obj_type(spWeight), m_spVelocity(SPNULL) {};

	/// DTOR
	~H1EnergyComponentSpace() {};


	/// \copydoc IComponentSpace<TGridFunction>::norm
	/// \copydoc IComponentSpace<TGridFunction>::distance
	using IComponentSpace<TGridFunction>::m_quadorder;
	using IComponentSpace<TGridFunction>::norm;
	using IComponentSpace<TGridFunction>::distance;

	/// \copydoc IComponentSpace<TGridFunction>::norm2
	double norm2(TGridFunction& uFine)
	{
		if (m_spVelocity.valid()) {
			const char* subsets = NULL;
			UserDataIntegrandSq<MathVector<TGridFunction::dim>, TGridFunction> integrand2(m_spVelocity, &uFine, 0.0);
			return IntegrateSubsets(integrand2, uFine, subsets, m_quadorder);
		} else {
			return H1EnergyNorm2<TGridFunction>(uFine, base_type::m_fctNames.c_str(), base_type::m_quadorder, NULL, weighted_obj_type::m_spWeight);
		}
	}

	/// \copydoc IComponentSpace<TGridFunction>::distance2
	double distance2(TGridFunction& uFine, TGridFunction& uCoarse)
	{ return H1EnergyDistance2<TGridFunction>(uFine, base_type::m_fctNames.c_str(), uCoarse, base_type::m_fctNames.c_str(), base_type::m_quadorder, m_spWeight); }

	/// for weighted norms
	using weighted_obj_type::set_weight;
	using weighted_obj_type::get_weight;
	using weighted_obj_type::m_spWeight;

	void set_velocity(SmartPtr<velocity_type> spVelocity)
	{ m_spVelocity = spVelocity;}

protected:
		SmartPtr<velocity_type> m_spVelocity;

};

/** Defines a H1 space for single component */
template <typename TGridFunction>
class H1ComponentSpace :
		public IComponentSpace<TGridFunction>
{
public:
	typedef IComponentSpace<TGridFunction> base_type;

	H1ComponentSpace(const char *fctNames) : base_type(fctNames) {};
	H1ComponentSpace(const char *fctNames, int order) : base_type(fctNames, order) {};
	H1ComponentSpace(const char *fctNames,  const char* ssNames, int order) : base_type(fctNames, ssNames, order) {};
	~H1ComponentSpace() {};

	using IComponentSpace<TGridFunction>::norm;
	using IComponentSpace<TGridFunction>::distance;

	/// \copydoc IComponentSpace<TGridFunction>::norm
	double norm2(TGridFunction& uFine)
	{ return H1Norm2<TGridFunction>(uFine, base_type::m_fctNames.c_str(), base_type::m_quadorder); }

	/// \copydoc IComponentSpace<TGridFunction>::norm
	double distance2(TGridFunction& uFine, TGridFunction& uCoarse)
	{ return H1Distance2<TGridFunction>(uFine, base_type::m_fctNames.c_str(), uCoarse, base_type::m_fctNames.c_str(), base_type::m_quadorder); }

};



//! Defines a composite space, (i.e., additive composition from other spaces)
/*! Employs a l2-type extension: \| u \|^2 := \sum_I \sigma_I \| u \|_I^2
 *  TODO: Merge with ApproximationSpace?
 * */
template <typename TGridFunction>
class CompositeSpace
		: public IGridFunctionSpace<TGridFunction>
{
public:
	typedef IGridFunctionSpace<TGridFunction> base_type;
	typedef IComponentSpace<TGridFunction> obj_type;
	typedef TimeDependentSpace<TGridFunction> time_dependent_obj_type;
	typedef std::pair<SmartPtr<obj_type>, number> weighted_obj_type;

	CompositeSpace() {};
	// virtual ~CompositeSpace() {};

	using base_type::norm;
	using base_type::distance;

	/// \copydoc IComponentSpace<TGridFunction>::norm
	double norm(TGridFunction& uFine)
	{ return(sqrt(norm2(uFine))); }

	/// \copydoc IComponentSpace<TGridFunction>::norm2
	double norm2(TGridFunction& uFine)
	{
		number unorm2 = 0.0;
		for (typename std::vector<weighted_obj_type>::iterator it = m_spWeightedSubspaces.begin();
				it!= m_spWeightedSubspaces.end(); ++it)
		{
			double snorm2 = it->first->norm2(uFine);
			unorm2 += it->second * snorm2;
			// std::cerr << "norm2:" << snorm2*snorm << std::endl;
		}
		return unorm2;
	}

	/// \copydoc IComponentSpace<TGridFunction>::distance2
	double distance2(TGridFunction& uFine, TGridFunction& uCoarse)
	{
		number unorm2 = 0.0;
		for (typename std::vector<weighted_obj_type>::iterator it = m_spWeightedSubspaces.begin();
					it!= m_spWeightedSubspaces.end(); ++it)
		{
			double sdist2 = it->first->distance2(uFine, uCoarse);
			unorm2 += it->second * sdist2;
			//	std::cerr << "distance:" << sdist2 << std::endl;
		}
		return unorm2;
	}

	/// \copydoc IComponentSpace<TGridFunction>::distance2
	double distance(TGridFunction& uFine, TGridFunction& uCoarse)
	{ return sqrt(distance2(uFine, uCoarse)); }

	/// add space to composite (with weight 1.0)
	void add(SmartPtr<obj_type> spSubSpace)
	{ m_spWeightedSubspaces.push_back(std::make_pair(spSubSpace, 1.0)); }

	/// add space to composite (with variable weight)
	void add(SmartPtr<obj_type> spSubSpace, number sigma)
	{ m_spWeightedSubspaces.push_back(std::make_pair(spSubSpace, sigma)); }

	/// print config string                                                                                                                                                                              
	std::string config_string() const
	{
	    std::stringstream ss;
	    ss << "CompositeSpace:" << std::endl;

	    for (typename std::vector<weighted_obj_type>::const_iterator it = m_spWeightedSubspaces.begin();
	    		it!= m_spWeightedSubspaces.end(); ++it)
	    { ss << it->first->config_string(); }

	    return ss.str();
	}


	//! Forward update to all members
	void update_time_data(number t)
	{
		for (typename std::vector<weighted_obj_type>::iterator it = m_spWeightedSubspaces.begin();
			it!= m_spWeightedSubspaces.end(); ++it)
		{
			SmartPtr<time_dependent_obj_type> spSpaceT = it->first.template cast_dynamic<time_dependent_obj_type> ();
			if (spSpaceT.valid()) spSpaceT->update_time_data(t);
		}
	}

	//! Check, if any object is time-dependent.
	bool is_time_dependent() const
	{
		for (typename std::vector<weighted_obj_type>::const_iterator it = m_spWeightedSubspaces.begin();
				it!= m_spWeightedSubspaces.end(); ++it)
		{
			SmartPtr<time_dependent_obj_type> spSpaceT = it->first.template cast_dynamic<time_dependent_obj_type>();
			if (spSpaceT.valid()) return true;
		}
		return false;
	}

	const std::vector<weighted_obj_type> &get_subspaces() const
	{ return m_spWeightedSubspaces; }

protected:
	std::vector<weighted_obj_type> m_spWeightedSubspaces;

};

/*
template <typename TGridFunction>
class TimeDependentCompositeSpace
	: public CompositeSpace<TGridFunction>, public ITimeDependentSpace<TGridFunction>
{
protected:
	using CompositeSpace<TGridFunction>::m_spSubspaces;

public:
	typedef typename IComponentSpace<TGridFunction> obj_type;


	/// add spaces to composite
	void add(SmartPtr<obj_type> spSubSpace)
	{ m_spSubspaces.push_back(make_sp(new ITimeDependentSpace(spSubSpace))); }

	//! Forward update to all members
	void update_time_data(number t)
	{
		for (typename std::vector<SmartPtr<obj_type> >::iterator it = m_spSubspaces.begin();
			it!= m_spSubspaces.end(); ++it)
		{ if ((*it)->is_time_dependent()) (*it)->update_time_data(t); }
	}

	//! Check, if all objects are time-dependent.
	bool is_time_dependent() const
	{
		for (typename std::vector<SmartPtr<obj_type> >::const_iterator it = m_spSubspaces.begin();
				it!= m_spSubspaces.end(); ++it)
		{ if ((*it)->is_time_dependent()) return true; }
		return false;
	}
};

 */

} // namespace ug
#endif
