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

namespace ug {
/// Abstract base class for (algebraic) vectors
template<typename TVector>
class IBanachSpace
{

public:
	virtual ~IBanachSpace() {}

	/// euclidean norm (default)
	virtual double norm(TVector &x) const
	{ return x.norm(); }

	virtual double distance(TVector& x, TVector& y) const
	{
		SmartPtr<TVector> delta = x.clone();
		*delta -= y;
		return norm(*delta);
	}
};



/// Abstract base class for grid functions
template<typename TGridFunction>
class IGridFunctionSpace  : public IBanachSpace<typename TGridFunction::vector_type>
{

public:
	typedef typename TGridFunction::vector_type vector_type;
	typedef TGridFunction grid_function_type;

	/// DTOR
	virtual ~IGridFunctionSpace() {}

	/// norm (for grid functions)
	virtual double norm(TGridFunction& x) const
	{ return x.norm(); }

	/// distance (for grid functions)
	virtual double distance(TGridFunction& x, TGridFunction& y) const
	{
		SmartPtr<TGridFunction> delta = x.clone();
		*delta -= y;
		return norm(*delta);
	}

	/// norm (for vectors)
	virtual double norm(vector_type &x) const
	{
		TGridFunction* gfX=dynamic_cast< TGridFunction*>(&x);
		UG_ASSERT(gfX!=NULL, "Huhh: GridFunction required!");
		return norm(*gfX);
	}

	/// distance (for vectors)
	virtual double distance(vector_type &x, vector_type &y) const
	{ return distance(static_cast<TGridFunction &>(x), static_cast<TGridFunction &>(y)); }




};



//! Estimate the error (based on the difference between two grid functions)
template <typename TGridFunction>
class IComponentSpace : public IGridFunctionSpace<TGridFunction>
{
protected:
	std::string m_fctNames;
	const char* m_ssNames;
	int m_quadorder;
	number m_scale;

public:
	typedef IGridFunctionSpace<TGridFunction> base_type;

	IComponentSpace(const IComponentSpace<TGridFunction> &val)
	{
		this->m_fctNames = val.m_fctNames;
		this->m_ssNames = val.m_ssNames;
		this->m_quadorder = val.m_quadorder;
		this->m_scale = val.m_scale;
	}

	IComponentSpace(const char *fctNames) :
		m_fctNames(fctNames), m_ssNames(NULL), m_quadorder(3), m_scale(1.0) {}

	IComponentSpace(const char *fctNames, int order) :
		m_fctNames(fctNames), m_ssNames(NULL), m_quadorder(order), m_scale(1.0) {}

	IComponentSpace(const char *fctNames, int order, number scale) :
		m_fctNames(fctNames), m_ssNames(NULL), m_quadorder(order), m_scale(scale) {}

	IComponentSpace(const char *fctNames, const char* ssNames, int order, number scale) :
		m_fctNames(fctNames), m_ssNames(ssNames), m_quadorder(order), m_scale(scale) {}

	virtual ~IComponentSpace() {};


	// suppress warnings
	using IGridFunctionSpace<TGridFunction>::norm;
	using IGridFunctionSpace<TGridFunction>::distance;

	/// norm (for SmartPtr)
	virtual double norm(SmartPtr<TGridFunction> u) const
	{ return norm(*u); }

	/// distance (for SmartPtr)
	double distance(SmartPtr<TGridFunction> uFine, SmartPtr<TGridFunction> uCoarse) const
	{ return distance(*uFine, *uCoarse); }

	/// print config string
	std::string config_string() const
	{
		std::stringstream ss;

		if (this->m_ssNames)
		ss << this->m_fctNames << ", " << this->m_ssNames << ", " << this->m_quadorder
			<< ", type=" <<", scale=" << this->m_scale << std::endl;

		else
			ss << this->m_fctNames << ",  (no name), " << this->m_quadorder
						<< ", type=" <<", scale=" << this->m_scale << std::endl;
		return ss.str();
	}

};



/** Evaluates difference between two grid functions in L2 norm */
template <typename TGridFunction>
class L2ComponentSpace : public IComponentSpace<TGridFunction>
{
public:
	typedef IComponentSpace<TGridFunction> base_type;

	L2ComponentSpace(const char *fctNames) : base_type(fctNames) {};
	L2ComponentSpace(const char *fctNames, int order) : base_type(fctNames, order) {};
	L2ComponentSpace(const char *fctNames, int order, number scale) : base_type(fctNames, order, scale) {};
	L2ComponentSpace(const char *fctNames, const char* ssNames, int order, number scale)
		: base_type(fctNames, ssNames, order, scale) {};
	~L2ComponentSpace() {};

	using IComponentSpace<TGridFunction>::norm;
	using IComponentSpace<TGridFunction>::distance;

	/// \copydoc IComponentSpace<TGridFunction>::norm
	double norm(TGridFunction& uFine) const
	{ return L2Norm(uFine, base_type::m_fctNames.c_str(), base_type::m_quadorder, base_type::m_ssNames); }

	/// \copydoc IComponentSpace<TGridFunction>::distance
	double distance(TGridFunction& uFine, TGridFunction& uCoarse) const
	{ return L2Error(uFine, base_type::m_fctNames.c_str(), uCoarse, base_type::m_fctNames.c_str(),
		base_type::m_quadorder, base_type::m_ssNames);}
};


/** Evaluates distance between two grid functions in H1 semi-norm */
template <typename TGridFunction>
class H1SemiComponentSpace : public IComponentSpace<TGridFunction>
{
public:
	typedef IComponentSpace<TGridFunction> base_type;
	typedef typename H1SemiDistIntegrand<TGridFunction>::weight_type weight_type;

	H1SemiComponentSpace(const char *fctNames) : base_type(fctNames) {};
	H1SemiComponentSpace(const char *fctNames, int order) : base_type(fctNames, order) {};
	H1SemiComponentSpace(const char *fctNames, int order, number scale) : base_type(fctNames, order, scale) {};
	H1SemiComponentSpace(const char *fctNames, int order, number scale, SmartPtr<weight_type> spWeight)
	: base_type(fctNames, order, scale),  m_spWeight(spWeight) {};

	~H1SemiComponentSpace() {};

	using IComponentSpace<TGridFunction>::norm;
	using IComponentSpace<TGridFunction>::distance;

	/// \copydoc IComponentSpace<TGridFunction>::norm
	double norm(TGridFunction& uFine) const
	{ return H1SemiNorm<TGridFunction>(uFine, base_type::m_fctNames.c_str(), base_type::m_quadorder, NULL, m_spWeight); }

	/// \copydoc IComponentSpace<TGridFunction>::norm
	double distance(TGridFunction& uFine, TGridFunction& uCoarse) const
	{ return H1SemiDistance<TGridFunction>(uFine, base_type::m_fctNames.c_str(), uCoarse, base_type::m_fctNames.c_str(), base_type::m_quadorder, m_spWeight); }


	void set_weight(SmartPtr<weight_type> spWeight)
	{ m_spWeight = spWeight; }

	SmartPtr<weight_type> get_weight()
	{ return m_spWeight; }

protected:
	SmartPtr<weight_type> m_spWeight;
};

/** Evaluates difference between two grid functions in H1 semi-norm */
template <typename TGridFunction>
class H1ComponentSpace :
		public IComponentSpace<TGridFunction>
{
public:
	typedef IComponentSpace<TGridFunction> base_type;

	H1ComponentSpace(const char *fctNames) : base_type(fctNames) {};
	H1ComponentSpace(const char *fctNames, int order) : base_type(fctNames, order) {};
	H1ComponentSpace(const char *fctNames, int order, number scale) : base_type(fctNames, order, scale) {};
	~H1ComponentSpace() {};

	using IComponentSpace<TGridFunction>::norm;
	using IComponentSpace<TGridFunction>::distance;

	/// \copydoc IComponentSpace<TGridFunction>::norm
	double norm(SmartPtr<TGridFunction> uFine) const
	{ return H1Norm<TGridFunction>(uFine, base_type::m_fctNames.c_str(), base_type::m_quadorder); }

	/// \copydoc IComponentSpace<TGridFunction>::norm
	double distance(SmartPtr<TGridFunction> uFine, SmartPtr<TGridFunction> uCoarse) const
	{ return H1Error<TGridFunction>(uFine, base_type::m_fctNames.c_str(), uCoarse, base_type::m_fctNames.c_str(), base_type::m_quadorder); }

};

} // namespace ug
#endif
