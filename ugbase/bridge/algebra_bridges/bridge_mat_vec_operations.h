/*
 * Copyright (c) 2014-2015:  G-CSC, Goethe University Frankfurt
 * Author: Martin Rupp
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


// include bridge
#include "bridge/bridge.h"
#include "bridge/util.h"
#include "bridge/util_algebra_dependent.h"

#ifndef BRIDGE_MAT_VEC_OPERATIONS_H_
#define BRIDGE_MAT_VEC_OPERATIONS_H_

namespace ug{

namespace bridge{
namespace AlgebraCommon{

template<typename TAlgebra>
class VecScaleAddClass
{
	typedef typename TAlgebra::vector_type vector_type;
public:
	VecScaleAddClass(double scale1, SmartPtr<vector_type> v1)
	{
		scaling.push_back(scale1);
		vecs.push_back(v1);
	}
	VecScaleAddClass(double scale1, SmartPtr<vector_type> v1, double scale2, SmartPtr<vector_type> v2)
	{
		scaling.push_back(scale1);
		vecs.push_back(v1);
		scaling.push_back(scale2);
		vecs.push_back(v2);
	}

	VecScaleAddClass(double scale, SmartPtr<VecScaleAddClass<TAlgebra > > vsac, double scale1, SmartPtr<vector_type> v1)
	{
		scaling.push_back(scale1);
		vecs.push_back(v1);
		for(size_t i=0; i<vsac->size(); i++)
		{
			scaling.push_back(scale*vsac->scaling[i]);
			vecs.push_back(vsac->vecs[i]);
		}
	}

	VecScaleAddClass(double scale1, SmartPtr<vector_type> v1, double scale, SmartPtr<VecScaleAddClass<TAlgebra > > vsac)
	{
		scaling.push_back(scale1);
		vecs.push_back(v1);
		for(size_t i=0; i<vsac->size(); i++)
		{
			scaling.push_back(scale*vsac->scaling[i]);
			vecs.push_back(vsac->vecs[i]);
		}
	}

	VecScaleAddClass(double scale, SmartPtr<VecScaleAddClass<TAlgebra > > vsac)
	{
		for(size_t i=0; i<vsac->size(); i++)
		{
			scaling.push_back(scale*vsac->scaling[i]);
			vecs.push_back(vsac->vecs[i]);
		}
	}

	VecScaleAddClass() {}

	VecScaleAddClass(const VecScaleAddClass<TAlgebra> &parent) :
		scaling(parent.scaling), vecs(parent.vecs)
	{	}


	size_t size() const { return scaling.size(); }

	SmartPtr<vector_type> eval()
	{
		SmartPtr<vector_type> p = vecs[0]->clone();
		assign(p);
		return p;
	}

	void assign(SmartPtr<vector_type> p)
	{
#ifdef UG_PARALLEL
		for(size_t i=0; i<size(); i++)
			vecs[i]->change_storage_type(PST_CONSISTENT);
#endif
		if(size() == 1)
			VecScaleAssign(*p, scaling[0], *vecs[0]);
		else if(size() == 2)
			VecScaleAdd(*p, scaling[0], *vecs[0], scaling[1], *vecs[1]);
		else if(size() == 3)
			VecScaleAdd(*p, scaling[0], *vecs[0], scaling[1], *vecs[1], scaling[2], *vecs[2]);
		else
		{
			UG_THROW("not supported ATM.");
		}
	}
private:
	std::vector<double> scaling;
	std::vector<SmartPtr<vector_type> > vecs;
};

template<typename TAlgebra>
void Assign(SmartPtr<typename TAlgebra::vector_type> p, SmartPtr<VecScaleAddClass<TAlgebra> > vsac)
{
	vsac->assign(p);
}

template<typename TAlgebra>
SmartPtr<typename TAlgebra::vector_type> Eval(SmartPtr<VecScaleAddClass<TAlgebra> > vsac)
{
	return vsac->eval();
}


template<typename TAlgebra>
double VecProd2(SmartPtr<typename TAlgebra::vector_type> v1, SmartPtr<typename TAlgebra::vector_type> v2)
{
	return v1->dotprod(*v2);
}

template<typename T>
SmartPtr<T> VecProdOp(SmartPtr< ILinearOperator<T, T> > op, SmartPtr<T> v)
{
	SmartPtr<T> v2 = v->clone();
#ifdef UG_PARALLEL
	v->change_storage_type(PST_CONSISTENT);
#endif
	op->apply(*v2, *v);
	return v2;
}


template<typename TAlgebra>
double VecScaleAddProd1(SmartPtr<typename TAlgebra::vector_type> v1, SmartPtr<VecScaleAddClass<TAlgebra> > vsac)
{
	SmartPtr<typename TAlgebra::vector_type> v2 = vsac->eval();
	return VecProd2<TAlgebra>(v1, v2);
}

template<typename TAlgebra>
double VecScaleAddProd2(SmartPtr<VecScaleAddClass<TAlgebra> > vsac, SmartPtr<typename TAlgebra::vector_type> v1)
{
	SmartPtr<typename TAlgebra::vector_type> v2 = vsac->eval();
	return VecProd2<TAlgebra>(v1, v2);
}

template<typename TAlgebra>
double VecNorm(SmartPtr<typename TAlgebra::vector_type> v)
{
	return v->norm();
}


template<typename TAlgebra>
double VecScaleAddNorm(SmartPtr<VecScaleAddClass<TAlgebra> > vsac)
{
	SmartPtr<typename TAlgebra::vector_type> v = vsac->eval();
	return v->norm();
}



//! calculates dest = alpha1*v1 + alpha2*v2
template<typename vector_t>
inline void VecScaleAdd2(vector_t &dest, double alpha1, const vector_t &v1, double alpha2, const vector_t &v2)
{
	VecScaleAdd(dest, alpha1, v1, alpha2, v2);
}

//! calculates dest = alpha1*v1 + alpha2*v2 + alpha3*v3
template<typename vector_t>
inline void VecScaleAdd3(vector_t &dest, double alpha1, const vector_t &v1, double alpha2, const vector_t &v2, double alpha3, const vector_t &v3)
{
	VecScaleAdd(dest, alpha1, v1, alpha2, v2, alpha3, v3);
}

}
}
}
#endif /* BRIDGE_MAT_VEC_OPERATIONS_H_ */
