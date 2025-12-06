/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
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

#include "data_export.h"
//ø #include "lib_disc/reference_element/reference_element_util.h"

namespace ug{

////////////////////////////////////////////////////////////////////////////////
// ValueDataExport
////////////////////////////////////////////////////////////////////////////////

template <int dim>
template <int refDim>
void ValueDataExport<dim>::eval_and_deriv(number vValue[],
                    const MathVector<dim> vGlobIP[],
                    number time, int si,
                    GridObject* elem,
                    const MathVector<dim> vCornerCoords[],
                    const MathVector<refDim> vLocIP[],
                    const size_t nip,
                    LocalVector* u,
                    bool bDeriv,
                    int s,
                    std::vector<std::vector<number> > vvvDeriv[],
                    const MathMatrix<refDim, dim>* vJT) const
{
//	abbreviation for component
	static constexpr int _C_ = 0;

//	reference object id
	const ReferenceObjectID roid = elem->reference_object_id();

//	local finite element id
	const LFEID& lfeID = this->function_group().local_finite_element_id(_C_);

//	access local vector by map
	u->access_by_map(this->map());

//	request for trial space
	try{
	const LocalShapeFunctionSet<refDim>& rTrialSpace
		 = LocalFiniteElementProvider::get<refDim>(roid, lfeID);

//	memory for shapes
	std::vector<number> vShape;

//	loop ips
	for(size_t ip = 0; ip < nip; ++ip)
	{
	//	evaluate at shapes at ip
		rTrialSpace.shapes(vShape, vLocIP[ip]);

	//	compute value at ip
		vValue[ip] = 0.0;
		for(size_t sh = 0; sh < vShape.size(); ++sh)
			vValue[ip] += (*u)(_C_, sh) * vShape[sh];

	//	store derivative
		if(bDeriv)
			for(size_t sh = 0; sh < vShape.size(); ++sh)
				vvvDeriv[ip][_C_][sh] = vShape[sh];
	}

	}
	UG_CATCH_THROW("ValueDataExport: Trial space missing, Reference Object: "
	               <<roid<<", Trial Space: "<<lfeID<<", refDim="<<refDim);
}

template <int dim>
void ValueDataExport<dim>::check_setup() const
{
	if((int)this->function_group().size() != 1)
		UG_THROW("ValueDataExport: Expected to work on "<<1<<" function, "
		         "but only set "<<this->function_group().size()<<" ("<<
		         this->function_group().names()<<").");
}

template <int dim>
bool ValueDataExport<dim>::continuous() const
{
	static constexpr int _C_ = 0;
	const LFEID& lfeID = this->function_group().local_finite_element_id(_C_);
	return LocalFiniteElementProvider::continuous(lfeID);
}

////////////////////////////////////////////////////////////////////////////////
// GradientDataExport
////////////////////////////////////////////////////////////////////////////////

template <int dim>
template <int refDim>
void GradientDataExport<dim>::eval_and_deriv(MathVector<dim> vValue[],
                    const MathVector<dim> vGlobIP[],
                    number time, int si,
                    GridObject* elem,
                    const MathVector<dim> vCornerCoords[],
                    const MathVector<refDim> vLocIP[],
                    const size_t nip,
                    LocalVector* u,
                    bool bDeriv,
                    int s,
                    std::vector<std::vector<MathVector<dim> > > vvvDeriv[],
                    const MathMatrix<refDim, dim>* vJT) const
{
//	abbreviation for component
	static constexpr int _C_ = 0;

//	reference object id
	const ReferenceObjectID roid = elem->reference_object_id();

//	local finite element id
	const LFEID& lfeID = this->function_group().local_finite_element_id(_C_);

//	access local vector by map
	u->access_by_map(this->map());

//	request for trial space
	try{
	const LocalShapeFunctionSet<refDim>& rTrialSpace
		 = LocalFiniteElementProvider::get<refDim>(roid, lfeID);

//	Reference Mapping
	MathMatrix<dim, refDim> JTInv;
	std::vector<MathMatrix<refDim, dim> > vJTtmp;
	if(!vJT){
		DimReferenceMapping<refDim, dim>& map
			= ReferenceMappingProvider::get<refDim, dim>(roid, vCornerCoords);

		vJTtmp.resize(nip);
		map.jacobian_transposed(&vJTtmp[0], vLocIP, nip);
		vJT = &vJTtmp[0];
	}

//	storage for shape function at ip
	std::vector<MathVector<refDim> > vLocGrad;
	MathVector<refDim> locGrad;

//	loop ips
	for(size_t ip = 0; ip < nip; ++ip)
	{
	//	evaluate at shapes at ip
		rTrialSpace.grads(vLocGrad, vLocIP[ip]);

	//	compute grad at ip
		VecSet(locGrad, 0.0);
		for(size_t sh = 0; sh < vLocGrad.size(); ++sh)
			VecScaleAppend(locGrad, (*u)(_C_, sh), vLocGrad[sh]);

		Inverse(JTInv, vJT[ip]);
		MatVecMult(vValue[ip], JTInv, locGrad);

	//	store derivative
		if(bDeriv)
			for(size_t sh = 0; sh < vLocGrad.size(); ++sh)
				MatVecMult(vvvDeriv[ip][_C_][sh], JTInv, vLocGrad[sh]);
	}

	}
	UG_CATCH_THROW("GradientDataExport: Trial space missing, Reference Object: "
				 <<roid<<", Trial Space: "<<lfeID<<", refDim="<<refDim);
}

template <int dim>
void GradientDataExport<dim>::check_setup() const
{
	if((int)this->function_group().size() != 1)
		UG_THROW("GradientDataExport: Expected to work on "<<1<<" function, "
		         "but only set "<<this->function_group().size()<<" ("<<
		         this->function_group().names()<<").");
}

////////////////////////////////////////////////////////////////////////////////
// VectorDataExport
////////////////////////////////////////////////////////////////////////////////

template <int dim>
template <int refDim>
void VectorDataExport<dim>::eval_and_deriv(MathVector<dim> vValue[],
                    const MathVector<dim> vGlobIP[],
                    number time, int si,
                    GridObject* elem,
                    const MathVector<dim> vCornerCoords[],
                    const MathVector<refDim> vLocIP[],
                    const size_t nip,
                    LocalVector* u,
                    bool bDeriv,
                    int s,
                    std::vector<std::vector<MathVector<dim> > > vvvDeriv[],
                    const MathMatrix<refDim, dim>* vJT) const
{
//	reference object id
	const ReferenceObjectID roid = elem->reference_object_id();

//	access local vector by map
	u->access_by_map(this->map());

	if(bDeriv)
		this->set_zero(vvvDeriv, nip);

	for(size_t d = 0; d < dim; ++d)
	{
	//	local finite element id
		const LFEID& lfeID = this->function_group().local_finite_element_id(d);

	//	request for trial space
		try{
		const LocalShapeFunctionSet<refDim>& rTrialSpace
			 = LocalFiniteElementProvider::get<refDim>(roid, lfeID);

	//	memory for shapes
		std::vector<number> vShape;

	//	loop ips
		for(size_t ip = 0; ip < nip; ++ip)
		{
		//	evaluate at shapes at ip
			rTrialSpace.shapes(vShape, vLocIP[ip]);

		//	compute value at ip
			vValue[ip] = 0.0;
			for(size_t sh = 0; sh < vShape.size(); ++sh)
				vValue[ip][d] += (*u)(d, sh) * vShape[sh];

		//	store derivative
			if(bDeriv)
				for(size_t sh = 0; sh < vShape.size(); ++sh)
					vvvDeriv[ip][d][sh][d] = vShape[sh];
		}

		}
		UG_CATCH_THROW("VectorDataExport: Trial space missing, Reference Object: "
					   <<roid<<", Trial Space: "<<lfeID<<", refDim="<<refDim);
	}
}

template <int dim>
void VectorDataExport<dim>::check_setup() const
{
	if((int)this->function_group().size() != dim)
		UG_THROW("VectorDataExport: Expected to work on "<<dim<<" functions, "
		         "but only set "<<this->function_group().size()<<" ("<<
		         this->function_group().names()<<").");
}

template <int dim>
bool VectorDataExport<dim>::continuous() const
{
	for(size_t d = 0; d < dim; ++d){
		const LFEID& lfeID = this->function_group().local_finite_element_id(d);
		if(!LocalFiniteElementProvider::continuous(lfeID)) return false;
	}
	return true;
}

////////////////////////////////////////////////////////////////////////////////
// explicit instantiations
////////////////////////////////////////////////////////////////////////////////

#ifdef UG_DIM_1
template class ValueDataExport<1>;
template class GradientDataExport<1>;
template class VectorDataExport<1>;
#endif
#ifdef UG_DIM_2
template class ValueDataExport<2>;
template class GradientDataExport<2>;
template class VectorDataExport<2>;
#endif
#ifdef UG_DIM_3
template class ValueDataExport<3>;
template class GradientDataExport<3>;
template class VectorDataExport<3>;
#endif

} // end namespace ug
