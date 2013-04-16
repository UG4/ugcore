/*
 * data_export_impl.h
 *
 *  Created on: 04.07.2012
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__DATA_EXPORT_IMPL__
#define __H__UG__LIB_DISC__SPATIAL_DISC__DATA_EXPORT_IMPL__

#include "data_export.h"
#include "lib_disc/reference_element/reference_element_util.h"

namespace ug{

////////////////////////////////////////////////////////////////////////////////
// ValueDataExport
////////////////////////////////////////////////////////////////////////////////

template <int dim>
template <int refDim>
void ValueDataExport<dim>::eval_and_deriv(number vValue[],
                    const MathVector<dim> vGlobIP[],
                    number time, int si,
                    GeometricObject* elem,
                    const MathVector<dim> vCornerCoords[],
                    const MathVector<refDim> vLocIP[],
                    const size_t nip,
                    LocalVector* u,
                    bool bDeriv,
                    int s,
                    std::vector<std::vector<number> > vvvDeriv[],
                    const MathMatrix<refDim, dim>* vJT) const
{
//	reference object id
	const ReferenceObjectID roid = elem->reference_object_id();

//	local finite element id
	const LFEID& lfeID = this->function_group().local_finite_element_id(_C_);

//	access local vector by map
	u->access_by_map(this->map());

//	request for trial space
	try{
	const LocalShapeFunctionSet<refDim>& rTrialSpace
		 = LocalShapeFunctionSetProvider::get<refDim>(roid, lfeID);

//	memory for shapes
	std::vector<number> vShape;

//	loop ips
	for(size_t ip = 0; ip < nip; ++ip)
	{
	//	evaluate at shapes at ip
		rTrialSpace.shapes(vShape, vLocIP[ip]);

	//	compute concentration at ip
		vValue[ip] = 0.0;
		for(size_t sh = 0; sh < vShape.size(); ++sh)
			vValue[ip] += (*u)(_C_, sh) * vShape[sh];
	}

	}
	UG_CATCH_THROW("ValueDataExport: Trial space missing, Reference Object: "
	               <<roid<<", Trial Space: "<<lfeID<<", refDim="<<refDim);

	if(bDeriv)
		UG_THROW("Not implemented.");
}

template <int dim>
bool ValueDataExport<dim>::continuous() const
{
	const LFEID& lfeID = this->function_group().local_finite_element_id(_C_);
	return LocalShapeFunctionSetProvider::continuous(lfeID);
}

////////////////////////////////////////////////////////////////////////////////
// GradientDataExport
////////////////////////////////////////////////////////////////////////////////

template <int dim>
template <int refDim>
void GradientDataExport<dim>::eval_and_deriv(MathVector<dim> vValue[],
                    const MathVector<dim> vGlobIP[],
                    number time, int si,
                    GeometricObject* elem,
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

//	local finite element id
	const LFEID& lfeID = this->function_group().local_finite_element_id(_C_);

//	access local vector by map
	u->access_by_map(this->map());

//	request for trial space
	try{
	const LocalShapeFunctionSet<refDim>& rTrialSpace
		 = LocalShapeFunctionSetProvider::get<refDim>(roid, lfeID);

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
	}

	}
	UG_CATCH_THROW("GradientDataExport: Trial space missing, Reference Object: "
				 <<roid<<", Trial Space: "<<lfeID<<", refDim="<<refDim);

	if(bDeriv)
		UG_THROW("Not implemented.");
}


////////////////////////////////////////////////////////////////////////////////
// DataExport
////////////////////////////////////////////////////////////////////////////////

template <typename TData, int dim>
DataExport<TData, dim>::DataExport(const char* functions)
: StdDependentUserData<DataExport<TData,dim>, TData, dim>(functions),
  m_id(ROID_UNKNOWN)
{
//	reset all evaluation functions
	clear_fct();
}

template <typename TData, int dim>
void DataExport<TData, dim>::clear_fct()
{
	for(size_t roid = 0; roid < NUM_REFERENCE_OBJECTS; ++roid){
		eval_fct<1>((ReferenceObjectID)roid).invalidate();
		eval_fct<2>((ReferenceObjectID)roid).invalidate();
		eval_fct<3>((ReferenceObjectID)roid).invalidate();
	}
}

template <typename TData, int dim>
template <typename TClass, int refDim>
void DataExport<TData, dim>::
set_fct(ReferenceObjectID id, TClass* obj,
        void (TClass::*func)(	TData vValue[],
								const MathVector<dim> vGlobIP[],
								number time, int si,
								const LocalVector& u,
								GeometricObject* elem,
								const MathVector<dim> vCornerCoords[],
								const MathVector<refDim> vLocIP[],
								const size_t nip,
								bool bDeriv,
								std::vector<std::vector<TData> > vvvDeriv[]))
{
	if(id >= NUM_REFERENCE_OBJECTS)
		UG_THROW("Reference Object id invalid: "<<id);

	eval_fct<refDim>(id) = Functor<refDim>(obj, func);
}

template <typename TData, int dim>
template <int refDim>
void DataExport<TData, dim>::
set_fct(ReferenceObjectID id,
        void (*func)(	TData vValue[],
        				const MathVector<dim> vGlobIP[],
						number time, int si,
						const LocalVector& u,
						GeometricObject* elem,
						const MathVector<dim> vCornerCoords[],
						const MathVector<refDim> vLocIP[],
						const size_t nip,
        				bool bDeriv,
        				std::vector<std::vector<TData> > vvvDeriv[]))
{
	if(id >= NUM_REFERENCE_OBJECTS)
		UG_THROW("Reference Object id invalid: "<<id);

	eval_fct<refDim>(id) = Functor<refDim>(func);
}

template <typename TData, int dim>
void DataExport<TData, dim>::set_roid(ReferenceObjectID id)
{
	if(id == ROID_UNKNOWN)
		UG_THROW("DataExport::set_roid: Setting unknown ReferenceObjectId.");

	m_id = id;
}

template <typename TData, int dim>
void DataExport<TData, dim>::check_setup() const
{
	if(m_id == ROID_UNKNOWN)
		UG_THROW("DataExport::check_setup: The reference element "
				"type has not been set for evaluation.");

	if(!eval_fct_set(m_id))
		UG_THROW("DataExport::check_setup: There is no evaluation "
				"function registered for data export and element type "<<m_id<<
				", but required. (world dim: "<<dim<<", ref dim: "<<
				ReferenceElementDimension(m_id)<<")");
}

template <typename TData, int dim>
bool DataExport<TData, dim>::eval_fct_set(ReferenceObjectID id) const
{
	const int d = ReferenceElementDimension(id);
	bool bRes = false;
	switch(d){
		case 1: if(eval_fct<1>(id).valid()) bRes = true; break;
		case 2: if(eval_fct<2>(id).valid()) bRes = true; break;
		case 3: if(eval_fct<3>(id).valid()) bRes = true; break;
		default: UG_THROW("DataExport: Dimension "<<d<<" not supported.");
	}
	return bRes;
}


template <typename TData, int dim>
void DataExport<TData, dim>::
add_needed_data(SmartPtr<ICplUserData<dim> > data)
{
	m_vDependData.push_back(data);
}

template <typename TData, int dim>
void DataExport<TData, dim>::
remove_needed_data(SmartPtr<ICplUserData<dim> > data)
{
	m_vDependData.erase(remove(m_vDependData.begin(),
	                           m_vDependData.end(),
	                           data),
	                           m_vDependData.end());
}

} // end namespace ug

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__DATA_EXPORT_IMPL__ */
