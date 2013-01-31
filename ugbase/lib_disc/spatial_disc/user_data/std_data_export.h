/*
 * std_data_export.h
 *
 *  Created on: 22.05.2012
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__STD_DATA_IMPORT__
#define __H__UG__LIB_DISC__SPATIAL_DISC__STD_DATA_IMPORT__

#include "lib_disc/reference_element/reference_mapping_provider.h"
#include "lib_disc/spatial_disc/user_data/data_export.h"

namespace ug{

////////////////////////////////////////////////////////////////////////////////
// ValueDataExport
////////////////////////////////////////////////////////////////////////////////

template <int dim>
class ValueDataExport
	: public StdDataExport<ValueDataExport<dim>,number,dim>
{
	public:
		ValueDataExport(const char* functions){this->set_symb_fct(functions);}

		template <int refDim>
		inline void evaluate(number vValue[],
		                     const MathVector<dim> vGlobIP[],
		                     number time, int si,
		                     LocalVector& u,
		                     GeometricObject* elem,
		                     const MathVector<dim> vCornerCoords[],
		                     const MathVector<refDim> vLocIP[],
		                     const size_t nip,
		                     const MathMatrix<refDim, dim>* vJT = NULL) const
		{
		//	reference object id
			const ReferenceObjectID roid = elem->reference_object_id();

		//	local finite element id
			const LFEID& lfeID = this->fct_grp().local_finite_element_id(_C_);

		//	access local vector by map
			u.access_by_map(this->fct_index_map());

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
					vValue[ip] += u(_C_, sh) * vShape[sh];
			}

			}
			UG_CATCH_THROW("ValueDataExport: Trial space missing, Reference Object: "
			               <<roid<<", Trial Space: "<<lfeID<<", refDim="<<refDim);
		}

	///	returns if provided data is continuous over geometric object boundaries
		virtual bool continuous() const
		{
			const LFEID& lfeID = this->fct_grp().local_finite_element_id(_C_);

			if(lfeID.type() == LFEID::LAGRANGE) return true;
			else return false;
		}

	protected:
	//	abbreviation for component
		static const int _C_ = 0;
};

////////////////////////////////////////////////////////////////////////////////
// GradientDataExport
////////////////////////////////////////////////////////////////////////////////

template <int dim>
class GradientDataExport
	: public StdDataExport<GradientDataExport<dim>, MathVector<dim>,dim>
{
	public:
		GradientDataExport(const char* functions){this->set_symb_fct(functions);}

		template <int refDim>
		inline void evaluate(MathVector<dim> vValue[],
		                       const MathVector<dim> vGlobIP[],
		                       number time, int si,
		                       LocalVector& u,
		                       GeometricObject* elem,
		                       const MathVector<dim> vCornerCoords[],
		                       const MathVector<refDim> vLocIP[],
		                       const size_t nip,
		                       const MathMatrix<refDim, dim>* vJT = NULL) const
		{
		//	reference object id
			const ReferenceObjectID roid = elem->reference_object_id();

		//	local finite element id
			const LFEID& lfeID = this->fct_grp().local_finite_element_id(_C_);

		//	access local vector by map
			u.access_by_map(this->fct_index_map());

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
					VecScaleAppend(locGrad, u(_C_, sh), vLocGrad[sh]);

				Inverse(JTInv, vJT[ip]);
				MatVecMult(vValue[ip], JTInv, locGrad);
			}

			}
			UG_CATCH_THROW("GradientDataExport: Trial space missing, Reference Object: "
						 <<roid<<", Trial Space: "<<lfeID<<", refDim="<<refDim);
		}

	///	returns if provided data is continuous over geometric object boundaries
		virtual bool continuous() const
		{
			return false;
		}

	protected:
	//	abbreviation for component
		static const int _C_ = 0;
};

} // end namespace ug

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__STD_DATA_IMPORT__ */
