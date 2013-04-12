/*
 * std_data_export.h
 *
 *  Created on: 22.05.2012
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__STD_DATA_IMPORT__
#define __H__UG__LIB_DISC__SPATIAL_DISC__STD_DATA_IMPORT__

#include "std/std_user_data.h"
#include "data_export.h"

#include "lib_disc/common/function_group.h"
#include "lib_disc/common/local_algebra.h"
#include "lib_disc/common/groups_util.h"
#include "common/util/string_util.h"

#include "lib_disc/reference_element/reference_mapping_provider.h"
#include "lib_disc/spatial_disc/user_data/data_export.h"
#include "lib_disc/local_finite_element/local_shape_function_set.h"

namespace ug{

///////////////////////////////////////////////////////////////////////////////
// Base class for Exports
///////////////////////////////////////////////////////////////////////////////

template <typename TImpl, typename TData, int dim>
class StdDataExport
	:  	public StdUserData<	StdDataExport<TImpl,TData,dim>,
	  						DataExport<TData,dim>,
	  						TData,dim>
{
	public:
		StdDataExport() : m_pFctPatt(NULL) {m_SymbFct.clear();}

		inline void evaluate (TData& value,
		                      const MathVector<dim>& globIP,
		                      number time, int si) const
		{
			UG_THROW("StdDataExport: Solution, element and local ips required "
					"for evaluation, but not passed. Cannot evaluate.");
		}

		inline void evaluate (TData vValue[],
		                      const MathVector<dim> vGlobIP[],
		                      number time, int si, const size_t nip) const
		{
			UG_THROW("StdDataExport: Solution, element and local ips required "
					"for evaluation, but not passed. Cannot evaluate.");
		}

		template <int refDim>
		inline void evaluate (TData& value,
		                      const MathVector<dim>& globIP,
		                      number time, int si,
		                      LocalVector& u,
		                      GeometricObject* elem,
		                      const MathVector<dim> vCornerCoords[],
		                      const MathVector<refDim>& locIP) const
		{
			getImpl().template evaluate<refDim>(&value,&globIP,time,si,u,elem,
			                                    vCornerCoords,&locIP,1,NULL);
		}

		template <int refDim>
		inline void evaluate(TData vValue[],
		                     const MathVector<dim> vGlobIP[],
		                     number time, int si,
		                     LocalVector& u,
		                     GeometricObject* elem,
		                     const MathVector<dim> vCornerCoords[],
		                     const MathVector<refDim> vLocIP[],
		                     const size_t nip,
		                     const MathMatrix<refDim, dim>* vJT = NULL) const
		{
			getImpl().template evaluate<refDim>(vValue,vGlobIP,time,si,u,elem,
			                                    vCornerCoords,vLocIP,nip, vJT);
		}

	///	returns that a grid function is needed for evaluation
		virtual bool requires_grid_fct() const {return true;}

	///	sets the associated function pattern
		virtual void set_function_pattern(const FunctionPattern& fctPatt)
		{
			m_pFctPatt = &fctPatt;
			extract_fct_grp();
		}

	///	sets the associated symbolic functions
		void set_symb_fct(const char* symbFct)
		{
			m_SymbFct = symbFct;
			extract_fct_grp();
		}

	protected:
	///	extracts the function group
		void extract_fct_grp()
		{
		//	if associated infos missing return
			if(m_pFctPatt == NULL || m_SymbFct.empty()) return;

		//	create function group of this elem disc
			try{
				m_FctGrp.set_function_pattern(*m_pFctPatt);
				m_FctGrp.add(TokenizeString(m_SymbFct));
			}UG_CATCH_THROW("StdDataExport: Cannot find  some symbolic function "
							"name in '"<<m_SymbFct<<"'.");

		//	create a mapping between all functions and the function group of this
		//	element disc.
			try{
				CreateFunctionIndexMapping(m_FctIndexMap, m_FctGrp, *m_pFctPatt);
			}UG_CATCH_THROW("StdDataExport: Cannot create Function Index Mapping"
							" for '"<<m_SymbFct<<"'.");
		}

	///	returns the function group
		const FunctionGroup& fct_grp() const {return m_FctGrp;}

	///	returns the function index mapping
		const FunctionIndexMapping& fct_index_map() const {return m_FctIndexMap;}

		protected:
	///	associated function pattern
		const FunctionPattern* m_pFctPatt;

	///	string of symbolic functions required
		std::string m_SymbFct;

	///	FunctionGroup corresponding to symb functions
		FunctionGroup m_FctGrp;

	///	associated function index mapping
		FunctionIndexMapping m_FctIndexMap;

	protected:
	///	access to implementation
		TImpl& getImpl() {return static_cast<TImpl&>(*this);}

	///	const access to implementation
		const TImpl& getImpl() const {return static_cast<const TImpl&>(*this);}
};

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
