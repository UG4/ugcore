/*
 * std_data_export.h
 *
 *  Created on: 22.05.2012
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__STD_DATA_IMPORT__
#define __H__UG__LIB_DISC__SPATIAL_DISC__STD_DATA_IMPORT__

#include "data_import_export.h"
#include "lib_disc/common/function_group.h"
#include "lib_disc/common/local_algebra.h"
#include "lib_disc/common/groups_util.h"

namespace ug{


template <typename TData, int dim, typename TImpl>
class StdDataExport
	: 	public DataExport<TData,dim>
{
	public:
		StdDataExport() : m_pFctPatt(NULL) {m_SymbFct.clear();}

		virtual void operator() (TData& value,
		                         const MathVector<dim>& globIP,
		                         number time, int si) const
		{
			UG_THROW("StdDataExport: Solution, element and local ips required "
					"for evaluation, but not passed. Cannot evaluate.");
		}

		virtual void operator() (TData value[],
		                         const MathVector<dim> globIP[],
		                         number time, int si, const size_t nip) const
		{
			UG_THROW("StdDataExport: Solution, element and local ips required "
					"for evaluation, but not passed. Cannot evaluate.");
		}

		////////////////
		// one value
		////////////////

		virtual void operator() (TData& value,
		                         const MathVector<dim>& globIP,
		                         number time, int si,
		                         LocalVector& u,
		                         GeometricObject* elem,
		                         const MathVector<dim> vCornerCoords[],
		                         const MathVector<1>& locIP) const
		{
			getImpl().template evaluate<1>(value,globIP,time,si,u,elem,vCornerCoords,locIP);
		}

		virtual void operator() (TData& value,
		                         const MathVector<dim>& globIP,
		                         number time, int si,
		                         LocalVector& u,
		                         GeometricObject* elem,
		                         const MathVector<dim> vCornerCoords[],
		                         const MathVector<2>& locIP) const
		{
			getImpl().template evaluate<2>(value,globIP,time,si,u,elem,vCornerCoords,locIP);
		}

		virtual void operator() (TData& value,
		                         const MathVector<dim>& globIP,
		                         number time, int si,
		                         LocalVector& u,
		                         GeometricObject* elem,
		                         const MathVector<dim> vCornerCoords[],
		                         const MathVector<3>& locIP) const
		{
			getImpl().template evaluate<3>(value,globIP,time,si,u,elem,vCornerCoords,locIP);
		}

		////////////////
		// vector of values
		////////////////

		virtual void operator()(TData vValue[],
		                        const MathVector<dim> vGlobIP[],
		                        number time, int si,
		                        LocalVector& u,
		                        GeometricObject* elem,
		                        const MathVector<dim> vCornerCoords[],
		                        const MathVector<1> vLocIP[],
		                        const size_t nip,
		                        const MathMatrix<1, dim>* vJT = NULL) const
		{
			getImpl().template evaluate<1>(vValue,vGlobIP,time,si,u,elem,
			                               vCornerCoords,vLocIP,nip, vJT);
		}

		virtual void operator()(TData vValue[],
		                        const MathVector<dim> vGlobIP[],
		                        number time, int si,
		                        LocalVector& u,
		                        GeometricObject* elem,
		                        const MathVector<dim> vCornerCoords[],
		                        const MathVector<2> vLocIP[],
		                        const size_t nip,
		                        const MathMatrix<2, dim>* vJT = NULL) const
		{
			getImpl().template evaluate<2>(vValue,vGlobIP,time,si,u,elem,
			                               vCornerCoords,vLocIP,nip, vJT);
		}

		virtual void operator()(TData vValue[],
		                        const MathVector<dim> vGlobIP[],
		                        number time, int si,
		                        LocalVector& u,
		                        GeometricObject* elem,
		                        const MathVector<dim> vCornerCoords[],
		                        const MathVector<3> vLocIP[],
		                        const size_t nip,
		                        const MathMatrix<3, dim>* vJT = NULL) const
		{
			getImpl().template evaluate<3>(vValue,vGlobIP,time,si,u,elem,
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
				ConvertStringToFunctionGroup(m_FctGrp, *m_pFctPatt, m_SymbFct.c_str());
			}UG_CATCH_THROW("StdDataExport: Cannot find  some symbolic function "
							"name in '"<<m_SymbFct<<"'.");

		//	get common fct grp
			FunctionGroup commonFctGroup(*m_pFctPatt);
			commonFctGroup.add_all();

		//	create a mapping between all functions and the function group of this
		//	element disc.
			try{
				CreateFunctionIndexMapping(m_FctIndexMap, m_FctGrp, commonFctGroup);
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
	: public StdDataExport<number,dim,ValueDataExport<dim> >
{
	public:
		ValueDataExport(const char* functions){this->set_symb_fct(functions);}

		template <int refDim>
		inline void evaluate (number& value,
		                        const MathVector<dim>& globIP,
		                        number time, int si,
		                        LocalVector& u,
		                        GeometricObject* elem,
							    const MathVector<dim> vCornerCoords[],
		                        const MathVector<refDim>& locIP) const
		{
			UG_ASSERT(this->fct_grp().size() == 1, "Exactly one fct needed.")

		//	reference object id
			const ReferenceObjectID roid = elem->reference_object_id();

		//	local finite element id
			const LFEID& lfeID = this->fct_grp().local_finite_element_id(_C_);

		//	access local vector by map
			u.access_by_map(this->fct_index_map());

		//	get trial space
			try{
			const DimLocalShapeFunctionSet<refDim>& rTrialSpace =
					LocalShapeFunctionSetProvider::get<refDim>(roid, lfeID);

		//	memory for shapes
			std::vector<number> vShape;

		//	evaluate at shapes at ip
			rTrialSpace.shapes(vShape, locIP);

		//	compute concentration at ip
			value = 0.0;
			for(size_t sh = 0; sh < vShape.size(); ++sh)
				value += u(_C_, sh) * vShape[sh];

			}catch(UG_ERROR_LocalShapeFunctionSetNotRegistered& ex){
				UG_THROW("ValueDataExport: "<< ex.get_msg()<<", Reference Object: "
				         <<roid<<", Trial Space: "<<lfeID<<", refDim="<<refDim);
			}
		}

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

		//	request for trial space
			try{
			const DimLocalShapeFunctionSet<refDim>& rTrialSpace
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

			}catch(UG_ERROR_LocalShapeFunctionSetNotRegistered& ex){
				UG_THROW("ValueDataExport: "<< ex.get_msg()<<", Reference Object: "
				         <<roid<<", Trial Space: "<<lfeID<<", refDim="<<refDim);
			}
		}

	///	returns if provided data is continuous over geometric object boundaries
		virtual bool is_continuous() const
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
	: public StdDataExport<MathVector<dim>,dim,GradientDataExport<dim> >
{
	public:
		GradientDataExport(const char* functions){this->set_symb_fct(functions);}

		template <int refDim>
		inline void evaluate (MathVector<dim>& value,
		                      const MathVector<dim>& globIP,
		                      number time, int si,
		                      LocalVector& u,
		                      GeometricObject* elem,
		                      const MathVector<dim> vCornerCoords[],
		                      const MathVector<refDim>& locIP) const
		{
		//	reference object id
			const ReferenceObjectID roid = elem->reference_object_id();

		//	local finite element id
			const LFEID& lfeID = this->fct_grp().local_finite_element_id(_C_);

			try{

		//	local shape function set
			const DimLocalShapeFunctionSet<refDim>& rTrialSpace
				 = LocalShapeFunctionSetProvider::get<refDim>(roid, lfeID);

		//	storage for shape function at ip
			std::vector<MathVector<refDim> > vLocGrad;
			MathVector<refDim> locGrad;

		//	evaluate at shapes at ip
			rTrialSpace.grads(vLocGrad, locIP);

		//	compute grad at ip
			VecSet(locGrad, 0.0);
			for(size_t sh = 0; sh < vLocGrad.size(); ++sh)
				VecScaleAppend(locGrad, u(_C_, sh), vLocGrad[sh]);

		//	reference mapping
			DimReferenceMapping<refDim, dim>& map
				= ReferenceMappingProvider::get<refDim, dim>(roid);

		//	update the mapping for the new corners
			map.update(vCornerCoords);

		// 	compute transformation inverse and determinate at ip
			MathMatrix<dim, refDim> JTInv;
			map.jacobian_transposed_inverse(JTInv, locIP);

		//	compute global gradient
			MatVecMult(value, JTInv, locGrad);

			}catch(UG_ERROR_LocalShapeFunctionSetNotRegistered& ex){
				UG_THROW("GradientDataExport: "<< ex.get_msg()<<", Reference Object: "
						 <<roid<<", Trial Space: "<<lfeID<<", refDim="<<refDim);
			}
			catch(UG_ERROR_ReferenceMappingMissing& ex){
				UG_THROW("GradientDataExport: " << ex.get_msg());
			}
		}

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
			//	must pass vJT
				UG_ASSERT(vJT != NULL, "Jacobian transposed needed.");

			//	reference object id
				const ReferenceObjectID roid = elem->reference_object_id();

			//	local finite element id
				const LFEID& lfeID = this->fct_grp().local_finite_element_id(_C_);

			//	request for trial space
				try{
				const DimLocalShapeFunctionSet<refDim>& rTrialSpace
					 = LocalShapeFunctionSetProvider::get<refDim>(roid, lfeID);

			//	storage for shape function at ip
				std::vector<MathVector<refDim> > vLocGrad;
				MathVector<refDim> locGrad;

			//	Reference Mapping
				MathMatrix<dim, refDim> JTInv;

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

				}catch(UG_ERROR_LocalShapeFunctionSetNotRegistered& ex){
					UG_THROW("GradientDataExport: "<< ex.get_msg()<<", Reference Object: "
					         <<roid<<", Trial Space: "<<lfeID<<", refDim="<<refDim);
				}
		}

	///	returns if provided data is continuous over geometric object boundaries
		virtual bool is_continuous() const
		{
			return false;
		}

	protected:
	//	abbreviation for component
		static const int _C_ = 0;
};

} // end namespace ug

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__STD_DATA_IMPORT__ */
