/*
 * Copyright (c) 2014-2015:  G-CSC, Goethe University Frankfurt
 * Author: Arne Nägel
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

#ifndef __H__UG__LIB_DISC__FUNCTION_SPACE__CONST_GRID_FUNCTION_USER_DATA__
#define __H__UG__LIB_DISC__FUNCTION_SPACE__CONST_GRID_FUNCTION_USER_DATA__


#include <map>
#include <string>

#include "common/common.h"

#include "lib_grid/tools/subset_group.h"

#include "lib_disc/common/function_group.h"
#include "lib_disc/common/groups_util.h"
#include "lib_disc/quadrature/quadrature.h"
#include "lib_disc/local_finite_element/local_finite_element_provider.h"
#include "lib_disc/spatial_disc/user_data/std_user_data.h"
#include "lib_disc/reference_element/reference_mapping_provider.h"


namespace ug{




/**
 * TImpl: CRTP base impl
 * TData: data type (number, ....)
 * TGridFunction
 */
template <typename TImpl, typename TData, typename TGridFunction>
class StdExplicitGridFunctionData
: public StdUserData<StdExplicitGridFunctionData<TImpl,TData, TGridFunction>, TData, TGridFunction::dim>
{

protected:
	///	access to implementation
		TImpl& getImpl() {return static_cast<TImpl&>(*this);}

	///	const access to implementation
		const TImpl& getImpl() const {return static_cast<const TImpl&>(*this);}
protected:
		SmartPtr<TGridFunction> m_spGridFct; ///< grid function
		size_t m_fct;                        ///< component of function
		LFEID m_lfeID;


		//	using base_type::get_dim;
	public:
		//	world dimension of grid function
		static const int dim = TGridFunction::dim;

		///	returns if data is constant
		virtual bool constant() const {return false;}

		///	returns if grid function is needed for evaluation
		/** (true, since local coordinates may not be sufficient)*/
		virtual bool requires_grid_fct() const {return true;}

		///	returns if provided data is continuous over geometric object boundaries
		//virtual bool continuous() const {return true;}
		virtual bool continuous() const {return getImpl().continuous(); }

		template <int refDim>
		void eval(LocalVector* u, GridObject* elem,
				const MathVector<dim> vCornerCoords[], bool bDeriv = false) {

			const int si = this->subset();


			for(size_t s = 0; s < this->num_series(); ++s)
			{
				getImpl().template evaluate<refDim>(this->values(s), this->ips(s), this->time(s), si,
						elem, vCornerCoords,
						this->template local_ips<refDim>(s), this->num_ip(s),
						u);
			}
		}

		template <int refDim>
		void eval(LocalVectorTimeSeries* u, GridObject* elem,
				const MathVector<dim> vCornerCoords[], bool bDeriv = false) {

			const int si = this->subset();


			for(size_t s = 0; s < this->num_series(); ++s)
			{
				getImpl().template evaluate<refDim>(this->values(s), this->ips(s), this->time(s), si,
						elem, vCornerCoords,
						this->template local_ips<refDim>(s), this->num_ip(s),
						&(u->solution(this->time_point(s))));
			}
		}

		virtual void compute(LocalVector* u, GridObject* elem,
				const MathVector<dim> vCornerCoords[], bool bDeriv = false){

			UG_ASSERT(elem->base_object_id() == this->dim_local_ips(),
					"local ip dimension and reference element dimension mismatch.");

			switch(this->dim_local_ips()){
			case 1: eval<1>(u,elem,vCornerCoords,bDeriv); break;
			case 2: eval<2>(u,elem,vCornerCoords,bDeriv); break;
			case 3: eval<3>(u,elem,vCornerCoords,bDeriv); break;
			default: UG_THROW("StdDependentUserData: Dimension not supported.");
			}
		}

		virtual void compute(LocalVectorTimeSeries* u, GridObject* elem,
				const MathVector<dim> vCornerCoords[], bool bDeriv = false){

			UG_ASSERT(elem->base_object_id() == this->dim_local_ips(),
					"local ip dimension and reference element dimension mismatch.");

			switch(this->dim_local_ips()){
			case 1: eval<1>(u,elem,vCornerCoords,bDeriv); break;
			case 2: eval<2>(u,elem,vCornerCoords,bDeriv); break;
			case 3: eval<3>(u,elem,vCornerCoords,bDeriv); break;
			default: UG_THROW("StdDependentUserData: Dimension not supported.");
			}
		}

		/// common constructor
		 StdExplicitGridFunctionData(SmartPtr<TGridFunction> spGridFct, const char* cmp)
		: m_spGridFct(spGridFct)
		{
			//	get function id of name
			m_fct = spGridFct->fct_id_by_name(cmp);

			//	check that function exists
			if(m_fct >= spGridFct->num_fct())
				UG_THROW("ExplicitGridFunctionValueData: Function space does not contain"
						" a function with name " << cmp << ".");

			//	local finite element id
			m_lfeID = spGridFct->local_finite_element_id(m_fct);
		};

	/*	virtual bool continuous() const
		{
			return LocalFiniteElementProvider::continuous(m_lfeID);
		}*/

		virtual void operator() (TData& value,
								const MathVector<dim>& globIP,
								number time, int si) const
		{
			UG_THROW("StdExplicitGridFunctionData: Solution, element and local ips required "
					 "for evaluation, but not passed. Cannot evaluate.");
		}

		virtual void operator()(TData vValue[],
								const MathVector<dim> vGlobIP[],
								number time, int si, const size_t nip) const
		{
			UG_THROW("StdExplicitGridFunctionData: Solution, element and local ips required "
					 "for evaluation, but not passed. Cannot evaluate.");
		}

		template <int refDim>
		inline void evaluate(TData vValue[],
							 const MathVector<dim> vGlobIP[],
							 number time, int si,
							 GridObject* elem,
							 const MathVector<dim> vCornerCoords[],
							 const MathVector<refDim> vLocIP[],
							 const size_t nip,
							 LocalVector* u,
							 const MathMatrix<refDim, dim>* vJT = NULL) const
		{
			// forward
			getImpl().template evaluate<refDim> (vValue, vGlobIP, time, si,
			elem, vCornerCoords, vLocIP, nip, u, vJT);
		}
};

/// This grid function can be used as an explicit input (w/o additional derivativ)
template<typename TGridFunction>
class ExplicitGridFunctionValue
: public StdExplicitGridFunctionData<ExplicitGridFunctionValue<TGridFunction>, number, TGridFunction>
{
public:
	//	world dimension of grid function
	static const int dim = TGridFunction::dim;
	typedef StdExplicitGridFunctionData<ExplicitGridFunctionValue<TGridFunction>, number, TGridFunction> base_type;

	using base_type::m_spGridFct;
	using base_type::m_lfeID;
	using base_type::m_fct;




	// constructor
	ExplicitGridFunctionValue(SmartPtr<TGridFunction> gf, const char* cmp)
	//: StdExplicitGridFunctionData<ExplicitGridFunctionValue<TGridFunction>, number, TGridFunction>(gf, cmp)
	: base_type(gf, cmp)
	{}


	/*template <int refDim>
				void eval_and_deriv(number vValue[],
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
				                    const MathMatrix<refDim, dim>* vJT = NULL) const
				*/

	template <int refDim>
	inline void evaluate(number vValue[],
					     const MathVector<dim> vGlobIP[],
					     number time, int si,
					     GridObject* elem,
					     const MathVector<dim> vCornerCoords[],
					     const MathVector<refDim> vLocIP[],
					     const size_t nip,
					     LocalVector* u,
					     const MathMatrix<refDim, dim>* vJT = NULL) const
	{
		//	reference object id
		const ReferenceObjectID roid = elem->reference_object_id();

		//	get trial space
		try{
			const LocalShapeFunctionSet<refDim>& rTrialSpace =
					LocalFiniteElementProvider::get<refDim>(roid, base_type::m_lfeID);

			//	memory for shapes
			std::vector<number> vShape;

			//	get multiindices of element
			std::vector<DoFIndex> ind;
			m_spGridFct->dof_indices(elem, m_fct, ind);

          
			//	loop ips
			for(size_t ip = 0; ip < nip; ++ip)
			{
				//	evaluate at shapes at ip
				rTrialSpace.shapes(vShape, vLocIP[ip]);

				// 	compute solution at integration point
				vValue[ip] = 0.0;
				for(size_t sh = 0; sh < vShape.size(); ++sh)
				{
					const number valSH = DoFRef(*m_spGridFct, ind[sh]);
					vValue[ip] += valSH * vShape[sh];
				}
			}

		}
		UG_CATCH_THROW("ExplicitGridFunctionValue: Shape Function Set missing for"
				" Reference Object: "<<roid<<", Trial Space: "
				<<m_lfeID<<", refDim="<<refDim);

	}

	bool continuous() const
	{return LocalFiniteElementProvider::continuous(m_lfeID);}

};



template <typename TGridFunction>
class ExplicitGridFunctionGradient
: public StdExplicitGridFunctionData<ExplicitGridFunctionGradient<TGridFunction>, MathVector<TGridFunction::dim>, TGridFunction >
{
	public:
	//	world dimension of grid function
	static const int dim = TGridFunction::dim;
	typedef StdExplicitGridFunctionData<ExplicitGridFunctionGradient<TGridFunction>, MathVector<TGridFunction::dim>, TGridFunction > base_type;
//	typedef CplUserData<MathVector<TGridFunction::dim>, TGridFunction::dim> user_data_type;

	using base_type::m_spGridFct;
	using base_type::m_lfeID;
	using base_type::m_fct;

	// specify final overriders


	/// constructor
	ExplicitGridFunctionGradient(SmartPtr<TGridFunction> spGridFct, const char* cmp)
	: base_type(spGridFct, cmp), m_diffCoeffMap()
	{};

	// evaluate gradient
	template <int refDim>
	void evaluate(MathVector<dim> vValue[],
		                    const MathVector<dim> vGlobIP[],
		                    number time, int si,
		                    GridObject* elem,
		                    const MathVector<dim> vCornerCoords[],
		                    const MathVector<refDim> vLocIP[],
		                    const size_t nip,
		                    LocalVector* u,
		                    const MathMatrix<refDim, dim>* vJT = NULL) const
		{
			//	reference object id
			const ReferenceObjectID roid = elem->reference_object_id();

			//	get reference element mapping by reference object id
			std::vector<MathMatrix<refDim, dim> > vJTTmp(nip);
			if(vJT == NULL){
				try{
					DimReferenceMapping<refDim, dim>& mapping
					= ReferenceMappingProvider::get<refDim, dim>(roid, vCornerCoords);

					//	compute transformation matrices
					mapping.jacobian_transposed(&(vJTTmp[0]), vLocIP, nip);

					//	store tmp Gradient
					vJT = &(vJTTmp[0]);
				}UG_CATCH_THROW("GridFunctionGradientData: failed.");
			}

			// scale with coeficient
			const int   subsetInd  = m_spGridFct->domain()->subset_handler()->get_subset_index(elem);
			const char* subsetName = m_spGridFct->domain()->subset_handler()->get_subset_name(subsetInd);

			double diffCoeff = get_subset_coeff(std::string(subsetName));

			//	get trial space
			try{
				const LocalShapeFunctionSet<refDim>& rTrialSpace =
						LocalFiniteElementProvider::get<refDim>(roid, m_lfeID);

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

					//	get multiindices of element
					std::vector<DoFIndex > ind;
					m_spGridFct->dof_indices(elem, m_fct, ind);

					//	compute grad at ip
					VecSet(locGrad, 0.0);
					for(size_t sh = 0; sh < vLocGrad.size(); ++sh)
					{
						const number valSH = DoFRef( *m_spGridFct, ind[sh]);
						VecScaleAppend(locGrad, valSH, vLocGrad[sh]);
					}

					Inverse(JTInv, vJT[ip]);
					MatVecMult(vValue[ip], JTInv, locGrad);

					// scale with diff coeff (TODO: matrix)
					VecScale(vValue[ip], vValue[ip], diffCoeff);
				}
			}
			UG_CATCH_THROW("ExplicitGridFunctionGradient: Shape Function Set missing for"
					" Reference Object: "<<roid<<", Trial Space: "
					<<m_lfeID<<", refDim="<<refDim);
		}

	bool continuous() const
	{ return false; }

	void add_subset_coeff(const std::string &key, double val)
	{
		m_diffCoeffMap.insert(std::pair<std::string, double>(key, val));
	}

	double get_subset_coeff(const std::string &key) const
	{
		std::map<std::string, double>::const_iterator it = m_diffCoeffMap.find(key);
		double val = (it != m_diffCoeffMap.end()) ? it->second : 1.0;
		return val;
	}

	protected:
		std::map<std::string, double> m_diffCoeffMap;
};


} // end namespace ug

#endif /* __H__UG__LIB_DISC__FUNCTION_SPACE__GRID_FUNCTION_USER_DATA__ */
