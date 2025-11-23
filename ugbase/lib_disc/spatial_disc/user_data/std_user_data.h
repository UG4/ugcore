/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__STD_USER_DATA__
#define __H__UG__LIB_DISC__SPATIAL_DISC__STD_USER_DATA__

#include "user_data.h"

namespace ug{

///////////////////////////////////////////////////////////////////////////////
// Wrapper Class for UserData
///////////////////////////////////////////////////////////////////////////////

/**
 * This class is used as a wrapper class for user data in order to ease the
 * implementation of the virtual evaluation operators through templated
 * methods evaluate.
 *
 * template \<int refDim\>
 * inline TRet evaluate(TData vValue[],
		                const MathVector\<dim\> vGlobIP[],
		                number time, int si,
		                GridObject* elem,
		                const MathVector\<dim\> vCornerCoords[],
		                const MathVector\<refDim\> vLocIP[],
		                const size_t nip,
		                LocalVector* u,
		                const MathMatrix\<refDim, dim\>* vJT = nullptr) const
 *
 */
template <typename TImpl, typename TData, int dim, typename TRet = void, typename TBase = CplUserData<TData, dim, TRet> >
class StdUserData : public TBase
{
	public:
	///	returns value for a global position
		virtual TRet operator () (TData& value,
								 const MathVector<dim>& globIP,
								 number time, int si) const = 0;

	///	returns value for global positions
		virtual void operator () (TData vValue[],
								const MathVector<dim> vGlobIP[],
								number time, int si, const size_t nip) const = 0;

	///	returns values for local and global positions
	///	\{
		virtual void operator () (TData vValue[],
		                        const MathVector<dim> vGlobIP[],
		                        number time, int si,
		                        GridObject* elem,
		                        const MathVector<dim> vCornerCoords[],
		                        const MathVector<1> vLocIP[],
		                        const size_t nip,
		                        LocalVector* u,
		                        const MathMatrix<1, dim>* vJT = nullptr) const
		{
			getImpl().template evaluate<1>(vValue,vGlobIP,time,si,elem,
										   vCornerCoords,vLocIP,nip,u,vJT);
		}

		virtual void operator () (TData vValue[],
		                        const MathVector<dim> vGlobIP[],
		                        number time, int si,
		                        GridObject* elem,
		                        const MathVector<dim> vCornerCoords[],
		                        const MathVector<2> vLocIP[],
		                        const size_t nip,
		                        LocalVector* u,
		                        const MathMatrix<2, dim>* vJT = nullptr) const
		{
			getImpl().template evaluate<2>(vValue,vGlobIP,time,si,elem,
										   vCornerCoords,vLocIP,nip,u,vJT);
		}

		virtual void operator () (TData vValue[],
		                        const MathVector<dim> vGlobIP[],
		                        number time, int si,
		                        GridObject* elem,
		                        const MathVector<dim> vCornerCoords[],
		                        const MathVector<3> vLocIP[],
		                        const size_t nip,
		                        LocalVector* u,
		                        const MathMatrix<3, dim>* vJT = nullptr) const
		{
			getImpl().template evaluate<3>(vValue,vGlobIP,time,si,elem,
										   vCornerCoords,vLocIP,nip,u,vJT);
		}

	///	\}
	protected:
	///	access to implementation
		TImpl& getImpl() {return static_cast<TImpl&>(*this);}

	///	const access to implementation
		const TImpl& getImpl() const {return static_cast<const TImpl&>(*this);}
};

template <typename TImpl, typename TData, int dim>
class StdDependentUserData
	: public StdUserData< StdDependentUserData<TImpl, TData, dim>,
	  	  	  	  	  	  TData, dim, void,
	  	  	  	  	  	  DependentUserData<TData, dim> >
{
	public:
		StdDependentUserData() = default;

		StdDependentUserData(const char* symbFct) {this->set_functions(symbFct);}
		StdDependentUserData(const std::string& symbFct) {this->set_functions(symbFct);}
		StdDependentUserData(const std::vector<std::string>& symbFct) {this->set_functions(symbFct);}

	public:
		void operator () (TData& value,
		                  const MathVector<dim>& globIP,
		                  number time, int si) const override {
			UG_THROW("StdDependentUserData: Solution, element and local ips required "
					"for evaluation, but not passed. Cannot evaluate.");
		}

		void operator () (TData vValue[],
		                  const MathVector<dim> vGlobIP[],
		                  number time, int si, const size_t nip) const override {
			UG_THROW("StdDependentUserData: Solution, element and local ips required "
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
							 const MathMatrix<refDim, dim>* vJT = nullptr) const
		{
			const_cast<TImpl*>(static_cast<const TImpl*>(this))->
					template eval_and_deriv<refDim>(vValue,vGlobIP,time,si,elem,
					                                vCornerCoords,vLocIP,nip,u,
					                                false,0,nullptr,vJT);
		}

		template <int refDim>
		void eval_deriv(LocalVector* u, GridObject* elem,
		                const MathVector<dim> vCornerCoords[], bool bDeriv = false) {

			const int si = this->subset();

			std::vector<std::vector<TData> >* vvvDeriv = nullptr;

			for(size_t s = 0; s < this->num_series(); ++s){
				
				if(bDeriv && this->m_vvvvDeriv[s].size() > 0)
					vvvDeriv = &this->m_vvvvDeriv[s][0];
				else
					vvvDeriv = nullptr;

				getImpl().template eval_and_deriv<refDim>(this->values(s), this->ips(s), this->time(s), si,
				                                 elem, vCornerCoords,
				                                 this->template local_ips<refDim>(s), this->num_ip(s),
				                                 u, bDeriv, s, vvvDeriv);
			}
		}

		template <int refDim>
		void eval_deriv(LocalVectorTimeSeries* u, GridObject* elem,
		                const MathVector<dim> vCornerCoords[], bool bDeriv = false) {

			const int si = this->subset();

			std::vector<std::vector<TData> >* vvvDeriv = nullptr;

			for(size_t s = 0; s < this->num_series(); ++s){
				
				bool bDoDeriv = bDeriv && this->at_current_time (s); // derivatives only for the 'current' time point!

				if(bDoDeriv && this->m_vvvvDeriv[s].size() > 0)
					vvvDeriv = &this->m_vvvvDeriv[s][0];
				else
					vvvDeriv = nullptr;

				getImpl().template eval_and_deriv<refDim>(this->values(s), this->ips(s), this->time(s), si,
				                                 elem, vCornerCoords,
				                                 this->template local_ips<refDim>(s), this->num_ip(s),
				                                 &(u->solution(this->time_point(s))), bDoDeriv, s, vvvDeriv);
			}
		}

		void compute(LocalVector* u, GridObject* elem,
							 const MathVector<dim> vCornerCoords[],
							 bool bDeriv = false) override
	{

			UG_ASSERT(elem->base_object_id() == this->dim_local_ips(),
			          "local ip dimension (" << this->dim_local_ips()
					  << ") and reference element dimension ("
					  << elem->base_object_id() << ") mismatch.");

			switch(this->dim_local_ips()){
				case 1: eval_deriv<1>(u,elem,vCornerCoords,bDeriv); break;
				case 2: eval_deriv<2>(u,elem,vCornerCoords,bDeriv); break;
				case 3: eval_deriv<3>(u,elem,vCornerCoords,bDeriv); break;
				default: UG_THROW("StdDependentUserData: Dimension not supported.");
			}
		}

		void compute(LocalVectorTimeSeries* u, GridObject* elem,
							 const MathVector<dim> vCornerCoords[],
							 bool bDeriv = false) override {

			UG_ASSERT(elem->base_object_id() == this->dim_local_ips(),
			          "local ip dimension and reference element dimension mismatch.");

			switch(this->dim_local_ips()){
				case 1: eval_deriv<1>(u,elem,vCornerCoords,bDeriv); break;
				case 2: eval_deriv<2>(u,elem,vCornerCoords,bDeriv); break;
				case 3: eval_deriv<3>(u,elem,vCornerCoords,bDeriv); break;
				default: UG_THROW("StdDependentUserData: Dimension not supported.");
			}
		}

	protected:
	///	access to implementation
		TImpl& getImpl() {return static_cast<TImpl&>(*this);}

	///	const access to implementation
		const TImpl& getImpl() const {return static_cast<const TImpl&>(*this);}

};

} // namespace ug

#endif