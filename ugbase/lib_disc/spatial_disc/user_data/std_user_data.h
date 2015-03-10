/*
 * std_user_data.h
 *
 *  Created on: 03.07.2012
 *      Author: andreasvogel
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
		                const MathMatrix\<refDim, dim\>* vJT = NULL) const
 *
 */
template <typename TImpl, typename TData, int dim, typename TRet = void, typename TBase = CplUserData<TData, dim, TRet> >
class StdUserData : public TBase
{
	public:
	///	returns value for a global position
		virtual TRet operator() (TData& value,
								 const MathVector<dim>& globIP,
								 number time, int si) const = 0;

	///	returns value for global positions
		virtual void operator()(TData vValue[],
								const MathVector<dim> vGlobIP[],
								number time, int si, const size_t nip) const = 0;

	///	returns values for local and global positions
	///	\{
		virtual void operator()(TData vValue[],
		                        const MathVector<dim> vGlobIP[],
		                        number time, int si,
		                        GridObject* elem,
		                        const MathVector<dim> vCornerCoords[],
		                        const MathVector<1> vLocIP[],
		                        const size_t nip,
		                        LocalVector* u,
		                        const MathMatrix<1, dim>* vJT = NULL) const
		{
			getImpl().template evaluate<1>(vValue,vGlobIP,time,si,elem,
										   vCornerCoords,vLocIP,nip,u,vJT);
		}

		virtual void operator()(TData vValue[],
		                        const MathVector<dim> vGlobIP[],
		                        number time, int si,
		                        GridObject* elem,
		                        const MathVector<dim> vCornerCoords[],
		                        const MathVector<2> vLocIP[],
		                        const size_t nip,
		                        LocalVector* u,
		                        const MathMatrix<2, dim>* vJT = NULL) const
		{
			getImpl().template evaluate<2>(vValue,vGlobIP,time,si,elem,
										   vCornerCoords,vLocIP,nip,u,vJT);
		}

		virtual void operator()(TData vValue[],
		                        const MathVector<dim> vGlobIP[],
		                        number time, int si,
		                        GridObject* elem,
		                        const MathVector<dim> vCornerCoords[],
		                        const MathVector<3> vLocIP[],
		                        const size_t nip,
		                        LocalVector* u,
		                        const MathMatrix<3, dim>* vJT = NULL) const
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
		StdDependentUserData(){}

		StdDependentUserData(const char* symbFct) {this->set_functions(symbFct);}
		StdDependentUserData(const std::string& symbFct) {this->set_functions(symbFct);}
		StdDependentUserData(const std::vector<std::string>& symbFct) {this->set_functions(symbFct);}

	public:
		virtual void operator() (TData& value,
								 const MathVector<dim>& globIP,
								 number time, int si) const
		{
			UG_THROW("StdDependentUserData: Solution, element and local ips required "
					"for evaluation, but not passed. Cannot evaluate.");
		}

		virtual void operator()(TData vValue[],
								const MathVector<dim> vGlobIP[],
								number time, int si, const size_t nip) const
		{
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
							 const MathMatrix<refDim, dim>* vJT = NULL) const
		{
			const_cast<TImpl*>(static_cast<const TImpl*>(this))->
					template eval_and_deriv<refDim>(vValue,vGlobIP,time,si,elem,
					                                vCornerCoords,vLocIP,nip,u,
					                                false,0,NULL,vJT);
		}

		template <int refDim>
		void eval_deriv(LocalVector* u, GridObject* elem,
		                const MathVector<dim> vCornerCoords[], bool bDeriv = false) {

			const int si = this->subset();

			std::vector<std::vector<TData> >* vvvDeriv = NULL;

			for(size_t s = 0; s < this->num_series(); ++s){
				
				if(bDeriv && this->m_vvvvDeriv[s].size() > 0)
					vvvDeriv = &this->m_vvvvDeriv[s][0];
				else
					vvvDeriv = NULL;

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

			std::vector<std::vector<TData> >* vvvDeriv = NULL;

			for(size_t s = 0; s < this->num_series(); ++s){
				
				bool bDoDeriv = bDeriv && this->at_current_time (s); // derivatives only for the 'current' time point!

				if(bDoDeriv && this->m_vvvvDeriv[s].size() > 0)
					vvvDeriv = &this->m_vvvvDeriv[s][0];
				else
					vvvDeriv = NULL;

				getImpl().template eval_and_deriv<refDim>(this->values(s), this->ips(s), this->time(s), si,
				                                 elem, vCornerCoords,
				                                 this->template local_ips<refDim>(s), this->num_ip(s),
				                                 &(u->solution(this->time_point(s))), bDoDeriv, s, vvvDeriv);
			}
		}

		virtual void compute(LocalVector* u, GridObject* elem,
							 const MathVector<dim> vCornerCoords[], bool bDeriv = false){

			UG_ASSERT(elem->base_object_id() == this->dim_local_ips(),
			          "local ip dimension and reference element dimension mismatch.");

			switch(this->dim_local_ips()){
				case 1: eval_deriv<1>(u,elem,vCornerCoords,bDeriv); break;
				case 2: eval_deriv<2>(u,elem,vCornerCoords,bDeriv); break;
				case 3: eval_deriv<3>(u,elem,vCornerCoords,bDeriv); break;
				default: UG_THROW("StdDependentUserData: Dimension not supported.");
			}
		}

		virtual void compute(LocalVectorTimeSeries* u, GridObject* elem,
							 const MathVector<dim> vCornerCoords[], bool bDeriv = false){

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

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__STD_USER_DATA__ */
