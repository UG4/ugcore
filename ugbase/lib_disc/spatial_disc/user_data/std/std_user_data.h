/*
 * std_user_data.h
 *
 *  Created on: 03.07.2012
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__STD_USER_DATA__
#define __H__UG__LIB_DISC__SPATIAL_DISC__STD_USER_DATA__

#include "../user_data.h"

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
		                LocalVector& u,
		                GeometricObject* elem,
		                const MathVector\<dim\> vCornerCoords[],
		                const MathVector\<refDim\> vLocIP[],
		                const size_t nip,
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

	///	returns value for local and global position
	///	\{
		virtual TRet operator() (TData& value,
		                         const MathVector<dim>& globIP,
		                         number time, int si,
		                         LocalVector& u,
		                         GeometricObject* elem,
		                         const MathVector<dim> vCornerCoords[],
		                         const MathVector<1>& locIP) const
		{
			return getImpl().template evaluate<1>(value,globIP,time,si,u,
			                                      elem,vCornerCoords,locIP);
		}

		virtual TRet operator() (TData& value,
		                         const MathVector<dim>& globIP,
		                         number time, int si,
		                         LocalVector& u,
		                         GeometricObject* elem,
		                         const MathVector<dim> vCornerCoords[],
		                         const MathVector<2>& locIP) const
		{
			return getImpl().template evaluate<2>(value,globIP,time,si,u,
			                                      elem,vCornerCoords,locIP);
		}

		virtual TRet operator() (TData& value,
		                         const MathVector<dim>& globIP,
		                         number time, int si,
		                         LocalVector& u,
		                         GeometricObject* elem,
		                         const MathVector<dim> vCornerCoords[],
		                         const MathVector<3>& locIP) const
		{
			return getImpl().template evaluate<3>(value,globIP,time,si,u,
			                                      elem,vCornerCoords,locIP);
		}
	///	\}

	///	returns value for global positions
		virtual void operator()(TData vValue[],
								const MathVector<dim> vGlobIP[],
								number time, int si, const size_t nip) const = 0;

	///	returns values for local and global positions
	///	\{
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
		StdDependentUserData() : m_pFctPatt(NULL){}

		StdDependentUserData(const char* functions) : m_pFctPatt(NULL){
			set_functions(functions);
		}

	public:
		virtual void compute(LocalVector* u, GeometricObject* elem,
							 const MathVector<dim> vCornerCoords[], bool bDeriv = false) = 0;

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

	public:
	///	returns if grid function is needed for evaluation
		virtual bool requires_grid_fct() const {return true;}

	///	sets the associated function pattern
		virtual void set_function_pattern(const FunctionPattern& fctPatt)
		{
			m_pFctPatt = &fctPatt;
			extract_fct_grp();
		}

	///	sets the associated symbolic functions
		void set_functions(const char* symbFct)
		{
			m_SymbFct = symbFct;
			extract_fct_grp();
		}

	protected:
	///	extracts the function group
		void extract_fct_grp()
		{
		//	if associated infos missing return
			if(m_pFctPatt == NULL) return;
			this->m_fctGrp.set_function_pattern(*m_pFctPatt);

			if(m_SymbFct.empty()){
				this->m_fctGrp.clear();
				return;
			}

		//	create function group of this elem disc
			try{
				this->m_fctGrp.clear();
				this->m_fctGrp.add(TokenizeString(m_SymbFct));
			}UG_CATCH_THROW("StdDependendDataExport: Cannot find  some symbolic function "
							"name in '"<<m_SymbFct<<"'.");

		//	create a mapping between all functions and the function group of this
		//	element disc.
			try{
				CreateFunctionIndexMapping(this->m_map, this->m_fctGrp, *m_pFctPatt);
			}UG_CATCH_THROW("StdDependendDataExport: Cannot create Function Index Mapping"
							" for '"<<m_SymbFct<<"'.");
		}

	protected:
	///	associated function pattern
		const FunctionPattern* m_pFctPatt;

	///	string of symbolic functions required
		std::string m_SymbFct;

	protected:
	///	access to implementation
		TImpl& getImpl() {return static_cast<TImpl&>(*this);}

	///	const access to implementation
		const TImpl& getImpl() const {return static_cast<const TImpl&>(*this);}

};

} // namespace ug

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__STD_USER_DATA__ */
