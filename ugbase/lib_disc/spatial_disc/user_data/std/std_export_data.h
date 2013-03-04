/*
 * std_export_data.h
 *
 *  Created on: 03.07.2012
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__STD_EXPORT_DATA__
#define __H__UG__LIB_DISC__SPATIAL_DISC__STD_EXPORT_DATA__

#include "std_user_data.h"
#include "../data_export.h"

#include "lib_disc/common/function_group.h"
#include "lib_disc/common/local_algebra.h"
#include "lib_disc/common/groups_util.h"
#include "common/util/string_util.h"

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
			getImpl().evaluate<refDim>(&value,&globIP,time,si,u,elem,
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

} // namespace ug

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__STD_EXPORT_DATA__ */
