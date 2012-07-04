/*
 * std_ip_data.h
 *
 *  Created on: 03.07.2012
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__STD_IPDATA__
#define __H__UG__LIB_DISC__SPATIAL_DISC__STD_IPDATA__

#include "lib_disc/common/function_group.h"
#include "lib_disc/common/local_algebra.h"
#include "lib_disc/common/groups_util.h"

#include "ip_data.h"
#include "data_import_export.h"

namespace ug{

///////////////////////////////////////////////////////////////////////////////
// Base class for Position-Time-IPData
///////////////////////////////////////////////////////////////////////////////

/**
 * This class is a base class for all position and time dependent user data.
 * The data thus does not on the a computed solution.
 * In order to use the interface, the deriving class must implement the method:
 *
 * inline TRet evaluate(TData& D, const MathVector<dim>& x, number time, int si) const
 *
 */
template <typename TImpl, typename TData, int dim, typename TRet = void>
class StdPositionIPData
	: 	public IPData<TData,dim,TRet>
{
	public:
		StdPositionIPData() {}

		virtual TRet operator() (TData& value,
								 const MathVector<dim>& globIP,
								 number time, int si) const
		{
			return getImpl().evaluate(value, globIP, time, si);
		}

		virtual void operator() (TData vValue[],
								 const MathVector<dim> vGlobIP[],
								 number time, int si, const size_t nip) const
		{
			for(size_t ip = 0; ip < nip; ++ip)
				getImpl().evaluate(vValue[ip], vGlobIP[ip], time, si);
		}

		////////////////
		// one value
		////////////////

		virtual TRet operator() (TData& value,
								 const MathVector<dim>& globIP,
								 number time, int si,
								 LocalVector& u,
								 GeometricObject* elem,
								 const MathVector<dim> vCornerCoords[],
								 const MathVector<1>& locIP) const
		{
			return getImpl().evaluate(value, globIP, time, si);
		}

		virtual TRet operator() (TData& value,
								 const MathVector<dim>& globIP,
								 number time, int si,
								 LocalVector& u,
								 GeometricObject* elem,
								 const MathVector<dim> vCornerCoords[],
								 const MathVector<2>& locIP) const
		{
			return getImpl().evaluate(value, globIP, time, si);
		}

		virtual TRet operator() (TData& value,
								 const MathVector<dim>& globIP,
								 number time, int si,
								 LocalVector& u,
								 GeometricObject* elem,
								 const MathVector<dim> vCornerCoords[],
								 const MathVector<3>& locIP) const
		{
			return getImpl().evaluate(value, globIP, time, si);
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
								const size_t nip) const
		{
			for(size_t ip = 0; ip < nip; ++ip)
				getImpl().evaluate(vValue[ip], vGlobIP[ip], time, si);
		}

		virtual void operator()(TData vValue[],
								const MathVector<dim> vGlobIP[],
								number time, int si,
								LocalVector& u,
								GeometricObject* elem,
								const MathVector<dim> vCornerCoords[],
								const MathVector<2> vLocIP[],
								const size_t nip) const
		{
			for(size_t ip = 0; ip < nip; ++ip)
				getImpl().evaluate(vValue[ip], vGlobIP[ip], time, si);
		}

		virtual void operator()(TData vValue[],
								const MathVector<dim> vGlobIP[],
								number time, int si,
								LocalVector& u,
								GeometricObject* elem,
								const MathVector<dim> vCornerCoords[],
								const MathVector<3> vLocIP[],
								const size_t nip) const
		{
			for(size_t ip = 0; ip < nip; ++ip)
				getImpl().evaluate(vValue[ip], vGlobIP[ip], time, si);
		}

	///	implement as a IPData
		virtual void compute(bool bDeriv = false)
		{
			const number t = this->time();
			const int si = this->subset();

			for(size_t s = 0; s < this->num_series(); ++s)
				for(size_t ip = 0; ip < this->num_ip(s); ++ip)
					getImpl().evaluate(this->value(s,ip), this->ip(s, ip), t, si);
		}

	///	callback, invoked when data storage changed
		virtual void value_storage_changed(const size_t seriesID)
		{
			const number t = this->time();
			const int si = this->subset();

			for(size_t ip = 0; ip < this->num_ip(seriesID); ++ip)
				getImpl().evaluate(this->value(seriesID,ip), this->ip(seriesID, ip), t, si);
		}

	///	returns if data is constant
		virtual bool constant_data() const {return false;}

	///	returns if grid function is needed for evaluation
		virtual bool requires_grid_fct() const {return false;}

	///	returns if provided data is continuous over geometric object boundaries
		virtual bool is_continuous() const {return true;}

	protected:
	///	access to implementation
		TImpl& getImpl() {return static_cast<TImpl&>(*this);}

	///	const access to implementation
		const TImpl& getImpl() const {return static_cast<const TImpl&>(*this);}
};



///////////////////////////////////////////////////////////////////////////////
// Base class for Export IPData
///////////////////////////////////////////////////////////////////////////////

/**
 * This class is a base class for export user data.
 * The data thus does depend on the a computed solution and space, time and subset.
 * In order to use the interface, the deriving class must implement the method:
 *
 * template <int refDim>
 * inline TRet evaluate(TData vValue[],
		                const MathVector<dim> vGlobIP[],
		                number time, int si,
		                LocalVector& u,
		                GeometricObject* elem,
		                const MathVector<dim> vCornerCoords[],
		                const MathVector<refDim> vLocIP[],
		                const size_t nip,
		                const MathMatrix<refDim, dim>* vJT = NULL) const
 *
 */
template <typename TImpl, typename TData, int dim>
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
			getImpl().template evaluate<1>(&value,&globIP,time,si,u,elem,
			                               vCornerCoords,&locIP,1,NULL);
		}

		virtual void operator() (TData& value,
								 const MathVector<dim>& globIP,
								 number time, int si,
								 LocalVector& u,
								 GeometricObject* elem,
								 const MathVector<dim> vCornerCoords[],
								 const MathVector<2>& locIP) const
		{
			getImpl().template evaluate<2>(&value,&globIP,time,si,u,elem,
			                               vCornerCoords,&locIP,1,NULL);
		}

		virtual void operator() (TData& value,
								 const MathVector<dim>& globIP,
								 number time, int si,
								 LocalVector& u,
								 GeometricObject* elem,
								 const MathVector<dim> vCornerCoords[],
								 const MathVector<3>& locIP) const
		{
			getImpl().template evaluate<3>(&value,&globIP,time,si,u,elem,
			                               vCornerCoords,&locIP,1,NULL);
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

} // namespace ug

#include "data_linker.h"

namespace ug{
///////////////////////////////////////////////////////////////////////////////
// Base class for Linker IPData
///////////////////////////////////////////////////////////////////////////////

/**
 * This class is a base class for linker user data.
 * In order to use the interface, the deriving class must implement the method:
 *
 * template <int refDim>
 * inline TRet evaluate(TData vValue[],
						const MathVector<dim> vGlobIP[],
						number time, int si,
						LocalVector& u,
						GeometricObject* elem,
						const MathVector<dim> vCornerCoords[],
						const MathVector<refDim> vLocIP[],
						const size_t nip,
						const MathMatrix<refDim, dim>* vJT = NULL) const
 *
 */
template <typename TImpl, typename TData, int dim>
class StdDataLinker
	: 	public DataLinker<TData,dim>
{
	public:
		StdDataLinker() {}

		virtual void operator() (TData& value,
								 const MathVector<dim>& globIP,
								 number time, int si) const
		{
			getImpl().evaluate(value,globIP,time,si);
		}

		virtual void operator() (TData vValue[],
								 const MathVector<dim> vGlobIP[],
								 number time, int si, const size_t nip) const
		{
			for(size_t ip = 0; ip < nip; ++ip)
				getImpl().evaluate(vValue[ip],vGlobIP[ip],time,si);
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
		virtual bool requires_grid_fct() const
		{
			for(size_t i = 0; i < this->m_vpIIPData.size(); ++i)
				if(this->m_vpIIPData[i]->requires_grid_fct())
					return true;
			return false;
		}

	///	returns if provided data is continuous over geometric object boundaries
		virtual bool is_continuous() const
		{
			bool bRet = true;
			for(size_t i = 0; i < this->m_vpIIPData.size(); ++i)
				bRet &= this->m_vpIIPData[i]->is_continuous();
			return bRet;
		}

	///	sets the associated function pattern
		virtual void set_function_pattern(const FunctionPattern& fctPatt)
		{
			for(size_t i = 0; i < this->m_vpIIPData.size(); ++i)
				this->m_vpIIPData[i]->set_function_pattern(fctPatt);
		}

	protected:
	///	access to implementation
		TImpl& getImpl() {return static_cast<TImpl&>(*this);}

	///	const access to implementation
		const TImpl& getImpl() const {return static_cast<const TImpl&>(*this);}
};


} // end namespace ug

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__STD_IPDATA__ */
