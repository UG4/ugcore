/*
 * data_export.h
 *
 *  Created on: 04.07.2012
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__DATA_EXPORT__
#define __H__UG__LIB_DISC__SPATIAL_DISC__DATA_EXPORT__

#include "lib_disc/common/function_group.h"
#include "lib_disc/common/local_algebra.h"
#include "lib_disc/common/groups_util.h"
#include "common/util/string_util.h"

#include "user_data.h"

namespace ug{

// predeclaration
class IElemDisc;

////////////////////////////////////////////////////////////////////////////////
// Data Export
////////////////////////////////////////////////////////////////////////////////

/// Data export
/**
 * A DataExport is user data produced by an element discretization.
 */
template <typename TData, int dim>
class DataExport : 	public DependentUserData<TData, dim>
{
	public:
	///	default constructor
		DataExport();

	///	implement compute() method of IUserData
		virtual void compute(LocalVector* u, GeometricObject* elem, bool bDeriv = false);

	///	sets the geometric object type
		virtual void set_roid(ReferenceObjectID id);

	///	register evaluation of export function
		template <typename TClass, int refDim>
		void set_fct(ReferenceObjectID id, TClass* obj,
		             void (TClass::*func)(
		            		 const LocalVector& u,
		            		 const MathVector<dim> vGlobIP[],
		            		 const MathVector<refDim> vLocIP[],
		            		 const size_t nip,
		            		 TData vValue[],
		            		 bool bDeriv,
		            		 std::vector<std::vector<TData> > vvvDeriv[]));

	///	register evaluation of export function
		template <int refDim>
		void set_fct(ReferenceObjectID id,
		             void (*func)(
		            		 const LocalVector& u,
		            		 const MathVector<dim> vGlobIP[],
		            		 const MathVector<refDim> vLocIP[],
		            		 const size_t nip,
		            		 TData vValue[],
		            		 bool bDeriv,
		            		 std::vector<std::vector<TData> > vvvDeriv[]));

	///	clears all export functions
		void clear_fct();

	///	returns if data depends on solution
		virtual bool zero_derivative() const {return false;}

	///	clear dependent data
		void clear() {m_vDependData.clear();}

	///	add data dependency
		void add_needed_data(SmartPtr<IUserData> data);

	///	remove needed data
		void remove_needed_data(SmartPtr<IUserData> data);

	///	number of other Data this data depends on
		virtual size_t num_needed_data() const {return m_vDependData.size();}

	///	return needed data
		virtual SmartPtr<IUserData> needed_data(size_t i) {return m_vDependData.at(i);}

	///	returns if the dependent data is ready for evaluation
		virtual void check_setup() const;

	protected:
		template <int refDim>
		struct traits{
			typedef boost::function<void (const LocalVector& u,
			                              const MathVector<dim> vGlobIP[],
			                              const MathVector<refDim> vLocIP[],
			                              const size_t nip,
			                              TData vValue[],
			                              bool bDeriv,
			                              std::vector<std::vector<TData> > vvvDeriv[])>
			EvalFunc;
		};

		template <int refDim>
		typename traits<refDim>::EvalFunc& eval_fct(ReferenceObjectID id) {return eval_fct(id, Int2Type<refDim>());}
		template <int refDim>
		const typename traits<refDim>::EvalFunc& eval_fct(ReferenceObjectID id) const {return eval_fct(id, Int2Type<refDim>());}

		typename traits<1>::EvalFunc& eval_fct(ReferenceObjectID id, Int2Type<1>) {return m_vCompFct1d[id];}
		typename traits<2>::EvalFunc& eval_fct(ReferenceObjectID id, Int2Type<2>) {return m_vCompFct2d[id];}
		typename traits<3>::EvalFunc& eval_fct(ReferenceObjectID id, Int2Type<3>) {return m_vCompFct3d[id];}
		const typename traits<1>::EvalFunc& eval_fct(ReferenceObjectID id, Int2Type<1>) const {return m_vCompFct1d[id];}
		const typename traits<2>::EvalFunc& eval_fct(ReferenceObjectID id, Int2Type<2>) const {return m_vCompFct2d[id];}
		const typename traits<3>::EvalFunc& eval_fct(ReferenceObjectID id, Int2Type<3>) const {return m_vCompFct3d[id];}
		typename traits<1>::EvalFunc m_vCompFct1d[NUM_REFERENCE_OBJECTS];
		typename traits<2>::EvalFunc m_vCompFct2d[NUM_REFERENCE_OBJECTS];
		typename traits<3>::EvalFunc m_vCompFct3d[NUM_REFERENCE_OBJECTS];

		bool eval_fct_set(ReferenceObjectID id) const;

		template <int refDim>
		void comp(const LocalVector& u, bool bDeriv);

	protected:
	/// current Geom Object
		ReferenceObjectID m_id;

	///	data the export depends on
		std::vector<SmartPtr<IUserData> > m_vDependData;
};


///////////////////////////////////////////////////////////////////////////////
// Base class for Export UserData
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

		////////////////
		// one value
		////////////////
		virtual void operator() (TData& value,
								 const MathVector<dim>& globIP,
								 number time, int si) const
		{
			UG_THROW("StdDataExport: Solution, element and local ips required "
					"for evaluation, but not passed. Cannot evaluate.");
		}

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
		virtual void operator() (TData value[],
								 const MathVector<dim> globIP[],
								 number time, int si, const size_t nip) const
		{
			UG_THROW("StdDataExport: Solution, element and local ips required "
					"for evaluation, but not passed. Cannot evaluate.");
		}

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

} // end namespace ug

#include "data_export_impl.h"

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__DATA_EXPORT__ */
