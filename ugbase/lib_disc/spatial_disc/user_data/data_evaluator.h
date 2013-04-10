/*
 * data_evaluator.h
 *
 *  Created on: 17.12.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__DATA_EVALUATOR__
#define __H__UG__LIB_DISC__SPATIAL_DISC__DATA_EVALUATOR__

#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h"
#include "lib_disc/spatial_disc/user_data/user_data.h"
#include "lib_disc/spatial_disc/user_data/data_import.h"

namespace ug{

enum ProcessType {PT_ALL=0, PT_STATIONARY, PT_INSTATIONARY, MAX_PROCESS};

/// helper class to evaluate data evaluation for local contributions during assembling
/**
 * This class is a helper class to facilitate the correct evaluation of data
 * during the assembling process. For a given set of IElemDisc this class
 * prepares the assembling, checks for DataImport/DataExports and schedules the
 * evaluation of the data in the correct order. In addition the coupling
 * contributions due to Import/Exports are computed by this class and can be
 * added to the local jacobian/defect.
 */
template <typename TDomain>
class DataEvaluator
{
	public:
	///	world dimension
		static const int dim = TDomain::dim;

	public:
	///	sets the elem discs to evaluate
		DataEvaluator(int discPart,
		              const std::vector<IElemDisc<TDomain>*>& vElemDisc,
		              const FunctionPattern& fctPat,
		              const int subset,
		              const bool bNonRegularGrid,
		              LocalVectorTimeSeries* locTimeSeries = NULL,
		              const std::vector<number>* vScaleMass = NULL,
		              const std::vector<number>* vScaleStiff = NULL);

	///	sets the time point for data evaluation
		void set_time_point(const size_t timePoint);

	///	returns if local time series is really needed for assembling
		bool time_series_needed() const {return m_bNeedLocTimeSeries;}

	///	returns if one of the element discs needs hanging dofs
		bool use_hanging() const {return m_bUseHanging;}

		////////////////////////////////////////////
		// Regular assembling
		////////////////////////////////////////////

	///	prepares the element for all time-dependent IElemDiscs
		template <typename TElem>
		void prepare_timestep_elem(TElem* elem, LocalVector& u);

	///	prepares the element loop for all IElemDiscs
		template <typename TElem>
		void prepare_elem_loop(int si);

	///	prepares the element for all IElemDiscs
		template <typename TElem>
		void prepare_elem(TElem* elem, LocalVector& u,
		                  const LocalIndices& ind,
		                  bool bDeriv = false);

	///	finishes the element for all time-dependent IElemDiscs
		template <typename TElem>
		void finish_timestep_elem(TElem* elem, const number time, LocalVector& u);

	///	compute local stiffness matrix for all IElemDiscs
		void add_JA_elem(LocalMatrix& A, LocalVector& u, GeometricObject* elem, ProcessType type = PT_ALL);

	///	compute local mass matrix for all IElemDiscs
		void add_JM_elem(LocalMatrix& M, LocalVector& u, GeometricObject* elem, ProcessType type = PT_ALL);

	///	compute local stiffness defect for all IElemDiscs
		void add_dA_elem(LocalVector& d, LocalVector& u, GeometricObject* elem, ProcessType type = PT_ALL);

	///	compute local stiffness defect for all IElemDiscs explicit
		void add_dA_elem_explicit(LocalVector& d, LocalVector& u, GeometricObject* elem, ProcessType type = PT_ALL);

	///	compute local mass defect for all IElemDiscs
		void add_dM_elem(LocalVector& d, LocalVector& u, GeometricObject* elem, ProcessType type = PT_ALL);

	///	compute local rhs for all IElemDiscs
		void add_rhs_elem(LocalVector& rhs, GeometricObject* elem, ProcessType type = PT_ALL);

	///	finishes the element loop for all IElemDiscs
		void finish_elem_loop();

	protected:
	///	computes all needed data on the element
		void compute_elem_data(LocalVector& u, GeometricObject* elem,
		                       const MathVector<dim> vCornerCoords[], bool bDeriv = false);

	///	adds the contribution due to coupling to local stiffness matrix
		void add_coupl_JA(LocalMatrix& J, LocalVector& u, ProcessType type = PT_ALL);

	///	adds the contribution due to coupling to local mass matrix
		void add_coupl_JM(LocalMatrix& J, LocalVector& u, ProcessType type = PT_ALL);

	///	clears imports and user data and mappings betweem commonFctGrp and local
		void clear_extracted_data_and_mappings();

	///	tries to add the last entry of vTryingToAdd to the eval data
		void add_data_to_eval_data(std::vector<SmartPtr<ICplUserData<dim> > >& vEvalData,
								   std::vector<SmartPtr<ICplUserData<dim> > >& vTryingToAdd);

	///	extracts imports and userdata from IElemDiscs
		void extract_imports_and_userdata(int discPart);

	///	clears all requested positions in user data
		void clear_positions_in_user_data();

	protected:
	///	disc part needed (MASS and/or STIFF and/or RHS)
		int m_discPart;

	///	struct to store data related to elem disc
		struct ElemDisc {
			IElemDisc<TDomain>* elemDisc;	//< element disc
			FunctionGroup fctGrp;			//< associated function group
			FunctionIndexMapping map;		//< mapping from fctPatt for the fct group
			bool needLocTimeSeries;			//< flag if time series needed
			ProcessType process;			//< type of process (STAT, INSTAT, ..)
		};

	///	elem disc data
		std::vector<ElemDisc> m_vElemDisc[MAX_PROCESS];

	///	underlying function pattern
		const FunctionPattern& m_fctPatt;

	///	flag if hanging nodes are used
		bool m_bUseHanging;

	///	flag indicating if any elem disc needs local time series
		bool m_bNeedLocTimeSeries;

	///	subset
		int m_subset;

	///	local time series (non-const since mapping may change)
		LocalVectorTimeSeries* m_pLocTimeSeries;

	///	current element
		GeometricObject* m_pElem;

	///	current Corner Coordinates
		const MathVector<dim>* m_vCornerCoords;

	////////////////////////////////
	// 	Data Import
	////////////////////////////////
	///	struct to store data related to an import
		struct Import{
			Import(IDataImport<dim>* _import,
			       FunctionIndexMapping& _map, FunctionIndexMapping& _connMap,
			       ProcessType _process) :
			import(_import), map(_map), connMap(_connMap), process(_process) {}

			IDataImport<dim>* import;		//< import
			FunctionIndexMapping map;		//< mapping for import fct group
			FunctionIndexMapping connMap;	//< mapping for data fct group
			ProcessType process;			//< type of process (STAT, INSTAT,...)
		};

		std::vector<Import> m_vImport[MAX_PROCESS][MAX_PART];

	////////////////////////////////
	// 	UserData
	////////////////////////////////
		std::vector<SmartPtr<ICplUserData<dim> > > m_vConstData;	    //< constant data
		std::vector<SmartPtr<ICplUserData<dim> > > m_vPosData;       //< position dependent data
		std::vector<SmartPtr<ICplUserData<dim> > > m_vDependentData; //< dependent data
		std::vector<FunctionIndexMapping> m_vDependentMap;
};

} // end namespace ug

#include "data_evaluator_impl.h"

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__DATA_EVALUATOR__ */
