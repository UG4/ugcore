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
		              ConstSmartPtr<FunctionPattern> fctPat,
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
		void prepare_timestep_elem(const number time, LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

	///	prepares the element loop for all IElemDiscs
		void prepare_elem_loop(const ReferenceObjectID id, int si);

	///	prepares the element for all IElemDiscs
		void prepare_elem(LocalVector& u, GridObject* elem, const ReferenceObjectID roid, const MathVector<dim> vCornerCoords[],
		                  const LocalIndices& ind, bool bDeriv = false);

	///	finishes the element loop for all IElemDiscs
		void finish_elem_loop();

	///	finishes the element for all time-dependent IElemDiscs
		void finish_timestep_elem(const number time, LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

	///	compute local stiffness matrix for all IElemDiscs
		void add_jac_A_elem(LocalMatrix& A, LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[], ProcessType type = PT_ALL);

	///	compute local mass matrix for all IElemDiscs
		void add_jac_M_elem(LocalMatrix& M, LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[], ProcessType type = PT_ALL);

	///	compute local stiffness defect for all IElemDiscs
		void add_def_A_elem(LocalVector& d, LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[], ProcessType type = PT_ALL);

	///	compute local stiffness defect for all IElemDiscs explicit
		void add_def_A_expl_elem(LocalVector& d, LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[], ProcessType type = PT_ALL);

	///	compute local mass defect for all IElemDiscs
		void add_def_M_elem(LocalVector& d, LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[], ProcessType type = PT_ALL);

	///	compute local rhs for all IElemDiscs
		void add_rhs_elem(LocalVector& rhs, GridObject* elem, const MathVector<dim> vCornerCoords[], ProcessType type = PT_ALL);
		
	///	prepares the element loop for all IElemDiscs for the computation of the error estimator
		void prepare_err_est_elem_loop(const ReferenceObjectID id, int si);

	///	prepares the element for all IElemDiscs
		void prepare_err_est_elem(LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[],
						  	  	  const LocalIndices& ind, bool bDeriv = false);

	///	compute contributions of the local error indicators in one element for all IElemDiscs
		void compute_err_est_A_elem(LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[], const LocalIndices& ind,
								  const number scaleMass = (number) 1.0, const number scaleStiff = (number) 1.0);

	///	compute contributions of the local error indicators in one element for all IElemDiscs
		void compute_err_est_M_elem(LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[], const LocalIndices& ind,
								  const number scaleMass = (number) 1.0, const number scaleStiff = (number) 1.0);

	///	compute contributions of the local error indicators in one element for all IElemDiscs
		void compute_err_est_rhs_elem(GridObject* elem, const MathVector<dim> vCornerCoords[], const LocalIndices& ind,
								  const number scaleMass = (number) 1.0, const number scaleStiff = (number) 1.0);

	///	finishes the error estimator element loop for all IElemDiscs
		void finish_err_est_elem_loop();

	protected:
	///	clears imports and user data and mappings betweem commonFctGrp and local
		void clear_extracted_data_and_mappings();

	///	tries to add the last entry of vTryingToAdd to the eval data
		void add_data_to_eval_data(std::vector<SmartPtr<ICplUserData<dim> > >& vEvalData,
								   std::vector<SmartPtr<ICplUserData<dim> > >& vTryingToAdd);

	///	extracts imports and userdata from IElemDiscs
		void extract_imports_and_userdata(int subsetIndex, int discPart);

	///	clears all requested positions in user data
		void clear_positions_in_user_data();

	protected:
	///	disc part needed (MASS and/or STIFF and/or RHS)
		int m_discPart;

	///	elem disc data
		std::vector<IElemDisc<TDomain>*> m_vElemDisc[MAX_PROCESS];

	///	underlying function pattern
		ConstSmartPtr<FunctionPattern> m_spFctPattern;

	///	flag if hanging nodes are used
		bool m_bUseHanging;

	///	flag indicating if any elem disc needs local time series
		bool m_bNeedLocTimeSeries;

	///	local time series (non-const since mapping may change)
		LocalVectorTimeSeries* m_pLocTimeSeries;

	///	imports that must be evaluated
		std::vector<IDataImport<dim>*> m_vImport[MAX_PROCESS][MAX_PART];

	///	user data that must be evaluated
	/// \{
		std::vector<SmartPtr<ICplUserData<dim> > > m_vConstData;	 //< constant data
		std::vector<SmartPtr<ICplUserData<dim> > > m_vPosData;       //< position dependent data
		std::vector<SmartPtr<ICplUserData<dim> > > m_vDependentData; //< dependent data
	/// \}
};

} // end namespace ug

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__DATA_EVALUATOR__ */
