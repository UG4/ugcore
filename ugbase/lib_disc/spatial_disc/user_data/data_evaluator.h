/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
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

	///	prepares all time-dependent IElemDiscs
		void prepare_timestep(number future_time, const number time, VectorProxyBase* u, size_t algebra_id);

	///	prepares the elements of all time-dependent IElemDiscs
		void prepare_timestep_elem(const number time, LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

	///	prepares the element loop for all IElemDiscs
		void prepare_elem_loop(const ReferenceObjectID id, int si);

	///	prepares the element for all IElemDiscs
		void prepare_elem(LocalVector& u, GridObject* elem, const ReferenceObjectID roid, const MathVector<dim> vCornerCoords[],
		                  const LocalIndices& ind, bool bDeriv = false);

	///	finishes the element loop for all IElemDiscs
		void finish_elem_loop();

	/// finishes all time-dependent IElemDiscs
		void finish_timestep(const number time, VectorProxyBase* u, size_t algebra_id);

	///	finishes the elements of all time-dependent IElemDiscs
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
