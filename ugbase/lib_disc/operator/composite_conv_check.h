/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
 * Author: Markus Breit
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

#ifndef __H__LIB_DISC__OPERATOR__COMPOSITE_CONVERGENCE_CHECK__
#define __H__LIB_DISC__OPERATOR__COMPOSITE_CONVERGENCE_CHECK__

#include <ostream>
#include <sstream>
#include <string>
#include <limits>
#include <algorithm>
#include <iomanip>
#include <vector>
#include <list>
#include <cmath>

#include "common/common.h"
#include "common/stopwatch.h"
#include "lib_algebra/operator/convergence_check.h"
#include "lib_disc/dof_manager/dof_distribution.h"
#include "lib_disc/function_spaces/approximation_space.h"
#include "lib_disc/common/function_group.h"

namespace ug{

////////////////////////////////////////////////////////////////////////////////
// Composite convergence check
////////////////////////////////////////////////////////////////////////////////

/** CompositeConvCheck
 *
 * This is an implementation of the convergence check interface,
 * that makes it possible to define required defect reductions on
 * the individual functions constituting the overall solution.
 */
template <typename TVector, typename TDomain>
class CompositeConvCheck : public IConvergenceCheck<TVector>
{
	public:
	/// constructors
	/// \{
		explicit CompositeConvCheck(SmartPtr<ApproximationSpace<TDomain> > approx);
		CompositeConvCheck(SmartPtr<ApproximationSpace<TDomain> > spApproxSpace,
		                   int maxSteps, number minDefect, number relReduction);
	/// \}

	/// destructor
		~CompositeConvCheck() override = default;

	/// set level of grid, where defect vectors come from
		void set_level(int level);

	/// sets maximum number of iteration steps
		void set_maximum_steps(int maxSteps) {m_maxSteps = maxSteps;}

	///	sets default values for non-explicitly specified cmps
		void set_rest_check(number minDefect, number relReduction){
			m_bCheckRest = true;
			m_restMinDefect = minDefect; m_restRelReduction = relReduction;
			update_rest_check();
		}

	///	disables rest check
		void disable_rest_check(){m_bCheckRest = false; update_rest_check();}

	///	sets check for single component
	/// \{
		void set_component_check(const std::string& vFctName,
		                         const std::vector<number>& vMinDefect,
		                         const std::vector<number>& vRelReduction);

		void set_component_check(const std::vector<std::string>& vFctName,
		                         const std::vector<number>& vMinDefect,
		                         const std::vector<number>& vRelReduction);

		void set_component_check(const std::vector<std::string>& vFctName,
		                         number minDefect,
								 number relReduction);

		void set_component_check(const std::string& fctName,
		                         number minDefect,
		                         number relReduction);
	/// \}

	///	sets check for all components in approximation space
		void set_all_component_check(number minDefect,
		                             number relReduction);

	///	sets check for group of components
	/// \{
		void set_group_check(const std::vector<std::string>& vFctName,
		                     number minDefect,
							 number relReduction);

		void set_group_check(const std::string& fctNames,
		                     number minDefect,
							 number relReduction);
	/// \}

	/// defect control
		void start_defect(number initialDefect) override;
		void start(const TVector& d) override;
		void update_defect(number newDefect) override;
		void update(const TVector& d) override;
		bool iteration_ended() override;
		bool post() override;

	/// information about current status
		[[nodiscard]] int step() const override {return m_currentStep;}
		[[nodiscard]] number defect() const override {return defect_all();};
		[[nodiscard]] number reduction() const override {return defect_all()/initial_defect_all();};
		[[nodiscard]] number rate() const override {return defect_all()/last_defect_all();}
		[[nodiscard]] number avg_rate() const override {return std::pow(defect_all()/initial_defect_all(), 1.0/m_currentStep);}

	/// output
		[[nodiscard]] int get_offset() const override {return m_offset;};
		void set_offset(int offset) override {m_offset = offset;};
		void set_symbol(char symbol) override {m_symbol = symbol;};
		void set_name(std::string name) override {m_name = name;};
		void set_info(std::string info) override {m_info = info;};

	///	sets if verbose
		void set_verbose(bool level) {m_verbose = level;};

	/// set whether always to report success when max iter is reached (useful for LIMEX)
		void set_supress_unsuccessful(bool bsupress) {m_supress_unsuccessful = bsupress;}

	///	enables time measurement
		void set_time_measurement(bool yesOrNo) {m_bTimeMeas = yesOrNo;};

	/// whether the underlying approximation space is adaptive
		void set_adaptive(bool adapt) {m_bAdaptive = adapt;}

	///	clones this instance
		SmartPtr<IConvergenceCheck<TVector> > clone() override;

	///	prints a line using prefixes
		void print_line(std::string line) override;

		[[nodiscard]] std::string config_string() const override
		{
			std::stringstream ss;
			ss << "CompositeConvCheck( max steps = " << m_maxSteps << ")";
			ss << " Components:\n";
			for(size_t i=0; i<m_CmpInfo.size(); i++)
				ss << " | " << m_CmpInfo[i].config_string() << "\n";
			return ss.str();
		}

	protected:
		void print_offset() const;
		bool is_valid_number(number value) const;
		const std::string& fctName(size_t i) {return m_CmpInfo[i].name;};

	///	extracts multi-indices for a fct-comp on a element type
		template <typename TBaseElem>
		void extract_dof_indices(ConstSmartPtr<DoFDistribution> dd);

	///	extracts multi-indices from dof distribution
		void extract_dof_indices(ConstSmartPtr<DoFDistribution> dd);

	/// calculates the 2-norm of the entries of the vector vec specified by index
		number norm(const TVector& vec, const std::vector<DoFIndex>& index);

	protected:
	///	ApproxSpace
		SmartPtr<ApproximationSpace<TDomain> > m_spApprox;

		struct NativCmpInfo{
			std::string name; 	///< Name of components

			number initDefect;	///< Initial Defect of component
			number currDefect;	///< Current Defect of component
			number lastDefect;	///< Last Defect if component

			std::vector<DoFIndex> vMultiIndex; ///< associated indices
		};

	///	info on natural components
		std::vector<NativCmpInfo> m_vNativCmpInfo;
		size_t m_numAllDoFs;

	///	returns defect for all components
		[[nodiscard]] number defect_all() const {
			number defect = 0.0;
			for(size_t fct = 0; fct < m_vNativCmpInfo.size(); ++fct)
				defect += pow(m_vNativCmpInfo[fct].currDefect, 2);
			return sqrt(defect);
		}

	///	returns last defect for all components
		[[nodiscard]] number last_defect_all() const {
			number defect = 0.0;
			for(size_t fct = 0; fct < m_vNativCmpInfo.size(); ++fct)
				defect += pow(m_vNativCmpInfo[fct].lastDefect, 2);
			return sqrt(defect);
		}

	///	returns initial defect for all components
		[[nodiscard]] number initial_defect_all() const {
			number defect = 0.0;
			for(size_t fct = 0; fct < m_vNativCmpInfo.size(); ++fct)
				defect += pow(m_vNativCmpInfo[fct].initDefect, 2);
			return sqrt(defect);
		}

	protected:
		struct CmpInfo
		{
			CmpInfo() : isRest(false) {}

			std::vector<int> vFct;	///< Index of components
			std::string name; 		///< Name of components

			number initDefect;	///< Initial Defect of component
			number currDefect;	///< Current Defect of component
			number lastDefect;	///< Last Defect if component

			number minDefect;	///< Minimal required Defect of component
			number relReduction;///< Relative reduction required for component

			bool isRest; 	///< Shows, that this is group of remaining cmps

			[[nodiscard]] std::string config_string() const
			{
				std::stringstream ss;
				if(isRest) ss << "[Remaining Components]";
				else ss << "Component " << name;
				ss << ": minDefect = " << minDefect << ", relReduction = " << relReduction;
				return ss.str();
			}

		};

	///	infos for each component
		std::vector<CmpInfo> m_CmpInfo;

	///	default Values
		bool m_bCheckRest;
		number m_restMinDefect;
		number m_restRelReduction;
		void update_rest_check();

	/// current step
		int m_currentStep;

	/// maximum number of steps to be performed
		int m_maxSteps;

	protected:
		/// verbose level
		bool m_verbose;

		/// whether to always report success once max iter is reached
		bool m_supress_unsuccessful;

		/// number of spaces inserted before output
		int m_offset;

		/// symbol for output appearance
		char m_symbol;

		/// name of iteration
		std::string m_name;

		/// info for iteration (e.g. preconditioner type)
		std::string m_info;

	protected:
		/// enables time measurement
		bool m_bTimeMeas;

		/// a stopwatch
		Stopwatch m_stopwatch;

		// TODO: maybe solve this using an adaption/distribution event listener
		/// adaptivity flag
		bool m_bAdaptive;
};

} // end namespace ug

#include "composite_conv_check_impl.h"

#endif