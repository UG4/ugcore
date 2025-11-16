/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__LIB_DISC__FUNCTION_SPACE__APPROXIMATION_SPACE__
#define __H__UG__LIB_DISC__FUNCTION_SPACE__APPROXIMATION_SPACE__

#include "lib_disc/common/revision_counter.h"
#include "lib_disc/dof_manager/dof_distribution_info.h"
#include "lib_disc/dof_manager/dof_distribution.h"
#include "lib_grid/tools/surface_view.h"
#include "lib_algebra/algebra_type.h"

namespace ug{

/// describes the ansatz spaces on a domain
/**
 * This class provides grid function spaces on a domain.
 *
 * The Domain defines a partition of the Grid/Multigrid in terms of subsets.
 * The user can add discrete functions on this subsets or unions of them.
 *
 * Once finalized, this function pattern is fixed. Internally DoF indices are
 * created. Using this Approximation Space the user can create GridFunctions
 * of the following types:
 *
 * - surface grid function = grid function representing the space on the surface grid
 * - level grid function = grid function representing the space on a level grid
 *  (NOTE: 	For a fully refined Multigrid a level grid covers the whole domain.
 *  		However for a locally/adaptively refined MultiGrid the level grid
 *  		solution is only living on a part of the domain)
 */
class IApproximationSpace : public DoFDistributionInfoProvider
{
	public:
	///	Type of Subset Handler
	using subset_handler_type = MGSubsetHandler;

	///	Type of Subset Handler
	using grid_type = MultiGrid;

	public:
	///	Constructor
		IApproximationSpace(SmartPtr<subset_handler_type> spMGSH,
		                    SmartPtr<grid_type> spMG);

	///	Constructor setting the grouping flag
		IApproximationSpace(SmartPtr<subset_handler_type> spMGSH,
		                    SmartPtr<grid_type>,
		                    const AlgebraType& algebraType);

	protected:
	///	initializing
		void init(SmartPtr<subset_handler_type> spMGSH,
		          SmartPtr<grid_type> spMG, const AlgebraType& algebraType);

	public:
	///	Destructor
		~IApproximationSpace();

	///	clears functions
		void clear() {m_spDoFDistributionInfo->clear();}

	public:
	/// add single solutions of LocalShapeFunctionSetID to the entire domain
	/**
	 * \param[in] 	name		name(s) of single solution (comma separated)
	 * \param[in] 	id			Shape Function set id
	 */
		void add(const std::vector<std::string>& vName, LFEID id){
			m_spDoFDistributionInfo->add(vName, id);
		}

	///	adds function using string to indicate finite element type
		void add(const std::vector<std::string>& vName, const char* type, int order);

	///	adds function using string to indicate finite element type
		void add(const std::vector<std::string>& vName, const char* type);

	///	adds function using string to indicate finite element type
		void add(const char* name, const char* type, int order);

	///	adds function using string to indicate finite element type
		void add(const char* name, const char* type);

	public:
	/// add single solutions of LocalShapeFunctionSetID to selected subsets
	/**
	 * \param[in] name			name(s) of single solution (comma separated)
	 * \param[in] id			Shape Function set id
	 * \param[in] subsets		Subsets separated by ','
	 */
		void add(const std::vector<std::string>& vName, LFEID id,
				 const std::vector<std::string>& vSubset){
			m_spDoFDistributionInfo->add(vName, id, vSubset);
		}

	///	adds function using string to indicate finite element type
		void add(const std::vector<std::string>& vName, const char* type, int order,
				 const std::vector<std::string>& vSubsets);

	///	adds function using string to indicate finite element type
	/**
	 * \param[in] name			name(s) of single solution (comma separated)
	 * \param[in] type			type of local finite element space
	 * \param[in] order			order of local finite element space
	 * \param[in] subsets		Subsets separated by ','
	 */
		void add(const char* name, const char* type, int order, const char* subsets);

	///	adds function using string to indicate finite element type
		void add(const std::vector<std::string>& vName, const char* type,
				 const std::vector<std::string>& vSubsets);

	///	adds function using string to indicate finite element type
	/**
	 * \param[in] name			name(s) of single solution (comma separated)
	 * \param[in] type			type of local finite element space
	 * \param[in] subsets		Subsets separated by ','
	 */
		void add(const char* name, const char* type, const char* subsets);

	public:
	/// get underlying subset handler
		ConstSmartPtr<MGSubsetHandler> subset_handler() const {return m_spMGSH;}

	///	returns the number of level
		size_t num_levels() const {return m_spMGSH->num_levels();}

	///	returns the approximation space
		ConstSmartPtr<SurfaceView> surface_view() const {return m_spSurfaceView;}

	///	returns if dofs are grouped
		bool grouped() const {return m_bGrouped;}

	///	returns if ghosts might be present on a level
		bool might_contain_ghosts(int lvl) const;

	///	returns if ghosts might be present on any level
		bool might_contain_ghosts() const;

	///	returns dof distribution for a grid level
	/// \{
		SmartPtr<DoFDistribution> dof_distribution(const GridLevel& gl, bool bCreate = true);
		SmartPtr<DoFDistribution> dd(const GridLevel& gl, bool bCreate = true);
	/// \}

	///	returns dof distribution for a grid level
	/// \{
		ConstSmartPtr<DoFDistribution> dof_distribution(const GridLevel& gl, bool bCreate = true) const;
		ConstSmartPtr<DoFDistribution> dd(const GridLevel& gl, bool bCreate = true) const;
	/// \}

	///	returns all currently created dof distributions
		std::vector<SmartPtr<DoFDistribution> > dof_distributions() const;

	///	returns dof distribution info
	/// \{
		ConstSmartPtr<DoFDistributionInfo> ddinfo() const {return m_spDoFDistributionInfo;}
		ConstSmartPtr<DoFDistributionInfo> dof_distribution_info() const {return ddinfo();}
	///	\}

	///	prints statistic about DoF Distribution
		void print_statistic(std::string flags) const;

	///	prints statistic about DoF Distribution
		void print_statistic() const;

	///	prints statistic on layouts
		void print_layout_statistic() const;


	///	initializes all level dof distributions
		void init_levels();

	///	initializes all surface dof distributions
		void init_surfaces();

	///	initializes all top surface dof distributions
		void init_top_surface();

	///	returns the current revision
		const RevisionCounter& revision() const {return m_RevCnt;}

	protected:
	///	creates a dof distribution
		void create_dof_distribution(const GridLevel& gl);

	///	creates surface SurfaceView if needed
		void surface_view_required();

	///	create dof distribution info
		void dof_distribution_info_required();

	protected:
	///	reinits all data after grid adaption
		void reinit();

	///	message hub id
		MessageHub::SPCallbackId m_spGridAdaptionCallbackID;
		MessageHub::SPCallbackId m_spGridDistributionCallbackID;
		bool m_bAdaptionIsActive;

	///	registers at message hub for grid adaption
		void register_at_adaption_msg_hub();

	/**	this callback is called by the message hub, when a grid change has been
	 * performed. It will call all necessary actions in order to keep the grid
	 * correct for computations.*/
		void grid_changed_callback(const GridMessage_Adaption& msg);

	///	called during parallel redistribution
		void grid_distribution_callback(const GridMessage_Distribution& msg);

	protected:
	/// multigrid, where elements are stored
		SmartPtr<MultiGrid> m_spMG;

	/// subsethandler, where elements are stored
		SmartPtr<MGSubsetHandler> m_spMGSH;

	///	Surface View
		SmartPtr<SurfaceView> m_spSurfaceView;

	///	suitable algebra type for the index distribution pattern
		AlgebraType m_algebraType;

	///	flag if DoFs should be grouped
		bool m_bGrouped;

	///	DofDistributionInfo
		SmartPtr<DoFDistributionInfo> m_spDoFDistributionInfo;

	///	revision counter
		RevisionCounter m_RevCnt;

	protected:
	///	MG Level DoF Distribution
		std::vector<SmartPtr<DoFDistribution> > m_vDD;

	///	Index Storage for Level (ghost / noghost)
	///	\{
		SmartPtr<DoFIndexStorage> m_spDoFIndexStrgForLevelNoGhost;
		SmartPtr<DoFIndexStorage> m_spDoFIndexStrgForLevelWithGhost;
	/// \}
};

/// base class for approximation spaces without type of algebra or dof distribution
template <typename TDomain>
class ApproximationSpace : public IApproximationSpace
{
	public:
	///	Domain type
	using domain_type = TDomain;

	///	World Dimension
		static constexpr int dim = domain_type::dim;

	///	Subset Handler type
	using subset_handler_type = typename domain_type::subset_handler_type;

	public:
	/// constructor
		ApproximationSpace(SmartPtr<TDomain> domain);

	/// constructor passing requested algebra type
		ApproximationSpace(SmartPtr<TDomain> domain, const AlgebraType& algebraType);

	/// Return the domain
		ConstSmartPtr<TDomain> domain() const {return m_spDomain;}

	///	Return the domain
		SmartPtr<TDomain> domain() {return m_spDomain;}

		int get_dim() const { return dim; }

	protected:
	///	Domain, where solution lives
		SmartPtr<TDomain> m_spDomain;
};


} // end namespace ug

#endif /* __H__UG__LIB_DISC__FUNCTION_SPACE__APPROXIMATION_SPACE__ */
