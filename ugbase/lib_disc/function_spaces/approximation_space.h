/*
 * approximation_space.h
 *
 *  Created on: 19.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__FUNCTION_SPACE__APPROXIMATION_SPACE__
#define __H__UG__LIB_DISC__FUNCTION_SPACE__APPROXIMATION_SPACE__

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
		typedef MGSubsetHandler subset_handler_type;

	///	Type of Subset Handler
		typedef MultiGrid grid_type;

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
		void print_statistic(int verboseLev = 1) const;

	///	prints statistic about DoF Distribution
		void print_statistic() const {print_statistic(1);}

	///	prints statistic on layouts
		void print_layout_statistic() const;


	///	initializes all level dof distributions
		void init_levels();

	///	initializes all surface dof distributions
		void init_surfaces();

	///	initializes all top surface dof distributions
		void init_top_surface();

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
		typedef TDomain domain_type;

	///	World Dimension
		static const int dim = domain_type::dim;

	///	Subset Handler type
		typedef typename domain_type::subset_handler_type subset_handler_type;

	public:
	/// constructor
		ApproximationSpace(SmartPtr<TDomain> domain);

	/// constructor passing requested algebra type
		ApproximationSpace(SmartPtr<TDomain> domain, const AlgebraType& algebraType);

	/// Return the domain
		ConstSmartPtr<TDomain> domain() const {return m_spDomain;}

	///	Return the domain
		SmartPtr<TDomain> domain() {return m_spDomain;}

	protected:
	///	Domain, where solution lives
		SmartPtr<TDomain> m_spDomain;
};


} // end namespace ug

#endif /* __H__UG__LIB_DISC__FUNCTION_SPACE__APPROXIMATION_SPACE__ */
