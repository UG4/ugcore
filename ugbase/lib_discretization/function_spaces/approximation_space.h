/*
 * approximation_space.h
 *
 *  Created on: 19.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__FUNCTION_SPACE__APPROXIMATION_SPACE__
#define __H__LIBDISCRETIZATION__FUNCTION_SPACE__APPROXIMATION_SPACE__

#ifdef UG_PARALLEL
	#include "lib_discretization/parallelization/parallelization.h"
#endif

#include "lib_discretization/dof_manager/mg_dof_manager.h"
#include "./grid_function.h"

namespace ug{

/// base class for approximation spaces without type of algebra or dof distribution
template <typename TDomain>
class IApproximationSpace : public FunctionPattern
{
	public:
	///	Domain type
		typedef TDomain domain_type;

	///	World Dimension
		static const int dim = domain_type::dim;

	///	Subset Handler type
		typedef typename domain_type::subset_handler_type subset_handler_type;

	public:
	/// Default constructor
		IApproximationSpace() : m_pDomain(NULL), m_pMGSH(NULL) {};

	///	Assign domain
		void assign_domain(domain_type& domain)
		{
			m_pDomain = &domain;
			m_pMGSH = &domain.get_subset_handler();

			return this->set_subset_handler(*m_pMGSH);
		}

	/// Return the domain
		const domain_type& get_domain() const {return *m_pDomain;}

	///	Return the domain
		domain_type& get_domain() {return *m_pDomain;}

	///	virtual destructor
		virtual ~IApproximationSpace()	{}

	protected:
	///	Domain, where solution lives
		domain_type* m_pDomain;

	/// grid or multigrid or subsethandler, where elements are stored
		subset_handler_type* m_pMGSH;
};

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
 *
 * The user is responsible to free the memory of grid functions.
 *
 * \tparam		TDomain				Domain type
 * \tparam		TDoFDistribution	DoF Distribution type
 * \tparam		TAlgebra			Algebra type
 */
template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
class ApproximationSpace : public IApproximationSpace<TDomain>{
	private:
	//	 to make it more readable
		typedef ApproximationSpace<TDomain, TDoFDistribution, TAlgebra> this_type;

	public:
	///	Type of Domain, where DoFs are defined
		typedef TDomain domain_type;

	///	Type of Grid, where DoFs are defined
		typedef typename domain_type::grid_type grid_type;

	///	Type of Subset Handler, where DoFs are defined
		typedef typename domain_type::subset_handler_type subset_handler_type;

	///	Type of Algebra used
		typedef TAlgebra algebra_type;

	///	Type of DoF Distribution
		typedef IDoFDistribution<TDoFDistribution> dof_distribution_type;

	///	Type of DoF Manager used
		#ifdef UG_PARALLEL
			typedef ParallelMGDoFManager<MGDoFManager<TDoFDistribution> > dof_manager_type;
		#else
			typedef MGDoFManager<TDoFDistribution> dof_manager_type;
		#endif

	///	Type of Grid function used
		#ifdef UG_PARALLEL
			typedef ParallelGridFunction<GridFunction<TDomain, TDoFDistribution, TAlgebra> > function_type;
		#else
			typedef GridFunction<TDomain, TDoFDistribution, TAlgebra> function_type;
		#endif

	public:
	///	Constructor
		ApproximationSpace() : m_bInit(false) {};

	///	Destructor
		~ApproximationSpace(){}

	///	initializes the Approximation Space
		bool init();

	///	returns if ansatz space is supported
		virtual bool supports_trial_space(LFEID& id) const
			{return TDoFDistribution::supports_trial_space(id);}

	///	sets the grouping of indices on the same object
		void set_grouping(bool bGrouped)
		{
			if(m_bInit)
				throw(UGFatalError("ApproximationSpace: cannot change grouping "
						"strategy after initialization."));
			m_MGDoFManager.set_grouping(bGrouped);
		}

	///	prints statistic about DoF Distribution
		void print_statistic(int verboseLev) const
			{m_MGDoFManager.print_statistic(verboseLev);}

	///	prints statistic about DoF Distribution
		void print_statistic() const {print_statistic(1);}

	///	prints statistic on layouts
		void print_layout_statistic(int verboseLev = 1) const
			{m_MGDoFManager.print_layout_statistic(verboseLev);}

	///	prints statistic on layouts
		void print_layout_statistic() const {print_layout_statistic(1);}

	///	prints statistic on local dof distribution
		void print_local_dof_statistic(int verboseLev = 1) const
			{m_MGDoFManager.print_local_dof_statistic(verboseLev);}

	///	prints statistic on local dof distribution
		void print_local_dof_statistic() const {print_local_dof_statistic(1);}

	///	defragments the index set of the DoF Distribution
		void defragment() {m_MGDoFManager.defragment();}

	/// create a new grid function of this approximation space
		function_type* create_level_function(size_t level);

	/// create a new grid function of this approximation space
		function_type* create_surface_function();

	///	returns the surface dof distribution
		dof_distribution_type& get_surface_dof_distribution();

	///	returns the surface dof distribution
		const dof_distribution_type& get_surface_dof_distribution() const;

	///	returns the surface view
		const SurfaceView* get_surface_view() const {return m_MGDoFManager.get_surface_view();}

	///	returns the level dof distributions
		std::vector<const dof_distribution_type*> get_level_dof_distributions() const
				{return m_MGDoFManager.get_level_dof_distributions();}

	///	returns the level dof distribution
		dof_distribution_type& get_level_dof_distribution(size_t level);

	///	returns the level dof distribution
		const dof_distribution_type& get_level_dof_distribution(size_t level) const;

	///	returns the number of level
		size_t num_levels() const {return m_MGDoFManager.num_levels();}

	protected:
	///	Init flag
		bool m_bInit;

	/// Dof manager used for this Approximation Space
		dof_manager_type m_MGDoFManager;
};

} // end namespace ug

// include implementation
#include "approximation_space_impl.h"

#endif /* __H__LIBDISCRETIZATION__FUNCTION_SPACE__APPROXIMATION_SPACE__ */
