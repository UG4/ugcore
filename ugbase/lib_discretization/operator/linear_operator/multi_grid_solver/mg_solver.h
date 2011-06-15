/*
 * mg_solver.h
 *
 *  Created on: 07.12.2009
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__MULTI_GRID_SOLVER__MG_SOLVER__
#define __H__LIB_DISCRETIZATION__MULTI_GRID_SOLVER__MG_SOLVER__

// extern includes
#include <vector>
#include <iostream>

// other ug4 modules
#include "common/common.h"
#include "lib_grid/lg_base.h"
#ifdef UG_PARALLEL
	#include "lib_grid/parallelization/distributed_grid.h"
#endif
#include "lib_algebra/operator/operator_iterator_interface.h"
#include "lib_algebra/operator/operator_inverse_interface.h"
#include "lib_algebra/operator/operator_interface.h"

// library intern headers
#include "lib_discretization/function_spaces/grid_function_util.h"
#include "lib_discretization/operator/linear_operator/assembled_linear_operator.h"

namespace ug{

/// geometric multi grid preconditioner
/**
 * This class implements one step of the geometric multi grid as a
 * preconditioner for linear iteration schemes such as linear iteration, cg
 * or bicgstab.
 *
 * The coarse grid spaces are build up according to the Approximation Space
 * that is set from outside. In addition an Assembling routine must be
 * specified that is used to assemble the coarse grid matrices.
 *
 * \tparam		TApproximationSpace		Type of Approximation Space
 * \tparam		TAlgebra				Type of Algebra
 */
template <typename TApproximationSpace, typename TAlgebra>
class AssembledMultiGridCycle :
	virtual public ILinearIterator<	typename TAlgebra::vector_type,
									typename TAlgebra::vector_type>
{
	public:
	///	Approximation Space
		typedef TApproximationSpace approximation_space_type;

	///	DoFDistribution Type
		typedef typename TApproximationSpace::dof_distribution_type dof_distribution_type;

	/// Implementation type of DoFDistribution
		typedef typename dof_distribution_type::implementation_type dof_distribution_impl_type;

	///	Algebra type
		typedef TAlgebra algebra_type;

	///	type of assembling
		typedef IAssemble<dof_distribution_impl_type, algebra_type> assemble_type;

	///	Vector type
		typedef typename algebra_type::vector_type vector_type;

	///	Matrix type
		typedef typename algebra_type::matrix_type matrix_type;

	///	Level Operator Type
		typedef AssembledLinearOperator<dof_distribution_impl_type, algebra_type> operator_type;

	///	Prolongation Operator
		typedef IProlongationOperator<vector_type, vector_type> prolongation_operator_type;

	///	Projection Operator
		typedef IProjectionOperator<vector_type, vector_type> projection_operator_type;

	///	Base Solver type
		typedef ILinearOperatorInverse<vector_type, vector_type> base_solver_type;

	///	Smoother type
		typedef ILinearIterator<vector_type, vector_type> smoother_type;

	///	Own type
		typedef ILinearIterator<vector_type, vector_type> base_type;

	public:
	/// default Constructor
		AssembledMultiGridCycle() :
			m_pAss(NULL), m_pApproxSpace(NULL),
			m_topLev(0), m_baseLev(0), m_bBaseParallel(true), m_cycleType(1),
			m_numPreSmooth(1), m_numPostSmooth(1),
			m_bFullRefined(false),
			m_pSmootherPrototype(NULL),
			m_pProjectionPrototype(NULL), m_pProlongationPrototype(NULL),
			m_pBaseSolver(NULL),
			m_pDebugWriter(NULL), m_dbgIterCnt(0)
		{};

	///////////////////////////////////////////////////////////////////////////
	//	Setup
	///////////////////////////////////////////////////////////////////////////

	/// sets the assembling procedure that is used to compute coarse grid matrices
		void set_discretization(assemble_type& ass)
			{m_pAss = &ass;}

	///	sets the approximation space that is used to build up the grid hierarchy
		void set_approximation_space(approximation_space_type& approxSpace)
			{m_pApproxSpace = &approxSpace;}

	///	sets the level where exact solving is performed in the mg cycle
		void set_base_level(int baseLevel) {m_baseLev = baseLevel;}

	///	sets the base solver that is used
		void set_base_solver(base_solver_type& baseSolver)
			{m_pBaseSolver = &baseSolver;}

	///	sets if the base solver is applied in parallel
		void set_parallel_base_solver(bool bParallel) {m_bBaseParallel = bParallel;}

	///	sets the cycle type (1 = V-cycle, 2 = W-cycle, ...)
		void set_cycle_type(int type) {m_cycleType = type;}

	///	sets the number of pre-smoothing steps to be performed
		void set_num_presmooth(int num) {m_numPreSmooth = num;}

	///	sets the number of post-smoothing steps to be performed
		void set_num_postsmooth(int num) {m_numPostSmooth = num;}

	///	sets the smoother that is used
		void set_smoother(smoother_type& smoother)
			{m_pSmootherPrototype = & smoother;}

	///	sets the prolongation operator
		void set_prolongation_operator(IProlongationOperator<vector_type, vector_type>& P)
			{m_pProlongationPrototype = &P;}

	///	sets the projection operator
		void set_projection_operator(IProjectionOperator<vector_type, vector_type>& P)
			{m_pProjectionPrototype = &P;}

	///////////////////////////////////////////////////////////////////////////
	//	Linear Solver interface methods
	///////////////////////////////////////////////////////////////////////////

	///	name
		virtual const char* name() const {return "Geometric MultiGrid";}

	/// Prepare for Operator J(u) and linearization point u (current solution)
		virtual bool init(ILinearOperator<vector_type, vector_type>& J, const vector_type& u);

	///	Prepare for Linear Operartor L
		virtual bool init(ILinearOperator<vector_type, vector_type>& L);

	///	Compute new correction c = B*d
		virtual bool apply(vector_type& c, const vector_type& d);

	///	Compute new correction c = B*d and return new defect d := d - A*c
		virtual bool apply_update_defect(vector_type& c, vector_type& d);

	///	Clone
		base_type* clone();

	///	Destructor
		~AssembledMultiGridCycle() {};

 	protected:
 	/// smooth on level l, restrict defect, call lmgc (..., l-1) and interpolate correction
		bool lmgc(vector_type& c, vector_type& d, size_t lev);

	/// performs smoothing on level l, nu times
		bool smooth(vector_type& c, vector_type& d, vector_type& t,
		            IMatrixOperator<vector_type, vector_type, matrix_type>& A,
		            smoother_type& S, size_t lev, int nu);

	///	returns the number of allocated levels
		size_t num_levels() const {return m_vLevData.size();}

	///	allocates the memory
		bool top_level_required(size_t topLevel);

	///	initializes common part
		bool init_common(bool nonlinear);

	///	initializes the smoother and base solver
		bool init_smoother();

	///	initializes the coarse grid matrices for non-linear case
		bool init_non_linear_level_operator();

	///	initializes the coarse grid matrices for linear case
		bool init_linear_level_operator();

	///	initializes the smoother and base solver
		bool init_base_solver();

	///	initializes the prolongation
		bool init_prolongation();

	///	initializes the prolongation
		bool init_projection();

	///	projects a grid function from the surface to the levels
		bool project_surface_to_level(std::vector<vector_type*> vLevelFunc,
		                              const vector_type& surfFunc);

	///	projects a grid function from the levels to the surface
		bool project_level_to_surface(vector_type& surfFunc,
		                              std::vector<const vector_type*> vLevelFunc);

	///	assembles the missing matrix part on the coarse level, that must be
	///	added if the correction has been computed to ensure a correctly updated
	///	defect. (i.e. assembles A^c, with d^f -= A^c * c^c)
		bool init_missing_coarse_grid_coupling(const vector_type* u);

	protected:
	/// operator to invert (surface grid)
		operator_type* m_pSurfaceOp;

	///	assembling routine for coarse grid matrices
		assemble_type* m_pAss;

	///	approximation space for level and surface grid
		approximation_space_type* m_pApproxSpace;

	///	top level (i.e. highest level in hierarchy. This is the surface level
	///	in case of non-adaptive refinement)
		size_t m_topLev;

	///	base level (where exact inverse is computed)
		size_t m_baseLev;

	///	flag, if to solve base problem in parallel
		bool m_bBaseParallel;

	///	cylce type (1 = V-cycle, 2 = W-cylcle, ...)
		int m_cycleType;

	///	number of Presmooth steps
		int m_numPreSmooth;

	///	number of Postsmooth steps
		int m_numPostSmooth;

	///	flag indicating if grid is full refined
		bool m_bFullRefined;

	///	mapping from surface to top level (only valid in case of full refinement)
		std::vector<size_t> m_vSurfToTopMap;

	///	prototype for smoother
		smoother_type* m_pSmootherPrototype;

	///	prototype for projection operator
		projection_operator_type* m_pProjectionPrototype;

	///	prototype for prolongation operator
		prolongation_operator_type* m_pProlongationPrototype;

		////////////////////////////////////
		// Storage for each grid level
		////////////////////////////////////

		struct LevData
		{
			LevData() : pLevDD(0), A(0), Smoother(0), Projection(0), Prolongation(0),
						u(0), c(0), d(0), t(0), CoarseGridContribution(0),
						SmoothMat(0),
						su(0), sc(0), sd(0), st(0),
#ifdef UG_PARALLEL
						masterLayout(0), slaveLayout(0),
#endif
						sel(0)
			{};

			void allocate(size_t lev,
						  approximation_space_type& approxSpace,
			              assemble_type& ass,
			              smoother_type& smoother,
			              projection_operator_type& projection,
			              prolongation_operator_type& prolongation);

			bool has_ghosts() const
			{
				#ifdef UG_PARALLEL
					return !pLevDD->get_vertical_master_layout().empty();
				#else
					return false;
				#endif
			}

			void free();

			~LevData()
			{}

		//	level DoF Distribution
			dof_distribution_type* pLevDD;

		//	operator
			operator_type* A;

		//	smoother
			smoother_type* Smoother;

		//	projection operator
			projection_operator_type* Projection;

		//	prolongation operator
			prolongation_operator_type* Prolongation;

		//	vectors needed
			vector_type *u, *c, *d, *t;

		//	missing coarse grid correction
			matrix_type *CoarseGridContribution;

		//	smaller matrix for smoothing
			PureMatrixOperator<vector_type, vector_type, matrix_type> *SmoothMat;

		//	vectors needed for smoothing
			vector_type *su, *sc, *sd, *st;

#ifdef UG_PARALLEL
		//	interfaces needed for smoothing
			IndexLayout *masterLayout, *slaveLayout;
#endif

		//	map for smoothing
			std::vector<size_t> vMap;
			std::vector<int> vMapMat;

		//	selector for smoothing elements
			Selector *sel;
		};

	///	storage for all level
		std::vector<LevData> m_vLevData;

	///	base solver for the coarse problem
		base_solver_type* m_pBaseSolver;

	///	operator for base solver
		operator_type m_BaseOperator;

		std::vector<vector_type*> level_defects()
		{
			std::vector<vector_type*> vVec;
			for(size_t i = 0; i < m_vLevData.size(); ++i)
				vVec.push_back(m_vLevData[i].d);
			return vVec;
		}

		std::vector<const vector_type*> const_level_defects() const
		{
			std::vector<const vector_type*> vVec;
			for(size_t i = 0; i < m_vLevData.size(); ++i)
				vVec.push_back(m_vLevData[i].d);
			return vVec;
		}

		std::vector<vector_type*> level_corrections()
		{
			std::vector<vector_type*> vVec;
			for(size_t i = 0; i < m_vLevData.size(); ++i)
				vVec.push_back(m_vLevData[i].c);
			return vVec;
		}

		std::vector<const vector_type*> const_level_corrections() const
		{
			std::vector<const vector_type*> vVec;
			for(size_t i = 0; i < m_vLevData.size(); ++i)
				vVec.push_back(m_vLevData[i].c);
			return vVec;
		}

		std::vector<vector_type*> level_solutions()
		{
			std::vector<vector_type*> vVec;
			for(size_t i = 0; i < m_vLevData.size(); ++i)
				vVec.push_back(m_vLevData[i].u);
			return vVec;
		}

#ifdef UG_PARALLEL
	/**
	 *	gathers the vector using vertical interfaces, returns if this proc
	 *  will still have dofs on the next level (iff has no vertical slaves)
	 */
		bool gather_vertical(vector_type& d);

	/**
	 *	broadcasts the vector using vertical interfaces.
	 */
		void broadcast_vertical(vector_type& t);

	/// communicator
		pcl::ParallelCommunicator<IndexLayout> m_Com;
#endif

	public:
	///	set debug output
	/**
	 * If a DebugWriter is passed by this method, the multi grid cycle writes
	 * the level/surface vectors and matrices for debug purposes.
	 *
	 * \param[in]	debugWriter		Debug Writer to use
	 */
		void set_debug(IDebugWriter<algebra_type>* debugWriter)
		{
			m_pDebugWriter = debugWriter;
		}

	protected:
	///	writes debug output for a level vector
	/**
	 * This method writes the level vector to a debug file, if a debug writer
	 * has been set.
	 *
	 * \param[in]		vec			Level Vector to write for debug purpose
	 * \param[in]		filename	Filename
	 * \param[in]		level		grid level corresponding to the vector
	 */
		bool write_level_debug(const vector_type& vec, const char* filename, size_t lev);

	///	writes debug output for a level matrix
	/**
	 * This method writes the level matrix to a debug file, if a debug writer
	 * has been set.
	 *
	 * \param[in]		mat			Level Matrix to write for debug purpose
	 * \param[in]		filename	Filename
	 * \param[in]		level		grid level corresponding to the matrix
	 */
		bool write_level_debug(const matrix_type& mat, const char* filename, size_t lev);

	///	writes debug output for a surface vector
	/**
	 * This method writes the surface vector to a debug file, if a debug writer
	 * has been set.
	 *
	 * \param[in]		vec			Level Vector to write for debug purpose
	 * \param[in]		filename	Filename
	 */
		bool write_surface_debug(const vector_type& vec, const char* filename);

	///	writes debug output for a surface matrix
	/**
	 * This method writes the surface matrix to a debug file, if a debug writer
	 * has been set.
	 *
	 * \param[in]		mat			Level Matrix to write for debug purpose
	 * \param[in]		filename	Filename
	 */
		bool write_surface_debug(const matrix_type& mat, const char* filename);

	///	Debug Writer
		IDebugWriter<algebra_type>* m_pDebugWriter;

	///	counter for debug, to distinguish the iterations
		int m_dbgIterCnt;
};

} // end namespace ug

// include implementation
#include "mg_solver_impl.hpp"

#endif /* __H__LIB_DISCRETIZATION__MULTI_GRID_SOLVER__MG_SOLVER__ */
