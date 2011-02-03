/* TODO:
 * Solve Schur complement w.r.t. "Pi system"
 * Fuer 'apply_F()' eigene Klasse einfuehren, mit 'apply_sub()' und so? - Laenge von Vector 'lambda'?
 * If finished: Check for unnecessary "set zero" operations!
 * Man koennte, statt einzelner Layout levels, gleich immer den ganzen vector von layouts holen (so alle benoetigten dort drin sind)!
 *
 * Outline of the algorithm (after Klawonn, A., Widlund, O., Dryja, M.:
 * "Dual-Primal FETI methods for three-dimensional elliptic problems with
 * heterogeneous coefficients", SIAM J. Numer. Anal., 40/2002a, 159--179.):
 * 1. eliminate the primal variables associated with the interior nodes,
 *    the vertex nodes designated as primal (as well as the Lagrange multipliers related to the optional constraints - if any)
 *    ==> reduced system with Schur complement \tilde{S} related to the subspace \tilde{W}_{\Delta},
 *        and a reduced right hand side \tilde{f}_{\Delta}
 * 2. Reformulate as a variational problem, with a vector of Lagrange multipliers
 *    ==> saddle point problem
 * 3. Eliminate subvectors u_{\Delta}
 *    ==> system for the dual variables: F \lambda = d
 * 4. Instead of solving ... action of the inverse of the Schur complement can be obtained
 *    by solving the entire linear system from which it originates after augmenting the given right hand side with zeros.
 *    Group all the interior and dual variables of each subdomain together
 *    and factor the resulting blocks in paralle across the subdomains using a good ordering algorithm.
 *
 *    The contributions of the remaining Schur complement, of the primal variables,
 *    can also be computed locally prior to subassembly and factorization of ths final, global part of the linear system of equations.
 */
/*
 * feti.h
 *
 *  Created on: 21.01.2011
 *      Author: iheppner, avogel
 */

#ifndef __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__FETI__
#define __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__FETI__

namespace ug{

#ifdef UG_PARALLEL

#include <iostream>
#include <sstream>
#include <string>
#include "lib_algebra/operator/operator_inverse_interface.h"
#include "lib_algebra/parallelization/parallelization.h"
#include "lib_algebra/operator/debug_writer.h"
#include "pcl/pcl.h"

/// Groups all Feti Layouts together and provides useful funktions on those
template <typename TAlgebra>
class FetiLayouts
{
	public:
	// 	Algebra type
		typedef TAlgebra algebra_type;

	// 	Vector type
		typedef typename TAlgebra::vector_type vector_type;

	// 	Matrix type
		typedef typename TAlgebra::matrix_type matrix_type;

	public:
		FetiLayouts() : m_pMasterStdLayout(NULL), m_pSlaveStdLayout(NULL) {}

	//	assigns standard layouts, creates the feti layouts
		void create_layouts(IndexLayout& stdMasterLayout,
		                    IndexLayout& stdSlaveLayout,
		                    pcl::ProcessCommunicator& stdProcessCom,
		                    size_t numIndices,
		                    pcl::IDomainDecompositionInfo& DDInfo,
		                    bool bDebug = false)
		{
		//	if no indices, nothing to do
			if(numIndices == 0) return;

		//	remember std layouts
			m_pMasterStdLayout = &stdMasterLayout;
			m_pSlaveStdLayout = &stdSlaveLayout;
			m_stdProcessCom = stdProcessCom;

		//	create FETI Layouts:
 // 	\todo: For some documentation info see mail by S. Reiter, 30. Januar 2011 16:10:52 MEZ
//             (if the information therein is not already outdated ...)
			BuildDomainDecompositionLayouts(m_masterDualLayout, m_slaveDualLayout,
							m_masterInnerLayout, m_slaveInnerLayout, m_masterDualNbrLayout,
							m_slaveDualNbrLayout, m_masterPrimalLayout, m_slavePrimalLayout,
							*m_pMasterStdLayout, *m_pSlaveStdLayout,
							(int)(numIndices - 1), DDInfo);

		//	create local feti block communicator
			int localSubdomID = DDInfo.map_proc_id_to_subdomain_id(pcl::GetProcRank());
			pcl::ProcessCommunicator worldComm;
			for(int i = 0; i < DDInfo.get_num_subdomains(); ++i){
				if(localSubdomID == i)
					m_localFetiBlockComm = worldComm.create_sub_communicator(true);
				else
					worldComm.create_sub_communicator(false);
			}

		//	test output
			if(bDebug)
			{
				pcl::ParallelCommunicator<IndexLayout> comTmp;
				PrintLayout(comTmp, m_masterInnerLayout, m_slaveInnerLayout);
			}
		}

	//	standard layouts
		IndexLayout& get_std_master_layout() {return *m_pMasterStdLayout;}
		IndexLayout& get_std_slave_layout() {return *m_pSlaveStdLayout;}
		pcl::ProcessCommunicator& get_std_process_communicator() {return m_stdProcessCom;}

	//	inner layouts
		IndexLayout& get_inner_master_layout() {return m_masterInnerLayout;}
		IndexLayout& get_inner_slave_layout() {return m_slaveInnerLayout;}
		pcl::ProcessCommunicator& get_inner_process_communicator() {return m_localFetiBlockComm;}

	//	dual layouts
		IndexLayout& get_dual_master_layout() {return m_masterDualLayout;}
		IndexLayout& get_dual_slave_layout() {return m_slaveDualLayout;}

	//	dual nbr layouts
		IndexLayout& get_dual_nbr_master_layout() {return m_masterDualNbrLayout;}
		IndexLayout& get_dual_nbr_slave_layout() {return m_slaveDualNbrLayout;}

	//	primal layouts
		IndexLayout& get_primal_master_layout() {return m_masterPrimalLayout;}
		IndexLayout& get_primal_slave_layout() {return m_slavePrimalLayout;}

	public:

	//	sets standard communication layouts and communicators
		void vec_use_std_communication(vector_type& vec)
		{
			vec.set_slave_layout(get_std_slave_layout());
			vec.set_master_layout(get_std_master_layout());
			vec.set_process_communicator(get_std_process_communicator());
		}

	//	sets inner communication layouts and communicators
		void vec_use_inner_communication(vector_type& vec)
		{
			vec.set_slave_layout(get_inner_slave_layout());
			vec.set_master_layout(get_inner_master_layout());
			vec.set_process_communicator(get_inner_process_communicator());
		}

	public:
		number vec_norm_on_dual(vector_type& vec)
		{
		//	forward to VecProc
			return sqrt(vec_prod_on_dual(vec, vec));
		}

		void vec_scale_add_on_dual(vector_type& vecDest,
		                           number alpha1, const vector_type& vecSrc1,
		                           number alpha2, const vector_type& vecSrc2)
		{
			VecScaleAddOnLayout(&vecDest, alpha1, &vecSrc1, alpha2, &vecSrc2, m_slaveDualLayout);
			VecScaleAddOnLayout(&vecDest, alpha1, &vecSrc1, alpha2, &vecSrc2, m_masterDualLayout);
		}

		void vec_scale_add_on_primal(vector_type& vecDest,
		                           number alpha1, const vector_type& vecSrc1,
		                           number alpha2, const vector_type& vecSrc2)
		{
			VecScaleAddOnLayout(&vecDest, alpha1, &vecSrc1, alpha2, &vecSrc2, m_slavePrimalLayout);
			VecScaleAddOnLayout(&vecDest, alpha1, &vecSrc1, alpha2, &vecSrc2, m_masterPrimalLayout);
		}

		void vec_scale_append_on_dual(vector_type& vecInOut,
		                              const vector_type& vecSrc1, number alpha1)
		{
			VecScaleAppendOnLayout(&vecInOut, &vecSrc1, alpha1, m_slaveDualLayout);
			VecScaleAppendOnLayout(&vecInOut, &vecSrc1, alpha1, m_masterDualLayout);
		}

		void vec_set_on_dual(vector_type& vecSrc, number alpha)
		{
			VecSetOnLayout(&vecSrc, alpha, m_slaveDualLayout);
			VecSetOnLayout(&vecSrc, alpha, m_masterDualLayout);
		}

		number vec_prod_on_dual(const vector_type& vecSrc1, const vector_type& vecSrc2)
		{
		//	reset result
			number prod = 0.0, prodTmp = 0.0;

		//	add prod on master
			VecProdOnLayout(prodTmp, &vecSrc1, &vecSrc2, m_masterDualLayout);
			prod += prodTmp;

		//	add prod on slave
			VecProdOnLayout(prodTmp, &vecSrc1, &vecSrc2, m_slaveDualLayout);
			prod += prodTmp;

		//	return result
			return prod;
		}

		void vec_set_on_primal(vector_type& vecSrc, number alpha)
		{
			VecSetOnLayout(&vecSrc, alpha, m_slavePrimalLayout);
			VecSetOnLayout(&vecSrc, alpha, m_masterPrimalLayout);
		}

		void vec_set_excl_primal(vector_type& vecInOut, number value)
		{
			VecSetExcludingLayout(&vecInOut, value, m_slavePrimalLayout);
			VecSetExcludingLayout(&vecInOut, value, m_masterPrimalLayout);
		}

		void vec_set_excl_dual(vector_type& vecInOut,number value)
		{
			VecSetExcludingLayout(&vecInOut, value, m_slaveDualLayout);
			VecSetExcludingLayout(&vecInOut, value, m_masterDualLayout);
		}

	public:
		void mat_set_dirichlet_on_dual(matrix_type& mat)
		{
			MatSetDirichletOnLayout(&mat, m_slaveDualLayout);
			MatSetDirichletOnLayout(&mat, m_masterDualLayout);
		}

		void mat_set_dirichlet_on_primal(matrix_type& mat)
		{
			MatSetDirichletOnLayout(&mat, m_slavePrimalLayout);
			MatSetDirichletOnLayout(&mat, m_masterPrimalLayout);
		}

	protected:
	//	Standard Layouts
		IndexLayout* m_pMasterStdLayout;
		IndexLayout* m_pSlaveStdLayout;
		pcl::ProcessCommunicator	m_stdProcessCom;

	//	Layouts and Communicator for Inner variables
		IndexLayout m_masterInnerLayout;
		IndexLayout m_slaveInnerLayout;
		pcl::ProcessCommunicator	m_localFetiBlockComm;

	//	Layouts for Dual variables
		IndexLayout m_masterDualLayout;
		IndexLayout m_slaveDualLayout;

	//	Layouts for Dual Neighbour variables
		IndexLayout m_masterDualNbrLayout;
		IndexLayout m_slaveDualNbrLayout;

	//	Layouts for Primal variables
		IndexLayout m_masterPrimalLayout;
		IndexLayout m_slavePrimalLayout;
};

///	Application of the "jump operator" \f$B_{\Delta}\f$:
/// 'ComputeDifferenceOnDelta()': Apply \f$B_{\Delta}\f$ to \f$u_{\Delta}\f$
/**
 * \f$B_{\Delta}\f$ computes the difference between the double-valued unknowns \f$u_{\Delta}\f$.
 * This computation is only unique up to the sign of the difference. Thus, we can
 * freely decide it, but have then to stay with the choice.
 *
 * \param[out]		diff				destination vector for computed differences on "Delta layout"
 * \param[in]		u					vector \f$u_{\Delta}\f$
 * \param[in]		masterLayoutIn		master layout to operate on (caller has to provide a "Delta layout")
 * \param[in]		slaveLayoutIn		slave  layout to operate on (caller has to provide a "Delta layout")
 * \param[in]		masterNbrLayoutIn	master layout to operate on (caller has to provide a "Delta Nbr layout")
 * \param[in]		slaveNbrLayoutIn	slave  layout to operate on (caller has to provide a "Delta Nbr layout")
 */
template <typename TVector>
void ComputeDifferenceOnDelta(TVector& diff, const TVector& u,
							  IndexLayout&    masterLayoutIn,
							  IndexLayout&     slaveLayoutIn,
							  IndexLayout& masterNbrLayoutIn,
							  IndexLayout& slaveNbrLayoutIn)
{
	// Reset all values
	diff.set(0.0);

	// All masters subtract the values of the slave, all slaves subtract the values
	// of the master (communication is performed) ...
	VecSubtractOnLayout(&diff, masterLayoutIn, slaveLayoutIn);

	// ... and slaves multiplies the result by '-1' (no communication is performed)
	VecScaleOnLayout(&diff, -1.0, slaveLayoutIn);

	// ... and copy master values to additional slaves living in "Dual neighbour layout" (communication is performed)
	VecCopyOnLayout(&diff, masterNbrLayoutIn, slaveNbrLayoutIn);

	return;
};

/// 'ComputeDifferenceOnDeltaTransposed()': Apply \f$B_{\Delta}^T\f$ to \f$\lambda\f$
/**
 * \f$B_{\Delta}\f$ computes the difference between the double-valued unknowns \f$u_{\Delta}\f$.
 * This computation is only unique up to the sign of the difference. Thus, we can
 * freely decide it, but have then to stay with the choice (no communication is performed).
 *
 * \param[out]		lambda			vector of Lagrange multipliers
 * \param[in]		f				vector \f$f_{\Delta}\f$ living on "Delta layout"
 * \param[in]		masterLayoutIn	master layout to operate on (caller has to provide a "Delta layout")
 * \param[in]		slaveLayoutIn	slave  layout to operate on (caller has to provide a "Delta layout")
 */
template <typename TVector>
void ComputeDifferenceOnDeltaTransposed(TVector& f, const TVector& lambda,
										IndexLayout& masterLayoutIn,
										IndexLayout& slaveLayoutIn,
										IndexLayout& slaveNbrLayoutIn)
{
	

	// (a) Reset all values
	f.set(0.0);
	//VecSetOnLayout(&f, 0.0, masterLayoutIn);
	//VecSetOnLayout(&f, 0.0, slaveLayoutIn);

	// (b) Copy values on \Delta
	// 1. All masters set their values equal to $\lambda$
	VecScaleAppendOnLayout(&f, &lambda,  1.0,   masterLayoutIn);
	// 2. All slaves set their values equal to $-\lambda$
	VecScaleAppendOnLayout(&f, &lambda, -1.0,    slaveLayoutIn);
	VecScaleAppendOnLayout(&f, &lambda, -1.0, slaveNbrLayoutIn);

	return;

};

/// operator implementation of the local Schur complement
/**
 * This operator is the application of the local Schur complement \f$S_{\Delta}^{i}\f$.
 * The underlying matrix must have at least two layouts. The first layout, layout level 0,
 * will be used to describe subdomain *internal* interfaces (i.e. "pure" processor interfaces),
 * all other layouts are used to identify the boundary,
 * layout level 1: \Delta (edges of subdomains),
 * layout level 2: \Pi (vertices of subdomains, a.k.a. "cross points") - will be constructed here,
 * and the Schur complement is build w.r.t. to these variables.
 */
template <typename TAlgebra>
class LocalSchurComplement
	: public ILinearOperator<	typename TAlgebra::vector_type,
	  	  	  	  	  	  	  	typename TAlgebra::vector_type>
{
	public:
	// 	Algebra type
		typedef TAlgebra algebra_type;

	// 	Vector type
		typedef typename TAlgebra::vector_type vector_type;

	// 	Matrix type
		typedef typename TAlgebra::matrix_type matrix_type;

	public:
	///	constructor
		LocalSchurComplement();

	///	name of solver
		virtual const char* name() const {return "Local Schur Complement Solver";}

	///	sets a sequential Dirichlet solver
	/**
	 * This method sets the Dirichlet Solver that is used to invert the
	 * inner matrix \f$A_{II}\f$
	 */
		void set_dirichlet_solver(ILinearOperatorInverse<vector_type, vector_type>& dirichletSolver)
		{
		//	remember the Dirichlet Solver
			m_pDirichletSolver = &dirichletSolver;
		}

	///	set debug output
		void set_debug(IDebugWriter<algebra_type>* debugWriter)
		{
			m_pDebugWriter = debugWriter;
		}

	///	set original matrix from which the local Schur complement is constructed
	/**
	 * Using this method, the original matrix A is set. Given the matrix in the
	 * form \f$ A = \sum\limits_{p=1}^{N} A^{(p)}\f$ in additive form, with
	 * \f{align*}
	 * A^{(p)}
	 * \begin{pmatrix}
	 * A_{II}^{(p)} & A_{I \Delta}^{(p)} & A_{I \Pi}^{(p)} \\
	 * A_{\Delta I}^{(p)} & A_{\Delta \Delta}^{(p)} & A_{\Delta \Pi}^{(p)} \\
	 * A_{\Pi I}^{(p)} & A_{\Pi \Delta}^{(p)} & A_{\Pi \Pi}^{(p)}
	 * \end{pmatrix}
	 * \f}
	 *
	 * the local Schur complement is the processwise application of the
	 * operator
	 *
	 * \f{align*}
	 * S_{\Delta \Delta} = A_{\Delta \Delta} -
	 * \begin{pmatrix} A_{I \Delta}^T & A_{\Pi \Delta}^T \end{pmatrix}
	 * \begin{pmatrix} A_{II} & A_{I \Pi} \\ A_{\Pi I} & A_{\Pi \Pi} \end{pmatrix}^{-1}
	 * \begin{pmatrix} A_{I \Delta} \\ A_{\Pi \Delta} \end{pmatrix}	 *
	 * \f}
	 */
		void set_matrix(IMatrixOperator<vector_type, vector_type, matrix_type>& A)
		{
		//	save current operator
			m_pOperator = &A;
		}

	///	sets the primal layouts
		void set_feti_layouts(FetiLayouts<algebra_type>& fetiLayouts)
		{
			m_pFetiLayouts = &fetiLayouts;
		}

	/// implementation of the operator for the solution dependent initialization.
		bool init(const vector_type& u) {return init();}

	///	initializes the solver for operator A
	/**
	 * This method must be called, before the apply() method can be invoked.
	 * It has to be called each time, when the matrix has been replaced. A deep
	 * copy of the matrix is then constructed and in this copy the rows belonging
	 * to the \f$\Delta\f$ and \f$\Pi\f$ unknowns are set to identity rows. This
	 * matrix is used in the solution of the local dirichlet problem.
	 */
		virtual bool init();

	///	applies the Schur complement built from matrix operator set via 'set_matrix()'
	/// to 'u' and returns the result 'f := S times u'
		virtual bool apply(vector_type& f, const vector_type& u);

	///	solves the system
		virtual bool apply_sub(vector_type& f, const vector_type& u);

		// destructor
		virtual ~LocalSchurComplement() {};

	protected:
		bool write_debug(const vector_type& vec, const char* filename)
		{
		//	if no debug writer set, we're done
			if(m_pDebugWriter == NULL) return true;

		//	write
			return m_pDebugWriter->write_vector(vec, filename);
		}

	protected:
	// 	Operator that is inverted by this Inverse Operator
		IMatrixOperator<vector_type,vector_type,matrix_type>* m_pOperator;

	// 	Parallel Matrix
		matrix_type* m_pMatrix;

	//	Feti Layouts
		FetiLayouts<algebra_type>* m_pFetiLayouts;

	//	Copy of matrix
		PureMatrixOperator<vector_type, vector_type, matrix_type> m_DirichletOperator;

	// 	Parallel Dirichlet Matrix
		matrix_type* m_pDirichletMatrix;

	// 	Linear Solver to invert the local Dirichlet problems
		ILinearOperatorInverse<vector_type,vector_type>* m_pDirichletSolver;

	//	Debug Writer
		IDebugWriter<algebra_type>* m_pDebugWriter;

}; /* end class LocalSchurComplement */

/* 1.7 Application of \f${\tilde{S}_{\Delta \Delta}}^{-1}\f$ */ 
/// operator implementation of the inverse of the Schur complement w.r.t. the "Delta unknowns"
/**
 * This operator provides the application of the inverse of the Schur complement w.r.t. the "Delta unknowns",
 * \f${\tilde{S}_{\Delta \Delta}}^{-1}\f$.
 */
template <typename TAlgebra>
class SchurComplementInverse
	: public ILinearOperatorInverse<	typename TAlgebra::vector_type,
	  	  	  	  	  	  	  			typename TAlgebra::vector_type>
{
	public:
	// 	Algebra type
		typedef TAlgebra algebra_type;

	// 	Vector type
		typedef typename TAlgebra::vector_type vector_type;

	// 	Matrix type
		typedef typename TAlgebra::matrix_type matrix_type;

	public:
	///	constructor
		SchurComplementInverse();

	///	name of class
		virtual const char* name() const {return "Schur Complement Inverse";}

	//	set debug output
		void set_debug(IDebugWriter<algebra_type>* debugWriter)
		{
			m_pDebugWriter = debugWriter;
		}

	///	sets the Neumann solver
		void set_neumann_solver(ILinearOperatorInverse<vector_type, vector_type>& neumannSolver)
		{
		//	remember the Dirichlet Solver
			m_pNeumannSolver = &neumannSolver;
		}

	///	sets the primal layouts
		void set_feti_layouts(FetiLayouts<algebra_type>& fetiLayouts)
		{
			m_pFetiLayouts = &fetiLayouts;
		}

	// 	Init for Linear Operator L
		virtual bool init(ILinearOperator<vector_type, vector_type>& L);


	// 	Init for Linear Operator J and Linearization point (current solution)
		virtual bool init(ILinearOperator<vector_type, vector_type>& J, const vector_type& u)
		{
			return init(J);
		}

	// 	Solve A*u = f, such that u = A^{-1} f
		virtual bool apply(vector_type& u, const vector_type& f);

	// 	Solve A*u = f, such that u = A^{-1} f
	// 	This is done by iterating: u := u + B(f - A*u)
	// 	In f the last defect f := f - A*u is returned
		virtual bool apply_return_defect(vector_type& u, vector_type& f);

	///	sets a convergence check
		void set_convergence_check(IConvergenceCheck& convCheck)
		{
			m_pConvCheck = &convCheck;
			m_pConvCheck->set_offset(3);
		}

	/// returns the convergence check
		IConvergenceCheck* get_convergence_check() {return m_pConvCheck;}

	//  destructor
		virtual ~SchurComplementInverse() {};

	protected:
		bool write_debug(const vector_type& vec, const char* filename)
		{
		//	if no debug writer set, we're done
			if(m_pDebugWriter == NULL) return true;

		//	write
			return m_pDebugWriter->write_vector(vec, filename);
		}

	protected:
	// 	Operator that is inverted by this Inverse Operator
		IMatrixOperator<vector_type,vector_type,matrix_type>* m_A;

	// 	Parallel Matrix to invert
		matrix_type* m_pMatrix;

	//	Feti Layouts
		FetiLayouts<algebra_type>* m_pFetiLayouts;

	//	Copy of matrix
		PureMatrixOperator<vector_type, vector_type, matrix_type> m_NeumannOperator;

	// 	Parallel Neumann Matrix
		matrix_type* m_pNeumannMatrix;

	//	Neumann Solver
		ILinearOperatorInverse<vector_type, vector_type>* m_pNeumannSolver;

	//	Process gathering the Schur Complement w.r.t. Primal unknowns
		int m_primalRootProc;

	//	Layouts for all to one communication
		IndexLayout m_masterAllToOneLayout;
		IndexLayout m_slaveAllToOneLayout;
		pcl::ProcessCommunicator m_allToOneProcessComm;

	//	Schur Complement operator for gathered matrix
		PureMatrixOperator<vector_type, vector_type, matrix_type> m_RootSchurComplementOp;

	//	Matrix for one proc schur complement
		matrix_type* m_pRootSchurComplementMatrix;

	// 	Convergence Check
		IConvergenceCheck* m_pConvCheck;

	//	Debug Writer
		IDebugWriter<algebra_type>* m_pDebugWriter;
}; /* end class 'SchurComplementInverse' */

/// TODO:
// * Use members for delta layouts in LogIndexLayoutOnAllProcs() etc.
// * members for pi layouts will become pointers if no longer built here by ExtractCrossPointLayouts()!

//*  Preliminary: the actual 'FETISolver' - if not derived from 'CGSolver'!
// In the moment more or less only stuff copied from 'DirichletDirichletSolver'!
/// operator implementation of the local Schur complement
/**
 * This operator implements a FETI-DP solver, see e.g.
 * "Domain Decomposition Methods -- Algorithms and Theory",
 * A. Toselli, O. Widlund, Springer 2004, sec. 1.3.5, p. 12ff.
 */
template <typename TAlgebra>
class FETISolver : public IMatrixOperatorInverse<	typename TAlgebra::vector_type,
													typename TAlgebra::vector_type,
													typename TAlgebra::matrix_type>
{
	public:
	// 	Algebra type
		typedef TAlgebra algebra_type;

	// 	Vector type
		typedef typename TAlgebra::vector_type vector_type;

	// 	Matrix type
		typedef typename TAlgebra::matrix_type matrix_type;

	public:
	///	constructor
		FETISolver();

	///	name of solver
		virtual const char* name() const {return "FETI Solver";}

	///	sets a convergence check
		void set_convergence_check(IConvergenceCheck& convCheck)
		{
			m_pConvCheck = &convCheck;
			m_pConvCheck->set_offset(3);
		}

	/// returns the convergence check
		IConvergenceCheck* get_convergence_check() {return m_pConvCheck;}

	///	sets the Dirichlet solver
		void set_dirichlet_solver(ILinearOperatorInverse<vector_type, vector_type>& dirichletSolver)
		{
		//	remember the Dirichlet Solver
			m_pDirichletSolver = &dirichletSolver;
		}

	///	sets the Neumann solver
		void set_neumann_solver(ILinearOperatorInverse<vector_type, vector_type>& neumannSolver)
		{
		//	remember the Dirichlet Solver
			m_pNeumannSolver = &neumannSolver;
		}

	///	sets the Dual Dirichlet solver
		void set_dual_dirichlet_solver(ILinearOperatorInverse<vector_type, vector_type>& dualDirichletSolver)
		{
		//	remember the Dirichlet Solver
			m_pDualDirichletSolver = &dualDirichletSolver;
		}

	//	set debug output
		void set_debug(IDebugWriter<algebra_type>* debugWriter)
		{
			m_pDebugWriter = debugWriter;
			m_LocalSchurComplement.set_debug(m_pDebugWriter);
			m_SchurComplementInverse.set_debug(m_pDebugWriter);
		}

	///	initializes the solver for operator A
		virtual bool init(IMatrixOperator<vector_type, vector_type, matrix_type>& A);

	///	function which applies matrix \f$F\f$ of the reduced system ("Delta system") to a vector \f$v\f$
	/**
	 * This function applies matrix \f$F := B_{\Delta} \tilde{S}^{-1} B_{\Delta}^T\f$
	 * to a vector \f$v\f$. \f$v\f$ can be:
	 * (a) unknown vector of Lagrange multipliers, lambda, and
	 * (b) search direction 'p' in cg method.
	 *
	 * \param[in]		v				vector \f$v\f$ living on "Delta layout"
	 * \param[out]		f				result of application of \f$F\f$
	 */
		bool apply_F(vector_type& f, const vector_type& v);

	///	function which computes right hand side vector 'd' of the dual unknowns
	/**
	 * This function computes \f$d := B_{\Delta} \tilde{S}_{\Delta \Delta}^{-1} \tilde{f}_{\Delta}\f$
	 * to a vector \f$v\f$. \f$v\f$ can be:
	 * (a) unknown vector of Lagrange multipliers, lambda, and
	 * (b) search direction 'p' in cg method.
	 *
	 * \param[in]		f				vector \f$\tilde{f}_{\Delta}\f$
	 * \param[out]		d				right hand side vector \f$d\f$ of reduced system
	 */
		bool compute_d(vector_type& d, const vector_type& f);

	///	function which computes the adapted right hand side vector '\tilde{f}_{\Delta}' of the reduced system ("Delta system")
	/**
	 * This function computes \f$ \tilde{f}_{\Delta} := f_{\Delta}
	 *				 - A_{\Delta \{I \Pi\}} (A_{\{I \Pi\} \{I \Pi\}})^{-1} f_{\{I \Pi\}}\f$
	 * to a vector \f$f\f$.
	 *
	 * \param[in]		f				vector \f$f\f$
	 * \param[Out]		tildeF			vector \f$\tilde{f}_{\Delta}\f$
	 */
		bool compute_tilde_f(vector_type& tildeF, const vector_type& f);

	///	function which applies diagonal scaling matrix \f$D_{\Delta}^{(i)}\f$ to a vector \f$v\f$
		bool Apply_ScalingMatrix(vector_type& s, const vector_type& v) // maybe restrict to layout
		{
			// scaling operator is identity
			s = v;
			// more general: m_Ddelta.apply(s,v), with additional member 'matrix_type m_Ddelta;'

			//	we're done
			return true;
		}

	///	function which applies matrix \f$M^{-1}^{(i)}\f$ to a vector \f$r
	/**
	 * This function applies matrix \f$M^{-1}^{(i)} := D_{\Delta}^{(i)} B_{\Delta}^{(i)} S_{\Delta}^{(i)} {B_{\Delta}^{(i)}}^T D_{\Delta}^{(i)}\f$
	 * to a vector \f$r\f$.
	 *
	 * \param[in]		r				vector \f$r\f$ living on "Delta layout"
	 * \param[out]		z				result of application of \f$$M^{-1}\f$
	 */
		bool apply_M_inverse(vector_type& z, const vector_type& r);

		bool apply_M_inverse_with_identity_scaling(vector_type& z, const vector_type& r);

	//	TODO:
	// 1. x = lambda; Compute 'd' before calling apply_return_defect()! Note: 'ApplyLinearSolver()' calls 'apply(), not 'apply_return_defect()'!
	// 2. With \f$\lambda\f$ found, back solve for \f$u_{\Delta}\f$:
	//    \f$u_{\Delta} = {\tilde{S}_{\Delta \Delta}}^{-1} ({\tilde{f}_{\Delta}} - B_{\Delta}^T \lambda).\f$
	// 3. Assemble this and the solutions for \f$u_{I}\f$ and \f$u_{\Pi}\f$ to the global solution vector

	///	solves the reduced system \f$F \lambda = d\f$ with preconditioned cg method
	///	and returns the last defect of iteration in rhs
	/// (derived from 'CGSolver::apply_return_defect()')
		virtual bool apply_return_defect(vector_type& lambda, vector_type& d);

	///	solves the system
		virtual bool apply(vector_type& x, const vector_type& b)
		{
		//	copy defect
			vector_type d; d.resize(b.size());
			d = b;

		//	solve on copy of defect
			return apply_return_defect(x, d);
		}

		// destructor
		virtual ~FETISolver() {};

	protected:
	//	Prepare the convergence check
		void prepare_conv_check()
		{
			m_pConvCheck->set_name(name());
			m_pConvCheck->set_symbol('%');
			m_pConvCheck->set_name(name());
	}

	protected:
		bool write_debug(const vector_type& vec, const char* filename)
		{
		//	add iter count to name
			std::string name(filename);
			char ext[20]; sprintf(ext, "_iter%03d", m_iterCnt);
			name.append(ext);

		//	if no debug writer set, we're done
			if(m_pDebugWriter == NULL) return true;

		//	write
			return m_pDebugWriter->write_vector(vec, name.c_str());
		}

		int m_iterCnt;

	protected:
	// 	Operator that is inverted by this Inverse Operator
		IMatrixOperator<vector_type,vector_type,matrix_type>* m_pOperator;

	// 	Parallel Matrix to invert
		matrix_type* m_pMatrix;

	//	Feti Layouts
		FetiLayouts<algebra_type> m_fetiLayouts;

	//	Local Schur complement for each feti subdomain
		LocalSchurComplement<algebra_type> m_LocalSchurComplement;

	//	Dirichlet solver for inverse of A_{II} in local schur complement
		ILinearOperatorInverse<vector_type, vector_type>* m_pDirichletSolver;

	//	SchurComplementInverse
		SchurComplementInverse<algebra_type> m_SchurComplementInverse;

	//	Neumann solver for inverse of A_{\{I,\Delta\}, \{I,\Delta\}} in the
	//	creation of the S_{\Pi \Pi} schur complement
		ILinearOperatorInverse<vector_type, vector_type>* m_pNeumannSolver;

	//	Copy of matrix
		PureMatrixOperator<vector_type, vector_type, matrix_type> m_DualDirichletOperator;

	// 	Parallel Neumann Matrix
		matrix_type* m_pDualDirichletMatrix;

	//	Neumann Solver
		ILinearOperatorInverse<vector_type, vector_type>* m_pDualDirichletSolver;

	// 	Convergence Check
		IConvergenceCheck* m_pConvCheck;

	//	Debug Writer
		IDebugWriter<algebra_type>* m_pDebugWriter;

	public:
		void set_domain_decomp_info(pcl::IDomainDecompositionInfo& ddInfo)
		{
			m_pDDInfo = &ddInfo;
		}

	private:
	//	pointer to Domain decomposition info object
		pcl::IDomainDecompositionInfo* m_pDDInfo;

}; /* end class FETISolver */
#endif /* UG_PARALLEL */

} // end namespace ug

#endif /* __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__FETI__ */
