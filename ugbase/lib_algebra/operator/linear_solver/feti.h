/* TODO:
 * If finished: Check for unnecessary "set zero" operations!
 *
 */
/*
 * feti.h
 *
 *  Created on: 21.01.2011
 *      Author: iheppner, avogel
 */

#ifndef __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__FETI__
#define __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__FETI__


#ifdef UG_PARALLEL

#include <iostream>
#include <sstream>
#include <string>
#include "lib_algebra/operator/interface/operator_inverse.h"
#include "lib_algebra/parallelization/parallelization.h"
#include "lib_algebra/operator/debug_writer.h"
#include "pcl/pcl.h"

/* 
 * Outline of FETI-DP algorithm (after Klawonn, A., Widlund, O., Dryja, M.:
 * "Dual-Primal FETI methods for three-dimensional elliptic problems with
 * heterogeneous coefficients", SIAM J. Numer. Anal., 40/2002a, 159--179.):
 *
 * 1. eliminate the primal variables associated with the interior nodes,
 *    the vertex nodes designated as primal (as well as the Lagrange
 *    multipliers related to the optional constraints - if any)
 *    ==> reduced system with Schur complement \tilde{S} related to the subspace
 *        \tilde{W}_{\Delta}, and a reduced right hand side \tilde{f}_{\Delta}.
 *        
 * 2. Reformulate as a variational problem, with a vector of Lagrange multipliers
 *    ==> saddle point problem.
 *
 * 3. Eliminate subvectors u_{\Delta}
 *    ==> system for the dual variables: F \lambda = d
 *
 * 4. Instead of solving ... action of the inverse of the Schur complement can be
 *    obtained by solving the entire linear system from which it originates after
 *    augmenting the given right hand side with zeros.
 *    Group all the interior and dual variables of each subdomain together
 *    and factor the resulting blocks in parallel across the subdomains using a
 *     good ordering algorithm.
 *
 *    The contributions of the remaining Schur complement, of the primal variables,
 *    can also be computed locally prior to subassembly and factorization of this
 *    final, global part of the linear system of equations.
 */

namespace ug{

/// Auxiliary class for handling of "FETI layouts"
/**
 * This class groups all layouts used in FETI method together and provides useful
 * functions on those.
 */
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
			BuildDomainDecompositionLayouts(
							m_masterDualLayout, m_slaveDualLayout,
							m_masterInnerLayout, m_slaveInnerLayout,
							m_masterDualNbrLayout, m_slaveDualNbrLayout,
							m_masterPrimalLayout, m_slavePrimalLayout,
							*m_pMasterStdLayout, *m_pSlaveStdLayout,
							(int)(numIndices - 1), DDInfo);
			UG_LOG("[BuildDomainDecompositionLayouts done]");

		//	create intra feti subdomain communicator
			int localSubdomID = DDInfo.map_proc_id_to_subdomain_id(pcl::GetProcRank());
			pcl::ProcessCommunicator worldComm;
			for(int i = 0; i < DDInfo.get_num_subdomains(); ++i){
				if(localSubdomID == i)
					m_intraFetiSubDomComm = worldComm.create_sub_communicator(true);
				else
					worldComm.create_sub_communicator(false);
			}
			UG_LOG("[intra feti sd comms created]");

		//	create set of unique indices (to avoid doubles due to several interfaces)
		//	collect from dual layouts
			CollectElements(m_vUniqueMasterDualIndex, m_masterDualLayout);
			sort_and_remove_doubles(m_vUniqueMasterDualIndex);

			CollectElements(m_vUniqueSlaveDualIndex, m_slaveDualLayout);
			sort_and_remove_doubles(m_vUniqueSlaveDualIndex);

			CollectElements(m_vUniqueMasterDualNbrIndex, m_masterDualNbrLayout);
			sort_and_remove_doubles(m_vUniqueMasterDualNbrIndex);

			CollectElements(m_vUniqueSlaveDualNbrIndex, m_slaveDualNbrLayout);
			sort_and_remove_doubles(m_vUniqueSlaveDualNbrIndex);
			UG_LOG("[CollectElements for duals done]");

		//	create union of all dual unknowns
			add_merge_sort_remove_doubles(m_vUniqueDualIndex, m_vUniqueMasterDualIndex);
			add_merge_sort_remove_doubles(m_vUniqueDualIndex, m_vUniqueSlaveDualIndex);
			add_merge_sort_remove_doubles(m_vUniqueDualIndex, m_vUniqueMasterDualNbrIndex);
			add_merge_sort_remove_doubles(m_vUniqueDualIndex, m_vUniqueSlaveDualNbrIndex);
			UG_LOG("[union of duals created]");

		//	collect from primal layouts
			CollectElements(m_vUniquePrimalMasterIndex, m_masterPrimalLayout);
			sort_and_remove_doubles(m_vUniquePrimalMasterIndex);

			CollectElements(m_vUniquePrimalSlaveIndex, m_slavePrimalLayout);
			sort_and_remove_doubles(m_vUniquePrimalSlaveIndex);

		//	create union of all primal unknowns
			add_merge_sort_remove_doubles(m_vUniquePrimalIndex, m_vUniquePrimalMasterIndex);
			add_merge_sort_remove_doubles(m_vUniquePrimalIndex, m_vUniquePrimalSlaveIndex);
			UG_LOG("[CollectElements for primals done] ");

		//	test output
			if(bDebug && false)
			{
				pcl::InterfaceCommunicator<IndexLayout> comTmp;
				UG_LOG("STANDARD LAYOUTS:\n");
				PrintLayout(m_stdProcessCom, comTmp, *m_pMasterStdLayout, *m_pSlaveStdLayout);
				UG_LOG("INNER LAYOUTS:\n");
				PrintLayout(m_stdProcessCom, comTmp, m_masterInnerLayout, m_slaveInnerLayout);
				UG_LOG("PRIMAL LAYOUTS:\n");
				PrintLayout(m_stdProcessCom, comTmp, m_masterPrimalLayout, m_slavePrimalLayout);
				UG_LOG("DUAL LAYOUTS:\n");
				PrintLayout(m_stdProcessCom, comTmp, m_masterDualLayout, m_slaveDualLayout);
				UG_LOG("DUAL NBR LAYOUTS:\n");
				PrintLayout(m_stdProcessCom, comTmp, m_masterDualNbrLayout, m_slaveDualNbrLayout);
			}
		}

	//	tests the previously created layouts:
		void test_layouts(bool bPrint)
		{
			// TODO: Besserer Check, dass die Layouts schon existieren? Das hier sollte aber reichen, da Konstruktor
			// 'm_pMasterStdLayout' mit NULL initialisiert!?
			if (m_pMasterStdLayout == NULL) {
				UG_LOG("LAYOUTS not yet created!\n");
				return;
			}
			pcl::InterfaceCommunicator<IndexLayout> comTmp;

			UG_LOG("TEST STANDARD LAYOUTS:\n");
			if (TestLayout(m_stdProcessCom, comTmp, *m_pMasterStdLayout, *m_pSlaveStdLayout, bPrint) != true) {
				UG_LOG("STANDARD LAYOUTS inconsistent!\n");
			} else {
				UG_LOG("STANDARD LAYOUTS are consistent!\n");
			}

			UG_LOG("TEST INNER LAYOUTS:\n");
			if (TestLayout(m_stdProcessCom, comTmp, m_masterInnerLayout, m_slaveInnerLayout, bPrint) != true) {
				UG_LOG("INNER LAYOUTS inconsistent!\n");
			} else {
				UG_LOG("INNER LAYOUTS are consistent!\n");
			}

			UG_LOG("TEST PRIMAL LAYOUTS:\n");
			if (TestLayout(m_stdProcessCom, comTmp, m_masterPrimalLayout, m_slavePrimalLayout, bPrint) != true) {
				UG_LOG("PRIMAL LAYOUTS inconsistent!\n");
			} else {
				UG_LOG("PRIMAL LAYOUTS are consistent!\n");
			}

			UG_LOG("TEST DUAL LAYOUTS:\n");
			if (TestLayout(m_stdProcessCom, comTmp, m_masterDualLayout, m_slaveDualLayout, bPrint) != true) {
				UG_LOG("DUAL LAYOUTS inconsistent!\n");
			} else {
				UG_LOG("DUAL LAYOUTS are consistent!\n");
			}

			UG_LOG("TEST DUAL NBR LAYOUTS:\n");
			if (TestLayout(m_stdProcessCom, comTmp, m_masterDualNbrLayout, m_slaveDualNbrLayout, bPrint) != true) {
				UG_LOG("DUAL NBR LAYOUTS inconsistent!\n");
			} else {
				UG_LOG("DUAL NBR LAYOUTS are consistent!\n");
			}
		}

	//	standard layouts
		IndexLayout& get_std_master_layout() {return *m_pMasterStdLayout;}
		IndexLayout& get_std_slave_layout() {return *m_pSlaveStdLayout;}
		pcl::ProcessCommunicator& get_std_process_communicator() {return m_stdProcessCom;}

	//	intra subdomain layouts
		IndexLayout& get_intra_sd_master_layout() {return m_masterInnerLayout;}
		IndexLayout& get_intra_sd_slave_layout() {return m_slaveInnerLayout;}
		pcl::ProcessCommunicator& get_intra_sd_process_communicator() {return m_intraFetiSubDomComm;}

	//	dual layouts
		IndexLayout& get_dual_master_layout() {return m_masterDualLayout;}
		IndexLayout& get_dual_slave_layout() {return m_slaveDualLayout;}
		const std::vector<IndexLayout::Element>& get_dual_master_indices() const {return m_vUniqueMasterDualIndex;}
		const std::vector<IndexLayout::Element>& get_dual_slave_indices() const {return m_vUniqueSlaveDualIndex;}

	//	dual nbr layouts
		IndexLayout& get_dual_nbr_master_layout() {return m_masterDualNbrLayout;}
		IndexLayout& get_dual_nbr_slave_layout() {return m_slaveDualNbrLayout;}
		const std::vector<IndexLayout::Element>& get_dual_nbr_master_indices() const {return m_vUniqueMasterDualNbrIndex;}
		const std::vector<IndexLayout::Element>& get_dual_nbr_slave_indices() const {return m_vUniqueSlaveDualNbrIndex;}

	//	primal layouts
		IndexLayout& get_primal_master_layout() {return m_masterPrimalLayout;}
		IndexLayout& get_primal_slave_layout() {return m_slavePrimalLayout;}
		const std::vector<IndexLayout::Element>& get_primal_master_indices() const {return m_vUniquePrimalMasterIndex;}
		const std::vector<IndexLayout::Element>& get_primal_slave_indices() const {return m_vUniquePrimalSlaveIndex;}

	//	union of indices
		const std::vector<IndexLayout::Element>& get_primal_indices() const {return m_vUniquePrimalIndex;}
		const std::vector<IndexLayout::Element>& get_dual_indices() const {return m_vUniqueDualIndex;}

	public:

	//	sets standard communication layouts and communicators
		void vec_use_std_communication(vector_type& vec)
		{
			vec.set_slave_layout(get_std_slave_layout());
			vec.set_master_layout(get_std_master_layout());
			vec.set_process_communicator(get_std_process_communicator());
		}

	//	sets intra subdomain communication layouts and communicators
		void vec_use_intra_sd_communication(vector_type& vec)
		{
			vec.set_slave_layout(get_intra_sd_slave_layout());
			vec.set_master_layout(get_intra_sd_master_layout());
			vec.set_process_communicator(get_intra_sd_process_communicator());
		}

	public:
/// calculates vecDest = alpha1*vecSource1 + alpha2*vecSource2 on dual
		void vec_scale_add_on_dual(vector_type& vecDest,
		                           number alpha1, const vector_type& vecSrc1,
		                           number alpha2, const vector_type& vecSrc2)
		{
			VecScaleAdd(vecDest, alpha1, vecSrc1, alpha2, vecSrc2, m_vUniqueDualIndex);
		}

/// calculates vecDest = alpha1*vecSource1 + alpha2*vecSource2 on primal
		void vec_scale_add_on_primal(vector_type& vecDest,
		                           number alpha1, const vector_type& vecSrc1,
		                           number alpha2, const vector_type& vecSrc2)
		{
			VecScaleAdd(vecDest, alpha1, vecSrc1, alpha2, vecSrc2, m_vUniquePrimalIndex);
		}

/// calculates vecInOut += alpha1*vecSrc1 on dual
		void vec_scale_append_on_dual(vector_type& vecInOut,
		                              const vector_type& vecSrc1, number alpha1)
		{
			VecScaleAdd(vecInOut, alpha1, vecSrc1, 1.0, vecInOut, m_vUniqueDualIndex);
		}

/// calculates vecInOut += alpha1*vecSrc1 on primal
		void vec_scale_append_on_primal(vector_type& vecInOut,
		                              const vector_type& vecSrc1, number alpha1)
		{
			VecScaleAdd(vecInOut, alpha1, vecSrc1, 1.0, vecInOut, m_vUniquePrimalIndex);
		}

/// sets vecDest = alpha on dual
		void vec_set_on_dual(vector_type& vecDest, number alpha)
		{
			VecSet(vecDest, alpha, m_vUniqueDualIndex);
		}

/// sets vecDest = alpha on primal
		void vec_set_on_primal(vector_type& vecDest, number alpha)
		{
			VecSet(vecDest, alpha, m_vUniquePrimalIndex);
		}

/// calculates vecDest = alpha*vecSrc on dual
		void vec_scale_assign_on_dual(vector_type& vecDest, const vector_type& vecSrc, number alpha)
		{
			VecScaleAssign(vecDest, alpha, vecSrc, m_vUniqueDualIndex);
		}

/// calculates vecDest = alpha*vecSrc on primal
		void vec_scale_assign_on_primal(vector_type& vecDest, const vector_type& vecSrc, number alpha)
		{
			VecScaleAssign(vecDest, alpha, vecSrc, m_vUniquePrimalIndex);
		}

/// calculates norm on dual
		number vec_norm_on_identified_dual(vector_type& vec)
		{
		//	forward to VecProc
			return sqrt(vec_prod_on_identified_dual(vec, vec));
		}

		number vec_prod_on_identified_dual(const vector_type& vecSrc1, const vector_type& vecSrc2)
		{
		//	compute norm on dual master indices
			double tProdLocal = VecProd(vecSrc1, vecSrc2, m_vUniqueMasterDualIndex);

		//	global value
			double tProdGlobal;

		//	sum up values
			m_stdProcessCom.allreduce(&tProdLocal, &tProdGlobal, 1,
											PCL_DT_DOUBLE, PCL_RO_SUM);

		//	return result
			return tProdGlobal;
		}

	public:
		void mat_set_dirichlet_on_dual(matrix_type& mat)
		{
			SetDirichletRow(mat, m_vUniqueDualIndex);
		}

		void mat_set_dirichlet_on_primal(matrix_type& mat)
		{
			SetDirichletRow(mat, m_vUniquePrimalIndex);
		}

	//	sets standard communication layouts and communicators
		void mat_use_std_communication(matrix_type& mat)
		{
			mat.set_slave_layout(get_std_slave_layout());
			mat.set_master_layout(get_std_master_layout());
			mat.set_process_communicator(get_std_process_communicator());
		}

	//	sets intra subdomain communication layouts and communicators
		void mat_use_intra_sd_communication(matrix_type& mat)
		{
			mat.set_slave_layout(get_intra_sd_slave_layout());
			mat.set_master_layout(get_intra_sd_master_layout());
			mat.set_process_communicator(get_intra_sd_process_communicator());
		}

	protected:
		void sort_and_remove_doubles(std::vector<IndexLayout::Element>& vIndex)
		{
		//	sort vector
			std::sort(vIndex.begin(), vIndex.end());

		//	remove doubles
			vIndex.erase(std::unique(vIndex.begin(), vIndex.end()),
			                         vIndex.end());
		}

		void add_merge_sort_remove_doubles(std::vector<IndexLayout::Element>& vOut,
		                                   const std::vector<IndexLayout::Element>& v)
		{
		//	on entry, vOut and v are supposed to be sorted

		//	original size of vOut
			size_t mid = vOut.size();

		//	add second vector
			vOut.insert(vOut.end(), v.begin(), v.end() );

		//	merge both
			std::inplace_merge(vOut.begin(), vOut.begin()+mid, vOut.end());

		//	remove doubles
			vOut.erase(std::unique(vOut.begin(), vOut.end()),
			                         vOut.end());
		}

	protected:
	//	Standard Layouts
		IndexLayout* m_pMasterStdLayout;
		IndexLayout* m_pSlaveStdLayout;
		pcl::ProcessCommunicator	m_stdProcessCom;

	//	Layouts and Communicator for Inner variables
		IndexLayout m_masterInnerLayout;
		IndexLayout m_slaveInnerLayout;
		pcl::ProcessCommunicator	m_intraFetiSubDomComm;

	//	Layouts for Dual variables
		IndexLayout m_masterDualLayout;
		IndexLayout m_slaveDualLayout;
		std::vector<IndexLayout::Element> m_vUniqueMasterDualIndex;
		std::vector<IndexLayout::Element> m_vUniqueSlaveDualIndex;

	//	Layouts for Dual Neighbour variables
		IndexLayout m_masterDualNbrLayout;
		IndexLayout m_slaveDualNbrLayout;
		std::vector<IndexLayout::Element> m_vUniqueMasterDualNbrIndex;
		std::vector<IndexLayout::Element> m_vUniqueSlaveDualNbrIndex;

	//	Layouts for Primal variables
		IndexLayout m_masterPrimalLayout;
		IndexLayout m_slavePrimalLayout;
		std::vector<IndexLayout::Element> m_vUniquePrimalMasterIndex;
		std::vector<IndexLayout::Element> m_vUniquePrimalSlaveIndex;

	//	unique indices for layouts
		std::vector<IndexLayout::Element> m_vUniqueDualIndex;
		std::vector<IndexLayout::Element> m_vUniquePrimalIndex;

}; /* end class 'FetiLayouts' */

///	Application of the "jump operator" \f$B_{\Delta}\f$
/**
 * This function applies \f$B_{\Delta}\f$ to \f$u_{\Delta}\f$. It computes the
 * difference (the "jump") between the double-valued unknowns \f$u_{\Delta}\f$
 * and communicate the results so that the difference vector is stored consistently
 * afterwards.
 *
 * Note: This computation is only unique up to the sign of the difference.
 * Thus, we can freely decide it, but have then to stay with the choice.
 *
 * \param[out]		diff					destination vector for differences
											computed on "Dual layout"
 * \param[in]		u						vector \f$u_{\Delta}\f$
 * \param[in]		dualMasterLayoutIn		dual master layout to operate on
 * \param[in]		dualSlaveLayoutIn		dual slave  layout to operate on
 * \param[in]		dualNbrMasterLayoutIn	dual master layout to operate on
 * \param[in]		dualNbrSlaveLayoutIn	dual slave  layout to operate on
 */
template <typename TVector>
void ComputeDifferenceOnDelta(TVector& diff, const TVector& u,
							  IndexLayout&    dualMasterLayoutIn,
							  IndexLayout&     dualSlaveLayoutIn,
							  IndexLayout& dualNbrMasterLayoutIn,
							  IndexLayout&  dualNbrSlaveLayoutIn)
{
	// Copy all values
	diff = u;

	// Communicate values:
	// (a) Slaves send their u values, every master subtracts the value of only
	//     one of his slaves ...
	VecSubtractOneSlaveFromMaster(&diff, dualMasterLayoutIn, dualSlaveLayoutIn);

	// (b) masters send their diff values to all slaves (but not vice versa)
	VecCopy(&diff, dualMasterLayoutIn, dualSlaveLayoutIn);

	// (c) and also to the additional slaves living in "Dual neighbour layout"
	VecCopy(&diff, dualNbrMasterLayoutIn, dualNbrSlaveLayoutIn);

	// now the vector 'diff' should be consistent!
	return;
};

/// 'ComputeDifferenceOnDeltaTransposed()': Apply \f$B_{\Delta}^T\f$
/**
 * This function applies \f$B_{\Delta}^T\f$ to a difference vector \f$d\f$,
 * (lying in the same space as \f$\lambda\f$).
 * \f$d\f$ is supposed to be stored consistently.
 *
 * For the application of \f$B_{\Delta}\f$ and the chosen scaling factors see
 * the documentation of 'ComputeDifferenceOnDelta()'.
 *
 * \param[out]		f						vector \f$f_{\Delta}\f$ living on "Dual layout"
 * \param[in]		diff					difference vector on "Dual layout"
 * \param[in]		vDualMasterIndex		dual  master   layout to operate on
 * \param[in]		vDualSlaveIndex		    dual slave     layout to operate on
 * \param[in]		vDualNbrSlaveIndex	    dual nbr slave layout to operate on
 */
template <typename TVector>
void ComputeDifferenceOnDeltaTransposed(TVector& f, const TVector& diff,
										const std::vector<IndexLayout::Element>& vDualMasterIndex,
										const std::vector<IndexLayout::Element>& vDualSlaveIndex,
										const std::vector<IndexLayout::Element>& vDualNbrSlaveIndex)
{
	// Copy values (no communication is performed):
	// (a) set f = +1 * d on masters ...
	VecScaleAssign(f, 1.0, diff, vDualMasterIndex);

	// (b) set f = -1 * d on slaves
	VecScaleAssign(f, -1.0, diff, vDualSlaveIndex);

	// (c) set f = +1 * d on slaves living in "Dual neighbour layout"
	//     ("+1" since in this context these dofs play the role of masters!)
	VecScaleAssign(f, +1.0, diff, vDualNbrSlaveIndex);
	return;

};

/// operator implementation of the local Schur complement
/**
 * This operator is the application of the local Schur complement \f$S_{\Delta}^{i}\f$.
 * The underlying matrix must have at least two layouts. The first layout, layout level 0,
 * will be used to describe subdomain *internal* interfaces (i.e. "pure" processor interfaces),
 * all other layouts are used to identify the boundary,
 * layout level 1: \f$ \Delta \f$ (edges of subdomains),
 * layout level 2: \f$ \Pi \f$ (vertices of subdomains, a.k.a. "cross points") - will be constructed here,
 * and the Schur complement is build w.r.t. to these variables.
 */
template <typename TAlgebra>
class LocalSchurComplement
	: public ILinearOperator<	typename TAlgebra::vector_type,
	  	  	  	  	  	  	  	typename TAlgebra::vector_type>,
	  public DebugWritingObject<TAlgebra>
{
	public:
	// 	Algebra type
		typedef TAlgebra algebra_type;

	// 	Vector type
		typedef typename TAlgebra::vector_type vector_type;

	// 	Matrix type
		typedef typename TAlgebra::matrix_type matrix_type;

	protected:
		using DebugWritingObject<TAlgebra>::write_debug;
		using DebugWritingObject<TAlgebra>::debug_writer;

	public:
	///	constructor
		LocalSchurComplement();

	///	name of solver
		virtual const char* name() const {return "Local Schur Complement Solver";}

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
		void set_matrix(SmartPtr<MatrixOperator<matrix_type, vector_type> > A)
		{
		//	save current operator
			m_spOperator = A;
		}

	///	sets a sequential Dirichlet solver
	/**
	 * This method sets the Dirichlet Solver that is used to invert the
	 * inner matrix \f$A_{II}\f$
	 */
		void set_dirichlet_solver(SmartPtr<ILinearOperatorInverse<vector_type> > dirichletSolver)
		{
		//	remember the Dirichlet Solver
			m_spDirichletSolver = dirichletSolver;
		}

	///	sets the primal layouts
		void set_feti_layouts(FetiLayouts<algebra_type>& fetiLayouts)
		{
			m_pFetiLayouts = &fetiLayouts;
		}

	/// implementation of the operator for the solution dependent initialization.
		void init(const vector_type& u) {init();}

	///	initializes the solver for operator A
	/**
	 * This method must be called, before the apply() method can be invoked.
	 * It has to be called each time, when the matrix has been replaced. A deep
	 * copy of the matrix is then constructed and in this copy the rows belonging
	 * to the \f$\Delta\f$ and \f$\Pi\f$ unknowns are set to identity rows. This
	 * matrix is used in the solution of the local dirichlet problem.
	 */
		virtual void init();

	///	applies the Schur complement built from matrix operator set via 'set_matrix()'
	/// to 'u' and returns the result 'f := S times u'
		virtual void apply(vector_type& f, const vector_type& u);

	///	solves the system
		virtual void apply_sub(vector_type& f, const vector_type& u);

	///	sets statistic slot where next iterate should be counted
		void set_statistic_type(std::string type) {m_statType = type;}

	///	prints some convergence statistic of inner solvers
		void print_statistic_of_inner_solver(bool bPrintOnlyAverages); // const;

		void clear_total_itercnt_of_inner_solvers() {m_totalIterCntOfInnerSolvers = 0;}
		int get_total_itercnt_of_inner_solvers() {return m_totalIterCntOfInnerSolvers;}

	// destructor
		virtual ~LocalSchurComplement() {};

	protected:
	// 	Operator that is inverted by this Inverse Operator
		SmartPtr<MatrixOperator<matrix_type,vector_type> > m_spOperator;

	// 	Parallel Matrix
		matrix_type* m_pMatrix;

	//	Feti Layouts
		FetiLayouts<algebra_type>* m_pFetiLayouts;

	//	Copy of matrix
		SmartPtr<MatrixOperator<matrix_type, vector_type> > m_spDirichletOperator;

	// 	Parallel Dirichlet Matrix
		matrix_type* m_pDirichletMatrix;

	// 	Linear Solver to invert the local Dirichlet problems
		SmartPtr<ILinearOperatorInverse<vector_type> > m_spDirichletSolver;

	//	Convergence history
		std::string m_statType;

		struct StepConv
		{
			int numIter3b;
			number lastDef3b;
		};

		std::map<std::string, std::vector<StepConv> > m_mvStepConv;

		int m_applyCnt;

		int m_totalIterCntOfInnerSolvers;

}; /* end class 'LocalSchurComplement' */

/* 1.7 Application of \f${\tilde{S}_{\Delta \Delta}}^{-1}\f$ */ 
/// operator implementation of the inverse of the Schur complement w.r.t. the "Delta unknowns"
/**
 * This operator provides the application of the inverse of the Schur complement w.r.t. the "Delta unknowns",
 * \f${\tilde{S}_{\Delta \Delta}}^{-1}\f$.
 */
template <typename TAlgebra>
class PrimalSubassembledMatrixInverse
	: public ILinearOperatorInverse<	typename TAlgebra::vector_type,
	  	  	  	  	  	  	  			typename TAlgebra::vector_type>,
	  public DebugWritingObject<TAlgebra>
{
	public:
	// 	Algebra type
		typedef TAlgebra algebra_type;

	// 	Vector type
		typedef typename TAlgebra::vector_type vector_type;

	// 	Matrix type
		typedef typename TAlgebra::matrix_type matrix_type;

	///	Base type
		typedef ILinearOperatorInverse<vector_type> base_type;

	protected:
		using base_type::convergence_check;
		using DebugWritingObject<TAlgebra>::write_debug;
		using DebugWritingObject<TAlgebra>::debug_writer;

	public:
	///	constructor
		PrimalSubassembledMatrixInverse();

	///	name of class
		virtual const char* name() const {return "Schur Complement Inverse";}

	///	sets the Neumann solver
		void set_neumann_solver(SmartPtr<ILinearOperatorInverse<vector_type> > neumannSolver)
		{
		//	remember the Dirichlet Solver
			m_spNeumannSolver = neumannSolver;
		}

	///	sets the coarse problem solver
		void set_coarse_problem_solver(SmartPtr<ILinearOperatorInverse<vector_type> > coarseProblemSolver)
		{
		//	remember the coarse problem Solver
			m_spCoarseProblemSolver = coarseProblemSolver;
		}

	///	sets the primal layouts
		void set_feti_layouts(FetiLayouts<algebra_type>& fetiLayouts)
		{
			m_pFetiLayouts = &fetiLayouts;
		}

	// 	Init for Linear Operator L
		virtual bool init(SmartPtr<ILinearOperator<vector_type> > L);


	// 	Init for Linear Operator J and Linearization point (current solution)
		virtual bool init(SmartPtr<ILinearOperator<vector_type> > J, const vector_type& u)
		{
			return init(J);
		}

	// 	Solve A*u = f, such that u = A^{-1} f
		virtual bool apply(vector_type& u, const vector_type& f);

	// 	Solve A*u = f, such that u = A^{-1} f
	// 	This is done by iterating: u := u + B(f - A*u)
	// 	In f the last defect f := f - A*u is returned
		virtual bool apply_return_defect(vector_type& u, vector_type& f);

	///	sets statistic slot where next iterate should be counted
		void set_statistic_type(std::string type) {m_statType = type;}

	///	prints some convergence statistic of inner solvers
		void print_statistic_of_inner_solver(bool bPrintOnlyAverages); // const;

		void clear_total_itercnt_of_inner_solvers() {m_totalIterCntOfInnerSolvers = 0;}
		int get_total_itercnt_of_inner_solvers() {return m_totalIterCntOfInnerSolvers;}

	///	set 'm_bTestOneToManyLayouts'
		void set_test_one_to_many_layouts(bool bTest) {m_bTestOneToManyLayouts = bTest;}

	//  destructor
		virtual ~PrimalSubassembledMatrixInverse() {};

	protected:
	// 	Operator that is inverted by this Inverse Operator ==> from which SC is built (05022011)
		SmartPtr<MatrixOperator<matrix_type,vector_type> > m_spOperator;

	// 	Parallel Matrix to invert ==> from which SC is built (05022011)
		matrix_type* m_pMatrix;

	//	Feti Layouts
		FetiLayouts<algebra_type>* m_pFetiLayouts;

	//	Copy of matrix
		SmartPtr<MatrixOperator<matrix_type, vector_type> > m_spNeumannOperator;

	// 	Parallel Neumann Matrix
		matrix_type* m_pNeumannMatrix;

	//	Neumann Solver
		SmartPtr<ILinearOperatorInverse<vector_type> > m_spNeumannSolver;

	//	Coarse problem Solver on root.
		SmartPtr<ILinearOperatorInverse<vector_type> > m_spCoarseProblemSolver;

	//	Process gathering the Schur Complement w.r.t. Primal unknowns
		int m_primalRootProc;

	//	Layouts for all to one communication
		IndexLayout m_masterAllToOneLayout;
		IndexLayout m_slaveAllToOneLayout;
		pcl::ProcessCommunicator m_allToOneProcessComm;

	//	Schur Complement operator for gathered matrix
		SmartPtr<MatrixOperator<matrix_type, vector_type> > m_spRootSchurComplementOp;

	//	Matrix for one proc Schur complement
		matrix_type* m_pRootSchurComplementMatrix;

	//	Convergence history
		std::string m_statType;

		struct StepConv
		{
			int numIter2a;
			int numIter7;
			int numIterSC;
			number lastDef2a;
			number lastDef7;
			number lastDefSC;
		};

		std::map<std::string, std::vector<StepConv> > m_mvStepConv;

	//	testing of 'one-to-many layouts'
		bool m_bTestOneToManyLayouts;

		int m_totalIterCntOfInnerSolvers;

}; /* end class 'PrimalSubassembledMatrixInverse' */

/// operator implementation of the FETI-DP solver
/**
 * This operator implements a FETI-DP solver, see e.g.
 * "Domain Decomposition Methods -- Algorithms and Theory",
 * A. Toselli, O. Widlund, Springer 2004, sec. 1.3.5, p. 12ff.
 */
template <typename TAlgebra>
class FETISolver : public IMatrixOperatorInverse<	typename TAlgebra::matrix_type,
													typename TAlgebra::vector_type>,
	public DebugWritingObject<TAlgebra>
{
	public:
	// 	Algebra type
		typedef TAlgebra algebra_type;

	// 	Vector type
		typedef typename TAlgebra::vector_type vector_type;

	// 	Matrix type
		typedef typename TAlgebra::matrix_type matrix_type;

	///	Base type
		typedef ILinearOperatorInverse<vector_type> base_type;

	protected:
		using base_type::convergence_check;
		using DebugWritingObject<TAlgebra>::write_debug;
		using DebugWritingObject<TAlgebra>::debug_writer;

	public:
	///	constructor
		FETISolver();

	///	name of solver
		virtual const char* name() const {return "FETI Solver";}

	///	sets the Dirichlet solver
		void set_dirichlet_solver(SmartPtr<ILinearOperatorInverse<vector_type> > dirichletSolver)
		{
		//	remember the Dirichlet Solver
			m_spDirichletSolver = dirichletSolver;
		}

	///	sets the Neumann solver
		void set_neumann_solver(SmartPtr<ILinearOperatorInverse<vector_type> > neumannSolver)
		{
		//	remember the Dirichlet Solver
			m_spNeumannSolver = neumannSolver;
		}

	///	sets the coarse problem solver
		void set_coarse_problem_solver(SmartPtr<ILinearOperatorInverse<vector_type> > coarseProblemSolver)
		{
		//	remember the coarse problem Solver
			m_spCoarseProblemSolver = coarseProblemSolver;
		}

	//	set debug output
		void set_debug(SmartPtr<IDebugWriter<algebra_type> > spDebugWriter)
		{
			m_LocalSchurComplement.set_debug(spDebugWriter);
			m_PrimalSubassembledMatrixInverse.set_debug(spDebugWriter);
			VectorDebugWritingObject<vector_type>::set_debug(spDebugWriter);
		}

	///	initializes the solver for operator A
		virtual bool init(SmartPtr<MatrixOperator<matrix_type, vector_type> > A);

	///	function which applies matrix \f$F\f$ of the reduced system ("Delta system") to a vector \f$v\f$
	/**
	 * This function applies matrix \f$F := B_{\Delta} \tilde{S}^{-1} B_{\Delta}^T\f$
	 * to a vector \f$v\f$. \f$v\f$ can be:
	 * (a) unknown vector of Lagrange multipliers, lambda, and
	 * (b) search direction 'p' in cg method.
	 *
	 * \param[in]		v				vector \f$v\f$ living on "Dual layout"
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

	///	function which applies diagonal scaling matrix \f$D_{\Delta}^{(i)}\f$ to a vector \f$v\f$
		bool apply_scaling_matrix(vector_type& s, const vector_type& v) // maybe restrict to layout
		{
		// 	copy values
			s = v;

		// 	scale by 1/2
			s *= 1./2.;

			// \todo: more general: m_Ddelta.apply(s,v), with additional member 'matrix_type m_Ddelta;'

		//	we're done
			return true;
		}

	///	function which applies matrix \f$ M^{-1}^{(i)} \f$ to a vector \f$ r \f$
	/**
	 * This function applies matrix \f$M^{-1}^{(i)} := D_{\Delta}^{(i)} B_{\Delta}^{(i)} S_{\Delta}^{(i)} {B_{\Delta}^{(i)}}^T D_{\Delta}^{(i)}\f$
	 * to a vector \f$r\f$.
	 *
	 * \param[in]		r				vector \f$r\f$ living on "Dual layout"
	 * \param[out]		z				result of application of \f$ M^{-1} \f$
	 */
		bool apply_M_inverse(vector_type& z, const vector_type& r);

	//	tests layouts:
		void test_layouts(bool bPrint);

	//	set member 'm_bTestOneToManyLayouts' of class 'PrimalSubassembledMatrixInverse':
		void set_test_one_to_many_layouts(bool bTest)
		{
			m_PrimalSubassembledMatrixInverse.set_test_one_to_many_layouts(bTest);
		}


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

	///	prints some convergence statistic of inner solvers
		void print_statistic_of_inner_solver(bool bPrintOnlyAverages) //const
		{
			// sum over all calls of inner solvers is done in the following print functions!
			m_PrimalSubassembledMatrixInverse.clear_total_itercnt_of_inner_solvers();
			m_PrimalSubassembledMatrixInverse.print_statistic_of_inner_solver(bPrintOnlyAverages);
			UG_LOG("\n");

			m_LocalSchurComplement.clear_total_itercnt_of_inner_solvers();
			m_LocalSchurComplement.print_statistic_of_inner_solver(bPrintOnlyAverages);

			int totalIterCntOfInnerSolvers = m_PrimalSubassembledMatrixInverse.get_total_itercnt_of_inner_solvers();
			totalIterCntOfInnerSolvers     += m_LocalSchurComplement.get_total_itercnt_of_inner_solvers();

			UG_LOG("Total number of calls of sub problem solvers: " << totalIterCntOfInnerSolvers);
			UG_LOG(" in " << std::setw(5) << m_iterCnt << " FETI iterations ");
			UG_LOG(" on " << std::setw(5) << m_pDDInfo->get_num_subdomains() << " FETI subdomains, "
				          << std::fixed << (double)totalIterCntOfInnerSolvers/m_iterCnt
				          << std::scientific << " per FETI iteration, "
				          << std::fixed << (double)totalIterCntOfInnerSolvers/m_pDDInfo->get_num_subdomains()
				          << std::scientific << " per FETI subdomain.\n");

		}

		// destructor
		virtual ~FETISolver() {};

	protected:
	//	Prepare the convergence check
		void prepare_conv_check()
		{
			convergence_check()->set_name(name());
			convergence_check()->set_symbol('%');

			//	set preconditioner string
			std::stringstream ss;
			ss << " (Inherent Preconditioner) ";
			convergence_check()->set_info(ss.str());
		}

	protected:
		virtual void write_debug(const vector_type& vec, const char* filename)
		{
		//	add iter count to name
			std::string name(filename);
			char ext[20]; sprintf(ext, "_iter%03d.vec", m_iterCnt);
			name.append(ext);

		//	if no debug writer set, we're done
			if(debug_writer() == NULL) return;

		//	write
			debug_writer()->write_vector(vec, name.c_str());
		}

		int m_iterCnt;

	protected:
	// 	Operator that is inverted by this Inverse Operator
		SmartPtr<MatrixOperator<matrix_type,vector_type> > m_spOperator;

	// 	Parallel Matrix to invert
		matrix_type* m_pMatrix;

	//	Feti Layouts
		FetiLayouts<algebra_type> m_fetiLayouts;

	//	Local Schur complement for each feti subdomain
		LocalSchurComplement<algebra_type> m_LocalSchurComplement;

	//	Dirichlet solver for inverse of \f$A_{II}\f$ in local Schur complement
		SmartPtr<ILinearOperatorInverse<vector_type> > m_spDirichletSolver;

	//	PrimalSubassembledMatrixInverse
		PrimalSubassembledMatrixInverse<algebra_type> m_PrimalSubassembledMatrixInverse;

	//	Neumann solver for inverse of \f$A_{\{I,\Delta\}, \{I,\Delta\}}\f$ in the
	//	creation of the S_{\Pi \Pi} Schur complement in PrimalSubassembledMatrixInverse
		SmartPtr<ILinearOperatorInverse<vector_type> > m_spNeumannSolver;

	//	Solver used in solving coarse problem on root.
	// 	It solves \f$S_{\Pi \Pi} u_{\Pi} = \tilde{f}_{\Pi}\f$ 
		SmartPtr<ILinearOperatorInverse<vector_type> > m_spCoarseProblemSolver;

	public:
		void set_domain_decomp_info(pcl::IDomainDecompositionInfo& ddInfo)
		{
			m_pDDInfo = &ddInfo;
		}

	private:
	//	pointer to Domain decomposition info object
		pcl::IDomainDecompositionInfo* m_pDDInfo;

}; /* end class 'FETISolver' */

} // end namespace ug

#endif /* UG_PARALLEL */

#endif /* __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__FETI__ */
