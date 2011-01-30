/* TODO:
 * Solve Schur complement w.r.t. "Pi system"
 * If finished: Check for unnecessary "set zero" operations!
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

///	algebra blocks are static.
template <> struct block_traits<int>
{
	enum{
		is_static = 1
	};
};

///	Extracts cross points from Delta-Layout
/**
 * This function extracts the cross points from a "Delta layout" and creates a new "Pi layout"
 * containing only these cross points which on the other hand are eliminated from the "Delta layout"
 * so that both index sets are disjunct.
 *
 * \param[in]		numIDs			number of dof's of actual processor
 * \param[in]		masterLayoutIn	master layout to extract
 * \param[in]		slaveLayoutIn	slave layout to extract
 * \param[in]		ddInfo			domain decomposition info
 * \param[out]		masterLayoutOut	master layout created
 * \param[out]		slaveLayoutOut	slave layout created
 */
inline void ExtractCrossPointLayouts(size_t numIDs,
									 IndexLayout& masterLayoutIn,
									 IndexLayout& slaveLayoutIn,
									 pcl::IDomainDecompositionInfo* ddInfoIn,
									 IndexLayout& masterCPLayoutOut, // CP: "Cross Point"
									 IndexLayout& slaveCPLayoutOut
	)
{
	int num_pi_dofs = 0;
	int num_pi_interfaces = 0;

	std::vector<int> vMultiplicity;
//	generate an id for each entry.
	vMultiplicity.clear();
	vMultiplicity.resize(numIDs, -1);

	int localProc   = pcl::GetProcRank();
	int localSubdom = ddInfoIn->map_proc_id_to_subdomain_id(localProc);

//	Add 1 for all master layouts an index is contained in
	for(IndexLayout::iterator interface_iter = masterLayoutIn.begin();
			interface_iter != masterLayoutIn.end(); ++interface_iter)
	{
		int targetProc   = masterLayoutIn.proc_id(interface_iter);
		int targetSubdom = ddInfoIn->map_proc_id_to_subdomain_id(targetProc);

	//	skip interface to own subdomain
		if(localSubdom == targetSubdom) continue;

	//	get interface
		IndexLayout::Interface& interface = masterLayoutIn.interface(interface_iter);
		
	//	loop over indices
		for( IndexLayout::Interface::iterator iter = interface.begin();
				iter != interface.end(); ++iter)
		{
		//  get index
			const size_t index = interface.get_element(iter);

		//	set value of vector to zero
			if(vMultiplicity[index] == -1)
				vMultiplicity[index] = targetSubdom;
			else if(vMultiplicity[index] != targetSubdom)
				vMultiplicity[index] = -2;
		}
	}

//	copy all ids from master to slave interfaces
	ComPol_VecCopy<std::vector<int> >	copyPol(&vMultiplicity);

	pcl::ParallelCommunicator<IndexLayout> communicator;
	communicator.send_data(masterLayoutIn, copyPol);
	communicator.receive_data(slaveLayoutIn, copyPol);
	communicator.communicate();

//	Select cross points for all master layouts an index is contained in
	for(IndexLayout::iterator interface_iter = masterLayoutIn.begin();
			interface_iter !=  masterLayoutIn.end(); ++interface_iter)
	{
	//	get neighbour proc id
		int targetProc = masterLayoutIn.proc_id(interface_iter);

	//	get interface
		IndexLayout::Interface& interface = masterLayoutIn.interface(interface_iter);

	//	loop over indices
		for(IndexLayout::Interface::iterator iter = interface.begin();
				iter != interface.end(); )
		{
		//  get index
			const size_t index = interface.get_element(iter);

			if(vMultiplicity[index] == -2)
			{
			//	create a cross point interface
				IndexLayout::Interface& indexInterface = masterCPLayoutOut.interface(targetProc);

				indexInterface.push_back(index);
				iter = interface.erase(iter);
				num_pi_dofs++;
				num_pi_interfaces++;
			}
			else
			{
				++iter;
			}
		}
	}

//UG_LOG_ALL_PROCS("'ExtractCrossPointLayouts()': Proc " << pcl::GetProcRank() <<  " has " << num_pi_dofs << " in layout 2, contained in " << num_pi_interfaces << " interfaces!\n");

//	Select cross points for all slave layouts an index is contained in
	for(IndexLayout::iterator interface_iter = slaveLayoutIn.begin();
			interface_iter !=  slaveLayoutIn.end(); ++interface_iter)
	{
	//	get neighbour proc id
		int targetProc = slaveLayoutIn.proc_id(interface_iter);

	//	get interface
		IndexLayout::Interface& interface = slaveLayoutIn.interface(interface_iter);

	//	loop over indices
		for(IndexLayout::Interface::iterator iter = interface.begin();
				iter != interface.end(); )
		{
		//  get index
			const size_t index = interface.get_element(iter);

			if(vMultiplicity[index] == -2)
			{
			//	create a cross point interface
				IndexLayout::Interface& indexInterface = slaveCPLayoutOut.interface(targetProc);

				indexInterface.push_back(index);
				iter = interface.erase(iter);
			}
			else
			{
				++iter;
			}
		}
	}
	
	return;
}

///	Application of the "jump operator" \f$B_{\Delta}\f$:
/// 'ComputeDifferenceOnDelta()': Apply \f$B_{\Delta}\f$ to \f$u_{\Delta}\f$
/**
 * \f$B_{\Delta}\f$ computes the difference between the double-valued unknowns \f$u_{\Delta}\f$.
 * This computation is only unique up to the sign of the difference. Thus, we can
 * freely decide it, but have then to stay with the choice.
 *
 * \param[out]		diff			destination vector for computed differences on "Delta layout"
 * \param[in]		u				vector \f$u_{\Delta}\f$
 * \param[in]		masterLayoutIn	master layout to operate on (caller has to provided a "Delta layout")
 * \param[in]		slaveLayoutIn	slave  layout to operate on (caller has to provided a "Delta layout")
 */
template <typename TVector>
void ComputeDifferenceOnDelta(TVector& diff, const TVector& u,
							  IndexLayout& masterLayoutIn,
							  IndexLayout& slaveLayoutIn)
{
	// Reset all values
	diff.set(0.0);

	// All masters subtract the values of the slave, all slaves subtract the values
	// of the master (communication is performed) ...
	VecSubtractOnLayout(&diff, masterLayoutIn, slaveLayoutIn);

	// ... and slaves multiplies the result by '-1' (no communication is performed)
	VecSetOnLayout(&diff, -1.0, slaveLayoutIn);

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
 * \param[in]		masterLayoutIn	master layout to operate on (caller has to provided a "Delta layout")
 * \param[in]		slaveLayoutIn	slave  layout to operate on (caller has to provided a "Delta layout")
 */
template <typename TVector>
void ComputeDifferenceOnDeltaTransposed(TVector& f, const TVector& lambda,
							  IndexLayout& masterLayoutIn,
							  IndexLayout& slaveLayoutIn)
{
	
	// 1. All masters set their values equal to $\lambda$
	// 2. All slaves set their values equal to $-\lambda$

	// (a) Reset all values
	f.set(0.0);
	//VecSetOnLayout(&f, 0.0, masterLayoutIn);
	//VecSetOnLayout(&f, 0.0, slaveLayoutIn);

	// (b) Copy values on \Delta
	VecScaleAddOnLayout(&f, &lambda, 1.0, masterLayoutIn);
	VecScaleAddOnLayout(&f, &lambda, -1.0, slaveLayoutIn);

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
		LocalSchurComplement()
		{}

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

	///	sets the cross point layouts
		void set_crosspoint_layouts(IndexLayout& slaveLayout, IndexLayout& masterLayout)
		{
			m_pSlaveCPLayout = &slaveLayout;
			m_pMasterCPLayout = &masterLayout;
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
		virtual bool init()
		{
		//	check that operator has been set
			if(m_pOperator == NULL)
			{
				UG_LOG("ERROR in 'LocalSchurComplement::init': No Operator A"
						" set.\n");
				return false;
			}

		//	check that Pi layouts have been set
			if(m_pSlaveCPLayout == NULL || m_pMasterCPLayout == NULL)
			{
				UG_LOG("ERROR in 'LocalSchurComplement::init': Master or Slave"
						" layout for cross points not set.\n");
				return false;
			}

		//	save matrix from which we build the Schur complement
			m_pMatrix = &m_pOperator->get_matrix();

			/* no longer necessary to check here - we have already pointers to "Pi layouts"
			   (set via 'set_crosspoint_layouts()')!:
		//	check that matrix has enough decomposition levels
			if(m_pMatrix->num_layouts() != 2)
			{
				UG_LOG("ERROR in 'LocalSchurComplement::init': The Operator must"
						" have at two layouts, but the current Operator has"
						" only " << m_pMatrix->num_layouts() << "\n");
				return false;
			}
			*/

		//	get matrix from dirichlet operator
			m_pDirichletMatrix = &m_DirichletOperator.get_matrix();

		//	Copy Matrix for Dirichlet Problem
			*m_pDirichletMatrix = *m_pMatrix;

		//	Set Dirichlet values on Pi
			MatSetDirichletOnLayout(m_pDirichletMatrix, *m_pSlaveCPLayout);
			MatSetDirichletOnLayout(m_pDirichletMatrix, *m_pMasterCPLayout);

		//	Set Dirichlet values on Delta
			MatSetDirichletOnLayout(m_pDirichletMatrix, m_pDirichletMatrix->get_slave_layout(1));
			MatSetDirichletOnLayout(m_pDirichletMatrix, m_pDirichletMatrix->get_master_layout(1));

		//	init sequential solver for Dirichlet problem
			if(m_pDirichletSolver != NULL)
				if(!m_pDirichletSolver->init(m_DirichletOperator))
				{
					UG_LOG("ERROR in 'LocalSchurComplement::init': Cannot init "
							"Sequential Dirichlet Solver for Operator A.\n");return false;
				}

		//	Debug output of matrices
			if(m_pDebugWriter != NULL)
			{
				m_pDebugWriter->write_matrix(m_DirichletOperator.get_matrix(),
				                             "FetiDirichletMatrix");
				m_pDebugWriter->write_matrix(m_pOperator->get_matrix(),
				                             "FetiOriginalMatrix");
			}

		//	we're done
			return true;
		} /* end 'LocalSchurComplement::init()' */

	///	applies the Schur complement built from matrix operator set via 'set_matrix()'
	/// to 'u' and returns the result 'f := S times u'
		virtual bool apply(vector_type& f, const vector_type& u)
		{
		//	check that matrix has been set
			if(m_pOperator == NULL)
			{
				UG_LOG("ERROR: In 'LocalSchurComplement::apply': "
								"Matrix A not set.\n");
				return false;
			}

		//	check dirichlet solver
			if(m_pDirichletSolver == NULL)
			{
				UG_LOG("ERROR: In 'LocalSchurComplement::apply':"
								" No sequential Dirichlet Solver set.\n");
				return false;
			}

		//	Check parallel storage type of matrix
			if(!m_pDirichletMatrix->has_storage_type(PST_ADDITIVE))
			{
				UG_LOG("ERROR: In 'LocalSchurComplement::apply': "
								"Inadequate storage format of matrix.\n");
				return false;
			}

		//	Check parallel storage type of vectors
			if (!u.has_storage_type(PST_CONSISTENT))
			{
				UG_LOG("ERROR: In 'LocalSchurComplement::apply': "
								"Inadequate storage format of Vector 'u' (should be consistent).\n");
				return false;
			}
			if(!f.has_storage_type(PST_ADDITIVE))
			{
				UG_LOG("ERROR: In 'LocalSchurComplement::apply': "
								"Inadequate storage format of Vector 'f' (should be additive).\n");
				return false;
			}

		//	Help vector
		//	\todo: it would be sufficient to copy only the layouts without copying the values
			vector_type uTmp; uTmp.create(u.size()); uTmp = u;

		//	1. Set values to zero on I and Pi
			// (a) Reset all values
			uTmp.set(0.0);

			// (b) Copy values on \Delta
			VecScaleAddOnLayout(&uTmp, &u, 1.0, u.get_slave_layout(1));
			VecScaleAddOnLayout(&uTmp, &u, 1.0, u.get_master_layout(1));

		//	2. Compute rhs f_{I} = A_{I \Delta} u_{\Delta}
			if(!m_DirichletOperator.apply(f, uTmp))
			{
				UG_LOG_ALL_PROCS("ERROR in 'LocalSchurComplement::apply': "
								 "Could not compute Rhs for Dirichlet problem on "
								 "Proc " << pcl::GetProcRank() << ".\n");
				return false;
			}
			// set values to zero on \Delta
			VecSetOnLayout(&f, 0.0, u.get_slave_layout(1));
			VecSetOnLayout(&f, 0.0, u.get_master_layout(1));

		//	3. Invert on inner unknowns u_{I} = A_{II}^{-1} f_{I}
			// (a) use the inner-FETI-block layouts
			uTmp.use_layout(0);
			f.use_layout(0);
			m_pDirichletMatrix->use_layout(0);

			// (b) invoke Dirichlet solver
			if(!m_pDirichletSolver->apply_return_defect(uTmp, f))
			{
				UG_LOG_ALL_PROCS("ERROR in 'LocalSchurComplement::apply': "
								 "Could not solve Dirichlet problem on Proc "
									<< pcl::GetProcRank() << ".\n");
				return false;
			}

		//	4. Compute result vector
			// (a) Scale u_{I} by -1
			uTmp *= -1.0;

			// (b) Add u_{\Delta} on \Delta
			VecScaleAddOnLayout(&uTmp, &u, 1.0, u.get_slave_layout(1));
			VecScaleAddOnLayout(&uTmp, &u, 1.0, u.get_master_layout(1));

			// (c) Multiply with full matrix
			if(!m_pOperator->apply(f, uTmp))
			{
				UG_LOG_ALL_PROCS("ERROR in 'LocalSchurComplement::apply': "
								 "Could not apply full matrix on "
								 "Proc " << pcl::GetProcRank() << ".\n");
				return false;
			}

		//	5. Reset all values for I, \Pi
			VecSetExcludingLayout(&f, 0.0, u.get_slave_layout(1));

			VecSetOnLayout(&f, 0.0, *m_pSlaveCPLayout);
			VecSetOnLayout(&f, 0.0, *m_pMasterCPLayout);

		//	we're done
			return true;
		} /* end 'LocalSchurComplement::apply()' */

	///	solves the system
		virtual bool apply_sub(vector_type& f, const vector_type& u)
		{
		//	create new rhs
			vector_type d; d.resize(f.size());

		//	solve
			if(!apply(d, u)) return false;

		//	subtract from vector
			f -= d;

		//	we're done
			return true;
		}

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

	//	Layouts for Pi
		IndexLayout* m_pMasterCPLayout;
		IndexLayout* m_pSlaveCPLayout;

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
		SchurComplementInverse()
		{}

	///	name of class
		virtual const char* name() const {return "Schur Complement Inverse";}

	///	sets the Neumann solver
		void set_neumann_solver(ILinearOperatorInverse<vector_type, vector_type>& neumannSolver)
		{
		//	remember the Dirichlet Solver
			m_pNeumannSolver = &neumannSolver;
		}

	///	sets the cross point layouts
		void set_crosspoint_layouts(IndexLayout& slaveLayout, IndexLayout& masterLayout)
		{
			m_pSlaveCPLayout = &slaveLayout;
			m_pMasterCPLayout = &masterLayout;
		}

	// 	Init for Linear Operator L
		virtual bool init(ILinearOperator<vector_type, vector_type>& L)
		{
		//	remember operator
			m_A = dynamic_cast<IMatrixOperator<vector_type, vector_type, matrix_type>*>(&L);

		//	check, that operator is correct
			if(m_A == NULL)
			{
				UG_LOG("ERROR in 'SchurComplementInverse::init':"
						" Wrong type of operator passed for init.\n");
				return false;
			}

		//	check that Pi layouts have been set
			if(m_pSlaveCPLayout == NULL || m_pMasterCPLayout == NULL)
			{
				UG_LOG("ERROR in 'SchurComplementInverse::init': Master or Slave"
						" layout for cross points not set.\n");
				return false;
			}

		//	save matrix from which we build the Schur complement
			m_pMatrix = &m_A->get_matrix();

			/* no longer necessary to check here - we have already pointers to "Pi layouts"!:
			   (set via 'set_crosspoint_layouts()')!:
		//	check that matrix has enough decomposition levels
			if(m_pMatrix->num_layouts() != 2)
			{
				UG_LOG("ERROR in 'SchurComplementInverse::init': The Operator must"
					   " have two layouts, but the current Operator has only "
					   << m_pMatrix->num_layouts() << "\n");
				return false;
			}
			*/

		//	get matrix from Neumann operator
			m_pNeumannMatrix = &m_NeumannOperator.get_matrix();

		//	Copy Matrix for Neumann Problem
			*m_pNeumannMatrix = *m_pMatrix;

		//	Set Dirichlet values on Pi
			MatSetDirichletOnLayout(m_pNeumannMatrix, *m_pSlaveCPLayout);
			MatSetDirichletOnLayout(m_pNeumannMatrix, *m_pMasterCPLayout);

		//	init sequential solver for Dirichlet problem
			if(m_pNeumannSolver != NULL)
				if(!m_pNeumannSolver->init(m_NeumannOperator))
				{
					UG_LOG("ERROR in 'SchurComplementInverse::init': Cannot init "
							"Sequential Neumann Solver for Operator A.\n");
					return false;
				}

		//	Debug output of matrices
			if(m_pDebugWriter != NULL)
			{
				m_pDebugWriter->write_matrix(m_NeumannOperator.get_matrix(),
				                             "FetiNeumannMatrix");
			}

		//	we're done
			return true;
		} /* end 'SchurComplementInverse::init()' */


	// 	Init for Linear Operator J and Linearization point (current solution)
		virtual bool init(ILinearOperator<vector_type, vector_type>& J, const vector_type& u)
		{
			return init(J);
		}

	// 	Solve A*u = f, such that u = A^{-1} f
		virtual bool apply(vector_type& u, const vector_type& f)
		{
			return false;
		}

	// 	Solve A*u = f, such that u = A^{-1} f
	// 	This is done by iterating: u := u + B(f - A*u)
	// 	In f the last defect f := f - A*u is returned
		virtual bool apply_return_defect(vector_type& u, vector_type& f)
		{
			return false;
		}

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

	//	Copy of matrix
		PureMatrixOperator<vector_type, vector_type, matrix_type> m_NeumannOperator;

	// 	Parallel Neumann Matrix
		matrix_type* m_pNeumannMatrix;

	//	Neumann Solver
		ILinearOperatorInverse<vector_type, vector_type>* m_pNeumannSolver;

	//	Layouts for Pi
		IndexLayout* m_pMasterCPLayout;
		IndexLayout* m_pSlaveCPLayout;

	// 	Convergence Check
		IConvergenceCheck* m_pConvCheck;

	//	Debug Writer
		IDebugWriter<algebra_type>* m_pDebugWriter;
};

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
		FETISolver() :
			m_A(NULL),
			m_pConvCheck(NULL),
			m_pDebugWriter(NULL)
		{}

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

	//	set debug output
		void set_debug(IDebugWriter<algebra_type>* debugWriter)
		{
			m_pDebugWriter = debugWriter;
		}

	///	initializes the solver for operator A
		virtual bool init(IMatrixOperator<vector_type, vector_type, matrix_type>& A)
		{
		//	remember A
			m_A = &A;

		//	get matrix
			m_pMatrix = &m_A->get_matrix();

		//	get layouts
		//	check that matrix has enough decomposition levels - have to be already
		//  created on script level (via 'BuildDomainDecompositionLayoutsTest2d')
			if(m_pMatrix->num_layouts() != 3)
			{
				UG_LOG("ERROR in 'FETISolver::init': The Operator must"
					   " have three layouts, but the current Operator has only "
					   << m_pMatrix->num_layouts() << "\n");
				return false;
			}
			m_pMasterDeltaLayout = &m_pMatrix->get_master_layout(1);
			m_pSlaveDeltaLayout  = &m_pMatrix->get_slave_layout(1);
			m_pMasterCPLayout    = &m_pMatrix->get_master_layout(2);
			m_pSlaveCPLayout     = &m_pMatrix->get_slave_layout(2);

		//	write layouts
			pcl::SynchronizeProcesses();
			UG_LOG("------------- BEFORE ------------\n")
			UG_LOG("------------- DELTA MASTER ------------\n")
			LogIndexLayoutOnAllProcs(*m_pMasterDeltaLayout, 1);
			pcl::SynchronizeProcesses();
			UG_LOG("------------- DELTA SLAVE  ------------\n")
			LogIndexLayoutOnAllProcs(*m_pSlaveDeltaLayout, 1);

			/* 
		//	create "PI layout" containing cross points by extracting them from "Delta layout" ...
			ExtractCrossPointLayouts(m_pMatrix->num_rows(), // number of dof's of current processor
									 m_pMasterDeltaLayout,   //m_pMatrix->get_master_layout(1),
									 m_pSlaveDeltaLayout,    //m_pMatrix->get_slave_layout(1),
									 m_pDDInfo,              // contains mapping "proc id" ==> "subdom id"
									 m_pMasterCPLayout,      //m_masterCPLayout,
									 m_pSlaveCPLayout);      //m_slaveCPLayout
			*/

		//	write layouts
			pcl::SynchronizeProcesses();
			UG_LOG("------------- AFTER ------------\n")
			UG_LOG("------------- DELTA MASTER ------------\n")
			LogIndexLayoutOnAllProcs(*m_pMasterDeltaLayout, 1);
			pcl::SynchronizeProcesses();
			UG_LOG("------------- DELTA SLAVE  ------------\n")
			LogIndexLayoutOnAllProcs(*m_pSlaveDeltaLayout, 1);

			pcl::SynchronizeProcesses();
			UG_LOG("------------- PI MASTER ------------\n")
			LogIndexLayoutOnAllProcs(*m_pMasterCPLayout, 1);
			pcl::SynchronizeProcesses();
			UG_LOG("------------- PI SLAVE  ------------\n")
			LogIndexLayoutOnAllProcs(*m_pSlaveCPLayout, 1);

		//  ----- INIT DIRICHLET SOLVER  ----- //

		//	check that dirichlet solver has been set
			if(m_pDirichletSolver == NULL)
			{
				UG_LOG("ERROR in FETISolver::init: No dirichlet solver set "
						" for inversion of A_{II} in Local Schur complement.\n");
				return false;
			}

		//	set cross point layouts
			m_LocalSchurComplement.set_crosspoint_layouts(*m_pSlaveCPLayout, *m_pMasterCPLayout);

		//	set dirichlet solver for local schur complement
			m_LocalSchurComplement.set_dirichlet_solver(*m_pDirichletSolver);

		//	set operator in local schur complement
			m_LocalSchurComplement.set_matrix(*m_A);

		//	init local Schur complement
			if(m_LocalSchurComplement.init() != true)
			{
				UG_LOG("ERROR in FETISolver::init: Can not init local Schur "
						"complement.\n");
				return false;
			}

		//  ----- INIT NEUMANN SOLVER  ----- //

		//	check that neumann solver has been set
			if(m_pNeumannSolver == NULL)
			{
				UG_LOG("ERROR in FETISolver::init: No neumann solver set "
						" for inversion of A_{I,Delta}{I,Delta} in Local Schur complement.\n");
				return false;
			}

		//	set cross point layouts
			m_SchurComplementInverse.set_crosspoint_layouts(*m_pSlaveCPLayout, *m_pMasterCPLayout);

		//	set neumann solver in SchurComplementInverse
			m_SchurComplementInverse.set_neumann_solver(*m_pNeumannSolver);

		//	init Schur complement inverse
			if(m_SchurComplementInverse.init(*m_A) != true)
			{
				UG_LOG("ERROR in FETISolver::init: Can not init Schur "
						"complement inverse.\n");
				return false;
			}

		//	we're done
			return true;
		} /* end 'FETISolver::init()' */

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
		template <typename TVector>
		bool apply_F(TVector& f, const TVector& v)
		{
			//	Help vector
			TVector fTmp; fTmp.create(v.size());

			//	0. Reset values of f, fTmp
			f.set(0.0); fTmp.set(0.0);

			//	1. Apply transposed jump operator: f = B_{\Delta}^T * v_{\Delta}:
			ComputeDifferenceOnDeltaTransposed(f, v, *m_pMasterDeltaLayout, *m_pSlaveDeltaLayout); //(v.get_master_layout(1), v.get_slave_layout(1));

			//  2. Apply SchurComplementInverse to f - TODO: implement 'm_SchurComplementInverse.apply()'!
			m_SchurComplementInverse.apply(fTmp, f);

			//	3. Apply jump operator to get the final 'f'
			ComputeDifferenceOnDelta(f, fTmp, *m_pMasterDeltaLayout, *m_pSlaveDeltaLayout); //(v.get_master_layout(1), v.get_slave_layout(1));

			//	we're done
			return true;
		}

	///	function which computes right hand side vector 'd' of the reduced system ("Delta system")
	/**
	 * This function computes \f$d := B_{\Delta} \tilde{S}^{-1} \tilde{f}_{\Delta}\f$
	 * to a vector \f$v\f$. \f$v\f$ can be:
	 * (a) unknown vector of Lagrange multipliers, lambda, and
	 * (b) search direction 'p' in cg method.
	 *
	 * \param[in]		f				vector \f$\tilde{f}_{\Delta}\f$
	 * \param[out]		d				right hand side vector \f$d\f$ of reduced system
	 */
		template <typename TVector>
		bool compute_d(TVector& d, const TVector& f)
		{
			//	Help vector
			TVector dTmp; dTmp.create(f.size());

			//	0. Reset values of f, fTmp
			d.set(0.0); dTmp.set(0.0);

			//  1. Apply SchurComplementInverse to 'f' - TODO: implement 'm_SchurComplementInverse.apply()'!
			m_SchurComplementInverse.apply(dTmp, f);

			//	2. Apply jump operator to get the final 'd'
			ComputeDifferenceOnDelta(d, dTmp, *m_pMasterDeltaLayout, *m_pSlaveDeltaLayout); //(f.get_master_layout(1), f.get_slave_layout(1));

			//	we're done
			return true;
		}


	///	function which applies diagonal scaling matrix \f$D_{\Delta}^{(i)}\f$ to a vector \f$v\f$
		template <typename TVector>
		bool Apply_ScalingMatrix(TVector& s, const TVector& v) // maybe restrict to layout
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
		template <typename TVector>
		bool apply_M_inverse(TVector& z, const TVector& r)
		{
			//	Help vector
			TVector zTmp; zTmp.create(r.size());

			//	0. Reset values of z, zTmp
			z.set(0.0); zTmp.set(0.0);

			//	1. Apply scaling: z := D_{\Delta}^{(i)} * r
			Apply_ScalingMatrix(z, r); // maybe restrict to layout

			//  2. Apply transposed jump operator: zTmp := B_{\Delta}^T * z
			ComputeDifferenceOnDeltaTransposed(zTmp, z, *m_pMasterDeltaLayout, *m_pSlaveDeltaLayout); //(r.get_master_layout(1), r.get_slave_layout(1));

			//	3. Apply local Schur complement: z := S_{\Delta}^{(i)} * zTmp
			m_LocalSchurComplement.apply(z, zTmp);

			//  4. Apply jump operator:  zTmp :=  B_{\Delta} * z
			ComputeDifferenceOnDelta(zTmp, z, *m_pMasterDeltaLayout, *m_pSlaveDeltaLayout); //(r.get_master_layout(1), r.get_slave_layout(1));

			//	5. Apply scaling: z := D_{\Delta}^{(i)} * zTmp to get the final 'z'
			Apply_ScalingMatrix(z, zTmp); // maybe restrict to layout

			//	we're done
			return true;
		}

		template <typename TVector>
		bool apply_M_inverse_with_identity_scaling(TVector& z, const TVector& r)
		{
			//	Help vector
			TVector zTmp; zTmp.create(r.size());

			//	0. Reset values of z, zTmp
			z.set(0.0); zTmp.set(0.0);

			//	1. Apply scaling: z := D_{\Delta}^{(i)} * r
			//Apply_ScalingMatrix(z, r); // maybe restrict to layout

			//  2. Apply transposed jump operator: zTmp := B_{\Delta}^T * z
			ComputeDifferenceOnDeltaTransposed(z, r, *m_pMasterDeltaLayout, *m_pSlaveDeltaLayout); //(r.get_master_layout(1), r.get_slave_layout(1));

			//	3. Apply local Schur complement: z := S_{\Delta}^{(i)} * zTmp
			m_LocalSchurComplement.apply(zTmp, z);

			//  4. Apply jump operator:  zTmp :=  B_{\Delta} * z
			ComputeDifferenceOnDelta(z, zTmp, *m_pMasterDeltaLayout, *m_pSlaveDeltaLayout); //(r.get_master_layout(1), r.get_slave_layout(1));

			//	5. Apply scaling: z := D_{\Delta}^{(i)} * zTmp to get the final 'z'
			//Apply_ScalingMatrix(z, zTmp); // maybe restrict to layout

			//	we're done
			return true;
		}

	//	TODO:
	// 1. x = lambda; Compute 'd' before calling apply_return_defect()! Note: 'ApplyLinearSolver()' calls 'apply(), not 'apply_return_defect()'!
	// 2. With \f$\lambda\f$ found, back solve for \f$u_{\Delta}\f$:
	//    \f$u_{\Delta} = {\tilde{S}_{\Delta \Delta}}^{-1} ({\tilde{f}_{\Delta}} - B_{\Delta}^T \lambda).\f$
	// 3. Assemble this and the solutions for \f$u_{I}\f$ and \f$u_{\Pi}\f$ to the global solution vector

	///	solves the reduced system \f$F \lambda = d\f$ with preconditioned cg method
	///	and returns the last defect of iteration in rhs
	/// (derived from 'CGSolver::apply_return_defect()')
		virtual bool apply_return_defect(vector_type& lambda, vector_type& d)
		{
			if(m_pConvCheck == NULL)
			{
				UG_LOG("ERROR: In 'FETISolver::apply': Convergence check not set.\n");
				return false;
			}

			#ifdef UG_PARALLEL
			if(!d.has_storage_type(PST_ADDITIVE) || !lambda.has_storage_type(PST_CONSISTENT))
				{
					UG_LOG("ERROR: In 'FETISolver::apply':Inadequate storage format of Vectors 'd' and 'lambda'.\n");
					return false;
				}
			#endif

		// 	rename r as d (for convenience)
			vector_type& r = d;

		// 	Build residuum:  r := d - F*lambda - TODO: mit apply_F(TVector& f, const TVector& v) machen!
			/*
 			if(!m_A->apply_sub(r, lambda)) 
				{UG_LOG("ERROR in 'FETISolver::apply': "
						"Unable to build defect. Aborting.\n");return false;}
			*/

		// 	create help vector (h will be consistent r)
		//	todo: 	It would be sufficient to copy only the pattern and
		//			without initializing, but in parallel we have to copy communicators
			vector_type t; t.create(r.size()); t = r;
			vector_type z; z.create(lambda.size()); z = lambda;
			vector_type p; p.create(lambda.size()); p = lambda;

		// 	Preconditioning: apply z = M^-1 * r
			if (!apply_M_inverse_with_identity_scaling(z, r))
			{
				UG_LOG("ERROR in 'FETISolver::apply': "
					   "Cannot apply preconditioner. Aborting.\n"); return false;
			}
			else z = r;

			// make z consistent
			if(!z.change_storage_type(PST_CONSISTENT))
			{
				UG_LOG("ERROR in 'FETISolver::apply': "
					   "Cannot convert z to consistent vector.\n"); return false;
			}

			prepare_conv_check();
			m_pConvCheck->start(r);

			number rho, rho_new, beta, alpha, alpha_denominator;
			rho = rho_new = beta = alpha = alpha_denominator = 0.0;

			// start search direction
			p = z;

			// start rho
			rho = VecProd(z, r);

			m_iterCnt = 0;
		// 	Iteration loop
			while(!m_pConvCheck->iteration_ended())
			{
				m_iterCnt++;
			// 	Build t = F*p (t is additive afterwards)
				if(!apply_F(t, p))
				{
					UG_LOG("ERROR in 'FETISolver::apply': Unable "
								"to build t = F*p. Aborting.\n"); return false;
				}

			// 	Compute alpha
				alpha_denominator = VecProd(t, p);
				alpha = rho/alpha_denominator;

			// 	Update lambda := lambda + alpha*p
				VecScaleAdd(lambda, 1.0, lambda, alpha, p);

			// 	Update r := r - alpha*t
				VecScaleAdd(r, 1.0, r, -alpha, t);

			// 	Check convergence
				m_pConvCheck->update(r);

			// 	Preconditioning: apply z = M^-1 * r
				if (!apply_M_inverse_with_identity_scaling(z, r))
				{
					UG_LOG("ERROR in 'FETISolver::apply': "
						   "Cannot apply preconditioner. Aborting.\n"); return false;
				}
				else z = r;

			// 	make z consistent
				if(!z.change_storage_type(PST_CONSISTENT))
				{
					UG_LOG("ERROR in 'FETISolver::apply': "
						   "Cannot convert z to consistent vector.\n"); return false;
				}

			// 	new rho
				rho_new = VecProd(z, r);

			// 	new beta
				beta = rho_new/rho;

			// 	new direction p:= beta*p + z
				VecScaleAdd(p, beta, p, 1.0, z);

			// 	update rho
				rho = rho_new;
			}

			return m_pConvCheck->post();
		} /* end 'FETISolver::apply_return_defect()' */

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
		IMatrixOperator<vector_type,vector_type,matrix_type>* m_A;

	// 	Parallel Matrix to invert
		matrix_type* m_pMatrix;

	//	Layouts for Delta
		IndexLayout* m_pMasterDeltaLayout;
		IndexLayout* m_pSlaveDeltaLayout;

	//	Layouts for Pi
		IndexLayout* m_pMasterCPLayout;
		IndexLayout* m_pSlaveCPLayout;
		//IndexLayout m_masterCPLayout;
		//IndexLayout m_slaveCPLayout;

	//	Local Schur complement for each feti subdomain
		LocalSchurComplement<algebra_type> m_LocalSchurComplement;

	//	Dirichlet solver for inverse of A_{II} in local schur complement
		ILinearOperatorInverse<vector_type, vector_type>* m_pDirichletSolver;

	//	SchurComplementInverse
		SchurComplementInverse<algebra_type> m_SchurComplementInverse;

	//	Neumann solver for inverse of A_{\{I,\Delta\}, \{I,\Delta\}} in the
	//	creation of the S_{\Pi \Pi} schur complement
		ILinearOperatorInverse<vector_type, vector_type>* m_pNeumannSolver;

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
