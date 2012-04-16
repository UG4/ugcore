/**
 * \file amg_base_impl.h
 *
 * \author Martin Rupp
 *
 * \date 01.12.2010
 *
 * implementation file for base amg functionality
 *
 * Goethe-Center for Scientific Computing 2010.
 */

#ifndef __H__UG__LIB_ALGEBRA__AMG_SOLVER__AMG_BASE_IMPL_H__
#define __H__UG__LIB_ALGEBRA__AMG_SOLVER__AMG_BASE_IMPL_H__

#include <fstream>
#include "stopwatch.h"
#include "amg_debug.h"
#ifdef UG_PARALLEL
#include "lib_algebra/parallelization/collect_matrix.h"
#include "lib_algebra/parallelization/parallel_vector.h"
#include "lib_algebra/parallelization/parallel_storage_type.h"
#include "lib_algebra/parallelization/consistency_check.h"
#endif

#include "lib_algebra/common/connection_viewer_output.h"

std::string GetProcFilename(std::string path, std::string name, std::string extension);

namespace ug{


#ifdef UG_PARALLEL
template<typename TMatrix, typename TVector>
void SetParallelVectorAsMatrix(TVector &v, TMatrix &m, ParallelStorageType t)
{
	v.set_communicator(m.communicator());
	v.set_process_communicator(m.process_communicator());
	v.set_layouts(m.master_layout(), m.slave_layout());
	v.set_storage_type(t);
}
#endif

template<typename TMatrix>
size_t GetMaxConnections(const TMatrix &A)
{
	size_t m=0;
	for(size_t i=0; i<A.num_rows(); i++)
	{
		//if(m < A.num_connections(i)) m = A.num_connections(i);

		size_t n=0;
		for(typename TMatrix::const_row_iterator it = A.begin_row(i); it != A.end_row(i); ++it)
			if(it.value() != 0.0) n++;
		if(m < n) m = n;
	}

	return m;
}

/// returns the number of non-zeroes (!= number of connections)
template<typename TMatrix>
size_t GetNNZs(const TMatrix &A)
{
	size_t m=0;
	for(size_t i=0; i<A.num_rows(); i++)
	{
		for(typename TMatrix::const_row_iterator it = A.begin_row(i); it != A.end_row(i); ++it)
			if(it.value() != 0.0) m++;
	}
	return m;
}

template<typename TAlgebra>
void AMGBase<TAlgebra>::calculate_level_information(size_t level, double createAMGlevelTiming)
{
	PROFILE_FUNC();
	AMGLevel &L = *levels[level];
	matrix_operator_type &A = *L.pA;

	LevelInformation &li = L.m_levelInformation;

	size_t nnz = GetNNZs(A); // A.total_num_connections();
	size_t maxConnections = GetMaxConnections(A);
	li.m_dCreationTimeMS = createAMGlevelTiming;
#ifdef UG_PARALLEL
	pcl::ProcessCommunicator &pc = A.process_communicator(); //!
	size_t N = A.num_rows() - A.slave_layout().num_interface_elements();
	li.set_nr_of_nodes(pc.allreduce(N, PCL_RO_MIN),	pc.allreduce(N, PCL_RO_MAX), pc.allreduce(N, PCL_RO_SUM));
	li.set_nnz(	pc.allreduce(nnz, PCL_RO_MIN), pc.allreduce(nnz, PCL_RO_MAX), pc.allreduce(nnz, PCL_RO_SUM));
	li.set_max_connections(pc.allreduce(maxConnections, PCL_RO_MAX));
	size_t localInterfaceElements = A.master_layout().num_interface_elements();
	li.m_iInterfaceElements = pc.allreduce(localInterfaceElements, PCL_RO_SUM);
#else
	size_t N = A.num_rows();
	li.set_nr_of_nodes(N, N, N);
	li.set_nnz(nnz, nnz, nnz);
	li.set_max_connections(maxConnections);
	li.m_iInterfaceElements = 0;
#endif
	if(level>0)
		li.set_coarsening_rate((double)li.get_nr_of_nodes()/(double)levels[level-1]->m_levelInformation.get_nr_of_nodes());
	else
		li.set_coarsening_rate(0.0);

	UG_DLOG(LIB_ALG_AMG, 1, li.tostring() << "\n");
	UG_DLOG(LIB_ALG_AMG, 1, " level took " << createAMGlevelTiming << " ms" << std::endl << std::endl);
}


#ifdef UG_PARALLEL
static void eh( MPI_Comm *comm, int *err, ... )
{
	UG_LOG("MPI ERROR!\n");

	char error_string[256];
	int length_of_error_string=256, error_class;

	MPI_Error_class(*err, &error_class);
	MPI_Error_string(error_class, error_string, &length_of_error_string);
	UG_ASSERT(0,"MPI Error: " << error_string << "\n" );
	return;
}
#endif

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename TAlgebra>
bool AMGBase<TAlgebra>::preprocess(matrix_operator_type& mat)
{
	AMG_PROFILE_FUNC();
	if(m_bOneInit && m_bInited)
		return true;

#ifdef UG_PARALLEL
    MPI_Errhandler newerr;

	MPI_Comm_create_errhandler( eh, &newerr );
	MPI_Comm_set_errhandler( MPI_COMM_WORLD, newerr );
	  //MPI_Comm_call_errhandler( MPI_COMM_WORLD, MPI_ERR_OTHER );
#endif

	cleanup();

	update_positions();

	// init m_amghelper for grid printing
	if(m_writeMatrices)
	{
		if(m_dbgPositions.size() > 0)
		{
			m_amghelper.positions.resize(1);
			m_amghelper.positions[0] = m_dbgPositions;
		}
		m_amghelper.parentIndex = &m_parentIndex;
		m_amghelper.dimension = m_dbgDimension;
	}

	UG_ASSERT(m_writeMatrices == false || m_dbgDimension != 0, "you need to provide position information for writing matrices");

	if(m_basesolver.invalid())
	{
		UG_LOG("amg_base::preprocess(): No base solver selected. Call set_base_solver(basesolver) to set a base solver.\n");
		return false;
	}
	if(m_presmoother.invalid())
	{
		UG_LOG("amg_base::preprocess(): No PreSmoother selected. Call set_presmoother(presmoother) to set a PreSmoother.\n");
		return false;
	}
	if(m_postsmoother.invalid())
	{
		UG_LOG("amg_base::preprocess(): No PostSmoother selected. Call set_postsmoother(postsmoother) to set a PostSmoother.\n");
		return false;
	}

	UG_DLOG(LIB_ALG_AMG, 1, "Starting AMG Setup." << std::endl << std::endl);

	stopwatch SWwhole;
	SWwhole.start();

	size_t level=0;

	mat.defragment();

	double createAMGlevelTiming=0;

#ifdef UG_PARALLEL
	m_agglomerateLevel = -1;
#endif

	levels.clear();
	AMGLevel *pL = new AMGLevel(0);
	pL->pA = &mat;
#ifdef UG_PARALLEL
	pL->com = mat.communicator();
	pL->bHasBeenMerged = false;
#endif
	levels.push_back(pL);

	for(; ; level++)
	{
		AMGLevel &L = *levels[level];

		calculate_level_information(level, createAMGlevelTiming);
		UG_LOG(get_level_information(level)->tostring() << "\n");

		if(L.m_levelInformation.get_nr_of_nodes() < m_maxNodesForBase
				|| L.m_levelInformation.get_fill_in() > m_dMaxFillBeforeBase
				|| level >= m_maxLevels-1)
		{
			AMG_PROFILE_BEGIN(amg_createDirectSolver);
			create_direct_solver(level);
			AMG_PROFILE_END();
			break;
		}


#ifdef UG_PARALLEL
		//PRINTPC(L.pA->process_communicator());
		agglomerate(level);
		precalc_level(level);
		if(m_agglomerateLevel == level)
		{
			UG_LOG("Processor " << pcl::GetProcRank() << " got agglomerated at level " << level << ".\n");
			break;
		}
#else
		precalc_level(level);
		L.pAgglomeratedA = L.pA;
#endif

		if(m_writeMatrices && m_amghelper.has_positions())
		{
			std::string name = m_writeMatrixPath + std::string("AMG_A_L") + ToString(level) + ".mat";
		#ifdef UG_PARALLEL
			name = GetParallelName(name, L.pAgglomeratedA->process_communicator(), true);
		#endif
			AMGWriteToFile(*L.pAgglomeratedA, level, level, name.c_str(), m_amghelper);
			//PRINTPC(L.pAgglomeratedA->get_process_communicator());
		}

		//smoothem_R[level].init(*m_A[level]);

		levels.push_back(new AMGLevel(level+1));
		AMGLevel &nextL = *levels[level+1];

		L.presmoother = m_presmoother->clone();

		//PRINTPC(L.pAgglomeratedA->get_process_communicator());
		L.presmoother->init(L.pAgglomeratedA);
		if(m_presmoother == m_postsmoother)
			L.postsmoother = L.presmoother;
		else
		{
			L.postsmoother = m_postsmoother->clone();
			L.postsmoother->init(L.pAgglomeratedA);
		}

		nextL.pA = new matrix_operator_type;
		matrix_operator_type &AH = *nextL.pA;
		prolongation_matrix_type &R = L.R;
		prolongation_matrix_type &P = L.P;

#ifdef UG_PARALLEL
		AH.set_storage_type(PST_ADDITIVE);
		R.set_storage_type(PST_ADDITIVE);
		P.set_storage_type(PST_CONSISTENT);

		nextL.com = L.com;

		AH.set_layouts(nextL.masterLayout, nextL.slaveLayout);
		AH.set_communicator(nextL.com);
		AH.set_process_communicator(L.pAgglomeratedA->process_communicator());
		P.set_process_communicator(L.pAgglomeratedA->process_communicator());
		R.set_process_communicator(L.pAgglomeratedA->process_communicator());
#endif

		stopwatch SW; SW.start();

		create_AMG_level(AH, R, *L.pAgglomeratedA, P, level);

		SW.stop();
		createAMGlevelTiming = SW.ms();

		// create vectors for AMG multigrid
		/////////////////////////////////////////

		L.corr.resize(L.pAgglomeratedA->num_rows());
		L.cH.resize(AH.num_rows());
		L.dH.resize(AH.num_rows());

	#ifdef UG_PARALLEL
	#ifdef UG_DEBUG
		ConsistencyCheck(levels[level]->is_fine, L.pAgglomeratedA->communicator(),
				L.pAgglomeratedA->process_communicator(),
				L.pAgglomeratedA->master_layout(), L.pAgglomeratedA->slave_layout(), "fine marks");
//			PrintLayoutIfBroken(AH.get_process_communicator(), nextL.com, nextL.masterLayout, nextL.slaveLayout);
	#endif
		SetParallelVectorAsMatrix(L.corr, *L.pAgglomeratedA, PST_CONSISTENT);
		SetParallelVectorAsMatrix(L.cH, AH, PST_CONSISTENT);
		SetParallelVectorAsMatrix(L.dH, AH, PST_ADDITIVE);
		/*UG_LOG("nextLevel Layouts:\n");
				PrintLayout(nextL.pA->process_communicator(), nextL.pA->communicator(), nextL.masterLayout, nextL.slaveLayout);*/
	#endif

		// todo: set size for variable sized blockvectors
		/*for(size_t i=0; i<N; i++)
			if(nodes[i].isCoarse())
			{
				int rows = GetRows(A.begin_row(i).value());
				UG_ASSERT(newIndex[i] >= 0, "");
				SetSize((*m_vec1[level+1])[newIndex[i]], rows);
				SetSize((*m_vec2[level+1])[newIndex[i]], rows);
			}*/
	}

	UG_LOG("finished. now init fsmoothing.\n");
	init_fsmoothing();

	m_usedLevels = level+1;
	UG_LOG("AMG Setup finished. Used Levels: " << m_usedLevels << ". ");
	m_dTimingWholeSetupMS = SWwhole.ms();
	UG_DLOG(LIB_ALG_AMG, 1, "AMG Setup took " << m_dTimingWholeSetupMS << " ms." << std::endl);

	// calc complexities
	double totalNNZs=0;
	double totalNrOfNodes=0;
	for(size_t i=0; i<m_usedLevels; i++)
	{
		totalNNZs += levels[i]->pA->total_num_connections();
		totalNrOfNodes += levels[i]->pA->num_rows();
	}

	m_dOperatorComplexity = totalNNZs / levels[0]->pA->total_num_connections();
	m_dGridComplexity = totalNrOfNodes / levels[0]->pA->num_rows();

	UG_DLOG(LIB_ALG_AMG, 1, "Operator Complexity: " << m_dOperatorComplexity << " grid complexity: "
			<< m_dGridComplexity << std::endl << std::endl);

	m_bInited = true;
	return true;
}




template<typename TAlgebra>
void AMGBase<TAlgebra>::create_direct_solver(size_t level)
{
	UG_ASSERT(block_traits< typename vector_type::value_type >::is_static, "dynamic not yet implemented");

	AMGLevel &L = *levels[level];
	matrix_operator_type &A = *L.pA;
	size_t static_nrUnknowns = block_traits< typename vector_type::value_type >::static_size;
	UG_DLOG(LIB_ALG_AMG, 1, "Creating level " << level << " (" << A.num_rows() << " nodes, total "
			<< A.num_rows()*static_nrUnknowns << " unknowns)" << std::endl << "Using Direct Solver on Matrix "
			<< A.num_rows()*static_nrUnknowns << "x" << A.num_rows()*static_nrUnknowns << ". ");
	(void) static_nrUnknowns;
	stopwatch SW; SW.start();

	UG_LOG("creating direct solver on level " << level << "\n");

#ifdef UG_PARALLEL
	pcl::ProcessCommunicator &pc = A.process_communicator(); //!
	if(pc.size()>1)
	{
		L.bHasBeenMerged = true;
		L.bLevelHasMergers=true;
		// create basesolver
		collect_matrix(A, *L.collectedA, L.agglomerateMasterLayout, agglomerateSlaveLayout);

		if(m_writeMatrices && m_amghelper.has_positions())
		{
			std::vector<MathVector<3> > vec = m_amghelper.positions[level];
			vec.resize(L.collectedA->num_rows());
			ComPol_VecCopy<std::vector<MathVector<3> > >	copyPol(&vec);
			pcl::InterfaceCommunicator<IndexLayout> &communicator = A.communicator();
			communicator.send_data(agglomerateSlaveLayout, copyPol);
			communicator.receive_data(L.agglomerateMasterLayout, copyPol);
			communicator.communicate();
			if(pcl::GetProcRank() == pc.get_proc_id(0))
			{
				std::string name = (std::string(m_writeMatrixPath) + "collectedA.mat");
				WriteMatrixToConnectionViewer(name.c_str(), *L.collectedA, &vec[0], m_amghelper.dimension);
			}
		}
		if(pcl::GetProcRank() == pc.get_proc_id(0))
		{
			L.agglomeratedPC = pc.create_sub_communicator(true);
			L.collectedA->set_master_layout(m_emptyLayout);
			L.collectedA->set_slave_layout(m_emptyLayout);
			L.collectedA->set_process_communicator(m_emptyPC);
			m_basesolver->init(L.collectedA);

			L.collectedA->set_layouts(m_emptyLayout, m_emptyLayout);
			L.collC.resize(L.collectedA->num_rows());
			L.collD.resize(L.collectedA->num_rows());
			SetParallelVectorAsMatrix(L.collC, *L.collectedA, PST_ADDITIVE);
			SetParallelVectorAsMatrix(L.collD, *L.collectedA, PST_ADDITIVE);
		}
		else
		{
			m_agglomerateLevel = level;
			L.agglomeratedPC = pc.create_sub_communicator(false);
		}
	}
	else
	{
		L.bLevelHasMergers=false;
		m_basesolver->init(L.pA);
	}
#else
	m_basesolver->init(&A);
#endif

	m_dTimingCoarseSolverSetupMS = SW.ms();
	UG_DLOG(LIB_ALG_AMG, 1, "Coarse Solver Setup took " << m_dTimingCoarseSolverSetupMS << "ms." << std::endl);
}

//!
//! amg constructor
template<typename TAlgebra>
AMGBase<TAlgebra>::AMGBase() :
	m_numPreSmooth(2),
	m_numPostSmooth(2),
	m_cycleType(1),

	m_maxLevels(100),
	m_usedLevels(0),

	m_presmoother(NULL),
	m_postsmoother(NULL),
	m_basesolver(NULL),

	m_pPositionProvider2d(NULL),
	m_pPositionProvider3d(NULL)
{
	m_vec4 = NULL;
	m_bInited = false;

	m_maxNodesForBase = 100;
	m_minNodesOnOneProcessor = 1;
	m_preferredNodesOnOneProcessor = 10000;
	m_dMaxFillBeforeBase = 0.5;

	m_writeMatrices = false;
	m_bFSmoothing = false;

	m_dbgDimension = 0;
	iteration_glboal=0;
	m_checkLevelPostIterations = 10;
	m_iNrOfPreiterationsCheck = 8;
	m_iYCycle = 0;
}


template<typename TAlgebra>
void AMGBase<TAlgebra>::cleanup()
{
	m_usedLevels = 0;
	m_bInited=false;

	if(m_vec4 != NULL) delete m_vec4; m_vec4 = NULL;

	for(size_t k=0; k<levels.size(); k++)
		delete levels[k];
}
//!
//! amg destructor
template<typename TAlgebra>
AMGBase<TAlgebra>::~AMGBase()
{
	UG_LOG("~AMGBase.\n");
	cleanup();
}

template<typename TAlgebra>
void AMGBase<TAlgebra>::init_fsmoothing()
{
	AMG_PROFILE_FUNC();
#ifdef UG_PARALLEL

	UG_LOG("\n\nInit FSmoothing...\n");
	// last level is basesolver level
	for(size_t k=0; k<levels.size()-1; k++)
	{
		//UG_LOG("---- level " << k << " ---- ");
		matrix_operator_type &mat = *levels[k]->pAgglomeratedA;

		ParallelVector<Vector< typename matrix_operator_type::value_type > > m_diag;
		size_t size = mat.num_rows();
		levels[k]->m_diagInv.resize(size);
		m_diag.resize(size);

		m_diag.set_layouts(mat.master_layout(), mat.slave_layout());
		m_diag.set_communicator(mat.communicator());

		// copy diagonal
		for(size_t i = 0; i < m_diag.size(); ++i)
			m_diag[i] = mat(i, i);

		//	make diagonal consistent
		m_diag.set_storage_type(PST_ADDITIVE);

		/*UG_LOG("______________________________________________\n");
		PrintLayout(mat.get_process_communicator(), mat.get_communicator(),
						mat.get_master_layout(), mat.get_slave_layout());*/
		m_diag.change_storage_type(PST_CONSISTENT);

		for(size_t i = 0; i < m_diag.size(); ++i)
			GetInverse(levels[k]->m_diagInv[i], m_diag[i]);
	}
#endif
}


template<typename TAlgebra>
bool AMGBase<TAlgebra>::f_smoothing(vector_type &corr, vector_type &d, size_t level)
{
	AMG_PROFILE_FUNC();
	AMGLevel &L = *levels[level];
	stdvector<bool> &is_fine = L.is_fine;
	matrix_operator_type &A = *L.pAgglomeratedA;
#ifdef UG_PARALLEL

	//UG_ASSERT(L.m_diagInv.size() == is_fine.size(), "fine markers do not match in size on level " << level);
	for(size_t i=0; i<L.m_diagInv.size(); i++)
	{
		if(is_fine[i])
			corr[i] = d[i] * L.m_diagInv[i];
		else
			corr[i] = 0.0;
	}
	corr.set_storage_type(PST_ADDITIVE);
	corr.change_storage_type(PST_CONSISTENT);
#else
	UG_ASSERT(A.num_rows() == is_fine.size(), "fine markers do not match in size on level " << level);
	for(size_t i=0; i<A.num_rows(); i++)
	{
		if(is_fine[i])
			corr[i] = d[i] / A(i,i);
		else
			corr[i] = 0.0;
	}
#endif
	A.matmul_minus(d, corr);
	return true;
}

template<typename TAlgebra>
bool AMGBase<TAlgebra>::solve_on_base(vector_type &c, vector_type &d, size_t level)
{
	AMG_PROFILE_FUNC();
	m_basesolver->apply_return_defect(c, d);
#if 0
	AMGLevel &L = *levels[level];
	matrix_operator_type &A = *L.pA;
#ifdef UG_PARALLEL
	if(A.process_communicator().size()>1)
	{
		vector_type collC;
		vector_type collD;
		pcl::InterfaceCommunicator<IndexLayout> &com = A.get_communicator();
		if(pcl::GetProcRank() == 0)
		{
			size_t N = L.collectedA->num_rows();
			collC.resize(N);
			collC.set(0.0);
			collC.set_master_layout(m_emptyLayout);
			collC.set_slave_layout(m_emptyLayout);
			collC.set_process_communicator(m_emptyPC);

			collD = d;
			collD.resize(N);
			for(size_t i=A.num_rows(); i<N; i++)
				collD[i] = 0.0;

			collD.set_master_layout(m_emptyLayout);
			collD.set_slave_layout(m_emptyLayout);
			collD.set_process_communicator(m_emptyPC);



			com.communicate();
		}
		// send d -> collD
		ComPol_VecAdd<vector_type > compolAdd(&collD, &d);
		pcl::InterfaceCommunicator<IndexLayout> &com = A.communicator();
		com.receive_data(L.agglomerateMasterLayout, compolAdd);
		com.send_data(agglomerateSlaveLayout, compolAdd);
		com.communicate();


		if(pcl::GetProcRank() == 0)
			m_basesolver->apply_return_defect(collC, collD);

		// send c -> collC
		ComPol_VecCopy<vector_type> compolCopy(&c, &collC);
		com.receive_data(L.agglomerateMasterLayout, compolPoly);
		com.send_data(agglomerateSlaveLayout, compolPoly);
		com.communicate();

		if(pcl::GetProcRank() == 0)
		{
			for(size_t i=0; i<A.num_rows(); i++)
				d[i] = collD[i];
			for(size_t i=0; i<A.num_rows(); i++)
				c[i] = collC[i];
		}
		else
			A.matmul_minus(d, c);
		d.set(0.0);
	}
	else
#endif
		m_basesolver->apply_return_defect(c, d);
#endif
	return true;
}



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// add_correction_and_update_defect:
//------------------------------------

template<typename TAlgebra>
bool AMGBase<TAlgebra>::add_correction_and_update_defect2(vector_type &c, vector_type &d,
		matrix_operator_type &A, size_t level)
{
	AMG_PROFILE_FUNC();
	//PRINTPC(A.get_process_communicator());
	// c = consistent, d = additiv

	if(level == m_usedLevels-1)
	{
		solve_on_base(c, d, level);
		return true;
	}

	AMGLevel &L = *levels[level];

	UG_ASSERT(c.size() == d.size() && c.size() == A.num_rows(),
			"c.size = " << c.size() << ", d.size = " << d.size() << ", A.size = " << A.num_rows() << ": not matching");



#ifdef UG_PARALLEL
	if(!d.has_storage_type(PST_ADDITIVE) || !c.has_storage_type(PST_CONSISTENT))
	{
		UG_LOG("ERROR: In 'amg::check':Inadequate storage format of Vectors.\n");
		return false;
	}
#endif
	vector_type &corr = levels[level]->corr;
	corr.set(0.0);

#ifdef UG_PARALLEL
	corr.set_storage_type(PST_CONSISTENT);
	c.set_storage_type(PST_CONSISTENT);
#endif


#if WRITEVEC_IN_SOLVER
	writevec((std::string("AMG_") + ToString(iteration_glboal++) + "aa_d_L").c_str(), d, level);
	writevec((std::string("AMG_") + ToString(iteration_glboal++) + "aa_c_L").c_str(), c, level);
#endif

	// presmooth
	for(size_t i=0; i < m_numPreSmooth; i++)
	{
		corr.set(0.0);
		L.presmoother->apply_update_defect(corr, d);
        c += corr;
	}
#ifdef UG_PARALLEL
	UG_ASSERT(c.has_storage_type(PST_CONSISTENT), "" );
#endif
	// pre f-smoothing
	if(m_bFSmoothing)
	{
		corr.set(0.0);
		f_smoothing(corr, d, level);
		c+=corr;
	}
#ifdef UG_PARALLEL
	UG_ASSERT(c.has_storage_type(PST_CONSISTENT), "" );
#endif

#if WRITEVEC_IN_SOLVER
	writevec((std::string("AMG_") + ToString(iteration_glboal++) + "pa_d_L").c_str(), d, level);
	writevec((std::string("AMG_") + ToString(iteration_glboal++) + "pa_c_L").c_str(), c, level);
#endif

	vector_type &cH = levels[level]->cH;
	vector_type &dH = levels[level]->dH;
	cH.set(0.0);

#ifdef UG_PARALLEL
	cH.set_storage_type(PST_CONSISTENT);
#endif

	dH.set(0.0);
	// restrict defect
	// R = additive, d additive -> consistent. dH is then additive.
	// dH = R*d;
	L.R.apply(dH, d);
#ifdef UG_PARALLEL
	UG_ASSERT(dH.has_storage_type(PST_ADDITIVE), "" );
#endif


#if WRITEVEC_IN_SOLVER
	writevec((std::string("AMG_") + ToString(iteration_glboal++) + "_dH_L").c_str(), dH, level+1);
#endif

	// apply lmgc on coarser nodes

	if(level+1 == m_usedLevels-1)
		add_correction_and_update_defect(cH, dH, level+1);
	else
	{
		if(m_iYCycle == 0)
		{
			for(int i=0; i< m_cycleType; i++)
				add_correction_and_update_defect(cH, dH, level+1);
		}
		else
		{
			double preHnorm = ConstTwoNorm(dH);
			for(size_t i=0; i<m_iYCycle; i++)
			{
				add_correction_and_update_defect(cH, dH, level+1);
				double nH = ConstTwoNorm(dH);
				if(nH/preHnorm < m_dYreduce || nH < m_dYabs)
					break;
			}
		}
	}

	corr.set(0.0);
	//cH.set(0.0);
	// interpolate correction. P = consistent, cH = consistent. -> corr consistent
	// corr = R*cH
	L.P.apply(corr, cH);
#ifdef UG_PARALLEL
	UG_ASSERT(corr.has_storage_type(PST_CONSISTENT), "" );
#endif

	// add coarse grid correction to level correction
	// c += corr;

	c += corr;

		// update defect
	// d = d - Ah*corr
	A.matmul_minus(d, corr);


#if WRITEVEC_IN_SOLVER
	writevec((std::string("AMG_") + ToString(iteration_glboal++) + "cc_d_L").c_str(), d, level);
		writevec((std::string("AMG_") + ToString(iteration_glboal++) + "cc_c_L" ).c_str(), c, level);
#endif


	// post f-smoothing
	if(m_bFSmoothing)
	{
		corr.set(0.0);
		f_smoothing(corr, d, level);
		c+=corr;
	}

	// postsmooth
	for(size_t i=0; i < m_numPostSmooth; i++)
	{
		corr.set(0.0);
		L.postsmoother->apply_update_defect(corr, d);
		c += corr;
	}

#if WRITEVEC_IN_SOLVER
	writevec((std::string("AMG_") + ToString(iteration_glboal++) + "pp_d_L").c_str(), d, level);
		writevec((std::string("AMG_") + ToString(iteration_glboal++) + "pp_c_L").c_str(), c, level);
#endif

	return true;
}



template<typename TAlgebra>
bool AMGBase<TAlgebra>::get_correction(vector_type &c, const vector_type &const_d)
{
	AMG_PROFILE_FUNC();
	UG_ASSERT(c.size() == const_d.size() && c.size() == levels[0]->pA->num_rows(),
				"c.size = " << c.size() << ", d.size = " << const_d.size() << ", A.size = " <<
					levels[0]->pA->num_rows() << ": not matching");

	if(m_vec4 == NULL)
		m_vec4 = new vector_type;

	m_vec4->resize(const_d.size());
#ifdef UG_PARALLEL
	m_vec4->set_layouts(c.master_layout(), c.slave_layout());
	m_vec4->set_storage_type(PST_ADDITIVE);
#endif
	vector_type &d = *m_vec4;


	d = const_d;

	/*if(m_usedLevels == 1)
	{
		solve_on_base(c, d, 0);
		return true;
	}
	else*/
	{
		c.set(0.0);
		return add_correction_and_update_defect(c, d);
	}
}



template<typename TAlgebra>
void AMGBase<TAlgebra>::tostring() const
{
	UG_LOG("AMGBase:\n");
	UG_LOG(" nr of pre smoothing steps (nu1)  = " << m_numPreSmooth<< std::endl);
	UG_LOG(" nr of post smoothing steps (nu2) = " << m_numPostSmooth << std::endl);
	UG_LOG("F-Smoothing is " << (m_bFSmoothing ? "[ON]" : "[OFF]") << std::endl);

	UG_LOG(" cycle type = ");
	if(m_cycleType == 1) { UG_LOG("V-cycle\n"); }
	else if(m_cycleType == 2) { UG_LOG("W-cycle\n"); }
	else { UG_LOG("gamma = " << m_cycleType << "\n"); }

	UG_LOG(" max levels = " << m_maxLevels << std::endl);

	if(m_presmoother.valid()) 	{UG_LOG(" presmoother is " << m_presmoother->name() << ".\n");}
	else						{UG_LOG(" no presmoother set!\n");}

	if(m_postsmoother.valid()) 	{UG_LOG(" postsmoother is " << m_postsmoother->name() << ".\n");}
	else						{UG_LOG(" no postsmoother set!\n");}

	if(m_basesolver.valid())	{UG_LOG(" basesolver set\n");}
	else				{UG_LOG(" no basesolver set!\n");}
	UG_LOG("Using base solver for max. " << m_maxNodesForBase << " nodes or " << m_dMaxFillBeforeBase*100 << "% fill in.\n")

	if(m_writeMatrices)
		UG_LOG("Writing Matrices to \"" << m_writeMatrixPath << "\"\n");

	if(m_bInited)
	{
		UG_LOG("AMG is initialized.\n");
		UG_LOG("Used Levels: " << m_usedLevels << "\n");
	}
}


template<typename TAlgebra>
void AMGBase<TAlgebra>::update_positions()
{
	AMG_PROFILE_FUNC();
	UG_ASSERT(m_dbgDimension != 0 || m_pPositionProvider2d == NULL
			|| m_pPositionProvider3d == NULL, "specify EITHER positions 2d or 3d");
	if(m_pPositionProvider2d == NULL && m_pPositionProvider3d == NULL)
		return;
	if(m_pPositionProvider2d)
	{
		std::vector<MathVector<2> > pos;
		m_pPositionProvider2d->get_positions(pos);
		set_positions(&pos[0], pos.size());
	}
	if(m_pPositionProvider3d)
	{
		std::vector<MathVector<3> > pos;
		m_pPositionProvider3d->get_positions(pos);
		set_positions(&pos[0], pos.size());
	}
}



template<typename matrix_type>
static std::string GetAMGFilename(const matrix_type &m, std::string name, int level, std::string extension,
		bool bWriteHeader=false)
{

#ifdef UG_PARALLEL
	return GetParallelName(std::string("AMG_") + name + "_L" + ToString(level) + extension,
			m.process_communicator(), bWriteHeader);
#else
	return std::string("AMG_") + name + "_L" + ToString(level) + extension;
#endif
}

template<typename TAlgebra>
template<typename TNodeType>
void AMGBase<TAlgebra>::write_debug_matrix_markers
	(size_t level, const TNodeType &nodes)
{
	matrix_type &A = *levels[level]->pA;
	// todo: replace some day with something like nodes.get_mark_count(), nodes.mark_name(i), nodes.is_marked(i)
	AMG_PROFILE_FUNC();
	std::fstream ffine((m_writeMatrixPath + GetAMGFilename(A, "fine", level, ".marks")).c_str(), std::ios::out);
	ffine << "1 0 0 1 0\n";
	std::fstream ffine2((m_writeMatrixPath + GetAMGFilename(A, "aggfine", level, ".marks")).c_str(), std::ios::out);
	ffine2 << "1 0.2 1 1 0\n";
	std::fstream fcoarse((m_writeMatrixPath + GetAMGFilename(A, "coarse", level, ".marks")).c_str(), std::ios::out);
	fcoarse << "0 0 1 1 2\n";
	std::fstream fother((m_writeMatrixPath + GetAMGFilename(A, "other", level, ".marks")).c_str(), std::ios::out);
	fother << "1 1 0 1 0\n";
	std::fstream fdirichlet((m_writeMatrixPath + GetAMGFilename(A, "dirichlet", level, ".marks")).c_str(), std::ios::out);
	fdirichlet << "0 1 1 1 0\n";
	for(size_t i=0; i < nodes.size(); i++)
	{
		//int o = m_amghelper.GetOriginalIndex(level, i);
		int o=i;
		if(nodes[i].is_aggressive_fine()) ffine2 << o << "\n";
		else if(nodes[i].is_fine_direct()) ffine << o << "\n";
		else if(nodes[i].is_coarse()) fcoarse << o << "\n";
		else if(nodes[i].is_dirichlet()) fdirichlet << o << "\n";
		else fother << o << "\n";
	}

}




template<typename TAlgebra>
template<typename TMatrix>
void AMGBase<TAlgebra>::
	write_debug_matrix(TMatrix &mat, size_t fromlevel, size_t tolevel, const char *name)
{
	if(m_writeMatrices)
	{
		AMG_PROFILE_FUNC();
		std::string filename = m_writeMatrixPath + name + ToString(fromlevel) + ".mat";
#ifdef UG_PARALLEL
		filename = GetParallelName(filename, mat.process_communicator(), true);
#endif

		int level = std::min(fromlevel, tolevel);
		AMGWriteToFile(mat, fromlevel, tolevel, filename.c_str(), m_amghelper);
		std::fstream f2(filename.c_str(), std::ios::out | std::ios::app);
		f2 << "c " << GetAMGFilename(mat, "fine", level, ".marks") << "\n";
		f2 << "c " << GetAMGFilename(mat, "aggfine", level, ".marks") << "\n";
		f2 << "c " << GetAMGFilename(mat, "coarse", level, ".marks") << "\n";
		f2 << "c " << GetAMGFilename(mat, "other", level, ".marks") << "\n";
		f2 << "c " << GetAMGFilename(mat, "dirichlet", level, ".marks") << "\n";
	}
}

template<typename TAlgebra>
void AMGBase<TAlgebra>::write_debug_matrices(matrix_type &AH, prolongation_matrix_type &R, const matrix_type &A,
		prolongation_matrix_type &P, size_t level)
{
	AMG_PROFILE_FUNC();
	UG_LOG("write matrices");

	//if(A.num_rows() < AMG_WRITE_MATRICES_MAX)
	{
		write_debug_matrix(A, level, level, "AMG_A_L");		UG_LOG(".");
		write_debug_matrix(P, level+1, level, "AMG_P_L");	UG_LOG(".");
		write_debug_matrix(R, level, level+1, "AMG_R_L");	UG_LOG(".");
	}

	//if(AH.num_rows() < AMG_WRITE_MATRICES_MAX)
	std::string name = m_writeMatrixPath + std::string("AMG_Au_L") + ToString(level+1) + ".mat";
#ifdef UG_PARALLEL
	name = GetParallelName(name, AH.process_communicator(), true);
#endif
	AMGWriteToFile(AH, level+1, level+1, name.c_str(), m_amghelper);

	UG_LOG(". done.\n");
}

template<typename TAlgebra>
bool AMGBase<TAlgebra>::injection(vector_type &vH, const vector_type &v, size_t level)
{
	size_t nH = levels[level]->P.num_cols(); (void) nH;
	size_t n = levels[level]->P.num_rows(); (void) n;
	UG_ASSERT(nH == vH.size() && n <= v.size(), nH << " != " << vH.size() << " or " << n << " > " << v.size())
	for(size_t j=0; j < nH; j++)
		vH[j] = v[m_parentIndex[level+1][j]];
	return true;
}


} // namespace ug

#include "amg_agglomeration.h"
#include "amg_checks.h"
#include "amg_communicate_prolongation.h"
#endif //  __H__UG__LIB_ALGEBRA__AMG_SOLVER__AMG_BASE_IMPL_H__

