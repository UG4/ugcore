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
#include "collect_matrix.h"
#endif

#include "lib_algebra/common/connection_viewer_output.h"

std::string GetProcFilename(std::string path, std::string name, std::string extension);

namespace ug{


#ifdef UG_PARALLEL
template<typename TMatrix, typename TVector>
void SetParallelVectorAsMatrix(TVector &v, TMatrix &m, ParallelStorageType t)
{
	v.set_communicator(m.get_communicator());
	v.set_process_communicator(m.get_process_communicator());
	v.set_layouts(m.get_master_layout(), m.get_slave_layout());
	v.set_storage_type(t);
}
#endif

template<typename TMatrix>
size_t GetMaxConnections(const TMatrix &A)
{
	size_t m=0;
	for(size_t i=0; i<A.num_rows(); i++)
		if(m < A.num_connections(i)) m = A.num_connections(i);
	return m;
}

template<typename TAlgebra>
void AMGBase<TAlgebra>::calculate_level_information(size_t level, double createAMGlevelTiming)
{
	AMGLevel &L = *levels[level];
	matrix_operator_type &A = *L.pA;

	LevelInformation &li = L.m_levelInformation;

	size_t nnz = A.total_num_connections();
	size_t maxConnections = GetMaxConnections(A);
#ifdef UG_PARALLEL
	size_t N = A.num_rows() - A.get_slave_layout().num_interface_elements();
	li.m_dCreationTimeMS = createAMGlevelTiming;
	li.set_nr_of_nodes(
			L.processCommunicator.allreduce(N, PCL_RO_MIN),
			L.processCommunicator.allreduce(N, PCL_RO_MAX),
			L.processCommunicator.allreduce(N, PCL_RO_SUM));
	li.set_nnz(
			L.processCommunicator.allreduce(nnz, PCL_RO_MIN),
			L.processCommunicator.allreduce(nnz, PCL_RO_MAX),
			L.processCommunicator.allreduce(nnz, PCL_RO_SUM));
	li.set_max_connections(L.processCommunicator.allreduce(maxConnections, PCL_RO_MAX));
	size_t localInterfaceElements = A.get_master_layout().num_interface_elements() + A.get_slave_layout().num_interface_elements();
	li.m_iInterfaceElements = L.processCommunicator.allreduce(localInterfaceElements, PCL_RO_SUM);
#else
	size_t N = A.num_rows();
	li.set_nr_of_nodes(N, N, N);
	li.set_nnz(nnz, nnz, nnz);
	li.set_max_connections(maxConnections);
	li.m_iInterfaceElements = 0;
#endif

	UG_LOG("nrOfCoarse: " << li.get_nr_of_nodes() << "\n");
	IF_DEBUG(LIB_ALG_AMG, 0)
	{
		UG_DLOG(LIB_ALG_AMG, 1, "AH: nnz: " << li.get_nnz() << " Density: " <<
				li.get_fill_in()*100.0 << "%, avg. nnz pre row: " <<
				li.get_avg_nnz_per_row() << ", max conns: " << li.get_max_connections() << std::endl);
		if(level>0)
		{
			UG_DLOG(LIB_ALG_AMG, 1, "Coarsening rate: " <<
				(100.0*li.get_nr_of_nodes())
				/ levels[level-1]->m_levelInformation.get_nr_of_nodes() <<
				"%" << std::endl);
			UG_DLOG(LIB_ALG_AMG, 1, " level took " << createAMGlevelTiming << " ms" << std::endl << std::endl);
		}
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//----------------
//! creates MG Hierachy for with matrix_operator_type A and temporary vectors for higher levels
//! @param A	matrix A.
template<typename TAlgebra>
bool AMGBase<TAlgebra>::preprocess(matrix_operator_type& mat)
{
	AMG_PROFILE_FUNC();

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

	if(m_basesolver==NULL)
	{
		UG_LOG("amg_base::preprocess(): No base solver selected. Call set_base_solver(basesolver) to set a base solver.\n");
		return false;
	}
	if(m_presmoother==NULL)
	{
		UG_LOG("amg_base::preprocess(): No PreSmoother selected. Call set_presmoother(presmoother) to set a PreSmoother.\n");
		return false;
	}
	if(m_postsmoother==NULL)
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
	AMGLevel *pL = new AMGLevel;
	pL->pA = &mat;
#ifdef UG_PARALLEL
	pL->processCommunicator = mat.get_process_communicator();
	pL->com = mat.get_communicator();
	pL->bHasBeenMerged = false;
#endif
	levels.push_back(pL);

	for(; ; level++)
	{
		AMGLevel &L = *levels[level];
		matrix_operator_type &A = *L.pA;

		calculate_level_information(level, createAMGlevelTiming);

		if(L.m_levelInformation.get_nr_of_nodes() < m_maxNodesForBase ||
				level >= m_maxLevels-1)
			// || A.total_num_connections()/(L*L) > m_dMaxFillBeforeBase)
		{
			AMG_PROFILE_BEGIN(amg_createDirectSolver);
			create_direct_solver(level);
			AMG_PROFILE_END();
			break;
		}
#ifdef UG_PARALLEL
		/*agglomerate(level);
		if(m_agglomerateLevel == level)
		{
			UG_LOG("got agglomerated.\n");
			break;
		}*/
#endif
		//smoothem_R[level].init(*m_A[level]);

		levels.push_back(new AMGLevel);
		AMGLevel &nextL = *levels[level+1];

		L.presmoother = m_presmoother->clone();
		L.presmoother->init(*L.pA);
		if(m_presmoother == m_postsmoother)
			L.postsmoother = L.presmoother;
		else
		{
			L.postsmoother = m_postsmoother->clone();
			L.postsmoother->init(*L.pA);
		}

		nextL.pA = new matrix_operator_type;
		matrix_operator_type &AH = *nextL.pA;
		prolongation_matrix_type &R = L.R;
		prolongation_matrix_type &P = L.P;

#ifdef UG_PARALLEL
		AH.set_storage_type(PST_ADDITIVE);
		R.set_storage_type(PST_ADDITIVE);
		P.set_storage_type(PST_ADDITIVE);

		nextL.com = L.com;
		nextL.processCommunicator = L.processCommunicator;

		AH.set_layouts(nextL.masterLayout, nextL.slaveLayout);
		AH.set_communicator(nextL.com);
		AH.set_process_communicator(nextL.processCommunicator);
#endif

		stopwatch SW; SW.start();

		create_AMG_level(AH, R, A, P, level);

		SW.stop();
		createAMGlevelTiming = SW.ms();

		// create vectors for AMG multigrid
		/////////////////////////////////////////

		L.corr.resize(A.num_rows());
		L.cH.resize(AH.num_rows());
		L.dH.resize(AH.num_rows());

	#ifdef UG_PARALLEL
		SetParallelVectorAsMatrix(L.corr, A, PST_CONSISTENT);
		SetParallelVectorAsMatrix(L.cH, AH, PST_ADDITIVE);
		SetParallelVectorAsMatrix(L.dH, AH, PST_ADDITIVE);
		/*UG_LOG("nextLevel Layouts:\n");
				PrintLayout(nextL.pA->get_communicator(), nextL.masterLayout, nextL.slaveLayout);*/
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
	if(A.get_process_communicator().size()>1)
	{
		L.bHasBeenMerged = true;
		// create basesolver
		matrix_operator_type collectedBaseA;
		collect_matrix(A, collectedBaseA, L.agglomerateMasterLayout, agglomerateSlaveLayout);

		if(m_writeMatrices && m_amghelper.has_positions())
		{
			std::vector<MathVector<3> > vec = m_amghelper.positions[level];
			vec.resize(collectedBaseA.num_rows());
			ComPol_VecCopy<std::vector<MathVector<3> > >	copyPol(&vec);
			pcl::ParallelCommunicator<IndexLayout> &communicator = A.get_communicator();
			communicator.send_data(agglomerateSlaveLayout, copyPol);
			communicator.receive_data(L.agglomerateMasterLayout, copyPol);
			communicator.communicate();
			if(pcl::GetProcRank() == L.processCommunicator.get_proc_id(0))
			{
				std::string name = (std::string(m_writeMatrixPath) + "collectedA.mat");
				WriteMatrixToConnectionViewer(name.c_str(), collectedBaseA, &vec[0], m_amghelper.dimension);
			}
		}
		if(pcl::GetProcRank() == L.processCommunicator.get_proc_id(0))
		{
			L.processCommunicator = L.processCommunicator.create_sub_communicator(true);
			collectedBaseA.set_master_layout(m_emptyLayout);
			collectedBaseA.set_slave_layout(m_emptyLayout);
			collectedBaseA.set_process_communicator(m_emptyPC);
			L.uncollectedA = A;
			A = collectedBaseA;
			m_basesolver->init(A);

			collectedBaseA.set_layouts(m_emptyLayout, m_emptyLayout);
			L.collC.resize(A.num_rows());
			L.collD.resize(A.num_rows());
			SetParallelVectorAsMatrix(L.collC, A, PST_ADDITIVE);
			SetParallelVectorAsMatrix(L.collD, A, PST_ADDITIVE);
		}
		else
		{
			m_agglomerateLevel = level;
			L.processCommunicator = L.processCommunicator.create_sub_communicator(false);
		}
	}
	else
		m_basesolver->init(A);
#else
	m_basesolver->init(A);
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
	m_minNodesOnOneProcessor = 100;
	m_dMaxFillBeforeBase = 0.5;

	m_writeMatrices = false;
	m_bFSmoothing = false;

	m_dbgDimension = 0;
	iteration_glboal=0;
}


template<typename TAlgebra>
void AMGBase<TAlgebra>::cleanup()
{
	for(size_t i=1; i < levels.size(); i++)
	{
		if(levels[i]->pA) delete levels[i]->pA;
		levels[i]->pA = NULL;
	}

	for(size_t i=0; i < levels.size(); i++)
	{
		if(levels[i]->postsmoother && levels[i]->postsmoother != levels[i]->presmoother)
			delete levels[i]->postsmoother;
		levels[i]->postsmoother = NULL;
		if(levels[i]->presmoother)
			delete levels[i]->presmoother;
	}

	m_usedLevels = 0;
	m_bInited=false;
}
//!
//! amg destructor
template<typename TAlgebra>
AMGBase<TAlgebra>::~AMGBase()
{
	cleanup();
}

template<typename TAlgebra>
void AMGBase<TAlgebra>::init_fsmoothing()
{
#ifdef UG_PARALLEL

	// last level is basesolver level
	for(size_t k=0; k<levels.size()-1; k++)
	{
		matrix_operator_type &mat = *levels[k]->pA;

		ParallelVector<Vector< typename matrix_operator_type::value_type > > m_diag;
		size_t size = mat.num_rows();
		levels[k]->m_diagInv.resize(size);
		m_diag.resize(size);

		m_diag.set_layouts(mat.get_master_layout(), mat.get_slave_layout());
		m_diag.set_communicator(mat.get_communicator());

		// copy diagonal
		for(size_t i = 0; i < m_diag.size(); ++i){
			m_diag[i] = mat(i, i);
		}

		//	make diagonal consistent
		m_diag.set_storage_type(PST_ADDITIVE);
		m_diag.change_storage_type(PST_CONSISTENT);

		for(size_t i = 0; i < m_diag.size(); ++i)
			GetInverse(levels[k]->m_diagInv[i], m_diag[i]);
	}
#endif
}


template<typename TAlgebra>
bool AMGBase<TAlgebra>::f_smoothing(vector_type &corr, vector_type &d, size_t level)
{
	AMGLevel &L = *levels[level];
	stdvector<bool> &is_fine = L.is_fine;
	matrix_operator_type &A = *L.pA;
#ifdef UG_PARALLEL

	UG_ASSERT(L.m_diagInv.size() == is_fine.size(), "fine markers do not match in size on level " << level);
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
	m_basesolver->apply_return_defect(c, d);
#if 0
	AMGLevel &L = *levels[level];
	matrix_operator_type &A = *L.pA;
#ifdef UG_PARALLEL
	if(A.get_process_communicator().size()>1)
	{
		vector_type collC;
		vector_type collD;
		if(pcl::GetProcRank() == 0)
		{
			size_t N = collectedBaseA.num_rows();
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
		}
		// send d -> collD
		ComPol_VecAdd<vector_type > compolAdd(&collD, &d);
		pcl::ParallelCommunicator<IndexLayout> &com = A.get_communicator();
		com.send_data(slaveColl, compolAdd);
		com.receive_data(masterColl, compolAdd);
		com.communicate();


		if(pcl::GetProcRank() == 0)
			m_basesolver->apply_return_defect(collC, collD);

		// send c -> collC
		ComPol_VecCopy<vector_type> compolCopy(&c, &collC);
		com.send_data(masterColl, compolCopy);
		com.receive_data(slaveColl, compolCopy);
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
bool AMGBase<TAlgebra>::add_correction_and_update_defect2(vector_type &c, vector_type &d, size_t level)
{
	if(level == m_usedLevels-1)
	{
		solve_on_base(c, d, level);
		return true;
	}

	AMGLevel &L = *levels[level];
	matrix_operator_type &A = *L.pA;

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
		UG_ASSERT(corr.has_storage_type(PST_CONSISTENT), "" );
		c += corr;
	}

	// pre f-smoothing
	if(m_bFSmoothing)
	{
		corr.set(0.0);
		f_smoothing(corr, d, level);
		UG_ASSERT(corr.has_storage_type(PST_CONSISTENT), "" );
		c+=corr;
	}


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
	// dH = R*d;
	L.R.apply(dH, d);

#if WRITEVEC_IN_SOLVER
	writevec((std::string("AMG_") + ToString(iteration_glboal++) + "_dH_L").c_str(), dH, level+1);
#endif

	// apply lmgc on coarser nodes

	if(level+1 == m_usedLevels-1)
		add_correction_and_update_defect(cH, dH, level+1);
	else
	{
		for(int i=0; i< m_cycleType; i++)
			add_correction_and_update_defect(cH, dH, level+1);
	}

	corr.set(0.0);
	//cH.set(0.0);
	// interpolate correction
	// corr = R*cH
	L.P.apply(corr, cH);

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
	UG_ASSERT(c.size() == const_d.size() && c.size() == levels[0]->pA->num_rows(),
				"c.size = " << c.size() << ", d.size = " << const_d.size() << ", A.size = " <<
					levels[0]->pA->num_rows() << ": not matching");

	if(m_vec4 == NULL)
		m_vec4 = new vector_type;

	m_vec4->resize(const_d.size());
#ifdef UG_PARALLEL
	m_vec4->set_layouts(c.get_master_layout(), c.get_slave_layout());
	m_vec4->set_storage_type(PST_ADDITIVE);
#endif
	vector_type &d = *m_vec4;


	d = const_d;

	if(m_usedLevels == 1)
	{
		solve_on_base(c, d, 0);
		return true;
	}
	else
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

	if(m_presmoother) 	{UG_LOG(" presmoother is " << m_presmoother->name() << ".\n");}
	else				{UG_LOG(" no presmoother set!\n");}

	if(m_postsmoother) 	{UG_LOG(" postsmoother is " << m_postsmoother->name() << ".\n");}
	else				{UG_LOG(" no postsmoother set!\n");}

	if(m_basesolver)	{UG_LOG(" basesolver set\n");}
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
	UG_ASSERT(m_pPositionProvider2d == NULL || m_pPositionProvider3d == NULL, "specify EITHER positions 2d or 3d");
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



template<typename TAlgebra>
template<typename TNodeType>
void AMGBase<TAlgebra>::write_debug_matrix_markers
	(size_t level, const TNodeType &nodes)
{
	// todo: replace some day with something like nodes.get_mark_count(), nodes.mark_name(i), nodes.is_marked(i)
	AMG_PROFILE_FUNC();
	std::fstream ffine(GetProcFilename(m_writeMatrixPath, std::string("AMG_fine_L") + ToString(level), ".marks").c_str(), std::ios::out);
	ffine << "1 0 0 1 0\n";
	std::fstream ffine2(GetProcFilename(m_writeMatrixPath, std::string("AMG_aggfine_L") + ToString(level), ".marks").c_str(), std::ios::out);
	ffine2 << "1 0.2 1 1 0\n";
	std::fstream fcoarse(GetProcFilename(m_writeMatrixPath, std::string("AMG_coarse_L") + ToString(level), ".marks").c_str(), std::ios::out);
	fcoarse << "0 0 1 1 2\n";
	std::fstream fother(GetProcFilename(m_writeMatrixPath, std::string("AMG_other_L") + ToString(level), ".marks").c_str(), std::ios::out);
	fother << "1 1 0 1 0\n";
	std::fstream fdirichlet(GetProcFilename(m_writeMatrixPath, std::string("AMG_dirichlet_L") + ToString(level), ".marks").c_str(), std::ios::out);
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


#ifdef UG_PARALLEL
inline void WritePMAT(std::string &path, std::string &name)
{
	std::string pmatname = path + name + ".pmat";
	std::fstream file(pmatname.c_str(), std::ios::out);
	file << pcl::GetNumProcesses() << "\n";
	for(int i=0; i<pcl::GetNumProcesses(); i++)
		file << name << "_" << i << ".mat\n";
}
#else
inline void WritePMAT(std::string &path, std::string &name)
{
}
#endif


template<typename TAlgebra>
template<typename TMatrix>
void AMGBase<TAlgebra>::
	write_debug_matrix(TMatrix &mat, size_t fromlevel, size_t tolevel, const char *name)
{
	if(m_writeMatrices)
	{
		AMG_PROFILE_FUNC();
		std::string filename = name + ToString(fromlevel);
		WritePMAT(m_writeMatrixPath, filename);

		int level = std::min(fromlevel, tolevel);
		filename = GetProcFilename(m_writeMatrixPath, filename,".mat");
		AMGWriteToFile(mat, fromlevel, tolevel, filename.c_str(), m_amghelper);
		std::fstream f2(filename.c_str(), std::ios::out | std::ios::app);
		f2 << "c " << GetProcFilename("", std::string("AMG_fine_L") + ToString(level), ".marks") << "\n";
		f2 << "c " << GetProcFilename("", std::string("AMG_aggfine_L") + ToString(level), ".marks") << "\n";
		f2 << "c " << GetProcFilename("", std::string("AMG_coarse_L") + ToString(level), ".marks") << "\n";
		f2 << "c " << GetProcFilename("", std::string("AMG_other_L") + ToString(level), ".marks") << "\n";
		f2 << "c " << GetProcFilename("", std::string("AMG_dirichlet_L") + ToString(level), ".marks") << "\n";
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
		AMGWriteToFile(AH, level+1, level+1, GetProcFilename(m_writeMatrixPath, ToString("AMG_A_L") + ToString(level+1),".mat").c_str(), m_amghelper);

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

