/**
 * \file amg_base.h
 *
 * \author Martin Rupp
 *
 * \date 01.12.2010
 *
 * implementation file for base amg functionality
 *
 * Goethe-Center for Scientific Computing 2010.
 */

#ifndef __H__LIB_DISCRETIZATION__AMG_SOLVER__AMG_BASE_IMPL_H__
#define __H__LIB_DISCRETIZATION__AMG_SOLVER__AMG_BASE_IMPL_H__

#include "stopwatch.h"
#include "amg_debug.h"

namespace ug{

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//----------------
//! creates MG Hierachy for with matrix_type A and temporary vectors for higher levels
//! @param A	matrix A.
template<typename TAlgebra>
bool amg_base<TAlgebra>::preprocess(matrix_type& mat)
{
	cleanup();
	A.resize(1);
	A[0] = &mat;
	m_bInited = false;
	return true;
}

template<typename TAlgebra>
bool amg_base<TAlgebra>::create_level_vectors(size_t level)
{
	size_t N = A[level]->num_rows();
	size_t iNrOfCoarse = A[level+1]->num_rows();
	// create vectors for AMG multigrid
	/////////////////////////////////////////

	vec3.resize(level+1);
	vec3[level] = new vector_type;
	vec3[level]->create(N);

	vec1.resize(level+2);
	vec1[level+1] = new vector_type;
	vec1[level+1]->create(iNrOfCoarse);

	vec2.resize(level+2);
	vec2[level+1] = new vector_type;
	vec2[level+1]->create(iNrOfCoarse);

	#ifdef UG_PARALLEL
	// todo: change this for later "right" parallel implementation
	//typedef pcl::SingleLevelLayout<pcl::OrderedInterface<size_t, std::vector> > IndexLayout


	// todo: implement "real" ParallelCommunicator
	vec3[level]->set_communicator(*com);
	vec1[level+1]->set_communicator(*com);
	vec2[level+1]->set_communicator(*com);

	vec3[level]->set_storage_type(PST_ADDITIVE);
	vec1[level+1]->set_storage_type(PST_ADDITIVE);
	vec2[level+1]->set_storage_type(PST_ADDITIVE);

	vec3[level]->set_layouts(pseudoLayout, pseudoLayout);
	vec1[level+1]->set_layouts(pseudoLayout, pseudoLayout);
	vec2[level+1]->set_layouts(pseudoLayout, pseudoLayout);

	#endif


	//UG_LOG(std::endl << "created vec3 on level " << level << ", vec1 and vec2 on level" << level +1);

	// todo: set size for variable sized blockvectors
	/*for(size_t i=0; i<N; i++)
		if(nodes[i].isCoarse())
		{
			int rows = GetRows(A.beginRow(i).value());
			UG_ASSERT(newIndex[i] >= 0, "");
			SetSize((*vec1[level+1])[newIndex[i]], rows);
			SetSize((*vec2[level+1])[newIndex[i]], rows);
		}*/
	return true;
}

template<typename TAlgebra>
bool amg_base<TAlgebra>::init()
{
	if(m_basesolver==NULL)
	{
		UG_LOG("amg_base::init(): No base solver selected. Call set_base_solver(basesolver) to set a base solver.\n");
		return false;
	}
	if(m_presmoother==NULL)
	{
		UG_LOG("amg_base::init(): No PreSmoother selected. Call set_presmoother(presmoother) to set a PreSmoother.\n");
		return false;
	}
	if(m_postsmoother==NULL)
	{
		UG_LOG("amg_base::init(): No PostSmoother selected. Call set_postsmoother(postsmoother) to set a PostSmoother.\n");
		return false;
	}

	// init amghelper for grid printing
	amghelper.positions = &dbg_positions[0];
	amghelper.size = A[0]->num_rows();
	amghelper.parentIndex = &parentIndex;
	amghelper.dimension = dbg_dimension;

	UG_LOG("Starting AMG Setup." << std::endl << std::endl);

	stopwatch SWwhole;
	SWwhole.start();

	UG_ASSERT(max_levels > 0, "max levels is " << max_levels << ", but has to be > 0");

	size_t level=0;
	for(; level< max_levels-1; level++)
	{
		SMO.resize(level+1);
		SMO[level].setmatrix(A[level]);

		m_presmoothers.resize(level+1);
		m_presmoothers[level] = m_presmoother->clone();
		m_presmoothers[level]->init(SMO[level]);

		m_postsmoothers.resize(level+1);
		if(m_presmoother == m_postsmoother)
			m_postsmoothers[level] = m_presmoothers[level];
		else
		{
			m_postsmoothers[level] = m_postsmoother->clone();
			m_postsmoothers[level]->init(SMO[level]);
		}

		double L = A[level]->num_rows();

		if(L < m_maxNodesForExact || A[level]->total_num_connections()/(L*L) > m_maxFillBeforeExact)
			break;
		//smoother[level].init(*A[level]);

		A.resize(level+2);
		A[level+1] = new matrix_type;

		P.resize(level+1);
		P[level] = new matrix_type;
		R.resize(level+1);
		R[level] = new matrix_type;

#ifdef UG_PARALLEL
		A[level+1]->set_storage_type(PST_ADDITIVE);
		P[level]->set_storage_type(PST_ADDITIVE);
		R[level]->set_storage_type(PST_ADDITIVE);
#endif


		stopwatch SWwhole; SWwhole.start();
		create_AMG_level(*A[level+1], *R[level], *A[level], *P[level], level);
		return true;
		SWwhole.stop();

		is_fine.resize(level+1);
		is_fine[level].resize(A[level]->num_rows(), false);
		for(size_t i=0; i<parentIndex.size(); i++)
			is_fine[level][parentIndex[level+1][i]] = true;

		// finish
		/////////////////////////////////////////

		size_t nnzCoarse = A[level+1]->total_num_connections();
		double nrOfFine = A[level]->num_rows();
		double nrOfCoarse = A[level+1]->num_rows();
		UG_LOG("AH: nnz: " << nnzCoarse << " Density: " <<
				nnzCoarse/(nrOfCoarse*nrOfCoarse)*100.0 << "%, avg. nnz pre row: " << nnzCoarse/nrOfCoarse << std::endl);

		UG_LOG("Coarsening rate: " << (100.0*nrOfCoarse)/nrOfFine << "%" << std::endl);

		UG_LOG(" level took " << SWwhole.ms() << " ms" << std::endl << std::endl);


	/*
	 if(m_writeMatrices && this->A[0]->num_rows() < AMG_WRITE_MATRICES_MAX)
	{
		UG_LOG("write matrices");
		AMGWriteToFile(*P[level], level+1, level, (m_writeMatrixPath + "AMG_P" + ToString(level) + ".mat").c_str(), amghelper);
		UG_LOG(".");
		AMGWriteToFile(*R[level], level, level+1, (m_writeMatrixPath + "AMG_R" + ToString(level) + ".mat").c_str(), amghelper);
		UG_LOG(".");
		AMGWriteToFile(*A[level+1], level+1, level+1, (m_writeMatrixPath + "AMG_A" + ToString(level+1) + ".mat").c_str(), amghelper);
		UG_LOG(". done.\n");
	}*/

		create_level_vectors(level);
	}

	UG_ASSERT(block_traits< typename vector_type::value_type >::is_static, "dynamic not yet implemented");
	size_t static_nrUnknowns = block_traits< typename vector_type::value_type >::static_size;

	UG_LOG("Creating level " << level << " (" << A[level]->num_rows() << " nodes, total "
			<< A[level]->num_rows()*static_nrUnknowns << " unknowns)" << std::endl << "Using Direct Solver on Matrix "
			<< A[level]->num_rows()*static_nrUnknowns << "x" << A[level]->num_rows()*static_nrUnknowns << ". ");

	stopwatch SW; SW.start();

	SMO.resize(level+1);
	SMO[level].setmatrix(A[level]);
	m_basesolver->init(SMO[level]);

	UG_LOG("Coarse Solver Setup took " << SW.ms() << "ms." << std::endl);

	used_levels = level+1;
	UG_LOG("AMG Setup finished. Used Levels: " << used_levels << ". ");
	UG_LOG("AMG Setup took " << SWwhole.ms() << " ms." << std::endl);

	// calc complexities
	double nnzs=0;
	double totallength=0;
	for(size_t i=0; i<used_levels; i++)
	{
		nnzs += A[i]->total_num_connections();
		totallength += A[i]->num_rows();
	}

	UG_LOG("Operator Complexity: " << nnzs/A[0]->total_num_connections() << " nodes complexity: "
			<< totallength/A[0]->num_rows() << std::endl << std::endl);

	return true;
}


//!
//! amg constructor
template<typename TAlgebra>
amg_base<TAlgebra>::amg_base() :
	m_numPreSmooth(2),
	m_numPostSmooth(2),
	m_cycleType(1),

	max_levels(100),
	used_levels(0),

	m_presmoother(NULL),
	m_postsmoother(NULL),
	m_basesolver(NULL)
{
	vec4 = NULL;
	m_bInited = false;

	m_maxNodesForExact = 100;
	m_maxFillBeforeExact = 0.5;

	m_writeMatrices = false;
	m_fDamp = 0.0;
}


template<typename TAlgebra>
void amg_base<TAlgebra>::cleanup()
{
	for(size_t i=1; i < A.size(); i++) {
		if(A[i]) delete A[i]; A[i] = NULL;
	}


	for(size_t i=0; i < vec1.size(); i++) { if(vec1[i]) delete vec1[i]; vec1[i] = NULL; }
	for(size_t i=0; i < vec2.size(); i++) { if(vec2[i]) delete vec2[i]; vec2[i] = NULL; }
	for(size_t i=0; i < vec3.size(); i++) { if(vec3[i]) delete vec3[i]; vec3[i] = NULL; }

	for(size_t i=0; i < m_postsmoothers.size(); i++) {
		if(m_postsmoothers[i] && (i < m_presmoothers.size() && m_presmoothers[i] == m_postsmoothers[i]) == false)
			delete m_postsmoothers[i];
		m_postsmoothers[i] = NULL;
	}

	for(size_t i=0; i < m_presmoothers.size(); i++) {
		if(m_presmoothers[i])
			delete m_presmoothers[i];
		m_presmoothers[i] = NULL;
	}

	used_levels = 0;
	m_bInited=false;
}
//!
//! amg destructor
template<typename TAlgebra>
amg_base<TAlgebra>::~amg_base()
{
	cleanup();
}


template<typename TAlgebra>
bool amg_base<TAlgebra>::do_f_smoothing(vector_type &corr, vector_type &d, size_t level)
{
	matrix_type &M = *A[level];
	for(size_t i=0; i<M.num_rows(); i++)
	{
		//if(is_fine[level][i])
			corr[i] = m_fDamp*d[i]/M(i,i);
		/*else
			corr[i] = 0.0;*/
	}

	MatMultAdd(d, 1.0, d, -1.0, M, corr);
	return true;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// get_correction_and_update_defect:
//------------------------------------

template<typename TAlgebra>
bool amg_base<TAlgebra>::get_correction_and_update_defect(vector_type &c, vector_type &d, size_t level)
{
	return true;
	UG_ASSERT(c.size() == d.size() && c.size() == A[level]->num_rows(),
			"c.size = " << c.size() << ", d.size = " << d.size() << ", A.size = " << A[level]->num_rows() << ": not matching");

	const matrix_type &Ah = *(A[level]);

	if(level == used_levels-1)
	{
		m_basesolver->apply_return_defect(c, d);
		return true;
	}

	vector_type &corr = *vec3[level];

	// presmooth
	// same as setting c.set(0.0).
	m_presmoothers[level]->apply_update_defect(c, d);
	if(m_fDamp != 0.0) { do_f_smoothing(corr, d, level); c+=corr; }
	for(size_t i=1; i < m_numPreSmooth; i++)
	{
		m_presmoothers[level]->apply_update_defect(corr, d);
		c += corr;
		if(m_fDamp != 0.0) { do_f_smoothing(corr, d, level); c+=corr; }
	}

	vector_type &cH = *vec1[level+1];
	vector_type &dH = *vec2[level+1];

#ifdef UG_PARALLEL
	cH.set_storage_type(PST_CONSISTENT);

	// restrict defect
	// dH = R[level]*d;
	dH.set_storage_type(PST_ADDITIVE);
#endif

	MatMult(dH, 1.0, *R[level], d);

	// apply lmgc on coarser nodes
	if(level+1 == used_levels-1)
		get_correction_and_update_defect(cH, dH, level+1);
	else
		for(int i=0; i< m_cycleType; i++)
			get_correction_and_update_defect(cH, dH, level+1);

	// interpolate correction
	// corr = P[level]*cH
#ifdef UG_PARALLEL
	corr.set_storage_type(PST_ADDITIVE);
#endif

	MatMult(corr, 1.0, *P[level], cH);

	// add coarse grid correction to level correction
	// c += corr;
#ifdef UG_PARALLEL
	corr.change_storage_type(PST_CONSISTENT);
#endif
	VecScaleAdd(c, 1.0, c, 1.0, corr);

	//update defect
	// d = d - Ah*corr
	MatMultAdd(d, 1.0, d, -1.0, Ah, corr);

	// postsmooth
	for(size_t i=0; i < m_numPostSmooth; i++)
	{
		if(m_fDamp != 0.0) { do_f_smoothing(corr, d, level); c+=corr; }
		m_postsmoothers[level]->apply_update_defect(corr, d);
		c += corr;
	}

	return true;
}



template<typename TAlgebra>
bool amg_base<TAlgebra>::get_correction(vector_type &c, const vector_type &const_d)
{
	return true;
	UG_ASSERT(c.size() == const_d.size() && c.size() == A[0]->num_rows(),
				"c.size = " << c.size() << ", d.size = " << const_d.size() << ", A.size = " << A[0]->num_rows() << ": not matching");


	if(vec4 == NULL)
	{
		vec4 = new vector_type;
		vec4->resize(const_d.size());
	}
	else if(vec4->size() != const_d.size())
		vec4->resize(const_d.size());

#ifdef UG_PARALLEL
	// todo: change this for later "right" parallel implementation
	UG_ASSERT(pcl::GetNumProcesses() == 1, "AMG currently only for 1 process");
	vec4->set_storage_type(PST_ADDITIVE);
#endif
	vector_type &d = *vec4;

	d = const_d;
	return get_correction_and_update_defect(c, d);
}

template<typename TAlgebra>
bool amg_base<TAlgebra>::check(const vector_type &const_c, const vector_type &const_d)
{
	vector_type d;
	d.resize(const_d.size());
	d = const_d;
	vector_type d_copy;

	for(size_t i=0; i<used_levels-1; i++)
	{
		UG_LOG("\nLEVEL " << i << "\n\n");
		vector_type c;
		c.resize(A[i]->num_rows());

		d_copy.resize(d.size());
		d_copy = d;

#if 0
		d = 0.0;
		for(size_t j=0; j<c.size(); j++)
		{
			if(A[i]->is_isolated(i))
				c[j] = 0.0;
			else
				c[j] = urand(-1.0, 1.0);
		}
#endif
		check_level(c, d, i);

		if(i+1 < used_levels)
		{
			d.resize(A[i+1]->num_rows());
			MatMult(d, 1.0, *R[i], d_copy);
		}
	}
	return true;
}

template<typename vector_type>
void writevec(cAMG_helper &amghelper, const char *filename, const vector_type &d, size_t level)
{
	fstream file(filename, ios::out);
	for(size_t i=0; i<d.size(); i++)
		file << amghelper.GetOriginalIndex(level, i) << " " << (d[i]) << std::endl;
}


template<typename TAlgebra>
bool amg_base<TAlgebra>::check_level(vector_type &c, vector_type &d, size_t level)
{
	const matrix_type &Ah = *(A[level]);
	vector_type corr; corr.resize(c.size());


	//UG_LOG("preprenorm: " << d.two_norm() << std::endl);
	/*for(size_t i=0; i<5; i++)
		get_correction_and_update_defect(corr, d, level);*/

	double prenorm = d.two_norm();
	UG_LOG("Prenorm = " << prenorm << "\n");

	// presmooth
	// same as setting c.set(0.0).

	double n1 = d.two_norm(), n2;

	m_presmoothers[level]->apply_update_defect(c, d);
	if(m_fDamp != 0.0) { do_f_smoothing(corr, d, level); c+=corr; }
	n2 = d.two_norm();	UG_LOG("presmoothing 1 " << ": " << n2/prenorm << "\t" <<n2/n1 << "\n");	n1 = n2;
	for(size_t i=1; i < m_numPreSmooth; i++)
	{
		m_presmoothers[level]->apply_update_defect(corr, d);
		c += corr;
		if(m_fDamp != 0.0) { do_f_smoothing(corr, d, level); c+=corr; }

		n2 = d.two_norm();	UG_LOG("presmoothing " << i+1 << ": " << n2/prenorm << "\t" <<n2/n1 << "\n");	n1 = n2;
	}

	if(m_writeMatrices) writevec(amghelper, (m_writeMatrixPath + "AMG_dp" + ToString(level) + ".values").c_str(), d, level);

	vector_type &cH = *vec1[level+1];
	vector_type &dH = *vec2[level+1];

#ifdef UG_PARALLEL
	cH.set_storage_type(PST_CONSISTENT);

	// restrict defect
	// dH = R[level]*d;
	dH.set_storage_type(PST_ADDITIVE);
#endif

	MatMult(dH, 1.0, *R[level], d);

	vector_type tc; tc.resize(cH.size());
	cH.set( 0.0);

	double nH1 = dH.norm();
	double preHnorm=nH1;
	size_t i;
	for(i=0; i<100; i++)
	{
		// apply lmgc on coarser nodes
		if(level+1 == used_levels-1)
		{
			get_correction_and_update_defect(tc, dH, level+1);
			cH += tc;
		}
		else
			for(int i=0; i< m_cycleType; i++)
			{
				get_correction_and_update_defect(tc, dH, level+1);
				cH += tc;
			}

		double nH2 = dH.norm();
		if(i < 6)
		{	UG_LOG("coarse correction (on coarse) " << i+1 << ": " << nH2/nH1 << "\n"); nH1 = nH2; }
		if(nH2/preHnorm < 0.01) { UG_LOG("coarse solver reduced by 0.01 in iteration " << i << std::endl); break; }

	}

	// interpolate correction
	// corr = P[level]*cH
#ifdef UG_PARALLEL
	corr.set_storage_type(PST_ADDITIVE);
#endif

	MatMult(corr, 1.0, *P[level], cH);

	// add coarse grid correction to level correction
	// c += corr;
#ifdef UG_PARALLEL
	corr.change_storage_type(PST_CONSISTENT);
#endif
	VecScaleAdd(c, 1.0, c, 1.0, corr);

	//update defect
	// d = d - Ah*corr
	MatMultAdd(d, 1.0, d, -1.0, Ah, corr);

	n2 = d.two_norm();	UG_LOG("complete coarse correction " << i+1 << ": " << n2/prenorm << "\t" <<n2/n1 << "\n");	n1 = n2;

	if(m_writeMatrices) writevec(amghelper, (m_writeMatrixPath + "AMG_dc" + ToString(level) + ".values").c_str(), d, level);

	// postsmooth
	for(size_t i=0; i < m_numPostSmooth; i++)
	{
		if(m_fDamp != 0.0) { do_f_smoothing(corr, d, level); c+=corr; }
		m_postsmoothers[level]->apply_update_defect(corr, d);
		c += corr;
		n2 = d.two_norm();	UG_LOG("postsmoothing " << i+1 << ": " << n2/prenorm << "\t" <<n2/n1 << "\n");	n1 = n2;
	}

	double postnorm = d.two_norm();
	UG_LOG("Level " << level << " reduction: " << postnorm/prenorm << std::endl);

	if(m_writeMatrices) writevec(amghelper, (m_writeMatrixPath + "AMG_d" + ToString(level) + ".values").c_str(), d, level);

	UG_LOG("\n\n");
	return true;
}



template<typename TAlgebra>
void amg_base<TAlgebra>::tostring() const
{
	UG_LOG("AMGBase:\n");
	UG_LOG(" nr of pre smoothing steps (nu1)  = " << m_numPreSmooth<< std::endl);
	UG_LOG(" nr of post smoothing steps (nu2) = " << m_numPostSmooth << std::endl);

	UG_LOG(" cycle type = ");
	if(m_cycleType == 1) { UG_LOG("V-cycle\n"); }
	else if(m_cycleType == 2) { UG_LOG("W-cycle\n"); }
	else { UG_LOG("gamma = " << m_cycleType << "\n"); }

	UG_LOG(" max levels = " << max_levels << std::endl);

	if(m_presmoother) 	{UG_LOG(" presmoother is " << m_presmoother->name() << ".\n");}
	else				{UG_LOG(" no presmoother set!\n");}

	if(m_postsmoother) 	{UG_LOG(" postsmoother is " << m_postsmoother->name() << ".\n");}
	else				{UG_LOG(" no postsmoother set!\n");}

	if(m_basesolver)	{UG_LOG(" basesolver set\n");}
	else				{UG_LOG(" no basesolver set!\n");}
}

} // namespace ug

#endif //  __H__LIB_DISCRETIZATION__AMG_SOLVER__AMG_BASE_IMPL_H__
