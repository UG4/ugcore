/**
 * \file amg_checks.h
 *
 * \author Martin Rupp
 *
 * \date 03.08.11
 *
 * Goethe-Center for Scientific Computing 2011
 *
 *
 * in this file are functions for checking the convergence behaviour on the AMG levels (for bugfixing etc.)
 */
#include <sstream>

#include <sys/stat.h>

namespace ug
{

#ifdef UG_PARALLEL
template<typename TAlgebra>
void AMGBase<TAlgebra>::write_interfaces()
{
	AMG_PROFILE_FUNC();
	if(!m_writeMatrices)
	{
		UG_LOG("no position data given.\n");
		return;
	}
	cAMG_helper &h = m_amghelper;
	for(size_t level=0; level<m_usedLevels; level++)
	{

		//	\TODO: Implement also Windows support. Comment in below afterwards
//		std::string path=std::string("/level") + ToString(level) + "/";
//		mkdir((std::string(m_writeMatrixPath) + path).c_str(), 0777);

		const char *filename = (/*std::string(m_writeMatrixPath) + path + */std::string("AMG_interface_L") + ToString(level) + "_" + ToString(pcl::GetProcRank()) + std::string(".mat")).c_str();
		std::fstream file(filename, std::ios::out);
		file << 1 << std::endl; // connection viewer version


		WritePositionsToStream(file, h.positions[level], h.dimension);

		file << 1 << std::endl;

		IndexLayout &masterLayout = levels[level]->pA->master_layout();
		IndexLayout &slaveLayout = levels[level]->pA->slave_layout();
		std::vector<std::string> v;
		v.resize(h.positions[level].size());


		//UG_LOG("==============================\n LEVEL "<< level << "\n============================\n")
		//PRINTLAYOUT(levels[level]->com, levels[level]->masterLayout, levels[level]->slaveLayout);
		//PRINTLAYOUT(levels[level]->com, masterLayout, slaveLayout);

		for(IndexLayout::iterator iter = masterLayout.begin(); iter != masterLayout.end(); ++iter)
		{
			IndexLayout::Interface &interface = masterLayout.interface(iter);
			int pid = masterLayout.proc_id(iter);
			for(IndexLayout::Interface::iterator iter = interface.begin(); iter != interface.end(); ++iter)
			{
				size_t i = interface.get_element(iter);
				v[i].append(std::string("M") + ToString(pid) + std::string(" "));
			}
		}

		for(IndexLayout::iterator iter = slaveLayout.begin(); iter != slaveLayout.end(); ++iter)
		{
			IndexLayout::Interface &interface = slaveLayout.interface(iter);
			int pid = slaveLayout.proc_id(iter);
			for(IndexLayout::Interface::iterator iter = interface.begin(); iter != interface.end(); ++iter)
			{
				size_t i = interface.get_element(iter);
				v[i].append(std::string("S") + ToString(pid) + std::string(" "));
			}
		}

		for(size_t i=0;i<v.size(); i++)
		{
			std::string &s = v[i];
			if(s.size() > 0)
				file << i << " " << i << " " << s << std::endl;
		}
	}
}
#else
template<typename TAlgebra>
void AMGBase<TAlgebra>::write_interfaces()
{
	return;
}
#endif



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// writevec
//--------------------------------
template<typename TAlgebra>
bool AMGBase<TAlgebra>::writevec(std::string filename, const vector_type &const_d, size_t level, const vector_type *solution)
{
	AMG_PROFILE_FUNC();
	UG_ASSERT(m_writeMatrices, "");
	vector_type solution2;

#ifdef UG_PARALLEL
	vector_type d;
	CloneVector(d, const_d);
	d = const_d;
	d.change_storage_type(PST_CONSISTENT);
	size_t pid = pcl::GetProcRank();

	if(solution)
	{
		CloneVector(solution2, *solution);
		solution2 = *solution;
		solution2.change_storage_type(PST_CONSISTENT);
	}
#else
	const vector_type &d = const_d;
	size_t pid = 0;
#endif

//	std::string path=std::string("/level") + ToString(level) + "/";
//	mkdir((std::string(m_writeMatrixPath) + path).c_str(), 0777);

	std::string name = (/*m_writeMatrixPath + path +*/ std::string("AMG_L") + ToString(level) + "_" + filename + ".mat");



#ifdef UG_PARALLEL
	name = GetParallelName(name, levels[level]->pA->process_communicator(), true);
	/*if(levels[level]->bHasBeenMerged && m_agglomerateLevel != level)
		AMGWriteToFile(levels[level]->uncollectedA, level, level, name.c_str(), m_amghelper);
	else*/
#endif
		AMGWriteToFile(*levels[level]->pA, level, level, name.c_str(), m_amghelper);
	std::fstream f(name.c_str(), std::ios::out | std::ios::app);
	//UG_LOG("writevec to " << name.c_str() << "\n");

	std::string name2 = (std::string("AMG_L") + ToString(level) + "_" + filename + "_p"+ ToString(pid) + ".values");
#ifdef UG_PARALLEL
	name2 = GetParallelName(name2, levels[level]->pA->process_communicator(), false);
#endif
	f << "v " << name2 << "\n";

	std::fstream file((/*std::string(m_writeMatrixPath)+ path +*/ name2).c_str(), std::ios::out);
	if(solution)
		for(size_t i=0; i<d.size(); i++)
			file << i << " " << d[i] - solution2[i] << std::endl;
	else
		for(size_t i=0; i<d.size(); i++)
			file << i << " " << d[i] << std::endl;
	return true;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// check
//--------------------------------
/**
 * this function checks convergence level for level by using \sa check_level
 *
 * \param 	const_c				c
 * \param	const_d				d
 */
template<typename TAlgebra>
bool AMGBase<TAlgebra>::check(const vector_type &const_c, const vector_type &const_d)
{
	AMG_PROFILE_FUNC();
	UG_LOG("\n");
	UG_LOG("         check     \n");
	UG_LOG("========================\n");

	/*if(m_usedLevels <= 1)
	{
		UG_LOG("No Multigrid hierachy.\n");
		return true;
	}*/
	vector_type c, d;
	CloneVector(c, const_c);
	CloneVector(d, const_d);
	c = const_c;
	d = const_d;

	//UG_LOG("\nLEVEL 0\n\n");
#ifdef UG_PARALLEL
	UG_LOG("agglomerate Level is " << m_agglomerateLevel << ", used levels is " << m_usedLevels << "\n");
#endif

	std::vector<checkResult> checkRes(m_usedLevels-1);
	for(size_t level=0; level<m_usedLevels; level++)
	{
		vector_type *pC, *pD;

		if(level==0)
		{
			pC = &c;
			pD = &d;
		}
		else
		{
			pC = &levels[level-1]->cH;
			pD = &levels[level-1]->dH;
			pC->set(0.0); // think about that
		}

		SmartPtr<matrix_operator_type> pA = levels[level]->pA;

#ifdef UG_PARALLEL
		if(levels[level]->bLevelHasMergers)
		{
			UG_LOG("level " << level << " has mergers\n");
			try
			{
			gather_vertical(*pD, levels[level]->collD, level, PST_ADDITIVE);

			gather_vertical(*pC, levels[level]->collC, level, PST_CONSISTENT);
			}
			catch(...)
			{
				UG_LOG("catch???");
			}
			if(isMergingSlave(level))
			{
				UG_LOG("breaking, because i'm slave on level " << level << "\n");
				break;
			}
			pD = &levels[level]->collD;
			pC = &levels[level]->collC;
			pA = levels[level]->collectedA;
		}
#endif
#ifdef UG_PARALLEL
		if(m_agglomerateLevel==(unsigned int)(-1) && level == m_usedLevels-1)
		{
			UG_LOG("breaking2 " << level << "\n");
			break;
		}
#endif

		vector_type d2;
		CloneVector(d2, *pD);
		d2 = *pD;

		UG_LOG("checking level " << level << "\n");
		check_level(*pC, *pD, *pA, level, checkRes[level]);


		levels[level]->R.apply(levels[level]->dH, d2);
		//levels[i]->R.apply(levels[i]->cH, c2); // other possiblity, think about that
	}

	UG_LOG("out\n");
#ifdef UG_PARALLEL
//	const_c.process_communicator().barrier();
#endif
	return true;
}




template<typename TVector>
double ConstTwoNorm(const TVector &const_d)
{
	TVector d;
	CloneVector(d, const_d);
	d = const_d;
	return d.norm();

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// check_level:
//--------------------------------
/**
 * this function checks convergence on a specified level.
 * it will log all convergence rates:
 * - presmoother
 * - fsmoother
 * - coarse correction
 * - f-postsmoother
 * - postsmoother
 *
 * as coarse correction we use the normal multigrid on a level below as often as necessary
 * to reduce the error on the coarse grid by 0.01 to distinguish between bad convergence on this level
 * and bad convergence on other levels
 *
 * \param 	c
 * \param	d
 * \param	level
 */
template<typename TAlgebra>
bool AMGBase<TAlgebra>::check_level(vector_type &c, vector_type &d, matrix_type &A,
		size_t level, checkResult &res, const vector_type *solution)
{
	AMG_PROFILE_FUNC();
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

	std::stringstream reason;
	//UG_LOG("preprenorm: " << d.norm() << std::endl);
	vector_type &corr = levels[level]->corr;

	double firstnorm = ConstTwoNorm(d), n2=firstnorm, n1=n2;
	//UG_LOG("firstnorm " <<  firstnorm << "\n");
	double prevrate=0;
	if(m_writeMatrices) writevec("cc0", c, level, solution);
	if(m_writeMatrices) writevec("c0", c, level);
	if(m_writeMatrices) writevec("d0", d, level);

	UG_LOG("Type    " << "  Iter#     Defect         Rate\n");
	UG_LOG("                 " << std::scientific << firstnorm << "\n");
	for(size_t i=0; i<m_iNrOfPreiterationsCheck; i++)
	{
		add_correction_and_update_defect2(c, d, *L.pAgglomeratedA, level);

		if(m_writeMatrices) writevec(std::string("cc")+ToString(i+1), c, level, solution);
		if(m_writeMatrices) writevec(std::string("c")+ToString(i+1), c, level);
		if(m_writeMatrices) writevec(std::string("d")+ToString(i+1), d, level);

		n1 = n2;
		n2 = ConstTwoNorm(d);

		if(n2 < m_dPreiterationsMimumDefect|| (prevrate/ (n2/n1) < 1.01 && prevrate/ (n2/n1) > 0.99) )
			break;

		prevrate = n2/n1;

		UG_LOG("preIter " << std::setw(4) << i+1 << ":    " << std::scientific << n2 <<
				"    " << std::scientific << n2/n1 << "\n");

	}




	corr.set(0.0);

	for(size_t iterations=0; iterations<1; iterations++)
	{

		double prenorm = ConstTwoNorm(d);
			UG_LOG("defect below is in relation to " << std::scientific << prenorm << ".\n");
			UG_LOG("Type    " << "  Iter#     Defect         Rate            RelDefect \n");
			//UG_LOG("Prenorm = " << prenorm << "\n");

			n1 = ConstTwoNorm(d);
			n2=n1;

			if(n1 < 1e-13) reason << "- error too small\n";

	//c.set(0.0);
#ifdef UG_PARALLEL
	corr.set_storage_type(PST_CONSISTENT);
#endif

	if(m_writeMatrices) writevec(ToString(iterations)+"S0_d_before", d, level);
	if(m_writeMatrices) writevec(ToString(iterations)+"S0_c_before", c, level, solution);

	// presmooth
	for(size_t i=0; i < m_numPreSmooth; i++)
	{
		L.presmoother->apply_update_defect(corr, d);
		c += corr;
		n2 = ConstTwoNorm(d);
		UG_LOG("preSm   " << std::setw(4) << i+1 << ":    " << std::scientific << n2 <<
						"    " << std::scientific << n2/n1 <<
						"    " << std::scientific << n2/prenorm << "\n");
		n1 = n2;
	}

	res.preSmoothing = n2/prenorm;

	if(m_writeMatrices) writevec(ToString(iterations)+"S1_d_presmoothed", d, level);
	if(m_writeMatrices) writevec(ToString(iterations)+"S1_c_presmoothed", c, level, solution);

	if(m_bFSmoothing)
	{
		f_smoothing(corr, d, level);
		c+=corr;
		n2 = ConstTwoNorm(d);
		UG_LOG("preFSm  " << std::setw(4) << 0 << ":    " << std::scientific << n2 <<
								"    " << std::scientific << n2/n1 <<
								"    " << std::scientific << n2/prenorm << "\n");
		n1 = n2;
	}

	res.preFSmoothing = n2/prenorm;

	if(m_writeMatrices) writevec(ToString(iterations)+"S1_d_fpresmoothed", d, level);
	if(m_writeMatrices) writevec(ToString(iterations)+"S1_c_fpresmoothed", c, level, solution);

	vector_type &cH = levels[level]->cH;
	vector_type &dH = levels[level]->dH;
	cH.set(0.0);
	dH.set(0.0);
#ifdef UG_PARALLEL
	cH.set_storage_type(PST_CONSISTENT);
#endif
	// restrict defect
	// dH = m_R[level]*d;
	L.R.apply(dH, d);

	if(m_writeMatrices) writevec(ToString(iterations)+"S2_dH", dH, level+1);
	double nH1 = ConstTwoNorm(dH);
	UG_LOG("norm of dH = " << nH1 << "\n");

	double preHnorm=nH1;
	size_t i=0;

#ifdef UG_PARALLEL
	if(level+1 == m_usedLevels-1 && level+1 != m_agglomerateLevel)
#else
	if(level+1 == m_usedLevels-1)
#endif
	{
#if 0
		// tried to check basesolver here

				AMGLevel &L = *levels[level+1];
				matrix_operator_type &A = *L.pA;
				pcl::InterfaceCommunicator<IndexLayout> &com = A.get_communicator();

				if(m_agglomerateLevel == level+1)
				{
					// send d -> collD
					ComPol_VecAdd<vector_type > compolAdd(&dH, &dH);
					com.send_data(agglomerateSlaveLayout, compolAdd);
					com.communicate();

					ComPol_VecCopy<vector_type> compolCopy(&cH, &cH);
					com.receive_data(agglomerateSlaveLayout, compolCopy);
					com.communicate();

					c.set_storage_type(PST_CONSISTENT);

					A.matmul_minus(dH, cH);
					return true;
				}
				else
				{
					L.collD.set(0.0);
					for(size_t i=0; i<dH.size(); i++)
						L.collD[i] = dH[i];

					L.collC.set_storage_type(PST_CONSISTENT);
					L.collD.set_storage_type(PST_ADDITIVE);
					// send d -> collD
					ComPol_VecAdd<vector_type > compolAdd(&L.collD, &dH);
					com.receive_data(L.agglomerateMasterLayout, compolAdd);
					com.communicate();

					L.collC.set(0.0);
					add_correction_and_update_defect2(L.collC, L.collD, level+1);

					// send collC -> c
					ComPol_VecCopy<vector_type> compolCopy(&cH, &L.collC);
					com.send_data(L.agglomerateMasterLayout, compolCopy);
					com.communicate();

					for(size_t i=0; i<cH.size(); i++)
						cH[i] = L.collC[i];
					cH.set_storage_type(PST_CONSISTENT);
					A.matmul_minus(dH, cH); // cannot use collD, because collD is not additive.
					return true;
				}

				double nH2 = ConstTwoNorm(dH);
				UG_LOG("base solver reduced by " << nH2/nH1 << " (on coarse)" << std::endl);
#else
		add_correction_and_update_defect(cH, dH, level+1);
#endif
		double nH2 = ConstTwoNorm(dH);
		UG_LOG("base solver reduced by " << nH2/nH1 << " (on coarse)" << std::endl);

		res.coarseDefect = nH2;
		res.lastCoarseReduction = nH2/nH1;
	}
	else
	{
		cH.set(0.0);
		double nH2=nH1;
		for(i=0; i<10; i++)
		{
			for(int j=0; j< m_cycleType; j++)
				add_correction_and_update_defect(cH, dH, level+1);
			nH1 = nH2;
			nH2 = ConstTwoNorm(dH);
			/*if(i == 0)
			{	UG_LOG("coarse correction (on coarse) " << i+1 << ": " << nH2/nH1 << " " << nH2/preHnorm << "\n"); }
			else UG_LOG(".");*/
			//if(m_writeMatrices) writevec(ToString(iterations)+(std::string("AMG_V")+ToString(i)+std::string("_c_L")).c_str(), cH, level+1);
			//if(m_writeMatrices) writevec(ToString(iterations)+(std::string("AMG_V")+ToString(i)+std::string("_d_L")).c_str(), dH, level+1);
			if(nH2/preHnorm < 1e-8 || nH2 < 1e-16)
				break;
		}
		res.iInnerIterations = i;
		res.coarseDefect = nH2;
		res.lastCoarseReduction = nH2/nH1;

		UG_LOG("\ncoarse correction (on coarse) " << i+1 << ": " << nH2/nH1 << " " << nH2 <<  "\n");


		/*cH.set(0.0);
		for(i=0; i<100; i++)
		{
			for(int j=0; j< m_cycleType; j++)
				add_correction_and_update_defect(cH, dH, level+1);

			double nH2 = ConstTwoNorm(dH);
			if(i < 6)
			{	UG_LOG("coarse correction (on coarse) " << i+1 << ": " << nH2/nH1 << "\n"); nH1 = nH2; }
			if(nH2/preHnorm < 0.01) { UG_LOG("coarse solver reduced by 0.01 in iteration " << i+1 << std::endl); break; }
		}*/
	}

	// interpolate correction
	// corr = m_P[level]*cH
	if(m_writeMatrices) writevec(ToString(iterations)+"S3_dH_solved", dH, level+1);
	if(m_writeMatrices) writevec(ToString(iterations)+"S3_cH_solved", cH, level+1);
	L.P.apply(corr, cH);
	if(m_writeMatrices) writevec(ToString(iterations)+"S3_corr_on_fine", corr, level);

#ifdef UG_PARALLEL
	cH.set_storage_type(PST_CONSISTENT);
#endif

	// add coarse grid correction to level correction

	c += corr;

	//update defect
	// d = d - Ah*corr
	A.matmul_minus(d, corr);

	n2 = ConstTwoNorm(d);
	UG_LOG("complete coarse correction " << i+1 << ": " << n2/prenorm << "\t" <<n2/n1 << "\n");	n1 = n2;
	res.coarseCorrection = n2/prenorm;

	if(m_writeMatrices) writevec(ToString(iterations)+"S4_d_coarseCorrected", d, level);
	if(m_writeMatrices) writevec(ToString(iterations)+"S4_c_coarseCorrected", c, level, solution);


	// post f-smoothing
	if(m_bFSmoothing)
	{
		f_smoothing(corr, d, level);
		c+=corr;
	}
	n2 = ConstTwoNorm(d);

	UG_LOG("postFSm " << std::setw(4) << 0 << ":    " << std::scientific << n2 <<
									"    " << std::scientific << n2/n1 <<
									"    " << std::scientific << n2/prenorm << "\n");
	res.postFSmoothing = n2/prenorm;

	if(m_writeMatrices) writevec(ToString(iterations)+"S5_d_fpost", d, level);
	if(m_writeMatrices) writevec(ToString(iterations)+"S5_c_fpost", c, level, solution);

	// postsmooth
	for(size_t i=0; i < m_numPostSmooth; i++)
	{
		L.postsmoother->apply_update_defect(corr, d);
		c += corr;
		n2 = ConstTwoNorm(d);
		UG_LOG("postSm  " << std::setw(4) << i+1 << ":    " << std::scientific << n2 <<
								"    " << std::scientific << n2/n1 <<
								"    " << std::scientific << n2/prenorm << "\n");
	}

	res.postSmoothing = n2/prenorm;


	if(res.postSmoothing > 0.1)
	{
		UG_LOG("\nERROR: Level " << level << " not working: total correction is " << res.postSmoothing << "\n");

		if(res.preSmoothing > 0.99)
		{	UG_LOG("- presmoother not working (reduction rate " << res.preSmoothing << ")\n"); }
		if(res.lastCoarseReduction > 0.2 || res.coarseDefect > 1e-10)
		{
			UG_LOG("- lower level not working (done " << res.iInnerIterations << " iterations):\n");
			if(res.lastCoarseReduction > 0.2) { UG_LOG(" - bad convergence rate (" << res.lastCoarseReduction << ")\n"); }
			if(res.coarseDefect > 1e-8) { UG_LOG(" - coarse reduction not sufficient (" << res.coarseDefect << ")\n"); }
		}
	}


	double postnorm = ConstTwoNorm(d);
	UG_LOG("Level " << level << " reduction: " << postnorm/prenorm << std::endl);

	if(m_writeMatrices) writevec(ToString(iterations)+"S5_d_postsmoothed", d, level);
	if(m_writeMatrices) writevec(ToString(iterations)+"S5_c_postsmoothed", c, level, solution);

	UG_LOG("Postnorm = " << postnorm << "\n\n");
	res.reduction = postnorm/prenorm;
	}


	//return postnorm/prenorm;
	return 0;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// check2
//--------------------------------
/**
 * for all i= 1 .. usedLevels-2
 * 	we check the convergence behaviour of the multigrid level if, on level i, we use an exact solver
 * \param 	const_c				c
 * \param	const_d				d
 */
template<typename TAlgebra>
bool AMGBase<TAlgebra>::check2(const vector_type &const_c, const vector_type &const_d)
{
	AMG_PROFILE_FUNC();
	UG_LOG("\n");
	UG_LOG("            check2\n");
	UG_LOG("==================================\n");
	UG_LOG("we check the convergence behaviour of the multigrid level if, on level i, we use an exact solver");

	if(m_usedLevels <= 1)
	{
		UG_LOG("No Multigrid hierachy.\n");
		return true;
	}

	vector_type c, d;
	CloneVector(c, const_c);
	CloneVector(d, const_d);

	for(size_t i=1; i<m_usedLevels-1; i++)
	{
		c = const_c;
		d = const_d;
		UG_LOG("Testing Level 0 to i=" << i << ":\n");
		double prenorm = ConstTwoNorm(d);
		add_correction_and_update_defect(c, d, 0, i);
		double postnorm = ConstTwoNorm(d);
		UG_LOG("Reduction is " << postnorm/prenorm << "\n\n");
	}

	UG_LOG("finished check2.\n");

	return true;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// add_correction_and_update_defect
//-----------------------------------
/**
 * \sa add_correction_and_update_defect
 * for all i= 1 .. usedLevels-2
 * 	we check the convergence behaviour of the multigrid level if, on level i, we use an exact solver
 * \param c
 * \param d
 * \param level
 * \param exactLevel	level on which we will use an exact solver, that is: doing MG on that level until norm(d) < 1e-12.
 */
template<typename TAlgebra>
bool AMGBase<TAlgebra>::add_correction_and_update_defect(vector_type &c, vector_type &d, size_t level, size_t exactLevel)
{
	AMG_PROFILE_FUNC();
	AMGLevel &L = *levels[level];
	const matrix_operator_type &A = *L.pA;
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
#ifdef UG_PARALLEL
	corr.set_storage_type(PST_CONSISTENT);
#endif

	// presmooth
	for(size_t i=0; i < m_numPreSmooth; i++)
	{
		L.presmoother->apply_update_defect(corr, d);
		c += corr;
	}

	// pre f-smoothing
	if(m_bFSmoothing)
	{
		f_smoothing(corr, d, level);
		c+=corr;
	}

	vector_type &cH = levels[level]->cH;
	vector_type &dH = levels[level]->dH;
	cH.set(0.0);
	dH.set(0.0);
#ifdef UG_PARALLEL
	cH.set_storage_type(PST_CONSISTENT);
#endif

	L.R.apply(dH, d);

	if(level+1 == m_usedLevels-1)
		solve_on_base(cH, dH, level+1);
	else
	{
		cH.set(0.0);
		if(level+1 < exactLevel)
		{
			for(int i=0; i< m_cycleType; i++)
				add_correction_and_update_defect(cH, dH, level+1, exactLevel);
		}
		else
		{
			UG_LOG("using normal MG cycle from level " << level+1 << "on: ");
			int k;
			for(k=0; k<20; k++)
			{
				for(int i=0; i< m_cycleType; i++)
					add_correction_and_update_defect(cH, dH, level+1);
				UG_LOG(".");
				if(ConstTwoNorm(dH) < 1e-12)
					break;
			}
			if(ConstTwoNorm(dH) > 1e-12) { UG_LOG("not converged after " << k << " steps\n"); }
			else { UG_LOG("took " << k << " steps\n") }
		}
	}

	L.P.apply(corr, cH);
	c += corr;
	A.matmul_minus(d, corr);

	// post f-smoothing
	if(m_bFSmoothing)
	{
		f_smoothing(corr, d, level);
		c+=corr;
	}

	// postsmooth
	for(size_t i=0; i < m_numPostSmooth; i++)
	{
		L.postsmoother->apply_update_defect(corr, d);
		c += corr;
	}

	return true;
}



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// check_fsmoothing
//--------------------------------
/**
 */
template<typename TAlgebra>
bool AMGBase<TAlgebra>::check_fsmoothing()
{
	AMG_PROFILE_FUNC();
	UG_LOG("\n");
	UG_LOG("            check_fsmoothing\n");
	UG_LOG("==========================================\n")
	UG_LOG("check the convergence behaviour of a pure f-smoother")

	if(m_usedLevels == 0) return false;

	for(size_t level=0; level<m_usedLevels-1; level++)
	{
		UG_LOG("LEVEL " << level << "\n");
		vector_type c, d;
		CloneVector(c, levels[level]->corr);
		CloneVector(d, levels[level]->corr);

		for(size_t j=0; j<c.size(); j++)
		{
			if(levels[level]->is_fine[j])
				c[j] = urand<double>(-1.0, 1.0);
			else
				c[j] = 0.0;
		}

#ifdef UG_PARALLEL
		c.set_storage_type(PST_ADDITIVE);
		c.change_storage_type(PST_CONSISTENT);
#endif
		// d = A*c
		levels[level]->pA->apply(d, c);
		c.set(0.0);
		double n=ConstTwoNorm(d);

		for(size_t j=0; j<10; j++)
		{
			f_smoothing(c, d, level);
			double n2 = ConstTwoNorm(d);
			UG_LOG(j << ": " << n2/n << "\n");
			n = n2;
			if(m_writeMatrices) writevec(std::string("fsmoothing_c_") + ToString(j), c, level);
			if(m_writeMatrices) writevec(std::string("fsmoothing_d_") + ToString(j), d, level);
		}
	}

	UG_LOG("finished check_fsmoothing.\n");

	return true;
}

}
