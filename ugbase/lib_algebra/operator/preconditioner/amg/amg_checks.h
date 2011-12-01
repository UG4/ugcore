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

namespace ug
{

#ifdef UG_PARALLEL
template<typename TAlgebra>
void AMGBase<TAlgebra>::write_interfaces()
{
	cAMG_helper &h = m_amghelper;
	for(size_t level=0; level<m_usedLevels; level++)
	{
		const char *filename = (std::string(m_writeMatrixPath) + std::string("AMG_interface_L") + ToString(level) + "_" + ToString(pcl::GetProcRank()) + std::string(".mat")).c_str();
		std::fstream file(filename, std::ios::out);
		file << 1 << std::endl; // connection viewer version


		WritePositionsToStream(file, h.positions[level], h.dimension);

		file << 1 << std::endl;

		IndexLayout &masterLayout = levels[level]->pA->get_master_layout();
		IndexLayout &slaveLayout = levels[level]->pA->get_slave_layout();
		std::vector<std::string> v;
		v.resize(h.positions[level].size());

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
/**
 * writes the vector d in a connection-viewer-vec format
 *
 * \param 	filename
 * \param	d
 * \param	level
 */
template<typename TAlgebra>
bool AMGBase<TAlgebra>::writevec(std::string filename, const vector_type &const_d, size_t level)
{
	UG_ASSERT(m_writeMatrices, "");
#ifdef UG_PARALLEL
	vector_type d;
	CloneVector(d, const_d);
	d = const_d;
	d.change_storage_type(PST_CONSISTENT);
	size_t pid = pcl::GetProcRank();
#else
	const vector_type &d = const_d;
	size_t pid = 0;
#endif
	std::string name = (std::string(m_writeMatrixPath) + filename + ToString(level) + "_" + ToString(pid) + ".mat");

#ifdef UG_PARALLEL
	if(levels[level]->bHasBeenMerged && m_agglomerateLevel != level)
		AMGWriteToFile(levels[level]->uncollectedA, level, level, name.c_str(), m_amghelper);
	else
#endif
		AMGWriteToFile(*levels[level]->pA, level, level, name.c_str(), m_amghelper);
	std::fstream f(name.c_str(), std::ios::out | std::ios::app);
	//UG_LOG("writevec to " << name.c_str() << "\n");

	std::string name2 = (std::string(m_writeMatrixPath) + filename + ToString(level) + "_" + ToString(pid)  + ".values");
	f << "v " << name2 << "\n";

	std::fstream file(name2.c_str(), std::ios::out);
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
	UG_LOG("\n");
	UG_LOG("         check     \n");
	UG_LOG("========================\n");

	if(m_usedLevels <= 1)
	{
		UG_LOG("No Multigrid hierachy.\n");
		return true;
	}
	vector_type c, d;
	CloneVector(c, const_c);
	CloneVector(d, const_d);
	c = const_c;
	d = const_d;

	UG_LOG("\nLEVEL 0\n\n");
	check_level(c, d, 0);

	if(m_usedLevels > 1)
	{
		//levels[0]->R.apply(levels[0]->dH, d);
		levels[0]->cH.set(0.0);
		levels[0]->R.apply(levels[0]->dH, const_d);
		//levels[0]->R.apply(levels[0]->cH, c);

		for(size_t i=0; i<m_usedLevels-2; i++)
		{
			UG_LOG("\nLEVEL " << i+1 << "\n\n");
			vector_type d2;
			CloneVector(d2, levels[i]->dH);
			d2 = levels[i]->dH;
			check_level(levels[i]->cH, levels[i]->dH, i+1);

			levels[i+1]->R.apply(levels[i+1]->dH, d2);
			levels[i+1]->cH.set(0.0);
		}
	}

	UG_LOG("finished check.\n");
	return true;
}

template<typename TVector>
double ConstTwoNorm(const TVector &const_d)
{
	TVector d;
	CloneVector(d, const_d);
	d = const_d;
	return d.two_norm();

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
bool AMGBase<TAlgebra>::check_level(vector_type &c, vector_type &d, size_t level)
{
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

	//UG_LOG("preprenorm: " << d.two_norm() << std::endl);
	/*for(size_t i=0; i<5; i++)
		add_correction_and_update_defect(corr, d, level);*/

	double prenorm = ConstTwoNorm(d);
	UG_LOG("Prenorm = " << prenorm << "\n");

	double n1 = ConstTwoNorm(d), n2;

	vector_type &corr = levels[level]->corr;
	corr.set(0.0);
	//c.set(0.0);
#ifdef UG_PARALLEL
	corr.set_storage_type(PST_CONSISTENT);
#endif

	if(m_writeMatrices) writevec("AMG_aa_d_L", d, level);
	if(m_writeMatrices) writevec("AMG_aa_c_L", c, level);

	// presmooth
	for(size_t i=0; i < m_numPreSmooth; i++)
	{
		L.presmoother->apply_update_defect(corr, d);
		c += corr;
		n2 = ConstTwoNorm(d);	UG_LOG("presmoothing " << i+1 << ": " << n2/prenorm << "\t" <<n2/n1 << "\n");	n1 = n2;
	}

	if(m_bFSmoothing)
	{
		f_smoothing(corr, d, level);
		c+=corr;
		n2 = ConstTwoNorm(d);	UG_LOG("pre f-smoothing: " << n2/prenorm << "\t" <<n2/n1 << "\n");	n1 = n2;
	}


	if(m_writeMatrices) writevec("AMG_pa_d_L", d, level);
	if(m_writeMatrices) writevec("AMG_pa_c_L", c, level);

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

	double nH1 = ConstTwoNorm(dH);

	double preHnorm=nH1;
	size_t i=0;

	if(level+1 == m_usedLevels-1)
	{
		add_correction_and_update_defect(cH, dH, level+1);
		double nH2 = ConstTwoNorm(dH);
		UG_LOG("base solver reduced by " << nH2/nH1 << " (on coarse)" << std::endl);
	}
	else
	{
		cH.set(0.0);
		for(int i=0; i<50; i++)
		{
			for(int j=0; j< m_cycleType; j++)
				add_correction_and_update_defect(cH, dH, level+1);
			double nH2 = ConstTwoNorm(dH);
			if(i < 6)
			{
				UG_LOG("coarse correction (on coarse) " << i+1 << ": " << nH2/nH1 << " " << nH2/preHnorm <<  "\n"); nH1 = nH2;
			}
			if(m_writeMatrices) writevec((std::string("AMG_V")+ToString(i)+std::string("_c_L")).c_str(), cH, level+1);
			if(m_writeMatrices) writevec((std::string("AMG_V")+ToString(i)+std::string("_d_L")).c_str(), dH, level+1);
			if(i > 4 && nH2/preHnorm < 0.01)
				break;
		}


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
	L.P.apply(corr, cH);

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

	if(m_writeMatrices) writevec("AMG_ap_d_L", d, level);
	if(m_writeMatrices) writevec("AMG_ap_c_L", c, level);


	// post f-smoothing
	if(m_bFSmoothing)
	{
		f_smoothing(corr, d, level);
		c+=corr;
	}
	n2 = ConstTwoNorm(d);	UG_LOG("post f-smoothing: " << n2/prenorm << "\t" <<n2/n1 << "\n");	n1 = n2;

	// postsmooth
	for(size_t i=0; i < m_numPostSmooth; i++)
	{
		L.postsmoother->apply_update_defect(corr, d);
		c += corr;
		n2 = ConstTwoNorm(d);	UG_LOG("postsmoothing " << i+1 << ": " << n2/prenorm << "\t" <<n2/n1 << "\n");	n1 = n2;
	}


	double postnorm = ConstTwoNorm(d);
	UG_LOG("Level " << level << " reduction: " << postnorm/prenorm << std::endl);

	if(m_writeMatrices) writevec("AMG_pp_d_L", d, level);
	if(m_writeMatrices) writevec("AMG_pp_c_L", c, level);

	UG_LOG("Postnorm = " << ConstTwoNorm(d) << "\n\n");
	return true;
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
			for(k=0; k<100; k++)
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
			writevec((std::string("AMG_fsmoothing_c_") + ToString(j) + std::string("d_L")).c_str(), c, level);
			writevec((std::string("AMG_fsmoothing_d_") + ToString(j) + std::string("d_L")).c_str(), d, level);
		}
	}

	UG_LOG("finished check_fsmoothing.\n");

	return true;
}

}
