namespace ug
{


template<typename TAlgebra>
bool AMGBase<TAlgebra>::check(const vector_type &const_c, const vector_type &const_d)
{

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
		levels[0]->R.apply(levels[0]->dH, d);
		levels[0]->cH.set(0.0);
		for(size_t i=0; i<m_usedLevels-1; i++)
		{
			UG_LOG("\nLEVEL " << i << "\n\n");
			check_level(levels[i-1]->cH, levels[i-1]->dH, i);

			if(i < m_usedLevels)
			{
				levels[i]->R.apply(levels[i]->dH, levels[i-1]->dH);
				levels[i]->cH.set(0.0);
			}
		}
	}
	return true;
}

template<typename TAlgebra>
bool AMGBase<TAlgebra>::writevec(std::string filename, const vector_type &d, size_t level)
{
	UG_ASSERT(m_writeMatrices, "");
#ifdef UG_PARALLEL
	size_t pid = pcl::GetProcRank();
#else
	size_t pid = 0;
#endif
	std::string name = (std::string(m_writeMatrixPath) + filename + ToString(level) + "_" + ToString(pid) + ".mat");
	AMGWriteToFile(*levels[level]->pA, level, level, name.c_str(), m_amghelper);
	std::fstream f(name.c_str(), std::ios::out | std::ios::app);

	std::string name2 = (std::string(m_writeMatrixPath) + filename + ToString(level) + "_" + ToString(pid)  + ".values");
	f << "v " << name2 << "\n";

	std::fstream file(name2.c_str(), std::ios::out);
	for(size_t i=0; i<d.size(); i++)
		file << m_amghelper.GetOriginalIndex(level, i) << " " << (d[i]) << std::endl;
	return true;
}


template<typename TAlgebra>
bool AMGBase<TAlgebra>::check_level(vector_type &c, vector_type &d, size_t level)
{
	AMGLevel &L = *levels[level];
	const matrix_operator_type &A = *L.pA;
	vector_type corr;


	//UG_LOG("preprenorm: " << d.two_norm() << std::endl);
	/*for(size_t i=0; i<5; i++)
		add_correction_and_update_defect(corr, d, level);*/

	double prenorm = d.two_norm();
	UG_LOG("Prenorm = " << prenorm << "\n");

	// presmooth
	// same as setting c.set(0.0).

	double n1 = d.two_norm(), n2;
	CloneVector(corr, c);
	L.presmoother->apply_update_defect(c, d);
	n2 = d.two_norm();	UG_LOG("presmoothing 1 " << ": " << n2/prenorm << "\t" <<n2/n1 << "\n");	n1 = n2;
	for(size_t i=1; i < m_numPreSmooth; i++)
	{
		L.presmoother->apply_update_defect(corr, d);
		c += corr;
		n2 = d.two_norm();	UG_LOG("presmoothing " << i+1 << ": " << n2/prenorm << "\t" <<n2/n1 << "\n");	n1 = n2;
	}

	if(m_bFSmoothing)
	{
		f_smoothing(corr, d, level);
		c+=corr;
		n2 = d.two_norm();	UG_LOG("pre f-smoothing: " << n2/prenorm << "\t" <<n2/n1 << "\n");	n1 = n2;
	}


	if(m_writeMatrices) writevec("AMG_dp_L", d, level);

	vector_type &cH = levels[level]->cH;
	vector_type &dH = levels[level]->dH;

	// restrict defect
	// dH = m_R[level]*d;
	L.R.apply(dH, d);

	double nH1 = dH.two_norm();

	double preHnorm=nH1;
	size_t i=0;

	if(level+1 == m_usedLevels-1)
	{
		solve_on_base(cH, dH, level+1);
		double nH2 = dH.two_norm();
		UG_LOG("base solver reduced by " << nH2/nH1 << " (on coarse)" << std::endl);
	}
	else
	{
		cH.set(0.0);
		for(i=0; i<100; i++)
		{
			for(int j=0; j< m_cycleType; j++)
				add_correction_and_update_defect(cH, dH, level+1);

			double nH2 = dH.two_norm();
			if(i < 6)
			{	UG_LOG("coarse correction (on coarse) " << i+1 << ": " << nH2/nH1 << "\n"); nH1 = nH2; }
			if(nH2/preHnorm < 0.01) { UG_LOG("coarse solver reduced by 0.01 in iteration " << i+1 << std::endl); break; }
		}
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

	n2 = d.two_norm();
	UG_LOG("complete coarse correction " << i+1 << ": " << n2/prenorm << "\t" <<n2/n1 << "\n");	n1 = n2;

	if(m_writeMatrices) writevec("AMG_dc_L", d, level);


	// post f-smoothing
	if(m_bFSmoothing)
	{
		f_smoothing(corr, d, level);
		c+=corr;
	}
	n2 = d.two_norm();	UG_LOG("post f-smoothing: " << n2/prenorm << "\t" <<n2/n1 << "\n");	n1 = n2;

	// postsmooth
	for(size_t i=0; i < m_numPostSmooth; i++)
	{
		L.postsmoother->apply_update_defect(corr, d);
		c += corr;
		n2 = d.two_norm();	UG_LOG("postsmoothing " << i+1 << ": " << n2/prenorm << "\t" <<n2/n1 << "\n");	n1 = n2;
	}


	double postnorm = d.two_norm();
	UG_LOG("Level " << level << " reduction: " << postnorm/prenorm << std::endl);

	if(m_writeMatrices) writevec("AMG_d_L", d, level);

	UG_LOG("\n\n");
	return true;
}


template<typename TAlgebra>
bool AMGBase<TAlgebra>::check2(const vector_type &const_c, const vector_type &const_d)
{
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
		UG_LOG("Testing Level 0 to " << i << ":\n");
		double prenorm = d.two_norm();
		add_correction_and_update_defect(c, d, 0, i);
		double postnorm = d.two_norm();
		UG_LOG("Reduction is " << postnorm/prenorm << "\n\n");
	}

	return true;
}


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
				if(dH.two_norm() < 1e-12)
					break;
			}
			if(dH.two_norm() > 1e-12) { UG_LOG("not converged after " << k << " steps\n"); }
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

}
