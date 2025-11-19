/*
 * Copyright (c) 2014-2015:  G-CSC, Goethe University Frankfurt
 * Authors: Ivo Muha, Martin Rupp
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

/*
 *      example usage (on unit square)
 *      	elemDisc = ConvectionDiffusion("c", "Inner", disc)
 *			elemDisc:set_diffusion(LognormalRandomField(100, 0, 1, 0.01))
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__LOGNORMAL_RANDOM_FIELD__
#define __H__UG__LIB_DISC__SPATIAL_DISC__LOGNORMAL_RANDOM_FIELD__

// extern headers
#include <vector>
#include "lib_disc/spatial_disc/user_data/std_glob_pos_data.h"

namespace ug{

template <typename TData, int dim, typename TRet = void>
class LognormalRandomField
		: public StdGlobPosData<LognormalRandomField<TData, dim, TRet>, TData, dim, TRet>
{

public:

	LognormalRandomField()
	{
		m_bNoExp=false;
		set_config(100, 0, 1, 0.1);
	};

	LognormalRandomField(size_t N, double mean_f, double sigma_f, double sigma)
	{
		m_bNoExp=false;
		set_config(N, mean_f, sigma_f, sigma);
	};


	virtual ~LognormalRandomField()
	{
	}
	inline TRet evaluate(TData& D, const MathVector<dim>& x, number time, int si) const;

	void set_no_exp() { m_bNoExp = true; }
	void set_config(size_t N, double mean_f, double sigma_f, double sigma);

	std::string config_string() const;

private:
	double gasdev();
    double undev();
    double eval_K(const MathVector<dim> &x)  const;

    MathVector<dim> m_sigma;

    double m_dMean_f;
    double m_dSigma_f;
    double m_N;
    bool m_bNoExp;
    double m_dSigma;

	std::vector<MathVector<dim> > m_vRandomQvec;
	std::vector<double> m_vRandomAlpha;

};

} // end ug

#include "lognormal_random_field_impl.h"

#endif