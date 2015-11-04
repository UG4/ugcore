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

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__LOGNORMAL_FIELD__ */
