// external libraries
#include <vector>
#include <cmath>

// ug libraries
#include "common/common.h"
#include "common/util/smart_pointer.h"

#include "lib_algebra/lib_algebra.h"

#include "lib_disc/time_disc/time_disc_interface.h"

namespace ug{

/*
 std::vector<size_t> steps {1, 2, 4};
 timex = new AitkenNevilleTimex(steps);

 //
 timex.set_solution(sol0, 0);  // single step
 timex.set_solution(sol1, 1);  // double step
 timex.set_solution(sol2, 2);  // four step

 //
 timex.apply()                // updating sol1 and sol2

 timex.get_error_estimate()

 */

namespace tools {
//! calculates dest = alpha1*v1 + alpha2*v2 + alpha3*v3. for doubles
inline void VecScaleAddWithNorm(double &dest, double alpha1, const double &v1, double alpha2, const double &v2, double alpha3, const double &v3, double &norm)
{
	dest = alpha1*v1 + alpha2*v2 + alpha3*v3;
	//norm += (alpha2*v2 + alpha3*v3)*(alpha2*v2 + alpha3*v3);
	norm = std::max(norm, fabs(alpha2*v2 + alpha3*v3));
}

//! calculates dest = alpha1*v1 + alpha2*v2 + alpha3*v3
template<typename vector_t, template <class T> class TE_VEC>
inline void VecScaleAddWithNorm(TE_VEC<vector_t> &dest, double alpha1, const TE_VEC<vector_t> &v1, double alpha2, const TE_VEC<vector_t> &v2, double alpha3, const TE_VEC<vector_t> &v3, double &norm)
{
	for(size_t i=0; i<dest.size(); i++)
		VecScaleAddWithNorm(dest[i], alpha1, v1[i], alpha2, v2[i], alpha3, v3[i], norm);
}

}

template <typename TVector>
class AitkenNevilleTimex
{
	public:
	///	vector type of solutions
		typedef TVector vector_type;

	public:

		/** Aitken Neville scheme with a given umber */
		AitkenNevilleTimex(std::vector<size_t> nsteps)
		: m_num_steps(nsteps),
		  m_solution(nsteps.size()),
		  m_subdiag_error_est(nsteps.size(), INFINITY) {};

		virtual ~AitkenNevilleTimex() {}

		void set_global_stepsize(number H) {m_stepsize=H;}
		number get_global_stepsize() {return m_stepsize;}


		 void set_solution(SmartPtr<vector_type> soli, int i)
		 { m_solution[i] = soli; }

		SmartPtr<vector_type> get_solution(size_t i)
		{ return m_solution[i]; }

		number get_error_estimate()
		{ return m_subdiag_error_est.back();}


		/**
		 * Triangular Aitken Neville extrapolation:
		 *
		 * T_11
		 * T_21 T_22
		 *
		 * T_N1        T_NN
		 *
		 * for T_ik
 		 * */
		void apply()
		{
			UG_ASSERT(  m_num_steps.size() == m_solution.size(), "Dimensions do not match");
			const size_t N = m_num_steps.size();

			//m_subdiag_error_est[0] = ;
			// process columns (left to right)
			for (size_t k=1; k<N; ++k)
			{

				// process rows (bottom up, allows recycling memory)
				for (size_t i=N-1; i>k; --i)
				{
					UG_ASSERT(m_solution[i].valid(), "Invalid SmarPtr!");
					UG_ASSERT(m_solution[i-1].valid(), "Invalid SmarPtr!");

					vector_type &solcoarse = *m_solution[i-1];
					vector_type &solfine = *m_solution[i];

					// (2^p -1)
					// m_solution[i] += (1.0/scal)*(m_solution[i]- m_solution[i-1]);
					number scaling = (m_num_steps[i]/m_num_steps[i-k]-1.0);
					VecScaleAdd(solfine, (1.0+1.0/scaling), solfine,
										-(1.0/scaling), solcoarse);

				}

				// subdiagonal error estimate
				{
					UG_ASSERT(m_solution[k].valid(), "Invalid SmarPtr!");
					UG_ASSERT(m_solution[k-1].valid(), "Invalid SmarPtr!");

					vector_type &solcoarse = *m_solution[k-1];
					vector_type &solfine = *m_solution[k];

					number subdiag_error_est=0.0;
					number scaling = (m_num_steps[k]/m_num_steps[k-1]-1.0);
					tools::VecScaleAddWithNorm(solfine, 1.0, solfine,
						(1.0/scaling), solfine,  -(1.0/scaling), solcoarse, subdiag_error_est);
					//m_subdiag_error_est[k] = sqrt(subdiag_error_est*scaling);
					m_subdiag_error_est[k]=subdiag_error_est*scaling;

					UG_LOG(" ErrorEst["<< k<<"]=" << m_subdiag_error_est[k] << ";" << std::endl);
				}

			}

		}

protected:
		number substep(size_t i) {return m_stepsize/m_num_steps[i];}

private:
		/** time step */
		number m_stepsize;
		static const int m_order=1;


		/** number of intermediate steps (per stage)*/
		std::vector<size_t> m_num_steps;

		/** vector of solutions (per stage)*/
		std::vector<SmartPtr<vector_type> > m_solution;

		/** sub-diagonal error estimate (per stage)*/
		std::vector<number> m_subdiag_error_est;


};



/** Estimate*/
/*
class TimeStepEstimator{

public:
	TimeStepEstimator(double tol) :
	m_rho(0.9), m_gamma(2.0), m_tol(tol)
	{

	}

	bool check(double eps, double &factor)
	{
		if (eps>=tol) return false;

		double val = (m_rho*m_tol)/eps;
		factor = pow(val, m_gamma);

	}


private:
	double m_rho;
	double m_gamma;
	double m_tol;
};
*/

}
