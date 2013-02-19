/*
 * domain_disc_base.h
 *
 *  Created on: 15.12.2011
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISC__SPATIAL_DISC__DOMAIN_DISC_BASE__
#define __H__LIB_DISC__SPATIAL_DISC__DOMAIN_DISC_BASE__

#include "domain_disc_interface.h"

namespace ug{

template <typename TDomain, typename TAlgebra>
class DomainDiscBase : public IDomainDiscretization<TAlgebra>
{
	public:
	///	Domain Type
		typedef TDomain domain_type;

	///	Algebra type
		typedef TAlgebra algebra_type;

	///	Type of algebra matrix
		typedef typename algebra_type::matrix_type matrix_type;

	///	Type of algebra vector
		typedef typename algebra_type::vector_type vector_type;

		using IDomainDiscretization<TAlgebra>::assemble_jacobian;
		using IDomainDiscretization<TAlgebra>::assemble_defect;
		using IDomainDiscretization<TAlgebra>::assemble_linear;
		using IDomainDiscretization<TAlgebra>::adjust_solution;
		using IDomainDiscretization<TAlgebra>::assemble_mass_matrix;
		using IDomainDiscretization<TAlgebra>::assemble_stiffness_matrix;
		using IDomainDiscretization<TAlgebra>::assemble_rhs;
		using IDomainDiscretization<TAlgebra>::prepare_timestep;
		using IDomainDiscretization<TAlgebra>::finish_timestep;

	public:
	///	sets the approximation space
		virtual void set_approximation_space(SmartPtr<ApproximationSpace<domain_type> > approxSpace)
		{
			m_spApproxSpace = approxSpace;
		}

	///	returns approximation space
		SmartPtr<ApproximationSpace<TDomain> > approximation_space()
		{
			return m_spApproxSpace;
		}

	///	returns approximation space
		ConstSmartPtr<ApproximationSpace<TDomain> > approximation_space() const
		{
			return m_spApproxSpace;
		}

		virtual void assemble_jacobian(matrix_type& J, const vector_type& u,
		                               GridLevel gl)
		{
			assemble_jacobian(J, u, dd(gl));
		}

		virtual void assemble_defect(vector_type& d, const vector_type& u,
		                           GridLevel gl)
		{
			assemble_defect(d, u, dd(gl));
		}

		virtual void assemble_linear(matrix_type& mat, vector_type& rhs,
		                           GridLevel gl)
		{
			assemble_linear(mat, rhs, dd(gl));
		}

		virtual void adjust_solution(vector_type& u, GridLevel gl)
		{
			adjust_solution(u, dd(gl));
		}

		virtual void assemble_mass_matrix(matrix_type& M, const vector_type& u,
		                                  GridLevel gl)
		{
			assemble_mass_matrix(M, u, dd(gl));
		}

		virtual void assemble_stiffness_matrix(matrix_type& A, const vector_type& u,
		                                       GridLevel gl)
		{
			assemble_stiffness_matrix(A, u, dd(gl));
		}

		virtual void assemble_rhs(vector_type& rhs, const vector_type& u,
		                          GridLevel gl)
		{
			assemble_rhs(rhs, u, dd(gl));
		}

		virtual void assemble_rhs(vector_type& rhs, GridLevel gl)
		{
			assemble_rhs(rhs, dd(gl));
		}

		// time dependent
		virtual	void prepare_timestep(ConstSmartPtr<VectorTimeSeries<vector_type> > vSol, GridLevel gl)
		{
			prepare_timestep(vSol, dd(gl));
		}

		virtual void assemble_jacobian(matrix_type& J,
		                               ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
		                               const number s_a, GridLevel gl)
		{
			assemble_jacobian(J, vSol, s_a, dd(gl));
		}

		virtual	void assemble_defect(vector_type& d,
		       	                     ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
		       	                     const std::vector<number>& vScaleMass,
		       	                     const std::vector<number>& vScaleStiff,
		       	                     GridLevel gl)
		{
			assemble_defect(d, vSol, vScaleMass, vScaleStiff, dd(gl));
		}

		virtual void assemble_linear(matrix_type& A, vector_type& b,
		                             ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
		                             const std::vector<number>& vScaleMass,
		                             const std::vector<number>& vScaleStiff,
		                             GridLevel gl)
		{
			assemble_linear(A, b, vSol, vScaleMass, vScaleStiff, dd(gl));
		}

		virtual void assemble_rhs(	 vector_type& b,
		                             ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
		                             const std::vector<number>& vScaleMass,
		                             const std::vector<number>& vScaleStiff,
		                             GridLevel gl)
		{
			assemble_rhs(b, vSol, vScaleMass, vScaleStiff, dd(gl));
		}

		virtual void adjust_solution(vector_type& u, number time, GridLevel gl)
		{
			adjust_solution(u, time, dd(gl));
		}

		virtual void finish_timestep(ConstSmartPtr<VectorTimeSeries<vector_type> > vSol, GridLevel gl)
		{
			finish_timestep(vSol, dd(gl));
		}

		virtual ~DomainDiscBase() {};

	protected:
	///	returns the level dof distribution
		ConstSmartPtr<DoFDistribution> dd(const GridLevel& gl) const
		{
			return m_spApproxSpace->dof_distribution(gl);
		}

	protected:
	///	Approximation Space
		SmartPtr<ApproximationSpace<domain_type> > m_spApproxSpace;
};


} // end namespace ug

#endif /* __H__LIB_DISC__SPATIAL_DISC__DOMAIN_DISC_BASE__ */
