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

template <typename TDomain, typename TAlgebra, typename TImpl>
class DomainDiscBase
	: 	public IDomainDiscretization<TAlgebra>
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

	public:
	///	sets the approximation space
		virtual void set_approximation_space(SmartPtr<ApproximationSpace<domain_type> > approxSpace)
		{
			m_spApproxSpace = approxSpace;
		}

		virtual void assemble_jacobian(matrix_type& J, const vector_type& u,
		                               GridLevel gl)
		{
			if(gl.type() == GridLevel::LEVEL)
				getImpl().template assemble_jacobian<LevelDoFDistribution>(J, u, lev_dd(gl));
			else if (gl.type() == GridLevel::SURFACE)
				getImpl().template assemble_jacobian<SurfaceDoFDistribution>(J, u, surf_dd(gl));
			else
				UG_THROW("Grid Level not recognized.");
		}

		virtual void assemble_defect(vector_type& d, const vector_type& u,
		                           GridLevel gl)
		{
			if(gl.type() == GridLevel::LEVEL)
				getImpl().template assemble_defect<LevelDoFDistribution>(d, u, lev_dd(gl));
			else if (gl.type() == GridLevel::SURFACE)
				getImpl().template assemble_defect<SurfaceDoFDistribution>(d, u, surf_dd(gl));
			else
				UG_THROW("Grid Level not recognized.");
		}

		virtual void assemble_linear(matrix_type& mat, vector_type& rhs,
		                           GridLevel gl)
		{
			if(gl.type() == GridLevel::LEVEL)
				getImpl().template assemble_linear<LevelDoFDistribution>(mat, rhs, lev_dd(gl));
			else if (gl.type() == GridLevel::SURFACE)
				getImpl().template assemble_linear<SurfaceDoFDistribution>(mat, rhs, surf_dd(gl));
			else
				UG_THROW("Grid Level not recognized.");
		}

		virtual void adjust_solution(vector_type& u, GridLevel gl)
		{
			if(gl.type() == GridLevel::LEVEL)
				getImpl().template adjust_solution<LevelDoFDistribution>(u, lev_dd(gl));
			else if (gl.type() == GridLevel::SURFACE)
				getImpl().template adjust_solution<SurfaceDoFDistribution>(u, surf_dd(gl));
			else
				UG_THROW("Grid Level not recognized.");
		}

		// mass- / stiffness - matrix, rhs
		virtual void assemble_mass_matrix(matrix_type& M, const vector_type& u,
		                               GridLevel gl)
		{
			if(gl.type() == GridLevel::LEVEL)
				getImpl().template assemble_mass_matrix<LevelDoFDistribution>(M, u, lev_dd(gl));
			else if (gl.type() == GridLevel::SURFACE)
				getImpl().template assemble_mass_matrix<SurfaceDoFDistribution>(M, u, surf_dd(gl));
			else
				UG_THROW("Grid Level not recognized.");
		}

		virtual void assemble_stiffness_matrix(matrix_type& A, const vector_type& u,
		                               GridLevel gl)
		{
			if(gl.type() == GridLevel::LEVEL)
				getImpl().template assemble_stiffness_matrix<LevelDoFDistribution>(A, u, lev_dd(gl));
			else if (gl.type() == GridLevel::SURFACE)
				getImpl().template assemble_stiffness_matrix<SurfaceDoFDistribution>(A, u, surf_dd(gl));
			else
				UG_THROW("Grid Level not recognized.");
		}

		virtual void assemble_rhs(vector_type& rhs, const vector_type& u,
		                           GridLevel gl)
		{
			if(gl.type() == GridLevel::LEVEL)
				getImpl().template assemble_rhs<LevelDoFDistribution>(rhs, u, lev_dd(gl));
			else if (gl.type() == GridLevel::SURFACE)
				getImpl().template assemble_rhs<SurfaceDoFDistribution>(rhs, u, surf_dd(gl));
			else
				UG_THROW("Grid Level not recognized.");
		}

		virtual void assemble_rhs(vector_type& rhs, GridLevel gl)
		{
			if(gl.type() == GridLevel::LEVEL)
				getImpl().template assemble_rhs<LevelDoFDistribution>(rhs, lev_dd(gl));
			else if (gl.type() == GridLevel::SURFACE)
				getImpl().template assemble_rhs<SurfaceDoFDistribution>(rhs, surf_dd(gl));
			else
				UG_THROW("Grid Level not recognized.");
		}


		// time dependent

		virtual	void prepare_timestep(ConstSmartPtr<VectorTimeSeries<vector_type> > vSol, GridLevel gl)
		{
			if(gl.type() == GridLevel::LEVEL)
				getImpl().template prepare_timestep<LevelDoFDistribution>(vSol, lev_dd(gl));
			else if (gl.type() == GridLevel::SURFACE)
				getImpl().template prepare_timestep<SurfaceDoFDistribution>(vSol, surf_dd(gl));
			else
				UG_THROW("Grid Level not recognized.");
		}

		virtual void assemble_jacobian(matrix_type& J,
		                               ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
		                               const number s_a, GridLevel gl)
		{
			if(gl.type() == GridLevel::LEVEL)
				getImpl().template assemble_jacobian<LevelDoFDistribution>(J, vSol, s_a, lev_dd(gl));
			else if (gl.type() == GridLevel::SURFACE)
				getImpl().template assemble_jacobian<SurfaceDoFDistribution>(J, vSol, s_a, surf_dd(gl));
			else
				UG_THROW("Grid Level not recognized.");
		}

		virtual	void assemble_defect(vector_type& d,
		       	                     ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
		       	                     const std::vector<number>& vScaleMass,
		       	                     const std::vector<number>& vScaleStiff,
		       	                     GridLevel gl)
		{
			if(gl.type() == GridLevel::LEVEL)
				getImpl().template assemble_defect<LevelDoFDistribution>(d, vSol, vScaleMass, vScaleStiff, lev_dd(gl));
			else if (gl.type() == GridLevel::SURFACE)
				getImpl().template assemble_defect<SurfaceDoFDistribution>(d, vSol, vScaleMass, vScaleStiff, surf_dd(gl));
			else
				UG_THROW("Grid Level not recognized.");
		}

		virtual void assemble_linear(matrix_type& A, vector_type& b,
		                             ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
		                             const std::vector<number>& vScaleMass,
		                             const std::vector<number>& vScaleStiff,
		                             GridLevel gl)
		{
			if(gl.type() == GridLevel::LEVEL)
				getImpl().template assemble_linear<LevelDoFDistribution>(A, b, vSol, vScaleMass, vScaleStiff, lev_dd(gl));
			else if (gl.type() == GridLevel::SURFACE)
				getImpl().template assemble_linear<SurfaceDoFDistribution>(A, b, vSol, vScaleMass, vScaleStiff, surf_dd(gl));
			else
				UG_THROW("Grid Level not recognized.");
		}

		virtual void assemble_rhs(	 vector_type& b,
		                             ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
		                             const std::vector<number>& vScaleMass,
		                             const std::vector<number>& vScaleStiff,
		                             GridLevel gl)
		{
			if(gl.type() == GridLevel::LEVEL)
				getImpl().template assemble_rhs<LevelDoFDistribution>(b, vSol, vScaleMass, vScaleStiff, lev_dd(gl));
			else if (gl.type() == GridLevel::SURFACE)
				getImpl().template assemble_rhs<SurfaceDoFDistribution>(b, vSol, vScaleMass, vScaleStiff, surf_dd(gl));
			else
				UG_THROW("Grid Level not recognized.");
		}

		virtual void adjust_solution(vector_type& u, number time, GridLevel gl)
		{
			if(gl.type() == GridLevel::LEVEL)
				getImpl().template adjust_solution<LevelDoFDistribution>(u, time, lev_dd(gl));
			else if (gl.type() == GridLevel::SURFACE)
				getImpl().template adjust_solution<SurfaceDoFDistribution>(u, time, surf_dd(gl));
			else
				UG_THROW("Grid Level not recognized.");
		}

		virtual	void finish_timestep(ConstSmartPtr<VectorTimeSeries<vector_type> > vSol, GridLevel gl)
		{
			if(gl.type() == GridLevel::LEVEL)
				getImpl().template finish_timestep<LevelDoFDistribution>(vSol, lev_dd(gl));
			else if (gl.type() == GridLevel::SURFACE)
				getImpl().template finish_timestep<SurfaceDoFDistribution>(vSol, surf_dd(gl));
			else
				UG_THROW("Grid Level not recognized.");
		}

	protected:
	///	returns the level dof distribution
		ConstSmartPtr<LevelDoFDistribution> lev_dd(GridLevel gl) const
		{
			return m_spApproxSpace->level_dof_distribution(gl.level());
		}

	///	returns the surface dof distribution
		ConstSmartPtr<SurfaceDoFDistribution> surf_dd(GridLevel gl) const
		{
			return m_spApproxSpace->surface_dof_distribution(gl.level());
		}

	protected:
	///	access to implementation
		TImpl& getImpl() {return static_cast<TImpl&>(*this);}

	///	const access to implementation
		const TImpl& getImpl() const {return static_cast<const TImpl&>(*this);}

	protected:
	///	Approximation Space
		SmartPtr<ApproximationSpace<domain_type> > m_spApproxSpace;
};


} // end namespace ug

#endif /* __H__LIB_DISC__SPATIAL_DISC__DOMAIN_DISC_BASE__ */
