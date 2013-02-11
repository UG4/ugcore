/*
 * constraint_base.h
 *
 *  Created on: 15.12.2011
 *      Author: andreasvogel
 */

#ifndef CONSTRAINT_BASE_H_
#define CONSTRAINT_BASE_H_

#include "constraint_interface.h"
#include "lib_disc/dof_manager/level_dof_distribution.h"
#include "lib_disc/dof_manager/surface_dof_distribution.h"
#include "lib_disc/spatial_disc/ass_adapter.h"

namespace ug{

template <typename TDomain, typename TAlgebra, typename TImpl>
class ConstraintBase
	: 	public IDomainConstraint<TDomain, TAlgebra>
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

	///	returns approximation space
		SmartPtr<ApproximationSpace<domain_type> > approximation_space()
		{
			return m_spApproxSpace;
		}

	///	sets the assemble index for index-wise assemble routine
		virtual void set_ass_index(){set_ass_index(0.0, false);}
		virtual void set_ass_index(size_t ind, bool index_set = true)
		{
			m_AssIndex.index_set = index_set;
			m_AssIndex.index = ind;
		}
		
	///	adapts jacobian to enforce constraints
		virtual void adjust_jacobian(matrix_type& J, const vector_type& u,
		                             GridLevel gl, number time = 0.0)
		{
			if(gl.type() == GridLevel::LEVEL)
				getImpl().template adjust_jacobian<LevelDoFDistribution>(J, u, lev_dd(gl), time);
			else if (gl.type() == GridLevel::SURFACE)
				getImpl().template adjust_jacobian<SurfaceDoFDistribution>(J, u, surf_dd(gl), time);
			else
				UG_THROW("Grid Level not recognized.");
		}

	///	adapts defect to enforce constraints
		virtual void adjust_defect(vector_type& d, const vector_type& u,
		                           GridLevel gl, number time = 0.0)
		{
			if(gl.type() == GridLevel::LEVEL)
				getImpl().template adjust_defect<LevelDoFDistribution>(d, u, lev_dd(gl), time);
			else if (gl.type() == GridLevel::SURFACE)
				getImpl().template adjust_defect<SurfaceDoFDistribution>(d, u, surf_dd(gl), time);
			else
				UG_THROW("Grid Level not recognized.");
		}

	///	adapts matrix and rhs (linear case) to enforce constraints
		virtual void adjust_linear(matrix_type& mat, vector_type& rhs,
		                           GridLevel gl, number time = 0.0)
		{
			if(gl.type() == GridLevel::LEVEL)
				getImpl().template adjust_linear<LevelDoFDistribution>(mat, rhs, lev_dd(gl), time);
			else if (gl.type() == GridLevel::SURFACE)
				getImpl().template adjust_linear<SurfaceDoFDistribution>(mat, rhs, surf_dd(gl), time);
			else
				UG_THROW("Grid Level not recognized.");
		}

	///	adapts a rhs to enforce constraints
		virtual void adjust_rhs(vector_type& rhs, const vector_type& u,
		                        GridLevel gl, number time = 0.0)
		{
			if(gl.type() == GridLevel::LEVEL)
				getImpl().template adjust_rhs<LevelDoFDistribution>(rhs, u, lev_dd(gl), time);
			else if (gl.type() == GridLevel::SURFACE)
				getImpl().template adjust_rhs<SurfaceDoFDistribution>(rhs, u, surf_dd(gl), time);
			else
				UG_THROW("Grid Level not recognized.");
		}

	///	sets the constraints in a solution vector
		virtual void adjust_solution(vector_type& u, GridLevel gl,
		                             number time = 0.0)
		{
			if(gl.type() == GridLevel::LEVEL)
				getImpl().template adjust_solution<LevelDoFDistribution>(u, lev_dd(gl), time);
			else if (gl.type() == GridLevel::SURFACE)
				getImpl().template adjust_solution<SurfaceDoFDistribution>(u, surf_dd(gl), time);
			else
				UG_THROW("Grid Level not recognized.");
		}

	///	returns the type of constraints
		virtual int type() const = 0;

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

	///	Assemble index
		AssIndex m_AssIndex;
};


} // end namespace ug

#endif /* CONSTRAINT_BASE_H_ */
