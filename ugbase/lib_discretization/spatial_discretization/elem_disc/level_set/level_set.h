/*
 * level_set_util.h
 *
 *  Created on: 01.07.2011
 *      Author: Christian Wehner  christian.wehner@gcsc.uni-frankfurt.de
 */

#ifndef LEVEL_SET_UTIL_H_
#define LEVEL_SET_UTIL_H_

#include <vector>
#include "lib_discretization/spatial_discretization/disc_util/finite_volume_geometry.h"

namespace ug{
	
template<typename TGridFunction>
class FV1LevelSetDisc
{
	///	domain type
		typedef typename TGridFunction::domain_type domain_type;
		
	///	algebra type
		typedef typename TGridFunction::algebra_type algebra_type;
	
	/// dof distribution type	
		typedef typename TGridFunction::dof_distribution_type dof_distribution_type;

		typedef typename TGridFunction::dof_distribution_type::implementation_type
				dof_distribution_impl_type;

	///	world dimension
		static const int dim = domain_type::dim;

		///for debug	type of scv-size attachment
			//	typedef typename Grid::VertexAttachmentAccessor<Attachment<number> > aaDiv;

	///	grid type
		typedef typename domain_type::grid_type grid_type;

	///	type of scv-size attachment
		typedef typename Grid::VertexAttachmentAccessor<Attachment<number> > aaSCV;

	///	type of gradient attachment
		typedef typename Grid::VertexAttachmentAccessor<Attachment<MathVector<dim> > > aaGrad;

		// 	Type of multi index vector
		typedef typename dof_distribution_type::multi_index_vector_type multi_index_vector_type;
		
		// edge iterator
		typedef geometry_traits<EdgeBase>::const_iterator EdgeBaseConstIterator;

		// vertex base iterator
		typedef geometry_traits<VertexBase>::const_iterator VertexBaseConstIterator;

	 public:
    ///	Constructor
      	    ///	Destructor
      	~FV1LevelSetDisc() {};

    ///	set time step
     	FV1LevelSetDisc() :
      		m_dt(0.0), m_time(0), m_gamma(1.0), m_delta(0.0),
      		m_reinit(false), m_analyticalSolution(true),
      		m_analyticalVelocity(true),m_externalVelocity(true),m_analyticalSource(true),
      		m_divFree(false), m_source(false), m_solutionNr(0),
      		m_velocityNr(0),m_sourceNr(0),m_nrOfSteps(1),m_bdryCondition(0),
      		m_static_values_type(0),m_maxCFL(0),m_print(false),m_timestep_nr(0),m_limiter(false)
      	{}

void set_dt(number deltaT){ UG_LOG("Set dt="<<deltaT<<"\n"); m_dt=deltaT; };

	/// set scale parameters for external velocity and velocity in normal direction
    	void set_vel_scale(number gamma,number delta){ };
    	void set_reinit(size_t n){ m_reinit=1;m_gamma=0;m_delta=1;m_nrOfSteps=n; };
		void set_div_bool(bool b){m_divFree=b;};
		void set_solution_nr(size_t n){m_solutionNr=n;};
		void set_velocity_nr(size_t n){m_velocityNr=n;};
		void set_source(size_t n){m_sourceNr=n;};
		void set_source_bool(bool b){m_source = b;};
		void set_static_values_type(size_t n){m_static_values_type=n;};
		void set_analytical_velocity_bool(bool b){m_analyticalVelocity=b;};
		void set_delta(number delta){m_delta =delta;}
		void set_time(double t){m_time = t;}
		void set_info(bool b){m_print=b;};
		void set_limiter(bool b){m_limiter=b;};
		void set_timestep_nr(size_t n){m_timestep_nr = n;};
		// set nr of time steps to perform in advect_lsf (default is 1)
		void set_nr_of_steps(size_t n){m_nrOfSteps = n;};
		///	adds a post process to be used when stepping the level set function
//		void add_post_process(IConstraint<dof_distribution_type, algebra_type>& pp){m_vPP.push_back(&pp);}
	    bool compute_error(TGridFunction& numsol);
		bool advect_lsf(TGridFunction& uNew,TGridFunction& u);
	    bool init_function(TGridFunction& u);
	    void set_neumann_boundary(const char* subsets){m_neumannSubsets = subsets;}
	///	adds a post process to be used when stepping the level set function
		void add_post_process(IConstraint<dof_distribution_impl_type, algebra_type>& pp) {m_vPP.push_back(&pp);}
		
	 protected:
	    number analytic_solution(number,MathVector<dim>);
		number analytic_source(number,MathVector<dim>);
	    bool analytic_velocity(MathVector<dim>&,number, MathVector<dim>);

	///	fills the scvVolume attachment for all element types
	    bool calculate_vertex_vol(TGridFunction& u, aaSCV& aaScvVolume);

//	    template <typename TElem>
//   	    bool calculate_vertex_grad_vol(grid_type& grid,TGridFunction& u, aaGrad& aaGradient, aaSCV& aaVolume );

	    bool calculate_vertex_grad_vol(TGridFunction& u, aaGrad& aaGradient, aaSCV& aaVolume );

		template <typename TElem>
		bool assemble_element(TElem& elem, DimFV1Geometry<dim>& geo, grid_type& grid,TGridFunction& uNew,const TGridFunction& uOld,aaGrad& aaGradient, aaSCV& aaVolume );

//fordebug		template <typename TElem>
//		bool assemble_divergence(TElem& elem,grid_type& grid,TGridFunction& uNew,aaDiv& aaDivergence,aaGrad& aaGradient );

		bool assign_dirichlet(TGridFunction&,int);
		bool limit_grad(TGridFunction& uOld, aaGrad& aaGradient);



		//bool limit_grad_alpha(TGridFunction& uOld,aaGrad& aaGradient,aaAlpha& aaAlpha);

	private:
	///	vector holding all scheduled post processes
		std::vector<IConstraint<dof_distribution_impl_type, algebra_type>*> m_vPP;
      	number m_dt;
		number m_time;
    	number m_gamma;
      	number m_delta;
    	bool m_reinit;
    	bool m_analyticalSolution;
	    bool m_analyticalVelocity;
	    bool m_externalVelocity;
		bool m_analyticalSource;
    	bool m_divFree;
		bool m_source;
        size_t m_solutionNr;
        size_t m_velocityNr;
		size_t m_sourceNr;
     	size_t m_nrOfSteps;
    	size_t m_bdryCondition;
		// keep values static in computation 0 nothing, 1 interface (reinitialization), 2 inside, 3 outside
		size_t m_static_values_type;
		number m_maxCFL;
		bool m_print;
		size_t m_timestep_nr;
		size_t m_limiter;
		std::string m_neumannSubsets;
};

} // end namespace ug

// include implementation
#include "level_set_impl.h"

#endif /* LEVEL_SET_UTIL_H_ */
