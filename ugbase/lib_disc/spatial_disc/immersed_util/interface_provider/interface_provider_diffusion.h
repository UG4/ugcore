/*
 * interface_handler_local.h
 *
 *  Created on: 15.01.2015
 *      Author: suze
 */

#ifndef INTERFACE_HANDLER_DIFFUSION_H_
#define INTERFACE_HANDLER_DIFFUSION_H_

#include "interface_provider_particle.h"


namespace ug{


template <int TWorldDim>
class DiffusionInterfaceProvider : public ParticleProviderSphere<TWorldDim>
{

	public:
///	world Dimension
	static const int dim = TWorldDim;

/// default constructor:
	DiffusionInterfaceProvider()
	{
		this->clear();
		UG_LOG("DiffusionInterfaceProvider constructor\n");
	};

/// destructor
	~DiffusionInterfaceProvider() {}

// setter methods
	void add_with_solution(number radius, const MathVector<dim>& center, number density, number solution)
	{
	// A. set given values:
 		this->m_vCenter.push_back(center);
 		this->m_vRadius.push_back(radius);
 		this->m_vDensity.push_back(density);

	// B. prepare data for free particle velocities, to be computed
	// B1: resize velocity vectors for TimeSeries-Data
		const size_t numTimePoints = 2;
		m_vvSolution.resize(numTimePoints);

	// B2: set default values if velocity of particle is free within fluid
 		for ( size_t i = 0; i < numTimePoints; ++i )
			m_vvSolution[i].push_back(solution);
 	}

	void add(number radius, const MathVector<dim>& center, number density)
	{
	// A. set given values:
 		this->m_vCenter.push_back(center);
		this->m_vRadius.push_back(radius);
		this->m_vDensity.push_back(density);

	// B. prepare data to set particle velocities
	// B1: resize velocity vectors for TimeSeries-Data
		const size_t numTimePoints = 2;
		m_vvSolution.resize(numTimePoints);

	// B2: set particle velocities
		for ( size_t i = 0; i < numTimePoints; ++i )
			m_vvSolution[i].push_back(0.0);

  	}

	number get_radius(int prtIndex)
	{ if ( (int)this->num_particles() > 1 )
		UG_THROW("DiffusionInterfaceProvider::num_particles(): number of given particles = " << this->num_particles() << " more than 1! ... not supposed to be!\n");

	  return this->m_vRadius[prtIndex]; }

	number get_density(int prtIndex)
	{ if ( (int)this->num_particles() > 1 )
		UG_THROW("DiffusionInterfaceProvider::num_particles(): number of given particles = " << this->num_particles() << " more than 1! ... not supposed to be!\n");
	  return this->m_vDensity[prtIndex]; }

	MathVector<dim> get_center(int prtIndex)
	{ if ( (int)this->num_particles() > 1 )
		UG_THROW("DiffusionInterfaceProvider::num_particles(): number of given particles = " << this->num_particles() << " more than 1! ... not supposed to be!\n");
	  return this->m_vCenter[prtIndex]; }

// called during ParticleMapper::modify_GlobalSol():
	void set_solution(number solution, size_t prtIndex, int timeSeriesInd)
	{ m_vvSolution[timeSeriesInd][prtIndex] = solution; }

/// get solution values
	number get_solution(size_t prtIndex, size_t timeSeriesInd)
		{ return m_vvSolution[timeSeriesInd][prtIndex]; }

	void print()
	{
		UG_LOG(" +++++++++++++++++++++++++ Particle Info ++++++++++++++++++++++++++++++ \n");
		UG_LOG("+++ num_particles = " << this->num_particles() << "\n");
		for (size_t p = 0; p < this->num_particles(); ++p)
		{
			UG_LOG("+++ center = " << this->get_center(p) << "\n");
			UG_LOG("+++ radius = " << this->get_radius(p) << "\n");
			UG_LOG("+++ density = " << this->get_density(p) << "\n");

			UG_LOG("+++ solution set to: " << get_solution(p, 0) << "\n\n");
		}
		UG_LOG(" ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ \n");
	}


	protected:

// 	m_vvSolution[i][j]: i := timeSeries-index; j := particle-index
	// ToDo switch indexing: [prtIndex][TimeSeriesIndes] !!
	std::vector<std::vector<number> > m_vvSolution;


};


}// end namespace ug


#endif /* INTERFACE_HANDLER_DIFFUSION_H_ */
