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
    virtual ~DiffusionInterfaceProvider() {}

// setter methods
	void add(number radius, const MathVector<dim>& center)
	{
	// set given values:
 		this->m_vCenter.push_back(center);
		this->m_vRadius.push_back(radius);
		this->m_vDensity.push_back(1.0);
  	}

// getter methods
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

// output methods
	void print()
	{
		UG_LOG(" +++++++++++++++++++++++++ Interface Info ++++++++++++++++++++++++++++++ \n");
		UG_LOG("+++ num_circles = " << this->num_particles() << "\n");
		for (size_t p = 0; p < this->num_particles(); ++p)
		{
			UG_LOG("+++ center = " << this->get_center(p) << "\n");
			UG_LOG("+++ radius = " << this->get_radius(p) << "\n");
			UG_LOG("+++ density = " << this->get_density(p) << "\n");
		}
		UG_LOG(" ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ \n");
	}


};


}// end namespace ug


#endif /* INTERFACE_HANDLER_DIFFUSION_H_ */
