/*
 * moving_particle.h
 *
 *  Created on: 20.01.2015
 *      Author: suze
 */

#ifndef INTERFACE_PROVIDER_PARTICLE_SPHERE_H_
#define INTERFACE_PROVIDER_PARTICLE_SPHERE_H_


namespace ug{

template<int TWorldDim>
void ParticleProviderSphere<TWorldDim>::
print_velocity(const MathVector<dim>& transSol, const MathVector<dim>& rotSol, const int prtIndex, const bool isTimedep, const number time, const char* filename)
{
    MathVector<dim>  center = this->get_center(prtIndex);
    const number radius = get_radius(prtIndex);
        
    number f1 = 0.0;
    number gravity = -9.81;
    if ( !isTimedep )
        f1 = transSol[0]*4.0/(radius*radius*gravity);
        
    std::string name(filename);
    
    char * cstr = new char [name.size()+1];
    strcpy (cstr, name.c_str());
    
    if ( !isTimedep )
    {
        FILE* print_velocity = fopen(name.c_str(), "a");
        fprintf(print_velocity,"%e \t %e \t ",radius, f1);
        for ( int d = 0; d < dim; ++d ) fprintf(print_velocity,"%e \t %e \t ", transSol[d], rotSol[d]);
        for ( int d = 0; d < dim; ++d ) fprintf(print_velocity,"%e \t ", center[d]);
        fprintf(print_velocity," # radius, f1, transVel[0], rotVel[0], transVel[1], rotVel[1], center_coords, (m_bSharedIP = ");
        fprintf(print_velocity, "0 )  \n ");
        fclose(print_velocity);
    }
    else
    {
        FILE* print_velocity = fopen(name.c_str(), "a");
        fprintf(print_velocity,"%e \t ",time);
        for ( int d = 0; d < dim; ++d )		fprintf(print_velocity,"%e \t %e \t ", transSol[d], rotSol[d]);
        for ( int d = 0; d < dim; ++d )		fprintf(print_velocity,"%e \t ", center[d]);
        fprintf(print_velocity," # time, transVel[0], rotVel[0], transVel[1], rotVel[1], center_coords[0], center[1]");
        fprintf(print_velocity, " \n ");
        fclose(print_velocity);
    }
        
}
    
template<int TWorldDim>
number ParticleProviderSphere<TWorldDim>::
Volume(int levIndex, size_t prtIndex)
{
    const number radius = get_radius(prtIndex);
    number volume;
            
    /*	if ( m_bVolNumeric ){
    volume = compute_volume(Index, prtIndex);
    }
    else
    {*/
    if ( dim == 2 )
        volume = 3.1415*radius*radius;
    if ( dim == 3 )
        volume = 4.0/3.0*3.1415*radius*radius*radius;
    //	}
    UG_LOG("Volume: " << volume << "\n");
            
    return volume;
}
        
template<int TWorldDim>
number ParticleProviderSphere<TWorldDim>::
Mass(const int levIndex, const int prtIndex, const number fluidDensity)
{
    const number prtDensity = this->get_density(prtIndex);
    const number effective_density = prtDensity - fluidDensity;
    number mass = Volume(levIndex, prtIndex) * effective_density;
    return mass;
}
        
template<int TWorldDim>
number ParticleProviderSphere<TWorldDim>::
Mass(const int levIndex, const int prtIndex, const number volume, const number fluidDensity)
{
    const number prtDensity = this->get_density(prtIndex);
    const number effective_density = prtDensity - fluidDensity;
    number mass = volume * effective_density;
            
    return mass;
}
        
template<int TWorldDim>
number ParticleProviderSphere<TWorldDim>::
MomOfInertia(const int levIndex, const int prtIndex, const number fluidDensity)
{
    const number radius = get_radius(prtIndex);
    number momOfInertia;
            
    if ( dim == 2 ) momOfInertia = 0.5*radius*radius*Mass(levIndex, prtIndex, fluidDensity);
    if ( dim == 3 ) momOfInertia = 0.4*radius*radius*Mass(levIndex, prtIndex, fluidDensity);
            
    return momOfInertia;
}
        
template<int TWorldDim>
number ParticleProviderSphere<TWorldDim>::
MomOfInertia(const int levIndex, const int prtIndex, const number volume, const number fluidDensity)
{
    const number radius = get_radius(prtIndex);
    number momOfInertia;
    
    if ( dim == 2 ) momOfInertia = 0.5*radius*radius*Mass(levIndex, prtIndex, volume, fluidDensity);
    if ( dim == 3 ) momOfInertia = 0.4*radius*radius*Mass(levIndex, prtIndex, volume, fluidDensity);
            
    return momOfInertia;
}
    

} // end namespace ug


#endif /* INTERFACE_PROVIDER_PARTICLE_SPHERE_H_ */
