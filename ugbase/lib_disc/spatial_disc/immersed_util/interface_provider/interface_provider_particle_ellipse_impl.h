/*
 * moving_particle.h
 *
 *  Created on: 20.01.2015
 *      Author: suze
 */

#ifndef INTERFACE_PROVIDER_PARTICLE_ELLIPSE_H_
#define INTERFACE_PROVIDER_PARTICLE_ELLIPSE_H_


namespace ug{
  
template<int TWorldDim>
void ParticleProviderEllipse<TWorldDim>::
print_velocity(const MathVector<dim>& transSol, const MathVector<dim>& rotSol, const int prtIndex, const bool isTimedep, const number time, const char* filename)
{
    MathVector<dim>  center = this->get_center(prtIndex);
    
    std::string name(filename);
            
    char * cstr = new char [name.size()+1];
    strcpy (cstr, name.c_str());
            
    if ( !isTimedep )
    {
        FILE* print_velocity = fopen(name.c_str(), "a");
        for ( int d = 0; d < dim; ++d ) fprintf(print_velocity,"%e \t %e \t ", transSol[d], rotSol[d]);
            for ( int d = 0; d < dim; ++d ) fprintf(print_velocity,"%e \t ", center[d]);
                fprintf(print_velocity," %e \t # transVel[0], rotVel[0], transVel[1], rotVel[1], center_coords, theta, (m_bSharedIP = ", get_theta(prtIndex));
                fprintf(print_velocity, "0 )  \n ");
                fclose(print_velocity);
    }
    else
    {
        FILE* print_velocity = fopen(name.c_str(), "a");
        fprintf(print_velocity,"%e \t ",time);
        for ( int d = 0; d < dim; ++d )		fprintf(print_velocity,"%e \t %e \t ", transSol[d], rotSol[d]);
            for ( int d = 0; d < dim; ++d )		fprintf(print_velocity,"%e \t ", center[d]);
                    fprintf(print_velocity,"  %e \t # time, transVel[0], rotVel[0], transVel[1], rotVel[1], center_coords[0], center[1], theta", get_theta(prtIndex));
        fprintf(print_velocity, " \n ");
        fclose(print_velocity);
    }
        
}
    
template<int TWorldDim>
number ParticleProviderEllipse<TWorldDim>::
Volume(int levIndex, size_t prtIndex)
{
    const number axis_x = get_elliptic_axis_x(prtIndex);
    const number axis_y = get_elliptic_axis_y(prtIndex);
    
    number volume;
    
    if ( dim == 2 )
        volume = 3.1415*axis_x*axis_y;
    if ( dim == 3 )
        UG_THROW("in 'ParticleProviderEllipse::Volume()': 3d-case (i.e. axis_z) not implemented!\n");
    //volume = 4.0/3.0*3.1415*axis_x*axis_y*axis_z;
    
    UG_LOG("Volume: " << volume << "\n");
    
    return volume;
}
        
template<int TWorldDim>
number ParticleProviderEllipse<TWorldDim>::
Mass(const int levIndex, const int prtIndex, const number fluidDensity)
{
    const number prtDensity = this->get_density(prtIndex);
    const number effective_density = prtDensity - fluidDensity;
    number mass = Volume(levIndex, prtIndex) * effective_density;
    return mass;
}
        
template<int TWorldDim>
number ParticleProviderEllipse<TWorldDim>::
Mass(const int levIndex, const int prtIndex, const number volume, const number fluidDensity)
{
    const number prtDensity = this->get_density(prtIndex);
    const number effective_density = prtDensity - fluidDensity;
    number mass = volume * effective_density;
    
    return mass;
}
        
template<int TWorldDim>
number ParticleProviderEllipse<TWorldDim>::
MomOfInertia(const int levIndex, const int prtIndex, const number fluidDensity)
{
    const number axis_x = get_elliptic_axis_x(prtIndex);
    const number axis_y = get_elliptic_axis_y(prtIndex);
    
    number momOfInertia;
    
    if ( dim == 2 ) momOfInertia = 0.25*axis_x*axis_y*Mass(levIndex, prtIndex, fluidDensity);
    if ( dim == 3 )
        UG_THROW("in 'ParticleProviderEllipse::MomOfInertia()': 3d-case not implemented!\n");
    
    return momOfInertia;
}
        
template<int TWorldDim>
number ParticleProviderEllipse<TWorldDim>::
MomOfInertia(const int levIndex, const int prtIndex, const number volume, const number fluidDensity)
{
    const number axis_x = get_elliptic_axis_x(prtIndex);
    const number axis_y = get_elliptic_axis_y(prtIndex);
    
    number momOfInertia;
    
    // I_x = 0.25 * pi * a * b^3 = 0.25 * b^2 * Area
    // I_y = 0.25 * pi * b * a^3 = 0.25 * a^2 * Area
    // => I_origin = I_x + I_y = 0.25 * (a^2+b^2) * Area
    // => identical to MomOfInertia of Circle for a = b  = r
    if ( dim == 2 ) momOfInertia = 0.25*(axis_x*axis_x + axis_y*axis_y)*Mass(levIndex, prtIndex, volume, fluidDensity);
    if ( dim == 3 )
        UG_THROW("in 'ParticleProviderEllipse::MomOfInertia()': 3d-case not implemented!\n");
    
    return momOfInertia;
}
    

} // end namespace ug


#endif /* INTERFACE_PROVIDER_PARTICLE_ELLIPSE_H_ */
