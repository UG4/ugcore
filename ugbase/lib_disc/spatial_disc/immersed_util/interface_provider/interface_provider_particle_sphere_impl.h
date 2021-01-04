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
add(number radius, const MathVector<dim>& center, number density)
{
    // A. set given values:
        this->m_vCenter.push_back(center);
        m_vRadius.push_back(radius);
        this->m_vDensity.push_back(density);
        std::vector<bool> bDummy(2,false);
        this->m_vGivenVelocity.push_back(bDummy);
        
    // B. prepare data for free particle velocities, to be computed
    // B1: resize velocity vectors for TimeSeries-Data
        const size_t numTimePoints = 2;
        this->m_vvLinearVelocity.resize(numTimePoints);
        this->m_vvAngularVelocity.resize(numTimePoints);
        
    // B2: set dummy values if velocity of particle is free within fluid
        const MathVector<dim> vDummy(0.0);
        for ( size_t i = 0; i < numTimePoints; ++i )
        {
            this->m_vvLinearVelocity[i].push_back(vDummy);
            this->m_vvAngularVelocity[i].push_back(vDummy);
        }
        
}
 
template<int TWorldDim>
void ParticleProviderSphere<TWorldDim>::
add_moving(number radius, const MathVector<dim>& center, number density,
            const MathVector<dim>& linearVel, const MathVector<dim>& angularVel)
{
    // A. set given values:
        this->m_vCenter.push_back(center);
        m_vRadius.push_back(radius);
        this->m_vDensity.push_back(density);
        std::vector<bool> dummy(2,true);
        this->m_vGivenVelocity.push_back(dummy);
        
    // B. prepare data to set particle velocities
    // B1: resize velocity vectors for TimeSeries-Data
        const size_t numTimePoints = 2;
        this->m_vvLinearVelocity.resize(numTimePoints);
        this->m_vvAngularVelocity.resize(numTimePoints);
        
    // B2: set particle velocities
        for ( size_t i = 0; i < numTimePoints; ++i )
        {
            this->m_vvLinearVelocity[i].push_back(linearVel);
            this->m_vvAngularVelocity[i].push_back(angularVel);
        }
}
    
template<int TWorldDim>
number ParticleProviderSphere<TWorldDim>::
get_LSvalue_byPosition(MathVector<dim> vrtPos)
{
    const int prtIndex = this->get_prtIndex();
    if ( prtIndex == -1 )
        UG_THROW("in 'ParticleProviderSphere:get_LSvalue_byPosition(vrtPos)': call of get_LSvalue_byPosition() with prtIndex = "
                     << prtIndex << " not allowed!\n");
        
    return get_LSvalue_byPosition(vrtPos, prtIndex);
}
    
template<int TWorldDim>
number ParticleProviderSphere<TWorldDim>::
get_LSvalue_byPosition(MathVector<dim> vrtPos, const int prtIndex)
{
    if ( prtIndex == -1 )
        UG_THROW("in 'ParticleProviderSphere:get_LSvalue_byPosition(vrtPos, prtIndex)': call of get_LSvalue_byPosition() with prtIndex = "
                     << prtIndex << " not allowed!\n");
    
    const number radius = get_radius(prtIndex);
    const MathVector<dim>& center = this->get_center(prtIndex);
    
    MathVector<dim> localCoords;
    VecSubtract(localCoords, vrtPos, center);
        
    number dist = VecDot(localCoords, localCoords);
    dist = sqrt(dist);
        
    return radius - dist;
    
}
    
template<int TWorldDim>
void ParticleProviderSphere<TWorldDim>::
update(number deltaT, const MathVector<dim> transSol, const MathVector<dim> rotSol, const int prtIndex)
{
    MathVector<dim> centerNew = this->get_center(prtIndex);
    VecScaleAdd(centerNew, 1.0, centerNew, deltaT, transSol);
    this->set_center(centerNew, prtIndex);
}
    

    
template<int TWorldDim>
number ParticleProviderSphere<TWorldDim>::
Volume(int levIndex, size_t prtIndex)
{
    const number radius = get_radius(prtIndex);
    number volume = 0.0;
    
    if ( dim == 2 )
        volume = 3.1415*radius*radius;
    if ( dim == 3 )
        volume = 4.0/3.0*3.1415*radius*radius*radius;
    
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
    number momOfInertia = 0.0;
            
    if ( dim == 2 ) momOfInertia = 0.5*radius*radius*Mass(levIndex, prtIndex, fluidDensity);
    if ( dim == 3 ) momOfInertia = 0.4*radius*radius*Mass(levIndex, prtIndex, fluidDensity);
            
    return momOfInertia;
}
        
template<int TWorldDim>
number ParticleProviderSphere<TWorldDim>::
MomOfInertia(const int levIndex, const int prtIndex, const number volume, const number fluidDensity)
{
    const number radius = get_radius(prtIndex);
    number momOfInertia = 0.0;
    
    if ( dim == 2 ) momOfInertia = 0.5*radius*radius*Mass(levIndex, prtIndex, volume, fluidDensity);
    if ( dim == 3 ) momOfInertia = 0.4*radius*radius*Mass(levIndex, prtIndex, volume, fluidDensity);
            
    return momOfInertia;
}
    
    
template<int TWorldDim>
bool ParticleProviderSphere<TWorldDim>::
get_LineCircle_Intersection(MathVector<dim>& Intersect, const MathVector<dim>& vrtPosOut,
                            const MathVector<dim>& vrtPosIn, int PrtIndex)
{
    std::vector<number> alphaOut;
    alphaOut.clear();
        
    return get_LineCircle_Intersection(Intersect, vrtPosOut, vrtPosIn, PrtIndex, alphaOut);
}
    
template<int TWorldDim>
bool ParticleProviderSphere<TWorldDim>::
get_LineCircle_Intersection(MathVector<dim>& Intersect, const MathVector<dim>& vrtPosOut,
                            const MathVector<dim>& vrtPosIn, int prtIndex, std::vector<number>& alphaOut)
{
        
    ///////////////////////////////////////////////////////////////////////////////////////
    //
    // 'vrtPosOut':= the starting point of the ray,
    // 'vrtPosIn' := the end point of the ray,
    // 'center'   := the center of sphere you're testing against
    // 'radius'   := the radius of that sphere
    // 'lineDir'  := direction vector of ray from start to end
    // 'rayDir'   := direction vector of ray from center to start
    //
    // 	Ansatz:
    //  (1) Intersect = vrtPosOut + alpha*lineDir
    //	(2) Intersect - center = radius
    //
    // 		=> (1)  Intersect[0] = vrtPosOut[0] + alpha*lineDir[0]
    // 			    Intersect[1] = vrtPosOut[1] + alpha*lineDir[1]
    // 		=> (2)  (Intersect[0] - center[0])^2 + (Intersect[1] - center[1])^2 = radius^2
    //
    // 	Plug (1) into (2) => ... => quadratic equation for alpha:
    //
    //	alpha^2 * <lineDir,lineDir> + alpha * 2*<lineDir, rayDir> + ( <rayDir, rayDir>-radius^2 ) = 0
    //
    //		-> a =  <lineDir,lineDir>, b =  <lineDir,linrayDireDir>, c =  <rayDir,rayDir> - radius^2
    //
    //
    // 	=> from (1): Intersect = vrtPosOut + alpha*(vrtPosIn - vrtPosOut)
    ///////////////////////////////////////////////////////////////////////////////////////
        
        number alpha;
        
        const number radius = get_radius(prtIndex);
        const MathVector<dim>& center = this->get_center(prtIndex);
        
        MathVector<dim> lineDir;
        MathVector<dim> rayDir;
        
        // lineDir = vrtPosIn - vrtPosOut;
        VecSubtract(lineDir, vrtPosIn, vrtPosOut);
        // rayDir = vrtPosOut - center;
        VecSubtract(rayDir, vrtPosOut, center);
        
        const number a = VecDot(lineDir, lineDir);
        const number b = 2.0 * VecDot(lineDir, rayDir);
        const number c = VecDot(vrtPosOut, vrtPosOut) + VecDot(center, center)
        - 2 * VecDot(vrtPosOut, center) - radius * radius;
        
        const number discriminant = b * b - 4 * a * c;
        
        // check that 'vrtPosOut' and 'vrtPosIn' really lie on different sides of the circle:
        if (discriminant < -1e-8)
            UG_THROW("Value of discriminant = " << discriminant << "\n");
        
        // discriminant = 0!
        const number alpha1 = (-b - sqrt(discriminant)) / (2.0 * a);
        const number alpha2 = (-b + sqrt(discriminant)) / (2.0 * a);
        
        if (alpha1 <= alpha2)
            alpha = alpha1;
        else
            alpha = alpha2;
        
        if (alpha < 0 || (alpha - 1.0) > 1e-8)
            UG_THROW("Error in 'get_LineCircle_Intersection()': alpha not valid;"
                     "should lie between 0 and 1: " << alpha << "\n");
        
        for (size_t d = 0; d < dim; ++d)
            Intersect[d] = vrtPosOut[d] + alpha * lineDir[d];
        
        alphaOut.clear();
        alphaOut.push_back(alpha);
        alphaOut.push_back(1.0 - alpha);
    
        return true;
}
  
    
template<int TWorldDim>
void ParticleProviderSphere<TWorldDim>::
print()
{
        UG_LOG(" +++++++++++++++++++++++++ ParticleProviderSphere Info ++++++++++++++++++++++++++++++ \n");
        UG_LOG("+++ num_particles = " << this->num_particles() << "\n");
        for (size_t p = 0; p < this->num_particles(); ++p)
        {
            UG_LOG(" +++++++++++++++++++++++++ particle " << p << " ++++++++++++++++++++++++++++++ \n\n");

            UG_LOG("+++ center = " << this->get_center(p) << "\n");
            UG_LOG("+++ density = " << this->get_density(p) << "\n");
            UG_LOG("+++ radius = " << get_radius(p) << "\n\n");
            
            UG_LOG("+++ linear velocity: " << this->get_linear_velocity(p, 0) << "\n");
            UG_LOG("+++ angular velocity: " << this->get_angular_velocity(p, 0) << "\n\n");
           
        }
        UG_LOG(" ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ \n");
}
    
template<int TWorldDim>
void ParticleProviderSphere<TWorldDim>::
print_velocity(const MathVector<dim>& transSol, const MathVector<dim>& rotSol, const int prtIndex,
                   const bool isTimedep, const number time, const char* filename)
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
    
} // end namespace ug


#endif /* INTERFACE_PROVIDER_PARTICLE_SPHERE_H_ */
