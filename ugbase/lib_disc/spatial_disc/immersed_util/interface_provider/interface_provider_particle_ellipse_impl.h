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
add(const MathVector<dim>& ellipticAxis, const MathVector<dim>& center, const number theta, number density)
{
    // A. set given values:
        this->m_vCenter.push_back(center);
        m_vEllipticAxis.push_back(ellipticAxis);
        m_vTheta.push_back(theta);
        
        this->m_vDensity.push_back(density);
        std::vector<bool> bDummy(2,false);
        this->m_vGivenVelocity.push_back(bDummy);
        
    // B. prepare data for free ellipse velocities, to be computed
    // B1: resize velocity vectors for TimeSeries-Data
        const size_t numTimePoints = 2;
        this->m_vvLinearVelocity.resize(numTimePoints);
        this->m_vvAngularVelocity.resize(numTimePoints);
        
    // B2: set dummy values if velocity of ellipse is free within fluid
        const MathVector<dim> vDummy(0.0);
        
        for ( size_t i = 0; i < numTimePoints; ++i )
        {
            this->m_vvLinearVelocity[i].push_back(vDummy);
            this->m_vvAngularVelocity[i].push_back(vDummy);
        }
        
}
    
template<int TWorldDim>
void ParticleProviderEllipse<TWorldDim>::
add_moving(const MathVector<dim>& ellipticAxis, const MathVector<dim>& center, const number theta, number density,
           const MathVector<dim>& linearVel, const MathVector<dim>& angularVel)
{
    // A. set given values:
        this->m_vCenter.push_back(center);
        m_vEllipticAxis.push_back(ellipticAxis);
        m_vTheta.push_back(theta);
        
        this->m_vDensity.push_back(density);
        std::vector<bool> dummy(2,true);
        this->m_vGivenVelocity.push_back(dummy);
        
    // B. prepare data to set ellipse velocities
    // B1: resize velocity vectors for TimeSeries-Data
        const size_t numTimePoints = 2;
        this->m_vvLinearVelocity.resize(numTimePoints);
        this->m_vvAngularVelocity.resize(numTimePoints);
        
    // B2: set ellipse velocities
        for ( size_t i = 0; i < numTimePoints; ++i )
        {
            this->m_vvLinearVelocity[i].push_back(linearVel);
            this->m_vvAngularVelocity[i].push_back(angularVel);
        }
}
    
    
// set vrtPos to origin and rotate by -theta
template<int TWorldDim>
void ParticleProviderEllipse<TWorldDim>::
rotate_vector(MathVector<dim>& vrtPos, const int prtIndex)
{
    const number theta = -1.0 * get_theta(prtIndex);
    const number radians = theta * PI / 180.0;
    const MathVector<dim> center = this->get_center(prtIndex);
        
// first move vtPos to origin:
    MathVector<dim> vrtPos_move;
    VecSubtract(vrtPos_move, vrtPos, center);
// buffer values in order to overwrite them:
    number vrtPos_x = vrtPos_move[0];
    number vrtPos_y = vrtPos_move[1];
        
// Rotate vrtPos by angle -theta:
    vrtPos_move[0] = vrtPos_x * cos(radians) - vrtPos_y * sin(radians);
    vrtPos_move[1] = vrtPos_x * sin(radians) + vrtPos_y * cos(radians);
        
// now move vtPos back to center:
    VecAdd(vrtPos, vrtPos_move, center);
}
    
// set vrtPos to center and rotate by theta
template<int TWorldDim>
void ParticleProviderEllipse<TWorldDim>::
rotate_vector_inverse(MathVector<dim>& vrtPos, const int prtIndex)
{
    const number theta = 1.0 * get_theta(prtIndex);
    const number radians = theta * PI / 180.0;
    const MathVector<dim> center = this->get_center(prtIndex);
        
// first move vtPos to origin:
    MathVector<dim> vrtPos_move;
    VecSubtract(vrtPos_move, vrtPos, center);
// buffer values in order to overwrite them:
    number vrtPos_x = vrtPos_move[0];
    number vrtPos_y = vrtPos_move[1];
        
// Rotate vrtPos by angle theta:
    vrtPos_move[0] = vrtPos_x * cos(radians) - vrtPos_y * sin(radians);
    vrtPos_move[1] = vrtPos_x * sin(radians) + vrtPos_y * cos(radians);
        
// now move vtPos to center:
    VecAdd(vrtPos, vrtPos_move, center);
        
}
    
    
template<int TWorldDim>
void ParticleProviderEllipse<TWorldDim>::
update(number deltaT, const MathVector<dim> transSol, const MathVector<dim> rotSol, const int prtIndex)
{
    //	update center of particle
        MathVector<dim> centerNew = this->get_center(prtIndex);
        VecScaleAdd(centerNew, 1.0, centerNew, deltaT, transSol);
        this->set_center(centerNew, prtIndex);
        
    //	update theta of particle
        number radianNew = get_theta(prtIndex)* PI / 180.0 + deltaT*rotSol[0];
        number thetaNew = radianNew * 180.0 / PI;
        set_theta(thetaNew, prtIndex);
        
}
    
template<int TWorldDim>
number ParticleProviderEllipse<TWorldDim>::
get_LSvalue_byPosition(MathVector<dim> vrtPos)
{
    const int prtIndex = this->get_prtIndex();
    if ( prtIndex == -1 )
        UG_THROW("in 'ParticleProviderEllipse:get_LSvalue_byPosition(vrtPos)': call of get_LSvalue_byPosition() with prtIndex = "
                    << prtIndex << " not allowed!\n");
        
    return get_LSvalue_byPosition(vrtPos, prtIndex);
        
}
    
template<int TWorldDim>
number ParticleProviderEllipse<TWorldDim>::
get_LSvalue_byPosition(MathVector<dim> vrtPos, const int prtIndex)
{
    if ( prtIndex == -1 )
        UG_THROW("in 'ParticleProviderEllipse:get_LSvalue_byPosition()': call of get_LSvalue_byPosition() with prtIndex = "
                    << prtIndex << " not allowed!\n");
// A) get center of ellipse as reference point
    const MathVector<dim> vrtPosIn = this->get_center(prtIndex);
    MathVector<dim> intersectionPnt;
    std::vector<number> alphaOut;
    alphaOut.resize(2,0.5);
        
// B) Rotate vrtPos by angle -theta and set center of ellipse to origin:
    rotate_vector(vrtPos, prtIndex);
        
// C) Compute LS_value for rotated vrtPos:
    number fake_LSvalue = IFF_LineEllipse_Intersection(intersectionPnt, vrtPos, vrtPosIn, prtIndex, alphaOut);
        
    return fake_LSvalue;
        
}

template<int TWorldDim>
bool ParticleProviderEllipse<TWorldDim>::
get_intersection_point(MathVector<dim>& Intersect, const MathVector<dim>& vrtPosOut,
                      const MathVector<dim>& vrtPosIn, const int PrtIndex)
{
// first rotate input points:
    MathVector<dim> vrtPosOut_rotated = vrtPosOut;
    MathVector<dim> vrtPosIn_rotated = vrtPosIn;
        
    rotate_vector(vrtPosOut_rotated, PrtIndex);
    rotate_vector(vrtPosIn_rotated, PrtIndex);
        
// returns intersection point already shifted back to real center!!
    get_LineEllipse_Intersection(Intersect, vrtPosOut_rotated, vrtPosIn_rotated, PrtIndex);
        
// NOW: rotate the intersection point BACK again!!
//      ---> (in 'get_LineEllipse_Intersection()' Intersect mit Ellipse through origin!
    rotate_vector_inverse(Intersect, PrtIndex);
        
    return true;
}
    

template<int TWorldDim>
bool ParticleProviderEllipse<TWorldDim>::
get_intersection_point(MathVector<dim>& Intersect, const MathVector<dim>& vrtPosOut,
                       const MathVector<dim>& vrtPosIn, const int PrtIndex, std::vector<number>& alphaOut)
{
// first rotate input points:
    MathVector<dim> vrtPosOut_rotated = vrtPosOut;
    MathVector<dim> vrtPosIn_rotated = vrtPosIn;
        
    rotate_vector(vrtPosOut_rotated, PrtIndex);
    rotate_vector(vrtPosIn_rotated, PrtIndex);
        
    // returns intersection point already shifted back to real center!!
    get_LineEllipse_Intersection(Intersect, vrtPosOut_rotated, vrtPosIn_rotated, PrtIndex, alphaOut);
        
    // NOW: rotate the intersection point BACK again!!
    //      ---> (in 'get_LineEllipse_Intersection()' Intersect mit Ellipse through origin!
    rotate_vector_inverse(Intersect, PrtIndex);
        
    return true;
}
 
template<int TWorldDim>
bool ParticleProviderEllipse<TWorldDim>::
get_LineEllipse_Intersection(MathVector<dim>& Intersect, const MathVector<dim>& vrtPosOut,
                             const MathVector<dim>& vrtPosIn, int PrtIndex)
{
    std::vector<number> alphaOut;
    alphaOut.clear();
    
    return get_LineEllipse_Intersection(Intersect, vrtPosOut, vrtPosIn, PrtIndex, alphaOut);
}
    
// Old version of 'get_LineEllipse_Intersection':
template<int TWorldDim>
bool ParticleProviderEllipse<TWorldDim>::
get_LineEllipse_Intersection(MathVector<dim>& Intersect, const MathVector<dim>& vrtPosOut,
                            const MathVector<dim>& vrtPosIn, int PrtIndex, std::vector<number>& alphaOut)
{
    ///////////////////////////////////////////////////////////////////////////////////////
    // The parameterized equation for a line segment through the points (x1, y1) and (x2, y2) is:
    //      x(t) = x1 + (x2-x1)t
    //      y(t) = y1 + (y2-y1)t
    // The following is the equation for an ellipse centered at the origin:
    //     (x/a)^2 + (y/b)^2 - 1 = 0
    // with a,b = elliptic axis
    
    // If you plug the equations for the line into the equation for the ellipse
    // AND after some computations, we get a quadratic equation of the form:
    //  At^2 + Bt + C = 0
    // with A = ..., B = ..., C = ... ---> see below in the code;)
    ///////////////////////////////////////////////////////////////////////////////////////
        
        number alpha;
        
        const number axis_x = get_elliptic_axis_x(PrtIndex);
        const number axis_y = get_elliptic_axis_y(PrtIndex);
        const MathVector<dim>& center = this->get_center(PrtIndex);
        
    // first move line by 'center'-vector, so that we can do the computations based on an ellipse with origin
        MathVector<dim> vrtPosIn_move;
        MathVector<dim> vrtPosOut_move;
        VecSubtract(vrtPosIn_move, vrtPosIn, center);
        VecSubtract(vrtPosOut_move, vrtPosOut, center);
        
    // lineDir = vrtPosIn_move - vrtPosOut_move;
        MathVector<dim> lineDir;
        VecSubtract(lineDir, vrtPosIn_move, vrtPosOut_move);
        
        const number a = (lineDir[0]*lineDir[0])/(axis_x*axis_x) + (lineDir[1]*lineDir[1])/(axis_y*axis_y);
        const number b = 2*vrtPosOut_move[0]*lineDir[0]/(axis_x*axis_x) + 2*vrtPosOut_move[1]*lineDir[1]/(axis_y*axis_y);
        const number c = (vrtPosOut_move[0]*vrtPosOut_move[0])/(axis_x*axis_x) + (vrtPosOut_move[1]*vrtPosOut_move[1])/(axis_y*axis_y) - 1.0;
        
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
            UG_THROW(
                     "Error in 'get_LineCircle_Intersection()': alpha not valid; should lie between 0 and 1: " << alpha << "\n");
        
        for (size_t d = 0; d < dim; ++d)
            Intersect[d] = vrtPosOut[d] + alpha * lineDir[d];
        
        alphaOut.clear();
        alphaOut.push_back(alpha);
        alphaOut.push_back(1.0 - alpha);
       
        return true;
        
}

    
    
// returns a value, which substitues the computation of an LSvalue, i.e.:
// return value > 0, if vrtPosOut on same side as vrtPosIn (= discriminant < 0)
// return value < 0, if vrtPosOut on other side as vrtPosIn (= discriminant > 0)
template<int TWorldDim>
number ParticleProviderEllipse<TWorldDim>::
IFF_LineEllipse_Intersection(MathVector<dim>& Intersect, const MathVector<dim>& vrtPosOut,
                            const MathVector<dim>& vrtPosIn, int PrtIndex, std::vector<number>& alphaOut)
{
    ///////////////////////////////////////////////////////////////////////////////////////
    // The parameterized equation for a line segment through the points (x1, y1) and (x2, y2) is:
    //      x(t) = x1 + (x2-x1)t
    //      y(t) = y1 + (y2-y1)t
    // The following is the equation for an ellipse centered at the origin:
    //     (x/axis_x)^2 + (y/axis_y)^2 - 1 = 0
    // with axis_x,axis_y = elliptic axis
    
    // If you plug the equations for the line into the equation for the ellipse
    // AND after some computations, we get a quadratic equation of the form:
    //  at^2 + bt + c = 0
    // with a = ..., b = ..., c = ... ---> see below in the code;)
    ///////////////////////////////////////////////////////////////////////////////////////
    
        const number axis_x = get_elliptic_axis_x(PrtIndex);
        const number axis_y = get_elliptic_axis_y(PrtIndex);
        const MathVector<dim>& center = this->get_center(PrtIndex);
        
    // first move line by 'center'-vector, so that we can do the computations based on an ellipse with origin
        MathVector<dim> vrtPosIn_move;
        MathVector<dim> vrtPosOut_move;
        MathVector<dim> center_move;
        VecSubtract(vrtPosIn_move, vrtPosIn, center);
        VecSubtract(vrtPosOut_move, vrtPosOut, center);
        VecSubtract(center_move, center, center);
        
        number factor1 = axis_x*axis_y;
        number factor2 = sqrt(axis_x*axis_x*vrtPosOut_move[1]*vrtPosOut_move[1] + axis_y*axis_y*vrtPosOut_move[0]*vrtPosOut_move[0]);
        number factor = factor1/factor2;
        
        Intersect[0] = factor * vrtPosOut_move[0];
        Intersect[1] = factor * vrtPosOut_move[1];
        
    // compute distances to center for all 3 points:
        number dist_vrtPosOut = VecDistance(vrtPosOut_move, center_move);
        number dist_Intersect = VecDistance(Intersect, center_move);
        
    // first check, that 'Intersect' does not lie ON the interface
    // ---> case relevant for check in 'CollectCorners_FlatTop_2d()
    //      ---> see 'interface_handler_particle_tools.h:CollectCorners_FlatTop_2d():
        /*
         // check for correct inersectionPnt
         if ( fabs(get_LSvalue_byPosition(intersectionPnt)) > 1e-6  )
         */
    
    const number threshold_max = -1e+6;
    if ( dist_Intersect < dist_vrtPosOut)
        return -1e+6;
    
    
    number dist1  = fabs(dist_Intersect - dist_vrtPosOut);
    number dist__ = VecDistance(vrtPosOut_move, Intersect);

    if ( fabs(dist1 - dist__) > 1e-6)
        UG_THROW("ohoh...in IFF_()\n");
    
    if ( fabs(dist__ - 0.001) < 1e-6)
        UG_LOG("stop!  IFF_()\n");

    return dist__;
    
    
    
    
    
        if ( fabs(dist_Intersect - dist_vrtPosOut) < 1e-6)
            return 0.0;
    
        if ( dist_Intersect > dist_vrtPosOut)
        {
            return -threshold_max;
        }
        
        return threshold_max;
        
}
    
    
template<int TWorldDim>
void ParticleProviderEllipse<TWorldDim>::
print()
{
        UG_LOG(" +++++++++++++++++++++++++ ParticleProviderEllipse Info ++++++++++++++++++++++++++++++ \n");
        UG_LOG("+++ num_particles = " << this->num_particles() << "\n");
        for (size_t p = 0; p < this->num_particles(); ++p)
        {
            UG_LOG(" +++++++++++++++++++++++++ particle " << p << " ++++++++++++++++++++++++++++++ \n\n");

            UG_LOG("+++ center = " << this->get_center(p) << "\n");
            UG_LOG("+++ density = " << this->get_density(p) << "\n");
            UG_LOG("+++ theta = " << get_theta(p) << "\n");
            UG_LOG("+++ elliptic_axis_x = " << get_elliptic_axis_x(p) << "\n");
            UG_LOG("+++ elliptic_axis_y = " << get_elliptic_axis_y(p) << "\n\n");
            
            UG_LOG("+++ linear velocity: " << this->get_linear_velocity(p, 0) << "\n");
            UG_LOG("+++ angular velocity: " << this->get_angular_velocity(p, 0) << "\n\n");
            
        }
        UG_LOG(" ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ \n");
}
    
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
