/*
 * moving_particle.h
 *
 *  Created on: 20.01.2015
 *      Author: suze
 */

#ifndef INTERFACE_PROVIDER_PARTICLE_H_
#define INTERFACE_PROVIDER_PARTICLE_H_


#ifdef UG_PARALLEL
 	#include "lib_grid/parallelization/load_balancer_util.h"
#endif


namespace ug{

class IInterfaceProvider
{
	public:

/// default constructor:
	IInterfaceProvider(){};

/// destructor

	~IInterfaceProvider() {}

};

template <int TWorldDim>
class ParticleProvider : public IInterfaceProvider
{

	public:
///	world Dimension
	static const int dim = TWorldDim;

/// default constructor:
	ParticleProvider()
	{
		clear();
		UG_LOG("ParticleProvider constructor\n");
	};

/// destructor
	~ParticleProvider() {}

	void set_density(number density, int prtIndex)
	{
		if ( (int)num_particles() < prtIndex )
			UG_THROW("ParticleProvider::set_density(): number of given particles = "
					<< num_particles() << " smaller than given prtIndex = " << prtIndex << "\n");
		m_vDensity[prtIndex] = density;
    }
	void set_center(MathVector<dim> center, int prtIndex)
	{
		if ( (int)num_particles() < prtIndex )
			UG_THROW("ParticleProvider::set_center(): number of given particles = "
					<< num_particles() << " smaller than given prtIndex = " << prtIndex << "\n");
		for ( size_t d = 0; d < dim; ++d )
			m_vCenter[prtIndex][d] = center[d];
	}

	void clear()
    {
        UG_LOG("clear particles!\n");
        m_vCenter.clear(); m_vDensity.clear(); m_vGivenVelocity.clear();
				  m_vvLinearVelocity.clear(); m_vvAngularVelocity.clear();
    }


// getter methods
	size_t num_particles() const{ return m_vCenter.size(); }

	number get_density(int prtIndex){ return m_vDensity[prtIndex]; }
	MathVector<dim> get_center(int prtIndex){ return m_vCenter[prtIndex]; }
	bool get_DoF_modus_linear(int prtIndex) { return m_vGivenVelocity[prtIndex][0]; }
	bool get_DoF_modus_angular(int prtIndex) { return m_vGivenVelocity[prtIndex][1]; }

// called during PartcleMapper::modify_GlobalSol():
	void set_linear_velocity(number solution, size_t prtIndex, int timeSeriesInd, int cmp)
	{ m_vvLinearVelocity[timeSeriesInd][prtIndex][cmp] = solution; }
	void set_linear_velocity(MathVector<dim> solution, size_t prtIndex, int timeSeriesInd)
	{ m_vvLinearVelocity[timeSeriesInd][prtIndex] = solution; }
	void set_angular_velocity(number solution, size_t prtIndex, int timeSeriesInd, int cmp)
	{ m_vvAngularVelocity[timeSeriesInd][prtIndex][cmp] = solution; }
	void set_angular_velocity(MathVector<dim> solution, size_t prtIndex, int timeSeriesInd)
	{ m_vvAngularVelocity[timeSeriesInd][prtIndex] = solution; }

/// get solution values
	MathVector<dim> get_linear_velocity(size_t prtIndex, size_t timeSeriesInd)
		{ return m_vvLinearVelocity[timeSeriesInd][prtIndex]; }
	MathVector<dim> get_angular_velocity(size_t prtIndex, size_t timeSeriesInd)
		{ return m_vvAngularVelocity[timeSeriesInd][prtIndex]; }

    virtual number get_LSvalue_byPosition(MathVector<dim> vrtPos, const int prtIndex)
    { UG_THROW("in 'ParticleProvider::get_LSvalue_byPosition()': needs to be implemented by derived class!\n");}

    virtual bool get_intersection_point(MathVector<dim>& Intersect, const MathVector<dim>& vrtPosOut, const MathVector<dim>& vrtPosIn, int PrtIndex)
    { UG_THROW("in 'ParticleProvider::get_intersection_point()': needs to be implemented by derived class!\n");}

    virtual bool get_intersection_point(MathVector<dim>& Intersect, const MathVector<dim>& vrtPosOut, const MathVector<dim>& vrtPosIn, int PrtIndex, std::vector<number>& alphaOut)
    { UG_THROW("in 'ParticleProvider::get_intersection_point()': needs to be implemented by derived class!\n");}


// careful! --> only dummy implementation, sice some methods in MovingParticle-class use it; in these cases it will be overwritten by
// the ParticleProviderSphere-class (hopefully) ;)
   virtual number get_radius(int prtIndex)
    { UG_THROW("in 'ParticleProvider::get_radius()': needs to be implemented by derived class!\n");}

    number get_theta(int prtIndex)
    { UG_THROW("in 'ParticleProvider::get_theta()': needs to be implemented by derived class!\n");}
    number set_theta(number theta, int prtIndex)
    { UG_THROW("in 'ParticleProvider::get_theta()': needs to be implemented by derived class!\n");}

    void print()
    { UG_THROW("in 'ParticleProvider::print()': needs to be implemented by derived class!\n");}

    /// methods called by local_to_global_mappe:
    virtual number Volume(int levIndex, size_t prtIndex)
    { UG_THROW("in 'ParticleProvider::Volume()': needs to be implemented by derived class!\n");}
    virtual number Mass(const int levIndex, const int prtIndex, const number fluidDensity)
    { UG_THROW("in 'ParticleProvider::Mass()': needs to be implemented by derived class!\n");}
    virtual number Mass(const int levIndex, const int prtIndex, const number volume, const number fluidDensity)
    { UG_THROW("in 'ParticleProvider::Mass()': needs to be implemented by derived class!\n");}
    virtual number MomOfInertia(const int levIndex, const int prtIndex, const number fluidDensity)
    { UG_THROW("in 'ParticleProvider::MomOfInertia()': needs to be implemented by derived class!\n");}
    virtual number MomOfInertia(const int levIndex, const int prtIndex, const number volume, const number fluidDensity)
    { UG_THROW("in 'ParticleProvider::MomOfInertia()': needs to be implemented by derived class!\n");}

    virtual void print_velocity(const MathVector<dim>& transSol, const MathVector<dim>& rotSol, const int prtIndex, const bool isTimedep, const number time, const char* filename)
    { UG_THROW("in 'ParticleProvider::print_velocity()': needs to be implemented by derived class!\n");}

    virtual void update(number deltaT, const MathVector<dim> transSol, const MathVector<dim> rotSol, const int prtIndex)
    { UG_THROW("in 'ParticleProvider::update()': needs to be implemented by derived class!\n");}

protected:
    std::vector<number> m_vDensity;
    std::vector<MathVector<dim> > m_vCenter;
    // 	m_vLinearVelocity[i][j]: i := timeSeries-index; j := particle-index
    // ToDo switch indexing: [prtIndex][TimeSeriesIndes] !!
    std::vector<std::vector<MathVector<dim> > > m_vvLinearVelocity;
    std::vector<std::vector<MathVector<dim> > > m_vvAngularVelocity;
    // [prtIndex][linear/angular]
    std::vector<std::vector<bool> > m_vGivenVelocity;
    
};
    
    
template <int TWorldDim>
class ParticleProviderSphere : public ParticleProvider<TWorldDim>
{
        
public:
    ///	world Dimension
    static const int dim = TWorldDim;
        
    /// default constructor:
    ParticleProviderSphere() : ParticleProvider<dim>()
    { UG_LOG("ParticleProviderSphere constructor\n") };
        
    /// destructor
    ~ParticleProviderSphere() {}
        
    // setter methods
    void add(number radius, const MathVector<dim>& center, number density)
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
        
    void add_moving(number radius, const MathVector<dim>& center, number density,
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

    void update(number deltaT, const MathVector<dim> transSol, const MathVector<dim> rotSol, const int prtIndex)
    {
        //	update center of particle
        MathVector<dim> centerNew = this->get_center(prtIndex);
        VecScaleAdd(centerNew, 1.0, centerNew, deltaT, transSol);
        this->set_center(centerNew, prtIndex);
    }
    
///////////////////////////////////////////////////////////////////////////////////////
// new setter methods
///////////////////////////////////////////////////////////////////////////////////////
    
    void set_radius(number radius, int prtIndex)
    {
        if ( (int)this->num_particles() < prtIndex )
            UG_THROW("ParticleProviderSphere::set_radius(): number of given particles = "
                        << this->num_particles() << " smaller than given prtIndex = " << prtIndex << "\n");
        m_vRadius[prtIndex] = radius;
    }
    
///////////////////////////////////////////////////////////////////////////////////////
// new getter methods
///////////////////////////////////////////////////////////////////////////////////////

    number get_radius(int prtIndex){ return m_vRadius[prtIndex]; }

    
///////////////////////////////////////////////////////////////////////////////////////
// methods from base class (necessary to be implemented!)
///////////////////////////////////////////////////////////////////////////////////////

    
    /// methods called by local_to_global_mappe:
    number Volume(int levIndex, size_t prtIndex);
    number Mass(const int levIndex, const int prtIndex, const number fluidDensity);
    number Mass(const int levIndex, const int prtIndex, const number volume, const number fluidDensity);
    number MomOfInertia(const int levIndex, const int prtIndex, const number fluidDensity);
    number MomOfInertia(const int levIndex, const int prtIndex, const number volume, const number fluidDensity);
  
    void print_velocity(const MathVector<dim>& transSol, const MathVector<dim>& rotSol, const int prtIndex, const bool isTimedep, const number time, const char* filename);

    number get_LSvalue_byPosition(MathVector<dim> vrtPos, const int prtIndex)
    {
      //  if ( prtIndex == 1 )
     //       UG_LOG("jihaa!...\n");
        
        const number radius = get_radius(prtIndex);
        const MathVector<dim>& center = this->get_center(prtIndex);
        
        MathVector<dim> localCoords;
        VecSubtract(localCoords, vrtPos, center);
        
        number dist = VecDot(localCoords, localCoords);
        dist = sqrt(dist);
        
        return radius - dist;
    }
    
    bool get_intersection_point(MathVector<dim>& Intersect, const MathVector<dim>& vrtPosOut, const MathVector<dim>& vrtPosIn, int PrtIndex)
    { return get_LineCircle_Intersection(Intersect, vrtPosOut, vrtPosIn, PrtIndex);}
    
    bool get_intersection_point(MathVector<dim>& Intersect, const MathVector<dim>& vrtPosOut, const MathVector<dim>& vrtPosIn, int PrtIndex, std::vector<number>& alphaOut)
    { return get_LineCircle_Intersection(Intersect, vrtPosOut, vrtPosIn, PrtIndex, alphaOut);}

    
 	bool get_LineCircle_Intersection(MathVector<dim>& Intersect, const MathVector<dim>& vrtPosOut,
							   const MathVector<dim>& vrtPosIn, int PrtIndex)
	{
	//	if ( dim == 3 ) UG_THROW("in 'get_LineCircle_Intersection()': not implemented for 3d case!\n");

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
		//			alpha^2 * <lineDir,lineDir> + alpha * 2*<lineDir, rayDir> + ( <rayDir, rayDir>-radius^2 ) = 0
		//
		// 			-> a =  <lineDir,lineDir>, b =  <lineDir,linrayDireDir>, c =  <rayDir,rayDir> - radius^2
		//
		//
		// 	=> from (1): Intersect = vrtPosOut + alpha*(vrtPosIn - vrtPosOut)
		///////////////////////////////////////////////////////////////////////////////////////

		number alpha;

		const number radius = get_radius(PrtIndex);
		const MathVector<dim>& center = this->get_center(PrtIndex);

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
			UG_THROW(
					"Error in 'get_LineCircle_Intersection()': alpha not valid; should lie between 0 and 1: " << alpha << "\n");

		for (size_t d = 0; d < dim; ++d)
			Intersect[d] = vrtPosOut[d] + alpha * lineDir[d];

		return true;
	}

  	bool get_LineCircle_Intersection(MathVector<dim>& Intersect, const MathVector<dim>& vrtPosOut,
							   const MathVector<dim>& vrtPosIn, int PrtIndex, std::vector<number>& alphaOut)
	{
		if ( (int)this->num_particles() > 1 )
			UG_THROW("get_LineCircle_Intersection(): number of given particles = " << this->num_particles() << " more than 1! ... not supposed to be!\n");

		//if ( dim == 3 ) UG_THROW("in 'get_LineCircle_Intersection()': not implemented for 3d case!\n");

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
		//			alpha^2 * <lineDir,lineDir> + alpha * 2*<lineDir, rayDir> + ( <rayDir, rayDir>-radius^2 ) = 0
		//
		// 			-> a =  <lineDir,lineDir>, b =  <lineDir,linrayDireDir>, c =  <rayDir,rayDir> - radius^2
		//
		//
		// 	=> from (1): Intersect = vrtPosOut + alpha*(vrtPosIn - vrtPosOut)
		///////////////////////////////////////////////////////////////////////////////////////

		number alpha;

		const number radius = get_radius(PrtIndex);
		const MathVector<dim>& center = this->get_center(PrtIndex);

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
			UG_THROW(
					"Error in 'get_LineCircle_Intersection()': alpha not valid; should lie between 0 and 1: " << alpha << "\n");

		for (size_t d = 0; d < dim; ++d)
			Intersect[d] = vrtPosOut[d] + alpha * lineDir[d];

		alphaOut.clear();
		alphaOut.push_back(alpha);
		alphaOut.push_back(1.0 - alpha);

		UG_LOG("alphaOut 0 = " << alphaOut[0] << "\n");
		UG_LOG("alphaOut 1 = " << alphaOut[1] << "\n");

		return true;
	}
    
///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////


    void print()
    {
        UG_LOG(" +++++++++++++++++++++++++ ParticleProviderSphere Info ++++++++++++++++++++++++++++++ \n");
        UG_LOG("+++ num_particles = " << this->num_particles() << "\n");
        for (size_t p = 0; p < this->num_particles(); ++p)
        {
            UG_LOG("+++ center = " << this->get_center(p) << "\n");
            UG_LOG("+++ density = " << this->get_density(p) << "\n");
            UG_LOG("+++ radius = " << get_radius(p) << "\n");

            if ( this->get_DoF_modus_linear(p) )
            {UG_LOG("+++ linear velocity set to: " << this->get_linear_velocity(p, 0) << "\n\n");}
            else
            {UG_LOG("+++ linear velocity FREE and set to: " << this->get_linear_velocity(p, 0) << "\n\n");}
            if ( this->get_DoF_modus_angular(p) )
            {UG_LOG("+++ angular velocity set to: " << this->get_angular_velocity(p, 0) << "\n\n");}
            else
            {UG_LOG("+++ lineangularar velocity FREE and set to: " << this->get_angular_velocity(p, 0) << "\n\n");}
            
        }
        UG_LOG(" ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ \n");
    }

	protected:
	std::vector<number> m_vRadius;


};

template <int TWorldDim>
class ParticleProviderEllipse : public ParticleProvider<TWorldDim>
{
        
public:
    ///	world Dimension
    static const int dim = TWorldDim;
        
    /// default constructor:
    ParticleProviderEllipse() : ParticleProvider<dim>()
    {
        this->clear();
        UG_LOG("ParticleProviderEllipse constructor\n");
        
    };
        
    /// destructor
    ~ParticleProviderEllipse() {}
        
    // setter methods
    void add(const MathVector<dim>& ellipticAxis, const MathVector<dim>& center, const number theta, number density)
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
        
    void add_moving(const MathVector<dim>& ellipticAxis, const MathVector<dim>& center, const number theta, number density,
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
 
    void update(number deltaT, const MathVector<dim> transSol, const MathVector<dim> rotSol, const int prtIndex)
    {
        //	update center of particle
        MathVector<dim> centerNew = this->get_center(prtIndex);
        VecScaleAdd(centerNew, 1.0, centerNew, deltaT, transSol);
        this->set_center(centerNew, prtIndex);
        
        UG_LOG("in update: theta vorher: " << get_theta(prtIndex) << "\n");
        
        //	update theta of particle
        UG_LOG("in update: thetaOld: " << get_theta(prtIndex) << "\n");
        UG_LOG("in update: radianOld: " << get_theta(prtIndex)* PI / 180.0 << "\n");

        number radianNew = get_theta(prtIndex)* PI / 180.0 + deltaT*rotSol[0];
        number thetaNew = radianNew * 180.0 / PI;
        set_theta(thetaNew, prtIndex);
        
        UG_LOG("in update: thetaNew: " << thetaNew << "\n");
        UG_LOG("in update: radianNew: " << radianNew << "\n");
        
        UG_LOG("in update: theta nachher: " << get_theta(prtIndex) << "\n");

    }
    
///////////////////////////////////////////////////////////////////////////////////////
// new setter methods
///////////////////////////////////////////////////////////////////////////////////////
 
    void set_elliptic_axis(MathVector<dim> ellitpicAxis, int prtIndex)
    {
        if ( (int)this->num_particles() < prtIndex )
            UG_THROW("EllipseProvider::set_center(): number of given particles = "
                         << this->num_particles() << " smaller than given prtIndex = " << prtIndex << "\n");
        for ( size_t d = 0; d < dim; ++d )
            m_vEllipticAxis[prtIndex][d] = ellitpicAxis[d];
    }

    void set_theta(number theta, int prtIndex)
    {
        if ( (int)this->num_particles() < prtIndex )
            UG_THROW("EllipseProvider::set_theta(): number of given particles = "
                         << this->num_particles() << " smaller than given prtIndex = " << prtIndex << "\n");
        m_vTheta[prtIndex] = theta;
    }

    void print_velocity(const MathVector<dim>& transSol, const MathVector<dim>& rotSol, const int prtIndex, const bool isTimedep, const number time, const char* filename);

///////////////////////////////////////////////////////////////////////////////////////
// new getter methods
///////////////////////////////////////////////////////////////////////////////////////

    MathVector<dim> get_elliptic_axis(int prtIndex){ return m_vEllipticAxis[prtIndex]; }
    number get_elliptic_axis_x(int prtIndex){ return m_vEllipticAxis[prtIndex][0]; }
    number get_elliptic_axis_y(int prtIndex){ return m_vEllipticAxis[prtIndex][1]; }
    number get_theta(int prtIndex){ return m_vTheta[prtIndex]; }
    
///////////////////////////////////////////////////////////////////////////////////////
// methods from base class (necessary to be implemented!)
///////////////////////////////////////////////////////////////////////////////////////

    /// methods called by local_to_global_mappe:
    number Volume(int levIndex, size_t prtIndex);
    number Mass(const int levIndex, const int prtIndex, const number fluidDensity);
    number Mass(const int levIndex, const int prtIndex, const number volume, const number fluidDensity);
    number MomOfInertia(const int levIndex, const int prtIndex, const number fluidDensity);
    number MomOfInertia(const int levIndex, const int prtIndex, const number volume, const number fluidDensity);
  
    
    // set vrtPos to origin and rotate by -theta
    void rotate_vector(MathVector<dim>& vrtPos, const int prtIndex)
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
    void rotate_vector_inverse(MathVector<dim>& vrtPos, const int prtIndex)
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
    
    number get_LSvalue_byPosition(MathVector<dim> vrtPos, const int prtIndex)
    {
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
    
    bool get_intersection_point(MathVector<dim>& Intersect, const MathVector<dim>& vrtPosOut, const MathVector<dim>& vrtPosIn, int PrtIndex)
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
    
    bool get_intersection_point(MathVector<dim>& Intersect, const MathVector<dim>& vrtPosOut, const MathVector<dim>& vrtPosIn, int PrtIndex, std::vector<number>& alphaOut)
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

    bool get_LineEllipse_Intersection(MathVector<dim>& Intersect, const MathVector<dim>& vrtPosOut,
                                      const MathVector<dim>& vrtPosIn, int PrtIndex)
    {
        std::vector<number> alphaOut; alphaOut.clear();
        get_LineEllipse_Intersection(Intersect, vrtPosOut, vrtPosIn, PrtIndex, alphaOut);
    }

// Old version of 'get_LineEllipse_Intersection':
    bool get_LineEllipse_Intersection(MathVector<dim>& Intersect, const MathVector<dim>& vrtPosOut,
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
        
    //    UG_LOG("alphaOut 0 = " << alphaOut[0] << "\n");
    //    UG_LOG("alphaOut 1 = " << alphaOut[1] << "\n");
        
        return true;
        
    }
    
    
  /*  bool get_LineEllipse_Intersection(MathVector<dim>& Intersect, const MathVector<dim>& vrtPosOut,
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
        
        MathVector<dim> lineDir;
        VecSubtract(lineDir, vrtPosIn, vrtPosOut);
        
    ///////////////////////////////////////////////////////////////////////////////////////
    //  Trick necessary, in order to apply computation of Intersect by unsing the above formular:
    //  The center needs to lie on the connecting line of vrtPosOut, vrtPosIn and Intersect!
    //      ==> projection of center onto lineDir necessary
    ///////////////////////////////////////////////////////////////////////////////////////

        // project center onto 'lineDir'
        MathVector<dim> center_project = this->get_center(PrtIndex);
        number alpha_project;
        
        // default: project in x-direction:
        size_t projCoord = 0;
        size_t keepCoord = 1;
        
        // BUT: check the following, in order to avoid, that projected center lies 'between' vrtPosOut and vrtPosIn:
        if ( vrtPosIn[keepCoord] < center_project[keepCoord] < vrtPosOut[keepCoord] )
        {
            projCoord = 1;
            keepCoord = 0;
        }
        else if ( vrtPosOut[keepCoord] < center_project[keepCoord] < vrtPosIn[keepCoord] )
        {
            projCoord = 1;
            keepCoord = 0;
        }
    
    // derive alpha_project and compute center_project:
        alpha_project = (vrtPosIn[keepCoord] - vrtPosOut[keepCoord])/lineDir[keepCoord];
        center_project[projCoord] = vrtPosOut[projCoord] + alpha_project * lineDir[projCoord];

    // tag: 'equationProj': == center_project[projCoord] = center[projCoord] + deltaCoord:
        number deltaCoord = center_project[projCoord] - center[projCoord];
        
    // apply deltaCoord-move to vrtPosIn and vrtPosOut:
    // see: 'equationProj'
        MathVector<dim> vrtPosIn_move;
        vrtPosIn_move[projCoord] = vrtPosIn[projCoord] + deltaCoord;
        vrtPosIn_move[keepCoord] = vrtPosIn[keepCoord];

        MathVector<dim> vrtPosOut_move;
        vrtPosOut_move[projCoord] = vrtPosOut[projCoord] + deltaCoord;
        vrtPosOut_move[keepCoord] = vrtPosOut[keepCoord];
        
    ///////////////////////////////////////////////////////////////////////////////////////
    // finally, the usual computations can start,
    //  using 'vrtPosOut_move', 'vrtPosIn_move' and 'center_project':
    ///////////////////////////////////////////////////////////////////////////////////////

        number factor1 = axis_x*axis_y;
        number factor2 = sqrt(axis_x*axis_x*vrtPosOut_move[1]*vrtPosOut_move[1] + axis_y*axis_y*vrtPosOut_move[0]*vrtPosOut_move[0]);
        number factor = factor1/factor2;
        
        Intersect[0] = factor * vrtPosOut_move[0];
        Intersect[1] = factor * vrtPosOut_move[1];
        
        // compute distances to center for all 3 points:
        number dist_vrtPosOut = VecDistance(vrtPosOut_move, center_project);
        number dist_Intersect = VecDistance(Intersect, center_project);
        
        UG_LOG("dist_Intersect = " << dist_Intersect << "\n");
        UG_LOG("dist_vrtPosOut = " << dist_vrtPosOut << "\n");
        
        // if both distances of intersection points (Intersect1, Intersect2) are bigger than the distance of vrtPosOut,
        // then vrtPos is inside circle:
        if ( dist_Intersect > dist_vrtPosOut)
            UG_THROW("get_Line_Ellipse_Intersection: should not be the case!\n");
        
        for (size_t d = 0; d < dim; ++d)
            alphaOut.push_back((Intersect[d] - vrtPosOut_move[d])/lineDir[d]);

        
        // re-map Intersect:
        for (size_t d = 0; d < dim; ++d)
            Intersect[d] += center[d];
        
        
      // check:
        if ( fabs(alphaOut[0] - alphaOut[1]) > 0.0000001)
        {
            UG_LOG("alphaOut[0] = " << alphaOut[0] << "\n");
            UG_LOG("alphaOut[1] = " << alphaOut[1] << "\n");
        
        }
        
        number dist_vrtPosIn = VecDistance(vrtPosIn_move, center_project);
        number dist_vrtPos_move = VecDistance(vrtPosIn_move, vrtPosOut_move);
        number dist_vrtPos = VecDistance(vrtPosIn, vrtPosOut);
        number lineDir_length = sqrt(VecDot(lineDir, lineDir));
 
        UG_LOG("dist_vrtPosIn = " << dist_vrtPosIn << "\n");

        
        number alpha = alphaOut[0];
        
        alphaOut.clear();
        alphaOut.push_back(alpha);
        alphaOut.push_back(1.0 - alpha);
        
        UG_LOG("alphaOut 0 = " << alphaOut[0] << "\n");
        UG_LOG("alphaOut 1 = " << alphaOut[1] << "\n");
        
        return true;
        
    }
    */
    // returns a number which substitues the computation of an LSvalue, i.e.:
    // return value > 0, if vrtPosOut on same side as vrtPosIn (= discriminant < 0)
    // return value < 0, if vrtPosOut on other side as vrtPosIn (= discriminant > 0)
    number IFF_LineEllipse_Intersection(MathVector<dim>& Intersect, const MathVector<dim>& vrtPosOut,
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
            
        if ( fabs(dist_Intersect - dist_vrtPosOut) < 1e-6)
            return 0.0;

        // if both distances of intersection points (Intersect1, Intersect2) are bigger than the distance of vrtPosOut,
        // then vrtPos is inside circle:
        if ( dist_Intersect > dist_vrtPosOut)
        {
            return 0.001;
        }

        return -0.001;
        
    }
    
    number IFF_LineCircle_Intersection(MathVector<dim>& Intersect, const MathVector<dim>& vrtPosOut,
                                     const MathVector<dim>& vrtPosIn, int PrtIndex, std::vector<number>& alphaOut)
    {
        if ( (int)this->num_particles() > 1 )
            UG_THROW("get_LineCircle_Intersection(): number of given particles = " << this->num_particles() << " more than 1! ... not supposed to be!\n");
        
        //if ( dim == 3 ) UG_THROW("in 'get_LineCircle_Intersection()': not implemented for 3d case!\n");
        
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
        //			alpha^2 * <lineDir,lineDir> + alpha * 2*<lineDir, rayDir> + ( <rayDir, rayDir>-radius^2 ) = 0
        //
        // 			-> a =  <lineDir,lineDir>, b =  <lineDir,linrayDireDir>, c =  <rayDir,rayDir> - radius^2
        //
        //
        // 	=> from (1): Intersect = vrtPosOut + alpha*(vrtPosIn - vrtPosOut)
        ///////////////////////////////////////////////////////////////////////////////////////
        
        number alpha;
        
        const number radius = 0.125;
        
        UG_THROW("in ParticleProviderEllipse::IFF_LineCircle_Intersection: be careful! The radius here is HARD CODED, since an ellipse has no radius.\n");
        const MathVector<dim>& center = this->get_center(PrtIndex);
        
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
        {
            UG_LOG("IFF_Line_Ellipse_Intersection: Value of discriminant = " << discriminant << "\n");
            return 0.001;
        }
        // discriminant = 0!
        const number alpha1 = (-b - sqrt(discriminant)) / (2.0 * a);
        const number alpha2 = (-b + sqrt(discriminant)) / (2.0 * a);
        
        if (alpha1 <= alpha2)
            alpha = alpha1;
        else
            alpha = alpha2;
        
        // if alpha < 0, both vrt are inside:
        if (alpha < 0 || (alpha - 1.0) > 1e-8)
        {
            UG_LOG("IFF_Line_Ellipse_Intersection: alpha < 0: " << alpha << "\n");
            return 0.001;
        }
        
        for (size_t d = 0; d < dim; ++d)
            Intersect[d] = vrtPosOut[d] + alpha * lineDir[d];
        
        alphaOut.clear();
        alphaOut.push_back(alpha);
        alphaOut.push_back(1.0 - alpha);
        
   //     UG_LOG("alphaOut 0 = " << alphaOut[0] << "\n");
   //     UG_LOG("alphaOut 1 = " << alphaOut[1] << "\n");
        
        return -0.001;
    }
    
///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
 

    
    void print()
    {
        UG_LOG(" +++++++++++++++++++++++++ ParticleProviderEllipse Info ++++++++++++++++++++++++++++++ \n");
        UG_LOG("+++ num_particles = " << this->num_particles() << "\n");
        for (size_t p = 0; p < this->num_particles(); ++p)
        {
            UG_LOG("+++ center = " << this->get_center(p) << "\n");
            UG_LOG("+++ density = " << this->get_density(p) << "\n");
            UG_LOG("+++ theta = " << get_theta(p) << "\n");
            UG_LOG("+++ elliptic_axis_x = " << get_elliptic_axis_x(p) << "\n");
            UG_LOG("+++ elliptic_axis_y = " << get_elliptic_axis_y(p) << "\n");
            
            if ( this->get_DoF_modus_linear(p) )
            {UG_LOG("+++ linear velocity set to: " << this->get_linear_velocity(p, 0) << "\n\n");}
            else
            {UG_LOG("+++ linear velocity FREE and set to: " << this->get_linear_velocity(p, 0) << "\n\n");}
            if ( this->get_DoF_modus_angular(p) )
            {UG_LOG("+++ angular velocity set to: " << this->get_angular_velocity(p, 0) << "\n\n");}
            else
            {UG_LOG("+++ lineangularar velocity FREE and set to: " << this->get_angular_velocity(p, 0) << "\n\n");}
            
        }
        UG_LOG(" ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ \n");
    }
        
protected:
    std::vector<MathVector<dim> > m_vEllipticAxis;
    std::vector<number> m_vTheta;
        
};
    

} // end namespace ug


#include "interface_provider_particle_sphere_impl.h"
#include "interface_provider_particle_ellipse_impl.h"


#endif /* INTERFACE_PROVIDER_PARTICLE_H_ */
