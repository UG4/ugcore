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

#include "interface_provider_base.h"

template <int TWorldDim>
class ParticleProvider : public IInterfaceProvider<TWorldDim>
{

	public:
///	world Dimension
	static const int dim = TWorldDim;

/// default constructor:
	ParticleProvider()
    :     m_prtIndex(-1), m_orientationInterface(1)
	{
		clear();
		UG_LOG("ParticleProvider constructor\n");
	};

/// destructor
	virtual ~ParticleProvider() {}

//////////////////////////////////////////////////////////
/// virtual base class methods, which need to be
///  implemented by derived class
//////////////////////////////////////////////////////////
    
/// methods called by CutElementHandler:
    virtual number get_LSvalue_byPosition(MathVector<dim> vrtPos)
    { UG_THROW("in 'ParticleProvider::get_LSvalue_byPosition(vrtPos)': needs to be implemented by derived class!\n");}
    
    virtual number get_LSvalue_byPosition(MathVector<dim> vrtPos, const int prtIndex)
    { UG_THROW("in 'ParticleProvider::get_LSvalue_byPosition(vrtPos, prtIndex)': needs to be implemented by derived class!\n");}
    
    virtual bool get_intersection_point(MathVector<dim>& Intersect, const MathVector<dim>& vrtPosOut,
                                        const MathVector<dim>& vrtPosIn, const int PrtIndex)
    { UG_THROW("in 'ParticleProvider::get_intersection_point()': needs to be implemented by derived class!\n");}
    
    virtual bool get_intersection_point(MathVector<dim>& Intersect, const MathVector<dim>& vrtPosOut,
                                        const MathVector<dim>& vrtPosIn, const int PrtIndex, std::vector<number>& alphaOut)
    { UG_THROW("in 'ParticleProvider::get_intersection_point()': needs to be implemented by derived class!\n");}
    
// updates the location of the interface for the next time step
    virtual void update(number deltaT, const MathVector<dim> transSol, const MathVector<dim> rotSol, const int prtIndex)
    { UG_THROW("in 'ParticleProvider::update()': needs to be implemented by derived class!\n");}

    
    virtual number get_radius(int prtIndex)
    { UG_THROW("in 'ParticleProvider::get_radius()': needs to be implemented by derived class!\n");}
    
    number get_theta(int prtIndex)
    { UG_THROW("in 'ParticleProvider::get_theta()': needs to be implemented by derived class!\n");}
    number set_theta(number theta, int prtIndex)
    { UG_THROW("in 'ParticleProvider::get_theta()': needs to be implemented by derived class!\n");}
    
    virtual void print()
    { UG_THROW("in 'ParticleProvider::print()': needs to be implemented by derived class!\n");}
    
/// methods called by LocalToGlobalMapper:
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
    

//////////////////////////////////////////////////////////
/// setter methods
//////////////////////////////////////////////////////////
    
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
        m_vDensity.clear(); m_vCenter.clear();
        m_vvLinearVelocity.clear(); m_vvAngularVelocity.clear();
        m_vGivenVelocity.clear();
    }

//////////////////////////////////////////////////////////
/// setter methods
//////////////////////////////////////////////////////////

	size_t num_particles() const{ return m_vCenter.size(); }

	number get_density(int prtIndex){ return m_vDensity[prtIndex]; }
	MathVector<dim> get_center(int prtIndex){ return m_vCenter[prtIndex]; }
	bool get_DoF_modus_linear(int prtIndex) { return m_vGivenVelocity[prtIndex][0]; }
	bool get_DoF_modus_angular(int prtIndex) { return m_vGivenVelocity[prtIndex][1]; }

// getter and setter for orientation of the interface
    const int get_prtIndex() const{ return m_prtIndex; }
    void set_prtIndex(const int prtIndex) { m_prtIndex = prtIndex; }
    
    
// getter and setter for orientation of the interface
    void set_orientation(const int orientation) { m_orientationInterface = orientation; }
    const int get_orientation() const{ return m_orientationInterface; }

// getter and setter for the velocity of the interface
//  --> called during PartcleMapper::modify_GlobalSol():
	void set_linear_velocity(number solution, size_t prtIndex, int timeSeriesInd, int cmp)
        { m_vvLinearVelocity[timeSeriesInd][prtIndex][cmp] = solution; }
	void set_linear_velocity(MathVector<dim> solution, size_t prtIndex, int timeSeriesInd)
        { m_vvLinearVelocity[timeSeriesInd][prtIndex] = solution; }
    void set_angular_velocity(number solution, size_t prtIndex, int timeSeriesInd, int cmp)
        { m_vvAngularVelocity[timeSeriesInd][prtIndex][cmp] = solution; }
	void set_angular_velocity(MathVector<dim> solution, size_t prtIndex, int timeSeriesInd)
        { m_vvAngularVelocity[timeSeriesInd][prtIndex] = solution; }

	MathVector<dim> get_linear_velocity(size_t prtIndex, size_t timeSeriesInd)
		{ return m_vvLinearVelocity[timeSeriesInd][prtIndex]; }
	MathVector<dim> get_angular_velocity(size_t prtIndex, size_t timeSeriesInd)
		{ return m_vvAngularVelocity[timeSeriesInd][prtIndex]; }


protected:
    std::vector<number> m_vDensity;
    std::vector<MathVector<dim> > m_vCenter;
    
// 	indexing: [TimeSeriesIndes][prtIndex]:
    std::vector<std::vector<MathVector<dim> > > m_vvLinearVelocity;
    std::vector<std::vector<MathVector<dim> > > m_vvAngularVelocity;
    
//  indexing: [prtIndex][linear/angular]
    std::vector<std::vector<bool> > m_vGivenVelocity;
  
    int m_prtIndex;
    int m_orientationInterface;
  
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
    virtual ~ParticleProviderSphere() {}
    
    
// initializing methods
    void add       (number radius, const MathVector<dim>& center, number density);
    void add_moving(number radius, const MathVector<dim>& center, number density,
                    const MathVector<dim>& linearVel, const MathVector<dim>& angularVel);
    
    
///////////////////////////////////////////////////////////////////////////////////////
// virtual base class methods (necessary to be implemented!)
///////////////////////////////////////////////////////////////////////////////////////

    number get_LSvalue_byPosition(MathVector<dim> vrtPos);
    number get_LSvalue_byPosition(MathVector<dim> vrtPos, const int prtIndex);

    bool get_intersection_point(MathVector<dim>& Intersect, const MathVector<dim>& vrtPosOut,
                                const MathVector<dim>& vrtPosIn, const int prtIndex)
    { return get_LineCircle_Intersection(Intersect, vrtPosOut, vrtPosIn, prtIndex);}
    
    bool get_intersection_point(MathVector<dim>& Intersect, const MathVector<dim>& vrtPosOut,
                                const MathVector<dim>& vrtPosIn, const int prtIndex, std::vector<number>& alphaOut)
    { return get_LineCircle_Intersection(Intersect, vrtPosOut, vrtPosIn, prtIndex, alphaOut);}

    bool get_LineCircle_Intersection(MathVector<dim>& Intersect, const MathVector<dim>& vrtPosOut,
                                     const MathVector<dim>& vrtPosIn, int prtIndex);
    bool get_LineCircle_Intersection(MathVector<dim>& Intersect, const MathVector<dim>& vrtPosOut,
                                     const MathVector<dim>& vrtPosIn, int prtIndex, std::vector<number>& alphaOut);

//	updates the center of particle
    void update(number deltaT, const MathVector<dim> transSol, const MathVector<dim> rotSol, const int prtIndex);

    void set_radius(number radius, int prtIndex)
    {
        if ( (int)this->num_particles() < prtIndex )
            UG_THROW("ParticleProviderSphere::set_radius(): number of given particles = "
                     << this->num_particles() << " smaller than given prtIndex = " << prtIndex << "\n");
        m_vRadius[prtIndex] = radius;
    }
    
    number get_radius(int prtIndex){ return m_vRadius[prtIndex]; }

/// methods called by class 'LocalToGlobalMapper':
    number Volume(int levIndex, size_t prtIndex);
    number Mass(const int levIndex, const int prtIndex, const number fluidDensity);
    number Mass(const int levIndex, const int prtIndex, const number volume, const number fluidDensity);
    number MomOfInertia(const int levIndex, const int prtIndex, const number fluidDensity);
    number MomOfInertia(const int levIndex, const int prtIndex, const number volume, const number fluidDensity);
  
    void print();
    void print_velocity(const MathVector<dim>& transSol, const MathVector<dim>& rotSol, const int prtIndex,
                        const bool isTimedep, const number time, const char* filename);
 	  
///////////////////////////////////////////////////////////////////////////////////////
/// class member
///////////////////////////////////////////////////////////////////////////////////////


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
        m_vEllipticAxis.clear();
        m_vTheta.clear();
        UG_LOG("ParticleProviderEllipse constructor\n");
    };
        
/// destructor
    virtual ~ParticleProviderEllipse() {}

// initializing methods
    void add       (const MathVector<dim>& ellipticAxis, const MathVector<dim>& center, const number theta, number density);
    void add_moving(const MathVector<dim>& ellipticAxis, const MathVector<dim>& center, const number theta, number density,
                    const MathVector<dim>& linearVel, const MathVector<dim>& angularVel);
 
    
///////////////////////////////////////////////////////////////////////////////////////
// virtual base class methods (necessary to be implemented!)
///////////////////////////////////////////////////////////////////////////////////////

    number get_LSvalue_byPosition(MathVector<dim> vrtPos);
    number get_LSvalue_byPosition(MathVector<dim> vrtPos, const int prtIndex);
    
    bool get_intersection_point(MathVector<dim>& Intersect, const MathVector<dim>& vrtPosOut,
                                const MathVector<dim>& vrtPosIn, const int PrtIndex);
    bool get_intersection_point(MathVector<dim>& Intersect, const MathVector<dim>& vrtPosOut,
                                const MathVector<dim>& vrtPosIn, const int PrtIndex, std::vector<number>& alphaOut);

    bool get_LineEllipse_Intersection(MathVector<dim>& Intersect, const MathVector<dim>& vrtPosOut,
                                      const MathVector<dim>& vrtPosIn, int PrtIndex);
    // Old version of 'get_LineEllipse_Intersection':
    bool get_LineEllipse_Intersection(MathVector<dim>& Intersect, const MathVector<dim>& vrtPosOut,
                                      const MathVector<dim>& vrtPosIn, int PrtIndex, std::vector<number>& alphaOut);
    
// returns a value, which substitutes the computation of an LSvalue, i.e.:
// return value > 0, if vrtPosOut on same side as vrtPosIn (= discriminant < 0)
// return value < 0, if vrtPosOut on other side as vrtPosIn (= discriminant > 0)
    number IFF_LineEllipse_Intersection(MathVector<dim>& Intersect, const MathVector<dim>& vrtPosOut,
                                        const MathVector<dim>& vrtPosIn, int PrtIndex, std::vector<number>& alphaOut);

/// methods called by class 'DiffusionInterfaceMapper':
    number Volume(int levIndex, size_t prtIndex);
    number Mass(const int levIndex, const int prtIndex, const number fluidDensity);
    number Mass(const int levIndex, const int prtIndex, const number volume, const number fluidDensity);
    number MomOfInertia(const int levIndex, const int prtIndex, const number fluidDensity);
    number MomOfInertia(const int levIndex, const int prtIndex, const number volume, const number fluidDensity);
    
    
    void print();
    void print_velocity(const MathVector<dim>& transSol, const MathVector<dim>& rotSol, const int prtIndex,
                        const bool isTimedep, const number time, const char* filename);
    
    void update(number deltaT, const MathVector<dim> transSol, const MathVector<dim> rotSol, const int prtIndex);
    
    void set_theta(number theta, int prtIndex)
    {
        if ( (int)this->num_particles() < prtIndex )
            UG_THROW("EllipseProvider::set_theta(): number of given particles = "
                         << this->num_particles() << " smaller than given prtIndex = " << prtIndex << "\n");
        m_vTheta[prtIndex] = theta;
    }

    number get_theta(int prtIndex){ return m_vTheta[prtIndex]; }

///////////////////////////////////////////////////////////////////////////////////////
// new getter methods
///////////////////////////////////////////////////////////////////////////////////////

    MathVector<dim> get_elliptic_axis(int prtIndex){ return m_vEllipticAxis[prtIndex]; }
    number get_elliptic_axis_x(int prtIndex){ return m_vEllipticAxis[prtIndex][0]; }
    number get_elliptic_axis_y(int prtIndex){ return m_vEllipticAxis[prtIndex][1]; }
    
    void set_elliptic_axis(MathVector<dim> ellitpicAxis, int prtIndex)
    {
        if ( (int)this->num_particles() < prtIndex )
            UG_THROW("EllipseProvider::set_center(): number of given particles = "
                     << this->num_particles() << " smaller than given prtIndex = " << prtIndex << "\n");
        for ( size_t d = 0; d < dim; ++d )
            m_vEllipticAxis[prtIndex][d] = ellitpicAxis[d];
    }
    
  
///////////////////////////////////////////////////////////////////////////////////////
// new methods
///////////////////////////////////////////////////////////////////////////////////////
 
// set vrtPos to origin and rotate by -theta
    void rotate_vector(MathVector<dim>& vrtPos, const int prtIndex);
    void rotate_vector_inverse(MathVector<dim>& vrtPos, const int prtIndex);


///////////////////////////////////////////////////////////////////////////////////////
// class member
///////////////////////////////////////////////////////////////////////////////////////
    
protected:
    std::vector<MathVector<dim> > m_vEllipticAxis;
    std::vector<number> m_vTheta;
        
};
    

} // end namespace ug


#include "interface_provider_particle_sphere_impl.h"
#include "interface_provider_particle_ellipse_impl.h"


#endif /* INTERFACE_PROVIDER_PARTICLE_H_ */
