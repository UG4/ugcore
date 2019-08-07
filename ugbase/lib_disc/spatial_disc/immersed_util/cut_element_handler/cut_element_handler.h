/*
 * interface_handler_local.h
 *
 *  Created on: 15.01.2015
 *      Author: suze
 */

#ifndef CUT_ELEMENT_HANDLER_H_
#define CUT_ELEMENT_HANDLER_H_

#include "lib_grid/multi_grid.h"
#include "lib_disc/function_spaces/grid_function.h"
#include "../interface_provider/interface_provider_diffusion.h"
#include "../interface_provider/interface_provider_particle.h"

#ifdef UG_PARALLEL
    #include "pcl/pcl_communication_structs.h"
    #include "pcl/pcl_process_communicator.h"
#endif

namespace ug{

enum ElementModus
{
	INSIDE_DOM = 0,
	OUTSIDE_DOM,
	CUT_BY_INTERFACE,
	CUT_BY_2_INTERFACE
};

enum VertexModus
{
	INSIDE = 0,
	OUTSIDE,
	ON_INTERFACE,
};


template <int TWorldDim>
class ICutElementHandler
{
    ///	world Dimension
    static const int dim = TWorldDim;
    
	public:

/// default constructor:
	ICutElementHandler(){};

/// destructor
	~ICutElementHandler() {}

	ElementModus get_element_modus(GridObject* elem)
	{ UG_THROW("in 'ICutElementHandler::get_element_modus()' : needs to be implemented by derived class!\n");}

	ElementModus compute_element_modus(GridObject* elem, const int interfaceOrientation)
	{ UG_THROW("in 'ICutElementHandler::compute_element_modus()': needs to be implemented by derived class!\n");}

    bool is_FTVertex(Vertex* vrt, size_t vrtIndex)
	{ UG_THROW("in 'ICutElementHandler::is_FTVertex()': needs to be implemented by derived class!\n");}

    bool is_OutsideVertex(Vertex* vrt, size_t vrtIndex)
	{ UG_THROW("in 'ICutElementHandler::is_OutsideVertex()': needs to be implemented by derived class!\n");}

    bool is_nearInterfaceVertex(Vertex* vrt, size_t vrtIndex)
	{ UG_THROW("in 'ICutElementHandler::is_nearInterfaceVertex()': needs to be implemented by derived class!\n");}

    virtual number get_LSvalue_byPosition(MathVector<dim> vrtPos)
    { UG_THROW("in 'ICutElementHandler::get_LSvalue_byPosition()': needs to be implemented by derived class!\n");}

    virtual number get_LSvalue(Vertex* vrt, const int prtIndex)
    { UG_THROW("in 'ICutElementHandler::get_LSvalue()': needs to be implemented by derived class!\n");}

};

        
template <int TWorldDim>
class CutElementHandlerFlatTop : public ICutElementHandler<TWorldDim>
{
        
public:
    ///	world Dimension
    static const int dim = TWorldDim;
        
    ///	Type of position coordinates
    typedef MathVector<dim> position_type;
    
    ///	Type of Position Attachment
    typedef Attachment<position_type> position_attachment_type;
    
    struct ParticleData {
        MathVector<dim> center;
        MathVector<dim> transVel;
        MathVector<dim> rotVel;
        bool valid;
    };
    
    //	struct CombinedParticleData {
    //		std::vector<ParticleData> data;
    //	};
    typedef std::vector<ParticleData> CombinedParticleData;
    
    
    ///	Type of Accessor to the Position Data Attachment
    typedef Grid::VertexAttachmentAccessor<position_attachment_type> position_accessor_type;
    
    typedef typename domain_traits<dim>::grid_base_object grid_base_object;
    
    /// list of element needed for looping cut OR all elements of a particle
    typedef std::vector<std::vector<grid_base_object*> > vvElemList;
    
    /// default constructor:
    CutElementHandlerFlatTop(SmartPtr<MultiGrid> mg, const char* fctNames,
                                SmartPtr<ParticleProviderSphere<dim> > interfaceProvider);
    CutElementHandlerFlatTop(SmartPtr<MultiGrid> mg, const char* fctNames,
                             SmartPtr<ParticleProviderEllipse<dim> > interfaceProvider);
    
    /// destructor
    ~CutElementHandlerFlatTop() {}
        
    /// calls 'copy_solution(topLevel) - update_for_multigrid(baseLev-topLev) - update_solution(topLevel)'
    template <typename TDomain>
    void init(ConstSmartPtr<DoFDistribution> dd, const int baseLevel, const int topLevel);
        
    void update_global_indices(ConstSmartPtr<DoFDistribution> dd, const int levIndex);
    void update_marker(ConstSmartPtr<DoFDistribution> dd, const int levIndex);
    void print_elem_lists(ConstSmartPtr<DoFDistribution> dd);
    void update_interface(const int topLevel, number deltaT);
    /// returns the levIndex of the corresponding gridLevel within 'm_Map'
    void update_multigrid_data(ConstSmartPtr<DoFDistribution> dd, const int levIndex);
    
    /// computing the element modus ONCE and marks according BoolMarker (called by 'InterfaceHandlerLocalParticle::update_marker()')
    ElementModus compute_element_modus(int prtIndex, GridObject* elem, const int interfaceOrientation);
    ElementModus compute_element_modus(int prtIndex, GridObject* elem)
    { return compute_element_modus(prtIndex, elem, 1); }
    ElementModus compute_element_modus(GridObject* elem, const int interfaceOrientation);
    ElementModus compute_element_modus(GridObject* elem)
    { return compute_element_modus(elem, 1); }
        
    /// derives the 'ElementModus' from the boolian of the BoolMarker
    ElementModus get_element_modus(GridObject* elem);
        
    /// checks if grid data is updated and returns 'levIndex'-pair for 'gridLevel' in 'm_Map'
    int get_Index(const GridLevel& gridLevel, ConstSmartPtr<DoFDistribution> dd);
    /// only call, when you are sure, that the data is available; if NOT: THROW error!
    int get_Index(const GridLevel& gridLevel);
    int get_Index_old(const GridLevel& gridLevel);
    
    /// resets all boolmarker of this class; called during update_marker(); necessary for time-dependent case
    void clear_bool_marker();
    
    /// update member 'm_spFlatTopVrtMarker' ...
    bool element_is_pyramid(grid_base_object* elem);

    bool is_outsideFluid_prtIndex(const int prtIndex, Vertex* vrt);

    bool is_outsideFluid(Vertex* vrt, const int interfaceOrientation);
    bool is_outsideFluid(int& PrtIndex, Vertex* vrt);
    bool is_outsideFluid(Vertex* vrt);
    bool is_outsideFluid_inverse(Vertex* vrt);
 
    
    /// returns true, if 'vrt' lies inside the 'prtIndex'-th particle OR:
    /// if the vrt lies outside the particle BUT near the interface,
    /// 	i.e.'m_pNearInterfaceVrtMarker->_is_marked(vrt) = true" && 'm_spFlatTopVrtMarker->is_marked(vrt) == true'
    bool set_nearInterface(Vertex* vrt, const number threshold);
    bool set_nearInterface(Vertex* vrt);

    /// returns true, if vertex lies INSIDE fluid, BUT near to interface => no DoF
    bool is_nearInterfaceVertex(Vertex* vrt, size_t vrtIndex)
    { return m_spNearInterfaceVrtMarker->is_marked(vrt); }
    
    /// returns true, if vertex lies OUTSIDE fluid, BUT near to interface => FT-vertex
    bool is_FTVertex(Vertex* vrt);
    bool is_FTVertex(Vertex* vrt, size_t vrtIndex);
    bool is_FTVertex( int& prtIndex, Vertex* vrt);
    
    /// returns true, if vertex lies OUTSIDE fluid => no DoF
    bool is_OutsideVertex(Vertex* vrt, size_t vrtIndex)
    { return m_spOutsideMarker->is_marked(vrt); }
    
    
    void set_threshold(size_t level, const number threshold)
    { m_vThresholdOnLevel[level] = threshold; }
    number get_threshold(Vertex* vrt)
    {
        const size_t lev = m_spMG->get_level(vrt);
        return m_vThresholdOnLevel[lev];
    }
    
    size_t num_particles() const			{ return m_spInterfaceProvider->num_particles();}
    number get_density(int prtIndex)		{ return m_spInterfaceProvider->get_density(prtIndex);}
    MathVector<dim> get_center(int prtIndex){ return m_spInterfaceProvider->get_center(prtIndex);}
    void set_center(MathVector<dim> center, int prtIndex){ m_spInterfaceProvider->set_center(center, prtIndex);}
    
    // careful! --> only dummy implementation, sice some methods in MovingParticle-class use it; in these cases it will be overwritten by
    // the ParticleProviderSphere-class (hopefully) ;)
    number get_radius(int prtIndex)         { return m_spInterfaceProvider->get_radius(prtIndex);}
    
    /// needed eg. for 'ParticleMapper::set_identity_mat_constraint()'
    bool is_extraDoF(DoFIndex dofIndex, int levIndex);
        
    MathMatrix<TWorldDim,TWorldDim> get_rotationMat(MathVector<TWorldDim> radialVector);
    
    bool get_DoF_modus_linear(int prtIndex) { return m_spInterfaceProvider->get_DoF_modus_linear(prtIndex);}
    bool get_DoF_modus_angular(int prtIndex){ return m_spInterfaceProvider->get_DoF_modus_angular(prtIndex);}
        
    // called during PartcleMapper::modify_GlobalSol():
    void set_extraSolTrans(number solution, size_t prtIndex, int timeSeriesInd, int cmp)
    { m_spInterfaceProvider->set_linear_velocity(solution, prtIndex, timeSeriesInd, cmp); }
    void set_extraSolTrans(MathVector<dim> solution, size_t prtIndex, int timeSeriesInd)
    { m_spInterfaceProvider->set_linear_velocity(solution, prtIndex, timeSeriesInd); }
    void set_extraSolRot(number solution, size_t prtIndex, int timeSeriesInd, int cmp)
    { m_spInterfaceProvider->set_angular_velocity(solution, prtIndex, timeSeriesInd, cmp); }
    void set_extraSolRot(MathVector<dim> solution, size_t prtIndex, int timeSeriesInd)
    { m_spInterfaceProvider->set_angular_velocity(solution, prtIndex, timeSeriesInd); }
    
    /// Access methods
    
    std::vector<DoFIndex> get_transInd(int levIndex, int prtIndex)
    {
        std::vector<DoFIndex> vTransInd;
        for ( size_t cmp = 0; cmp < dim; ++cmp )
            vTransInd.push_back(m_vvvGlobalIndices_linearVel[levIndex][prtIndex][cmp]);
        return vTransInd;
    }
    std::vector<DoFIndex> get_rotInd(int levIndex, int prtIndex)
    {
        std::vector<DoFIndex> vRotInd;
        for ( size_t cmp = 0; cmp < dim; ++cmp )
            vRotInd.push_back(m_vvvGlobalIndices_angularVel[levIndex][prtIndex][cmp]);
        return vRotInd;
    }
    DoFIndex get_transInd_Comp(int levIndex, int prtIndex, size_t cmp)
    { 	if ( cmp == dim ) UG_THROW("no acces to pressure index in transInd! EXIT...\n");
        return m_vvvGlobalIndices_linearVel[levIndex][prtIndex][cmp]; }
    DoFIndex get_rotInd_Comp(int levIndex, int prtIndex, size_t cmp)
    { 	if ( cmp == dim ) UG_THROW("no acces to pressure index in rotInd! EXIT...\n");
        return m_vvvGlobalIndices_angularVel[levIndex][prtIndex][cmp]; }
        
    
        
    void get_global_indices(std::vector<DoFIndex>& transInd, std::vector<DoFIndex>& rotInd,
                            const int levIndex, const int prtIndex)
    {
        transInd = get_transInd(levIndex, prtIndex);
        rotInd = get_rotInd(levIndex, prtIndex);
    }
    
    /// get solution values
    MathVector<dim> get_transSol(size_t prtIndex, size_t timeSeriesInd)
    { return m_spInterfaceProvider->get_linear_velocity(prtIndex, timeSeriesInd); }
    MathVector<dim> get_rotSol(size_t prtIndex, size_t timeSeriesInd)
    { return m_spInterfaceProvider->get_angular_velocity(prtIndex, timeSeriesInd); }
        
    // called during PartcleMapper::modify_GlobalSol():
    void set_solution(number solution, size_t prtIndex, int timeSeriesInd, int cmp)
    { m_spInterfaceProvider->set_solution(solution, prtIndex, timeSeriesInd); }
    
    /// get solution values
    number get_solution(size_t prtIndex, size_t timeSeriesInd)
    { return m_spInterfaceProvider->get_solution(prtIndex, timeSeriesInd); }
    
    int get_prtIndex() { UG_THROW("attention: member 'm_prtIndex' is only contained in 'InterfaceHandlerLocalParticle', NOT in 'CutElementHandler'!\n"); }
    int get_prtIndex(size_t dof)
    {
        UG_THROW("attention: You want to get the particle index of an array. Make sure that the data was set!\n");
        return m_vPrtIndex[dof];
    }
    
    //	bool valid_prt_information(int prtIndex) {
    //		return particleData[prtIndex].valid;
    //	}
    
    void set_mpi_routine(int val) {
        active_mpi_routine = val;
    }
    
#ifdef UG_PARALLEL
    void synchronize_particles(int levIndex) {
        
        bool verbose = false;
        
#ifdef UG_DEBUG
        verbose = true;
#endif
        
        pcl::ProcessCommunicator com;
        
        if (active_mpi_routine == 1) {
            ///////////////////////////////////////
            //	Synchronisation using Allreduce	 //
            ///////////////////////////////////////
            for (size_t p = 0; p < num_particles(); ++p)
            {
                if (verbose) {
                    UG_LOG("Synchronize particle " << p << ".\n");
                    UG_LOG("m_vvvElemListCut.size(): " << m_vvvElemListCut.size() << ".\n");
                    UG_LOG("m_vvvElemListCut[levIndex].size(): " << m_vvvElemListCut[levIndex].size() << ".\n");
                    UG_LOG("m_vvvElemListCut[levIndex][" << p << "].size(): " << m_vvvElemListCut[levIndex][p].size() << ".\n");
                }
                std::vector<grid_base_object*> ElemList = m_vvvElemListCut[levIndex][p];
                if (verbose)
                    UG_LOG("Process " << pcl::ProcRank() << ": 1_ synchronize_particles() ElemList.size(): " << ElemList.size() << "\n");
                if ( ElemList.size() == 0 ) {
                    // send zero vector as center to all other processes
                    std::vector<double> zero_values(3*dim,0.0);
                    std::vector<double> values(3*dim,0.0);
                    com.allreduce(&zero_values[0], &values[0], 3*dim, MPI_DOUBLE, PCL_RO_SUM);
                    MathVector<dim> center;
                    MathVector<dim> linearVelocity;
                    MathVector<dim> angularVelocity;
                    for (size_t i = 0; i < dim; ++i) {
                        center[i] = values[i];
                        linearVelocity[i] = values[dim+i];
                        angularVelocity[i] = values[2*dim+i];
                    }
                    if (verbose) {
                        UG_LOG("Process " << pcl::ProcRank() << ": setting center to " << center << "\n");
                        UG_LOG("Process " << pcl::ProcRank() << ": setting linear velocity to " << linearVelocity << "\n");
                        UG_LOG("Process " << pcl::ProcRank() << ": setting angular velocity to " << angularVelocity << ")\n");
                    }
                    m_spInterfaceProvider->set_center(center, p);
                    m_spInterfaceProvider->set_linear_velocity(linearVelocity, p, 0);
                    m_spInterfaceProvider->set_angular_velocity(angularVelocity, p, 0);
                } else {
                    // send correct center value to all other processes. Since all other processes should send zeros the allreduce with PCL_RO_SUM results in the correct center values on all processes.
                    MathVector<dim> center = m_spInterfaceProvider->get_center(p);
                    MathVector<dim> linearVelocity = m_spInterfaceProvider->get_linear_velocity(p, 0);
                    MathVector<dim> angularVelocity = m_spInterfaceProvider->get_angular_velocity(p, 0);
                    std::vector<double> values(3*dim,0.0);
                    for (size_t i = 0; i < dim; ++i) {
                        values[i] = center[i];
                        values[dim+i] = linearVelocity[i];
                        values[2*dim+i] = angularVelocity[i];
                    }
                    MathVector<3*dim> new_values;
                    com.allreduce(&values[0], &new_values[0], 3*dim, MPI_DOUBLE, PCL_RO_SUM);
                    MathVector<dim> new_center;
                    MathVector<dim> new_linearVelocity;
                    MathVector<dim> new_angularVelocity;
                    
                    for (size_t i = 0; i < dim; ++i) {
                        new_center[i] = new_values[i];
                        new_linearVelocity[i] = new_values[dim+i];
                        new_angularVelocity[i] = new_values[2*dim+i];
                    }
                    if (verbose)
                        UG_LOG("Process " << pcl::ProcRank() << ": setting center to " << new_center << "\n");
                    m_spInterfaceProvider->set_center(new_center, p);
                    if (verbose)
                        UG_LOG("Process " << pcl::ProcRank() << ": setting linear velocity to " << new_linearVelocity << "\n");
                    m_spInterfaceProvider->set_linear_velocity(new_linearVelocity, p, 0);
                    if (verbose)
                        UG_LOG("Process " << pcl::ProcRank() << ": setting angular velocity to " << new_angularVelocity << "\n");
                    m_spInterfaceProvider->set_angular_velocity(new_angularVelocity, p, 0);
                }
            }
        } else {
            
            ////////////////////////////////////////
            //	Synchronisation using Interfaces  //
            ////////////////////////////////////////
            // Create up to date particle data struct
            if (verbose)
                UG_LOG("Create data struct.\n");
            CombinedParticleData particleValues;
            ParticleData dataSet;
            for (size_t p = 0; p < num_particles(); ++p) {
                std::vector<grid_base_object*> ElemList = m_vvvElemListCut[levIndex][p];
                if ( ElemList.size() == 0 ) {
                    MathVector<dim> center;
                    if (verbose)
                        UG_LOG(pcl::ProcRank() << " sends zeros.\n");
                    for (std::size_t i = 0; i < dim; ++i){
                        center[i] = 0.0;
                        dataSet.transVel[i] = 0.0;
                        dataSet.rotVel[i] = 0.0;
                    }
                    dataSet.center = center;
                    dataSet.valid = false;
                } else {
                    dataSet.center = m_spInterfaceProvider->get_center(p);
                    dataSet.transVel = m_spInterfaceProvider->get_linear_velocity(p, 0);
                    dataSet.rotVel = m_spInterfaceProvider->get_angular_velocity(p, 0);
                    dataSet.valid = true;
                    if (verbose) {
                        UG_LOG(pcl::ProcRank() << " sends\n" );
                        UG_LOG("(" << dataSet.center[0] << ", " << dataSet.center[1] << ")\n");
                        UG_LOG("(" << dataSet.transVel[0] << ", " << dataSet.transVel[1] << ")\n");
                        UG_LOG("(" << dataSet.rotVel[0] << ", " << dataSet.rotVel[1] << ")\n");
                        UG_LOG(std::boolalpha << dataSet.valid);
                    }
                }
                particleValues.push_back(dataSet);
            }
            if (verbose)
                UG_LOG("Struct created.\n");
            com.barrier();
            
            ParticleData default_new_data;
            MathVector<dim> zero_vec;
            for (std::size_t i = 0; i < dim; ++i){
                zero_vec[i] = 0.0;
            }
            default_new_data.valid = false;
            default_new_data.center = zero_vec;
            default_new_data.transVel = zero_vec;
            default_new_data.rotVel = zero_vec;
            CombinedParticleData new_particleData(num_particles(),default_new_data);
            
            const GridLayoutMap& glm = m_spMG->distributed_grid_manager()->grid_layout_map();
            for (std::size_t proc = 0; proc < pcl::NumProcs(); ++proc) {
                if (verbose)
                    UG_LOG("proc " << proc << " sends to ");
                if (glm.has_layout<Vertex>(INT_H_MASTER))
                {
                    const typename GridLayoutMap::Types<Vertex>::Layout& vrt_hm_layout = glm.get_layout<Vertex>(INT_H_MASTER);
                    if (vrt_hm_layout.interface_exists(proc, levIndex))
                    {
                        for (GridLayoutMap::Types<Vertex>::Layout::const_iterator it = vrt_hm_layout.begin(levIndex); it != vrt_hm_layout.begin(levIndex); ++it) {
                            com.send_data(&particleValues, sizeof(particleValues), vrt_hm_layout.proc_id(it), proc);
                            //UG_LOG(vrt_hm_layout.proc_id(it) << " ");
                        }
                    }
                }
                if (glm.has_layout<Vertex>(INT_H_SLAVE))
                {
                    const typename GridLayoutMap::Types<Vertex>::Layout& vrt_hm_layout = glm.get_layout<Vertex>(INT_H_SLAVE);
                    if (vrt_hm_layout.interface_exists(proc, levIndex))
                    {
                        for (GridLayoutMap::Types<Vertex>::Layout::const_iterator it = vrt_hm_layout.begin(levIndex); it != vrt_hm_layout.begin(levIndex); ++it) {
                            com.receive_data(&new_particleData, sizeof(new_particleData), vrt_hm_layout.proc_id(it), proc);
                        }
                    }
                }
            }
            if (verbose)
                UG_LOG("\n")
                for (size_t p = 0; p < num_particles(); ++p) {
                    if (new_particleData[p].valid) {
                        m_spInterfaceProvider->set_center(new_particleData[p].center, p);
                        m_spInterfaceProvider->set_linear_velocity(new_particleData[p].transVel, p, 0);
                        m_spInterfaceProvider->set_angular_velocity(new_particleData[p].rotVel, p, 0);
                    }
                }
        }
    }
#endif
    
    /// methods called by local_to_global_mappe:
    number Volume(int levIndex, size_t prtIndex)
                        { return m_spInterfaceProvider->Volume(levIndex, prtIndex);}
    number Mass(const int levIndex, const int prtIndex, const number fluidDensity)
                        { return m_spInterfaceProvider->Mass(levIndex, prtIndex, fluidDensity);}
    number Mass(const int levIndex, const int prtIndex, const number volume, const number fluidDensity)
                        { return m_spInterfaceProvider->Mass(levIndex, prtIndex, volume, fluidDensity);}
    number MomOfInertia(const int levIndex, const int prtIndex, const number fluidDensity)
                        { return m_spInterfaceProvider->MomOfInertia(levIndex, prtIndex, fluidDensity);}
    number MomOfInertia(const int levIndex, const int prtIndex, const number volume, const number fluidDensity)
                        { return m_spInterfaceProvider->MomOfInertia(levIndex, prtIndex, volume, fluidDensity);}
    
    void update_prtCoords(const int topLevel, const number deltaT);
    
    number get_LSvalue_byPosition(MathVector<dim> vrtPos, const int prtIndex) { return m_spInterfaceProvider->get_LSvalue_byPosition(vrtPos, prtIndex);}
    number get_LSvalue_byPosition(MathVector<dim> vrtPos)
    {//ToDo
        UG_THROW("in CutElementHandler:getLSvalue_byPosition: no prtIndex-Parameter!\n");
        int prtIndex = 0;
        return get_LSvalue_byPosition(vrtPos, prtIndex);
    }
    number get_LSvalue(Vertex* vrt, const int prtIndex) { return get_LSvalue_byPosition(m_aaPos[vrt], prtIndex); }
    
    bool get_intersection_point(MathVector<dim>& Intersect, const MathVector<dim>& vrtPosOut, const MathVector<dim>& vrtPosIn, int PrtIndex)
    { return this->m_spInterfaceProvider->get_intersection_point(Intersect, vrtPosOut, vrtPosIn, PrtIndex); }
    bool get_intersection_point(MathVector<dim>& Intersect, const MathVector<dim>& vrtPosOut, const MathVector<dim>& vrtPosIn, int PrtIndex, std::vector<number>& alphaOut)
    { return this->m_spInterfaceProvider->get_intersection_point(Intersect, vrtPosOut, vrtPosIn, PrtIndex, alphaOut); }
    
    void print_velocity(const MathVector<dim>& transSol, const MathVector<dim>& rotSol, const int prtIndex, const bool isTimedep, const number time, const char* filename)
    {  m_spInterfaceProvider->print_velocity(transSol, rotSol, prtIndex, isTimedep, time, filename); }

    SmartPtr<MultiGrid> m_spMG;
    position_attachment_type m_aPos;		///<Position Attachment
    position_accessor_type m_aaPos;			///<Accessor
        
    const char* m_fctNames;
    //ToDo:brauche ich das??
    std::vector<int> m_vPrtIndex;
    
    /// global data updatet during 'update_for_time_stepping()'
    /// [levIndex][prtIndex][cmp]
    std::vector< std::vector< std::vector<DoFIndex> > > m_vvvGlobalIndices_linearVel;
    std::vector< std::vector< std::vector<DoFIndex> > > m_vvvGlobalIndices_angularVel;
        
    SmartPtr<BoolMarker> m_spFlatTopVrtMarker;  		// marks: vertices for 'is_outside()'
    SmartPtr<BoolMarker> m_spNearInterfaceVrtMarker;  	// marks: vertices for 'is_nearInterface()', i.e. potentially inside domain
    // necessary to distinguish between (FlatTop && outside) and (FlatTop && nearInterface)
    // Remark: outside AND nearInterface are FlatTopVrt, BUT nearInterface can lie inside the domain
    //			=> e.g during 'CollectCorners_FlatTop_2d()' call of 'get_intersection_point()'
    //				would fail if these cases are NOT distinguished!
        
    SmartPtr<BoolMarker> m_spCutMarker;  	// marks: elements, edges
    SmartPtr<BoolMarker> m_spOutsideMarker;	// marks: elements, edges, vertices (vrt, for which 'is_outsideFluid == true, i.e. outside OR nearInterface
        
    /// [levIndex][prtIndex][counter]
    std::vector<vvElemList> m_vvvElemList;
    std::vector<vvElemList> m_vvvElemListCut;
    std::vector<vvElemList> m_vvvElemListOutside;
        
    /// associotion of 'GridLevel' with an 'Index' value for access to m_vvv-lists:
    std::map<GridLevel, size_t> m_Map;
        
    /// contains radius, center and density of all given particles
    SmartPtr<ParticleProvider<dim> > m_spInterfaceProvider;
        
    std::vector<number> m_vThresholdOnLevel;
  
    //    CombinedParticleData particleData;
    std::size_t active_mpi_routine;

};
    


template <int TWorldDim>
class CutElementHandlerImmersed : public ICutElementHandler<TWorldDim>
{

	public:
	///	world Dimension
	static const int dim = TWorldDim;

    ///	Type of position coordinates
    typedef MathVector<dim> position_type;

    ///	Type of Position Attachment
    typedef Attachment<position_type> position_attachment_type;

    ///	Type of Accessor to the Position Data Attachment
    typedef Grid::VertexAttachmentAccessor<position_attachment_type> position_accessor_type;

	typedef typename domain_traits<dim>::grid_base_object grid_base_object;

/// list of element needed for looping cut OR all elements of a particle
    typedef std::vector<std::vector<grid_base_object*> > vvElemList;

/// default constructor:
    CutElementHandlerImmersed(SmartPtr<MultiGrid> mg, const char* fctNames,
			   SmartPtr<DiffusionInterfaceProvider<dim> > interfaceProvider);

/// destructor
	~CutElementHandlerImmersed() {}

 /// computing the element modus ONCE and marks according BoolMarker (called by 'InterfaceHandlerLocalParticle::update_marker()')
    ElementModus compute_element_modus(GridObject* elem)
    { return compute_element_modus(elem, 1); }

    ElementModus compute_element_modus(GridObject* elem, const int interfaceOrientation);

/// derives the 'ElementModus' from the boolian of the BoolMarker
    ElementModus get_element_modus(GridObject* elem)
    { return m_elementModus; }
    VertexModus get_vertex_modus(Vertex* vrt, size_t vrtIndex)
    { return m_vVertexMode[vrtIndex];}
    VertexModus get_vertex_modus(size_t vrtIndex)
    { return m_vVertexMode[vrtIndex];}
    VertexModus get_vertex_modus(size_t vrtIndex, const int interfaceOrientation)
    { return m_vvVertexMode[interfaceOrientation][vrtIndex];}

/// only call, when you are sure, that the data is available; if NOT: THROW error!
     int get_Index(const GridLevel& gridLevel, ConstSmartPtr<DoFDistribution> dd);
     int get_Index(const GridLevel& gridLevel);
     int get_Index_old(const GridLevel& gridLevel);

     bool is_outsideDomain(Vertex* vrt, const int interfaceOrientation)
     {
    	if (interfaceOrientation == 1 )
    		return is_outsideDomain(vrt);
    	else if (interfaceOrientation == -1 )
    		return is_outsideDomain_inverse(vrt);
    	else UG_THROW("in 'CutElementHandler::is_outsideDomain()': no valid orientation given!\n");
     }
     bool is_outsideDomain(Vertex* vrt);
     bool is_outsideDomain_inverse(Vertex* vrt);

/// returns true, if 'vrt' lies inside the 'prtIndex'-th particle OR:
/// if the vrt lies outside the particle BUT near the interface,
/// 	i.e.'m_pNearInterfaceVrtMarker->_is_marked(vrt) = true" && 'm_spFlatTopVrtMarker->is_marked(vrt) == true'
     bool is_nearInterface(Vertex* vrt, const number threshold);
     bool is_nearInterface(Vertex* vrt)
     {
         const number threshold = get_threshold(vrt);
         return is_nearInterface(vrt, threshold);
     }

     //ToDO: not needed?
     VertexModus compute_vertex_modus(Vertex* vrt, const int interfaceOrientation);

    size_t get_or_insert_vertex_near(const MathVector<dim>& vrtPos);
    const size_t get_numNearVerticesPos() const { return m_verticesNearPos.size(); }

    bool is_element_near_interface() { return m_bElementNearInterface;}

    /// returns true, if vertex lies INSIDE fluid, BUT near to interface => no DoF
    bool check_vertex_modus(VertexModus vrtModus, size_t vrtIndex, const int interfaceOrientation)
    {
        size_t localInd = -0.5*interfaceOrientation + 0.5;
        if ( get_vertex_modus(vrtIndex, localInd) == vrtModus )
            return true;
        return false;
    }
    bool check_vertex_modus(VertexModus vrtModus, size_t vrtIndex)
    {
        if ( get_vertex_modus(vrtIndex) == vrtModus )
            return true;
        return false;
    }
    /// returns true, if vertex lies INSIDE fluid, BUT near to interface => no DoF
    bool check_vertex_modus(VertexModus vrtModus, Vertex* vrt, size_t vrtIndex)
    {
        if ( get_vertex_modus(vrt, vrtIndex) == vrtModus )
            return true;
        return false;
    }

    /// returns true, if vertex lies OUTSIDE fluid, BUT near to interface => FT-vertex
     bool is_FTVertex(Vertex* vrt, size_t vrtIndex)
     { return !check_vertex_modus(INSIDE, vrt, vrtIndex); }

    /// returns true, if vertex lies OUTSIDE fluid => no DoF
     bool is_OutsideVertex(Vertex* vrt, size_t vrtIndex)
     { return check_vertex_modus(OUTSIDE, vrt, vrtIndex); }

     bool is_nearInterfaceVertex(Vertex* vrt, size_t vrtIndex)
     { return check_vertex_modus(ON_INTERFACE, vrt, vrtIndex); }

     void set_threshold(size_t level, const number threshold)
     { m_vThresholdOnLevel[level] = threshold; }
     number get_threshold(Vertex* vrt)
     {
        const size_t lev = m_spMG->get_level(vrt);
        return m_vThresholdOnLevel[lev];
     }

     size_t num_particles() const			 { return m_spInterfaceProvider->num_particles();}
 //    number get_radius(int prtIndex)		 { return m_spInterfaceProvider->get_radius(prtIndex);}
     number get_theta(int prtIndex)          { return m_spInterfaceProvider->get_theta(prtIndex);}
     number get_density(int prtIndex)		 { return m_spInterfaceProvider->get_density(prtIndex);}
     MathVector<dim> get_center(int prtIndex){ return m_spInterfaceProvider->get_center(prtIndex);}

    // called during PartcleMapper::modify_GlobalSol():
    void set_solution(number solution, size_t prtIndex, int timeSeriesInd, int cmp)
    { m_spInterfaceProvider->set_solution(solution, prtIndex, timeSeriesInd); }

    /// get solution values
    number get_solution(size_t prtIndex, size_t timeSeriesInd)
    { return m_spInterfaceProvider->get_solution(prtIndex, timeSeriesInd); }

    void print_velocity(const MathVector<dim>& transSol, const MathVector<dim>& rotSol, const int prtIndex, const bool isTimedep, const number time, const char* filename)
    {  m_spInterfaceProvider->print_velocity(transSol, rotSol, prtIndex, isTimedep, time, filename); }

    SmartPtr<MultiGrid> m_spMG;
    position_attachment_type m_aPos;		///<Position Attachment
    position_accessor_type m_aaPos;			///<Accessor

    const char* m_fctNames;
    //ToDo:brauche ich das??
    std::vector<int> m_vPrtIndex;
    /// associotion of 'GridLevel' with an 'Index' value for access to m_vvv-lists:
    std::map<GridLevel, size_t> m_Map;

    /// contains radius, center and density of all given particles
    SmartPtr<DiffusionInterfaceProvider<dim> > m_spInterfaceProvider;
 
    std::vector<number> m_vThresholdOnLevel;
    std::vector<VertexModus> m_vVertexMode;
    std::vector<std::vector<VertexModus> > m_vvVertexMode;
    
    ElementModus m_elementModus;
    bool m_bElementNearInterface;
    std::map<MathVector<dim>, size_t> m_MapNearVertices;
    std::vector<MathVector<dim> > m_verticesNearPos;


};


}// end namespace ug

#include "cut_element_handler_FT_impl.h"
#include "cut_element_handler_immersed_impl.h"

#endif /* CUT_ELEMENT_HANDLER_H_ */
