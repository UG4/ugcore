
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
#include "../interface_provider/interface_provider_base.h"
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

// computes the modus of the element as cut or non-cut (= inside or outside of the domain)
	ElementModus compute_element_modus(GridObject* elem)
	{ UG_THROW("in 'ICutElementHandler::compute_element_modus()': needs to be implemented by derived class!\n");}

// returns the computed element modus, which is either stored as element bool marker, or locally
// as property of the element
    ElementModus get_element_modus(GridObject* elem)
    { UG_THROW("in 'ICutElementHandler::get_element_modus()' : needs to be implemented by derived class!\n");}

    
// methods used to get the 'VertexModus' of a vertex
// REMARK: 2 arguments:
//          Vertex* vrt:     for global acces via BoolMarker
//          size_t vrtIndex: for local access via m_vvVertexModus

    bool is_onInterfaceVertex(Vertex* vrt, size_t vrtIndex)
	{ UG_THROW("in 'ICutElementHandler::is_onInterfaceVertex()': needs to be implemented by derived class!\n");}

    bool is_OutsideVertex(Vertex* vrt, size_t vrtIndex)
	{ UG_THROW("in 'ICutElementHandler::is_OutsideVertex()': needs to be implemented by derived class!\n");}

    bool is_nearInterfaceVertex(Vertex* vrt, size_t vrtIndex)
	{ UG_THROW("in 'ICutElementHandler::is_nearInterfaceVertex()': needs to be implemented by derived class!\n");}

    
     
};

template <int TWorldDim>
class CutElementHandlerBase : public ICutElementHandler<TWorldDim>
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
        
    /// constructor:

    // REMARK: due to call of the call of the methog 'm_spInterfaceProvider->get_LSvalue_byPosition()',
    //         the interfaceProvider needs to be specified already by the constructor!
    //         --> the call of 'get_LSvalue_byPosition()' will be forwarded by 'InterfaceProviderBase'
    //             to the derived class!
        CutElementHandlerBase(SmartPtr<MultiGrid> mg, SmartPtr<DiffusionInterfaceProvider<dim> > interfaceProvider);
        CutElementHandlerBase(SmartPtr<MultiGrid> mg, SmartPtr<ParticleProviderSphere<dim> > interfaceProvider);
        CutElementHandlerBase(SmartPtr<MultiGrid> mg, SmartPtr<ParticleProviderEllipse<dim> > interfaceProvider);
    
        
    /// destructor
        virtual ~CutElementHandlerBase() {}
    
    //////////////////////////////////////////////////////////
    /// init for time dependent update and initialisation:
    //////////////////////////////////////////////////////////
        
    /// calls 'copy_solution(topLevel) - update_for_multigrid(baseLev-topLev) - update_solution(topLevel)'
        template <typename TDomain>
        void init(ConstSmartPtr<DoFDistribution> dd, const int baseLevel, const int topLevel);

    /// sets all BoolMarker and fills the 'm_vvElemList...'
        virtual void update_interface_data(ConstSmartPtr<DoFDistribution> dd, const int levIndex);
    /// resets all boolmarker of this class; called during update_interface_data(); necessary for time-dependent case
        void clear_bool_marker();
    
    /// initializes crucial data for all cut element computations
    ///     especially'ElementModus' and 'VertexModus' categories
        virtual void update_multigrid_data(ConstSmartPtr<DoFDistribution> dd, const int levIndex)
            { update_interface_data(dd, levIndex); }
    
    //////////////////////////////////////////////////////////
    /// virtual base class methods:
    //////////////////////////////////////////////////////////
    
    // basic methods to compute the modus of an element and a vertex:
    
    /// returns boolian for 'does 'vrt' lie near interface?'
        virtual bool is_outsideDomain(Vertex* vrt);
    /// returns boolian for 'does 'vrt' lie outside the domain?'
        virtual bool is_nearInterface(Vertex* vrt);
 
    // methods needed during call of 'is_outsideDomain(), 'is_nearInterface()'
    // computes the 'level-set' value of a point w.r.t. the interface, in order to
    //  derive its 'VertexModus', i.e. INSIDE, OUTSIDE, ON_INTERFACE
    //  --> evaluate the distance between a vrt and the interface
        number get_LSvalue(Vertex* vrt)
        { return get_LSvalue_byPosition(m_aaPos[vrt]); }
   
        virtual number get_LSvalue_byPosition(MathVector<dim> vrtPos) = 0;
  
    //////////////////////////////////////////////////////////
    /// (1A) 'compute' and (1B) 'get' the 'ElementModus':
        
        
    // (1A) computes the modus of the element as cut or non-cut
    // (= inside or outside of the domain) and marks according BoolMarker
        virtual ElementModus compute_element_modus(GridObject* elem);
        
    // (1B) returns the computed element modus:
    //      (1) local access by member 'm_elementModus'
    //      (2) global acces via 'GridObject* elem' as parameter (can be stored using BoolMarker for elements!)
    
    // (1B)(1) LOCAL access method
        ElementModus get_element_modus() { return m_elementModus; }
        
    // (1B)(2) GLOBAL acces method
        ElementModus get_element_modus(GridObject* elem); // --> ToDo: siehe FT impl!!

        
    /// checks if cut element data is updated and returns 'levIndex'-pair for 'gridLevel' in 'm_Map'
        int get_Index(const GridLevel& gridLevel, ConstSmartPtr<DoFDistribution> dd);
    /// only call, when you are sure, that the data is available; if NOT: THROW error!
        int get_Index(const GridLevel& gridLevel);
        
        
    //////////////////////////////////////////////////////////////////////////////////
    /// (2A) compute, (2B) get, (C) check the 'VertexModus':
    //////////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////////
    //  (2A) computes the modus of the vertex as INSIDE, OUTSIDE, ON_INTERFACE

        VertexModus compute_vertex_modus(Vertex* vrt);

    //////////////////////////////////////////////////////////////////////////////////
    // (2B) returns the computed vertex modus:
    //      --> LOCAL acces (only valid for the currently assembled element!) + check via 'size_t vrtIndex'
        VertexModus get_vertex_modus(size_t vrtIndex, const int localInd)
            { return m_vvVertexMode[localInd][vrtIndex];}
    
        bool check_vertex_modus(VertexModus vrtModus, size_t vrtIndex, const int interfaceOrientation)
        {
            size_t localInd = -0.5*interfaceOrientation + 0.5;
            if ( get_vertex_modus(vrtIndex, localInd) == vrtModus )
                return true;
            return false;
        }
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// (C) base class methods 'is_onInterfaceVertex()', 'is_OutsideVertex()', 'is_nearInterfaceVertex()'
    //  (C1) global version via 'Vertex* vrt'    (i.e. valid during whole assembling process)
    //  (C2) local version via 'size_t vrtIndex' (i.e. valid only for assembling process on the current element)
    //  (C3) combined version
    
    //  (C1) global version via 'Vertex* vrt'
        bool is_onInterfaceVertex_globalAccess(Vertex* vrt, size_t vrtIndex)
        {  return m_spInterfaceVrtMarker->is_marked(vrt); }
    
        bool is_OutsideVertex_globalAccess(Vertex* vrt, size_t vrtIndex)
        {  return m_spOutsideMarker->is_marked(vrt); }
    
        bool is_nearInterfaceVertex_globalAccess(Vertex* vrt, size_t vrtIndex)
        {  return m_spNearInterfaceVrtMarker->is_marked(vrt); }
    
    //  (C2) local version via 'size_t vrtIndex'
        bool is_onInterfaceVertex_localAccess(Vertex* vrt, size_t vrtIndex)
        { return !check_vertex_modus(INSIDE, vrtIndex, get_orientation()); }
    
        bool is_OutsideVertex_localAccess(Vertex* vrt, size_t vrtIndex)
        { return check_vertex_modus(OUTSIDE, vrtIndex, get_orientation()); }
    
        bool is_nearInterfaceVertex_localAccess(Vertex* vrt, size_t vrtIndex)
        { return check_vertex_modus(ON_INTERFACE, vrtIndex, get_orientation()); }

    //  (C3) combined version
        bool is_onInterfaceVertex(Vertex* vrt, size_t vrtIndex)
        {
            if ( !m_bBoolMarkerInit ) return is_onInterfaceVertex_localAccess(vrt, vrtIndex);
            else                      return is_onInterfaceVertex_globalAccess(vrt, vrtIndex);
        }
        bool is_OutsideVertex(Vertex* vrt, size_t vrtIndex)
        {
            if ( !m_bBoolMarkerInit ) return is_OutsideVertex_localAccess(vrt, vrtIndex);
            else                      return is_OutsideVertex_globalAccess(vrt, vrtIndex);
        }
        bool is_nearInterfaceVertex(Vertex* vrt, size_t vrtIndex)
        {
            if ( !m_bBoolMarkerInit ) return is_nearInterfaceVertex_localAccess(vrt, vrtIndex);
            else                      return is_nearInterfaceVertex_globalAccess(vrt, vrtIndex);
        }

        
    //////////////////////////////////////////////////////////
    /// further base class methods:
    //////////////////////////////////////////////////////////
        
    // getter and setter for orientation of interface
        void set_orientation(const int orientation) { m_spInterfaceProvider->set_orientation(orientation); }
        const int get_orientation() const{ return m_spInterfaceProvider->get_orientation(); }
    
    // the 'threshold' defines the bandwidth around the immersed interface, in which a node
    //  counts as 'OUTSIDE' or 'ON_INTERFACE' during call of 'is_outside()', 'is_nearInterface()
        void set_threshold(size_t level, const number threshold)
        { m_vThresholdOnLevel[level] = threshold; }
       
        number get_threshold(Vertex* vrt)
        { return m_vThresholdOnLevel[m_spMG->get_level(vrt)]; }
    
    // returns the number of cut Elements on the level 'levIndex':
        size_t get_numCutElements(const int levIndex)
            {  return m_vvElemListCut[levIndex].size(); }

    
    ///////////////////////////////////////////////////////////////////////
    /// moving interface methods:
    ///////////////////////////////////////////////////////////////////////
    
    // updates the location of the interface in case of a moving interface
    // no base implementation provided; dependent on the specific application
        void update_interface(const int topLevel, const number deltaT)
        { UG_THROW("in 'CutElementHandlerBase::update_interface()': needs to be implemented by derived class!\n");}

    
    ///////////////////////////////////////////////////////////////////////
    /// class member
    ///////////////////////////////////////////////////////////////////////
    
        
    // data for geometric location access of vertices Vertex* vrt
        SmartPtr<MultiGrid> m_spMG;
        position_attachment_type m_aPos;		///<Position Attachment
        position_accessor_type m_aaPos;			///<Accessor
    
    /// associotion of 'GridLevel' with an 'Index' value for access to m_vvv-lists:
        std::map<GridLevel, size_t> m_Map;
    
    /// maximal distance, up to which a node will be couted as 'on interface'
        std::vector<number> m_vThresholdOnLevel;

    // data for LOCAL access
        ElementModus m_elementModus;
    
    // indexing: m_vVertexMode[...] = m_vvVertexMode[current Orientation][...]
        std::vector<VertexModus> m_vVertexMode;
    
    // indexing: [orientation][VertexModus]
    //  -> for both orientations of the interface, the cut element
    //      is one part of the element, inheriting according vertices
    //  -> their according VertexMode is stored in [][VertexModus]
    
        std::vector<std::vector<VertexModus> > m_vvVertexMode;
    
        bool m_bBoolMarkerInit;
        
    // data for GLOBAL access
    // REMARK: for explanation of the different types of vertices, elements see File!
        SmartPtr<BoolMarker> m_spInterfaceVrtMarker;  		// marks: vertices for 'is_outside() = true'
                                                            //  AND which lie on a cut element
        SmartPtr<BoolMarker> m_spNearInterfaceVrtMarker;  	// marks: vertices for 'is_nearInterface() = true'
                                                            //          i.e. potentially inside domain
        SmartPtr<BoolMarker> m_spCutMarker;                 // marks: elements, edges, vertices
        SmartPtr<BoolMarker> m_spOutsideMarker;             // marks: elements, edges, vertices
                                                            //    --> for vertices: 'true' for 'is_outside() = true'
                                                            //         OR 'is_nearInterface() = true'
        SmartPtr<BoolMarker> m_spInsideMarker;              // marks: elements, edges, vertices
                                                            //    --> for vertices: 'true' for 'is_outside() = false'
                                                            //          ( => is_nearInterface() = false'!)
        
        
    /// element lists computed during 'update_interface_data()' to enable looping over
    //      distinguished element types, i.e. all, cut elements or outside elements
    //  indexing: [levIndex][counter]
        std::vector<std::vector<grid_base_object* > > m_vvElemList;
        std::vector<std::vector<grid_base_object* > > m_vvElemListCut;
        std::vector<std::vector<grid_base_object* > > m_vvElemListOutside;
    
    
    /// contains interface infos (like radius, center ...)
        SmartPtr<IInterfaceProvider<dim> > m_spInterfaceProvider;
    
};
    
    
    
    
template <int TWorldDim>
class CutElementHandler_FlatTop : public CutElementHandlerBase<TWorldDim>
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
        CutElementHandler_FlatTop(SmartPtr<MultiGrid> mg, const char* fctNames,
                                 SmartPtr<ParticleProviderSphere<dim> > interfaceProvider);
        CutElementHandler_FlatTop(SmartPtr<MultiGrid> mg, const char* fctNames,
                                 SmartPtr<ParticleProviderEllipse<dim> > interfaceProvider);
        
        /// destructor
        virtual ~CutElementHandler_FlatTop() {}
        
    
    ////////////////////////////////////////////////////////////////////////////////
    /// methods called during initialisation of the interface:
    ////////////////////////////////////////////////////////////////////////////////

    /// initializes crucial data for all cut element computations
    ///     especially'ElementModus' and 'VertexModus' categories
        void update_multigrid_data(ConstSmartPtr<DoFDistribution> dd, const int levIndex);

    /// methods called during 'update_multigrid_data':
        
    /// sets all BoolMarker and fills the 'm_vvElemList...'
        void update_interface_data(ConstSmartPtr<DoFDistribution> dd, const int levIndex);
        
    /// computes the global indices of the nodes within a particle, where the DoFs of the
    //  particle velocities will be (algorithmically!) located, i.e. the solution is stored
    //    --> for the rigid body motion simulation, the DoFs of the particle will be located
    //          simply in any "unused" = "free" node inside the particle
        void update_global_indices(ConstSmartPtr<DoFDistribution> dd, const int levIndex);
        
    ////////////////////////////////////////////////////////////////////////////////
    /// virtual base class methods, needed for all basic computations
    ////////////////////////////////////////////////////////////////////////////////
 
    /// for explanations regarding the relation between near-interface, outside,
    /// on-interface: see the File
        
    /// returns boolian for 'does 'vrt' lie near interface?'
        bool is_nearInterface(Vertex* vrt)
        { UG_THROW("not implemented, since the alternative method" <<
                   "'is_nearInterface_with_given_index()' is used by this class!\n");}
    
    // returns boolian for 'is near interface?' for a given 'prtIndex'
        bool is_nearInterface_with_given_index(const int prtIndex, Vertex* vrt);

    /// returns boolian for 'does 'vrt' lie outside the domain?'
        bool is_outsideDomain(Vertex* vrt)
        { return is_outsideFluid(vrt); }
        
        bool is_outsideFluid(Vertex* vrt);
        
    // returns boolian for 'is outside?' for a given 'prtIndex'
        bool is_insideParticle_with_given_index(const int prtIndex, Vertex* vrt);
        
    // Remark: not only returns the boolian for 'is outise?', but also writes the
    //      according 'PrtIndex' into data, in which the 'vrt' is located
        bool is_outsideFluid(int& PrtIndex, Vertex* vrt);
        
    ////////////////////////////////////////////////////////////////////////////////
    // methods needed during call of 'is_outsideDomain(), 'is_nearInterface()':
    ////////////////////////////////////////////////////////////////////////////////
        
    // computes the 'level-set' value of a point w.r.t. the interface, in order to
    //  derive its 'VertexModus', i.e. INSIDE, OUTSIDE, ON_INTERFACE
    //  --> evaluate the distance between a vrt and the interface
        number get_LSvalue_byPosition(MathVector<dim> vrtPos)
            { return m_spInterfaceProvider->get_LSvalue_byPosition(vrtPos); }
        
    // new 'get_LSvalue()'-method needed for 'is_outsideDomain()' and 'is_nearInterface()':
    //  -->  with prtIndex as additional parameter compared to base class
        number get_LSvalue(Vertex* vrt, const int prtIndex)
            { return get_LSvalue_byPosition(this->m_aaPos[vrt], prtIndex); }
        
        number get_LSvalue_byPosition(MathVector<dim> vrtPos, const int prtIndex)
            { return m_spInterfaceProvider->get_LSvalue_byPosition(vrtPos, prtIndex);}
        
    ////////////////////////////////////////////////////////////////////////////////
    /// (A) 'compute' and (B) 'get' the 'ElementModus':
    ////////////////////////////////////////////////////////////////////////////////

    /// (A) computes the modus of the element as cut or non-cut
    /// (= inside or outside of the domain)  and marks according BoolMarker
        ElementModus compute_element_modus(int prtIndex, GridObject* elem);
        ElementModus compute_element_modus(GridObject* elem)
        { UG_THROW("Attention: the alternative implementation with additional parameter 'prtIndex' has to be calles!\n"); }
        
    /// derives the 'ElementModus' from the boolian of the BoolMarker
    //    => global access to information (not only locally, during assembling
    //          for that specific 'elem' (see base class for explanation LOCAL, GLOBAL)
        ElementModus get_element_modus(GridObject* elem);
        
    //////////////////////////////////////////////////////////////////////////////////
    /// checks for the 'VertexModus' of a vrt based on the
    //      --> GLOBAL acces, see comments in base class)
    ///     --> BoolMarker for global access set during 'update_interface_data()'
        
    /// returns true, if vertex lies OUTSIDE fluid AND near to interface
        bool is_onInterfaceVertex(Vertex* vrt, size_t vrtIndex = 0)
            { return this->m_spInterfaceVrtMarker->is_marked(vrt); }

    /// returns true, if vertex lies OUTSIDE fluid => no DoF
        bool is_OutsideVertex(Vertex* vrt, size_t vrtIndex = 0)
            { return this->m_spOutsideMarker->is_marked(vrt); }
 
    /// returns true, if vertex lies INSIDE fluid AND near to interface
        bool is_nearInterfaceVertex(Vertex* vrt, size_t vrtIndex = 0)
            { return this->m_spNearInterfaceVrtMarker->is_marked(vrt); }
        
        
    /////////////////////////////////////////////////////////////////////////////////////////////
    /// methods, which handle the DoFs (and associated global indices) for the particle
    ///   velocities, which are located in specified grid nodes inside each particle
    /////////////////////////////////////////////////////////////////////////////////////////////
        
    // checks weather node is transDoF OR rotDoF:
    /// returns true, if the local 'dofIndex' of a node of an element, which lies in the particle,
    //      is one of the selected DoFs for the particle velocities
    //  => outside domain, but a DoF!
        bool is_extraDoF(DoFIndex dofIndex, int levIndex);
      
    // flags being set during initialisation, which indicate, whether the velocities of a particle
    //  will be a DoF or not
    //  --> for simulations of a fixed particle (e.g. cylinder benchmark), the particle velocities are
    //      not a DoF;
        bool get_DoF_modus_linear(int prtIndex) { return m_spInterfaceProvider->get_DoF_modus_linear(prtIndex);}
        bool get_DoF_modus_angular(int prtIndex){ return m_spInterfaceProvider->get_DoF_modus_angular(prtIndex);}
        
    /// Access methods for the global indices of the DoFs of the particle velocities
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
        
     // returns single index for a given velocity component 'cmp'
        DoFIndex get_transInd_Comp(int levIndex, int prtIndex, size_t cmp)
            { 	if ( cmp == dim ) UG_THROW("no acces to pressure index in transInd! EXIT...\n");
                return m_vvvGlobalIndices_linearVel[levIndex][prtIndex][cmp]; }
        DoFIndex get_rotInd_Comp(int levIndex, int prtIndex, size_t cmp)
            { 	if ( cmp == dim ) UG_THROW("no acces to pressure index in rotInd! EXIT...\n");
                return m_vvvGlobalIndices_angularVel[levIndex][prtIndex][cmp]; }
       
    // returns all indices for transVel and rotVel
        void get_global_indices(std::vector<DoFIndex>& transInd, std::vector<DoFIndex>& rotInd,
                                const int levIndex, const int prtIndex)
        {
            transInd = get_transInd(levIndex, prtIndex);
            rotInd = get_rotInd(levIndex, prtIndex);
        }
        
    // writes the solution of the particle velocities into the data storage provided in the
    // 'ParticleProvider' called during PartcleMapper::modify_GlobalSol():
        void set_extraSolTrans(number solution, size_t prtIndex, int timeSeriesInd, int cmp)
            { m_spInterfaceProvider->set_linear_velocity(solution, prtIndex, timeSeriesInd, cmp); }
        void set_extraSolTrans(MathVector<dim> solution, size_t prtIndex, int timeSeriesInd)
            { m_spInterfaceProvider->set_linear_velocity(solution, prtIndex, timeSeriesInd); }
        void set_extraSolRot(number solution, size_t prtIndex, int timeSeriesInd, int cmp)
            { m_spInterfaceProvider->set_angular_velocity(solution, prtIndex, timeSeriesInd, cmp); }
        void set_extraSolRot(MathVector<dim> solution, size_t prtIndex, int timeSeriesInd)
            { m_spInterfaceProvider->set_angular_velocity(solution, prtIndex, timeSeriesInd); }
        
    /// get solution values
        MathVector<dim> get_transSol(size_t prtIndex, size_t timeSeriesInd)
            { return m_spInterfaceProvider->get_linear_velocity(prtIndex, timeSeriesInd); }
        MathVector<dim> get_rotSol(size_t prtIndex, size_t timeSeriesInd)
            { return m_spInterfaceProvider->get_angular_velocity(prtIndex, timeSeriesInd); }
        
    // returns the rotation matrix associated with the radial vector
        MathMatrix<TWorldDim,TWorldDim> get_rotationMat(MathVector<TWorldDim> radialVector);
        

    /////////////////////////////////////////////////////////////////////////////////////////////
    /// some functionality forwarded to the 'ParticleProvider' for according computations
    /////////////////////////////////////////////////////////////////////////////////////////////
        
        size_t num_particles() const			{ return m_spInterfaceProvider->num_particles();}
        number get_density(int prtIndex)		{ return m_spInterfaceProvider->get_density(prtIndex);}
        MathVector<dim> get_center(int prtIndex){ return m_spInterfaceProvider->get_center(prtIndex);}
        void set_center(MathVector<dim> center, int prtIndex){ m_spInterfaceProvider->set_center(center, prtIndex);}
        number get_radius(int prtIndex)         { return m_spInterfaceProvider->get_radius(prtIndex);}
        
    /// methods called by local-to-global-mapper 'ParticleMapper':
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
        

    //////////////////////////////////////////////////////////////////////////////////
    /// furhter functionality
    //////////////////////////////////////////////////////////////////////////////////
        
    /// for moving interfaces, the location needs to be updated
        void update_interface(const int topLevel, const number deltaT)
            { update_prtCoords(topLevel, deltaT); }
        
        void update_prtCoords(const int topLevel, const number deltaT);
        
    // computes the intersection point of the interface with an edge of the cut element
        bool get_intersection_point(MathVector<dim>& Intersect, const MathVector<dim>& vrtPosOut,
                                    const MathVector<dim>& vrtPosIn, const int PrtIndex)
            { return this->m_spInterfaceProvider->get_intersection_point(Intersect, vrtPosOut, vrtPosIn, PrtIndex); }
        bool get_intersection_point(MathVector<dim>& Intersect, const MathVector<dim>& vrtPosOut,
                                    const MathVector<dim>& vrtPosIn, const int PrtIndex, std::vector<number>& alphaOut)
            { return this->m_spInterfaceProvider->get_intersection_point(Intersect, vrtPosOut, vrtPosIn, PrtIndex, alphaOut); }
        
        
    // getter and setter for the index of the particle, cutting the current element
        void set_prtIndex(const int prtIndex) { m_spInterfaceProvider->set_prtIndex(prtIndex); }
        const int get_prtIndex() const        { return m_spInterfaceProvider->get_prtIndex(); }
        int get_prtIndex(size_t dof)
        {
            UG_THROW("attention: You want to get the particle index of an array. Make sure that the data was set!\n");
            return m_vPrtIndex[dof];
        }
        
    // called during update_elem() in InterfaceHandler
        void compute_and_set_prtIndex(GridObject* elem);
        
    /// returns true, if 'elem' is a pyramid; used for setting the BoolMarker appropriately
        bool element_is_pyramid(grid_base_object* elem);
        
        
    /////////////////////////////////////////////////////////////////////////////////////////////
    // output
    /////////////////////////////////////////////////////////////////////////////////////////////

    /// writes particle velocities into a file -> lua-call
        void print_velocity(const MathVector<dim>& transSol, const MathVector<dim>& rotSol, const int prtIndex,
                            const bool isTimedep, const number time, const char* filename)
            {  m_spInterfaceProvider->print_velocity(transSol, rotSol, prtIndex, isTimedep, time, filename); }
        
        void print_elem_lists(ConstSmartPtr<DoFDistribution> dd);

    // returns the number of cut Elements on the level 'levIndex':
        size_t get_numCutElements( const int levIndex, const size_t prtIndex)
        { return m_vvvElemListCut[levIndex][prtIndex].size(); }

    /////////////////////////////////////////////////////////////////////////////////////////////
    /// methods for parallel computations
    /////////////////////////////////////////////////////////////////////////////////////////////

        void set_mpi_routine(int val)
        { active_mpi_routine = val; }
        
#ifdef UG_PARALLEL
        void synchronize_particles(int levIndex);
#endif
/*
        bool valid_prt_information(int prtIndex)
            { return particleData[prtIndex].valid; }
*/
        
    //////////////////////////////////////////////////////////
    /// class member
    //////////////////////////////////////////////////////////
        
        const char* m_fctNames; 
       
    //ToDo:brauche ich das??
        std::vector<int> m_vPrtIndex;
        
    /// global indices of the DoFs for the particle velocities
    /// indexing: [levIndex][prtIndex][cmp]
        std::vector< std::vector< std::vector<DoFIndex> > > m_vvvGlobalIndices_linearVel;
        std::vector< std::vector< std::vector<DoFIndex> > > m_vvvGlobalIndices_angularVel;
      
    /// element lists computed during 'update_interface_data()' to enable looping over
    //      distinguished element types, i.e. all, cut elements or outside elements
    /// indexing: [levIndex][prtIndex][counter]
        std::vector<vvElemList> m_vvvElemList;
        std::vector<vvElemList> m_vvvElemListCut;
        std::vector<vvElemList> m_vvvElemListOutside;
        
    /// contains radius, center and density of all given particles
        SmartPtr<ParticleProvider<dim> > m_spInterfaceProvider;
        
        
    //    CombinedParticleData particleData;
        std::size_t active_mpi_routine;
        
        std::map<Vertex*, int> m_vPrtIndices;

};
    
    
    
    
template <int TWorldDim>
class CutElementHandler_TwoSided : public CutElementHandlerBase<TWorldDim>
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
        CutElementHandler_TwoSided(SmartPtr<MultiGrid> mg, const char* fctNames,
                                  SmartPtr<DiffusionInterfaceProvider<dim> > interfaceProvider);
        
        /// destructor
        virtual ~CutElementHandler_TwoSided() {}
        

    ////////////////////////////////////////////////////////////////////////////////
    /// virtual base class methods, needed for all basic computations
    ////////////////////////////////////////////////////////////////////////////////
        
    /// for explanations regarding the relation between near-interface, outside,
    /// on interface: see the 'Info File 1'

    /// returns boolian for 'does 'vrt' lie near interface?'
        bool is_nearInterface(Vertex* vrt);
        
    /// returns boolian for 'does 'vrt' lie outside the domain?'
        bool is_outsideDomain(Vertex* vrt);
        
    // Remark: not only returns the boolian for 'is outise?', but also writes the
    //      according 'PrtIndex' into data, in which the 'vrt' is located
        bool is_outsideDomain(int& PrtIndex, Vertex* vrt);
        
        
    ////////////////////////////////////////////////////////////////////////////////
    // methods needed during call of 'is_outsideDomain(), 'is_nearInterface()':
    ////////////////////////////////////////////////////////////////////////////////
        
    // computes the 'level-set' value of a point w.r.t. the interface, in order to
    //  derive its 'VertexModus', i.e. INSIDE, OUTSIDE, ON_INTERFACE
    //  --> evaluate the distance between a vrt and the interface
        number get_LSvalue_byPosition(MathVector<dim> vrtPos)
        { return m_spInterfaceProvider->get_LSvalue_byPosition(vrtPos); }
        
    // new 'get_LSvalue()'-method needed for 'is_outsideDomain()' and 'is_nearInterface()':
    //  -->  with prtIndex as parameter
        number get_LSvalue(Vertex* vrt, const int prtIndex)
        { return get_LSvalue_byPosition(this->m_aaPos[vrt], prtIndex); }
        
        number get_LSvalue_byPosition(MathVector<dim> vrtPos, const int prtIndex)
        { return m_spInterfaceProvider->get_LSvalue_byPosition(vrtPos, prtIndex);}
  
      
    ////////////////////////////////////////////////////////////////////////////////////////
    /// checks for the 'VertexModus' of a vrt
    ////////////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////////////
    /// REMARK: We can NOT use the BoolMarker for the methods 'is_onInterfaceVertex()',
    ///         'is_OutsideVertex(), and 'is_nearInterfaceVertex()', since for the
    ///         two-sided case, they will be overwritten after switch to the other
    ///         orientation of the interface
    //              => LOCAL acces, as specified in the base class, see comments there!
    ////////////////////////////////////////////////////////////////////////////////////////
    
    /// returns true, if vertex lies OUTSIDE fluid AND near to interface
        bool is_onInterfaceVertex(Vertex* vrt, size_t vrtIndex)
        {  return !check_vertex_modus(INSIDE, vrtIndex, this->get_orientation()); }
        
    /// returns true, if vertex lies OUTSIDE fluid => no DoF
        bool is_OutsideVertex(Vertex* vrt, size_t vrtIndex)
        { return check_vertex_modus(OUTSIDE, vrtIndex, this->get_orientation()); }
        
    /// returns true, if vertex lies INSIDE fluid AND near to interface
        bool is_nearInterfaceVertex(Vertex* vrt, size_t vrtIndex)
        { return check_vertex_modus(ON_INTERFACE, vrtIndex, this->get_orientation()); }
        
    //////////////////////////////////////////////////////////////////////////
    /// access methods for the 'VertexModus' WITHOUT 'Vertex'-Parameter:
    //      --> needed for access during loc_to_glob_mapper!
        
        VertexModus get_vertex_modus(size_t vrtIndex, const int localind) { return this->m_vvVertexMode[localind][vrtIndex];}
        
    /// returns true, if vertex lies INSIDE fluid and near to interface => no DoF
        bool check_vertex_modus(VertexModus vrtModus, size_t vrtIndex, const int interfaceOrientation)
        {
            size_t localInd = -0.5*interfaceOrientation + 0.5;
            if ( get_vertex_modus(vrtIndex, localInd) == vrtModus )
                return true;
            return false;
        }
        
        
    ///////////////////////////////////////////////////////////////////////////////////////
    ///  further getter and setter methods
    ///////////////////////////////////////////////////////////////////////////////////////
        
    // getter and setter for prtIndex of the particle, cutting the current element
        void set_prtIndex(const int prtIndex) { m_spInterfaceProvider->set_prtIndex(prtIndex); }
        const int get_prtIndex()         const{ return m_spInterfaceProvider->get_prtIndex(); }
        
        int get_prtIndex(size_t dof)
        {
            UG_THROW("attention: You want to get the particle index of an array. Make sure that the data was set!\n");
            return m_vPrtIndex[dof];
        }
        
     // called during update_elem() in InterfaceHandler
        void compute_and_set_prtIndex(GridObject* elem);
        
        size_t get_or_insert_vertex_near(const MathVector<dim>& vrtPos);
        const size_t get_numNearVerticesPos() const { return m_verticesNearPos.size(); }
        
        // ToDo: 'm_bElementNearInterface' und method brauche ich nicht mehr!
        bool is_element_near_interface() { return m_bElementNearInterface;}
        
    
    /////////////////////////////////////////////////////////////////////////////////////////////
    /// some functionality forwarded to the 'ParticleProvider' for according computations
    /////////////////////////////////////////////////////////////////////////////////////////////
        
        size_t num_particles() const			 { return m_spInterfaceProvider->num_particles();}
        number get_theta(int prtIndex)           { return m_spInterfaceProvider->get_theta(prtIndex);}
        number get_density(int prtIndex)		 { return m_spInterfaceProvider->get_density(prtIndex);}
        MathVector<dim> get_center(int prtIndex) { return m_spInterfaceProvider->get_center(prtIndex);}
        
    /////////////////////////////////////////////////////////////////////////////////////////////
    // output
    /////////////////////////////////////////////////////////////////////////////////////////////
    
    /// writes particle velocities into a file -> lua-call

    
        void print_velocity(const MathVector<dim>& transSol, const MathVector<dim>& rotSol, const int prtIndex,
                            const bool isTimedep, const number time, const char* filename)
            {  m_spInterfaceProvider->print_velocity(transSol, rotSol, prtIndex, isTimedep, time, filename); }
        
    //////////////////////////////////////////////////////////
    /// class member
    //////////////////////////////////////////////////////////
        
        //ToDo:brauche ich das??
        std::vector<int> m_vPrtIndex;

        std::map<Vertex*, int> m_vPrtIndices;

        
    /// element lists computed during 'update_interface_data()' to enable looping over
    //      distinguished element types, i.e. all, cut elements or outside elements
    /// indexing: [levIndex][prtIndex][counter]
        std::vector<vvElemList> m_vvvElemList;
        std::vector<vvElemList> m_vvvElemListCut;
        std::vector<vvElemList> m_vvvElemListOutside;
        
    /// contains radius, center and density of all given particles
        SmartPtr<DiffusionInterfaceProvider<dim> > m_spInterfaceProvider;
        
        
        bool m_bElementNearInterface;
        std::map<MathVector<dim>, size_t> m_MapNearVertices;
        std::vector<MathVector<dim> > m_verticesNearPos;
       
        
};

}// end namespace ug

#include "cut_element_handler_base_impl.h"
#include "cut_element_handler_flat_top_impl.h"
#include "cut_element_handler_two_sided_impl.h"


#endif /* CUT_ELEMENT_HANDLER_H_ */
