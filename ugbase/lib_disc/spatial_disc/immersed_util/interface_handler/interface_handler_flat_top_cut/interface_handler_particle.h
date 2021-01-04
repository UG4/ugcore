/*
 * interface_handler_local.h
 *
 *  Created on: 15.01.2015
 *      Author: suze
 */

#ifndef INTERFACE_HANDLER_LOCAL_PARTICLE_H_
#define INTERFACE_HANDLER_LOCAL_PARTICLE_H_

#include "lib_disc/spatial_disc/disc_util/geom_provider.h"
#include "lib_grid/multi_grid.h"

#include "../interface_handler_base.h"

#include <map>

namespace ug{

        
template <int TDim, int TWorldDim, class TInterfaceHandler>
class DimFV1FTGeometry;

template <int TWorldDim>
class InterfaceHandlerLocalParticle : public InterfaceHandlerLocalBase<TWorldDim>
//class FlatTopHandler : public IInterfaceHandlerLocal
{

	public:
	/// world dimension
		static const int dim = TWorldDim;

	///	max number of geometric objects in a dimension
	// 	(most objects in 1 dim, i.e. number of edges, but +1 for 1D)
		static const int maxMid = DimFV1FTGeometry<dim, dim, InterfaceHandlerLocalParticle<dim> >::maxNumSCVF + 1;

	///	used traits
		typedef fv1_dim_traits<dim, dim> traits;

	/// used boundary face type
		typedef typename DimFV1FTGeometry<dim, dim, InterfaceHandlerLocalParticle<dim> >::BF interfaceBF;

	///	Type of position coordinates
		typedef MathVector<dim> position_type;

	///	Type of Position Attachment
		typedef Attachment<position_type> position_attachment_type;

	///	Type of Accessor to the Position Data Attachment
	 	typedef Grid::VertexAttachmentAccessor<position_attachment_type> position_accessor_type;


  		InterfaceHandlerLocalParticle(SmartPtr<CutElementHandler_FlatTop<dim> > cutElementHandler,
  		 				number fluidDensity, number fluidKinVisc);


    	virtual ~InterfaceHandlerLocalParticle()	{}

    
    //////////////////////////////////////////////////////////////////////////
    /// virtual methods in 'IInterfaceHandlerLocal': need to be impl. here
    //////////////////////////////////////////////////////////////////////////

    // simple forwarding to the associated 'CutElementHandler' class:
    
        bool is_onInterfaceVertex(Vertex* vrt, size_t vrtIndex = 0)
            { return m_spCutElementHandler->is_onInterfaceVertex(vrt, vrtIndex);	}

        bool is_OutsideVertex(Vertex* vrt, size_t vrtIndex = 0)
            { return m_spCutElementHandler->is_OutsideVertex(vrt, vrtIndex); }
    
        bool is_nearInterfaceVertex(Vertex* vrt, size_t vrtIndex = 0)
            { return m_spCutElementHandler->is_nearInterfaceVertex(vrt, vrtIndex); }
    
    
    //////////////////////////////////////////////////////////////////////////
    /// virtual methods in 'InterfaceHandlerLocalBase': need to be impl. here
    //////////////////////////////////////////////////////////////////////////
    
    // forwarding to m_spCutElementHandler:
        number get_LSvalue_byPosition(MathVector<dim> vrtPos)
        { return m_spCutElementHandler->get_LSvalue_byPosition(vrtPos); }


    // forwarding to m_spCutElementHandler:
        bool get_intersection_point(MathVector<dim>& Intersect, Vertex* vrtOutsideCirc, Vertex* vrtInsideCirc)
        {
            const int orientation = this->get_orientation();

            const MathVector<dim>& vrtPosOut = this->m_aaPos[vrtOutsideCirc];
            const MathVector<dim>& vrtPosIn  = this->m_aaPos[vrtInsideCirc];
            
            const int prtIndex = get_prtIndex();
            if ( prtIndex == -1 ) UG_THROW("'get_intersection_point()': value of prtIndex not valid!\n");

            if ( orientation == 1 )
                return m_spCutElementHandler->get_intersection_point(Intersect, vrtPosOut, vrtPosIn, prtIndex);
        // inverse order of 'vrtPosOut' and 'vrtPosIn' for call of 'get_intersection_point()'
        // to avoid error for alpha < 0:
            else if ( orientation == -1 )
            {
                return m_spCutElementHandler->get_intersection_point(Intersect, vrtPosIn, vrtPosOut, prtIndex);
            }
            else
                UG_THROW("in InterfaceHandlerLocalDiffusion::get_intersection_point(): orientationInterface not set!\n");
        }
    
    // forwarding to m_spCutElementHandler:
        bool get_intersection_point(MathVector<dim>& Intersect, Vertex* vrtOutsideCirc, Vertex* vrtInsideCirc, std::vector<number>& alphaOut)
        {
            const int orientation = this->get_orientation();

            const MathVector<dim>& vrtPosOut = this->m_aaPos[vrtOutsideCirc];
            const MathVector<dim>& vrtPosIn  = this->m_aaPos[vrtInsideCirc];
        
            const int prtIndex = get_prtIndex();
            if ( prtIndex == -1 ) UG_THROW("'get_intersection_point()': value of prtIndex not valid!\n");

            if ( orientation == 1 )
                return m_spCutElementHandler->get_intersection_point(Intersect, vrtPosOut, vrtPosIn, prtIndex, alphaOut);
        // inverse order of 'vrtPosOut' and 'vrtPosIn' for call of 'get_intersection_point()'
        // to avoid error for alpha < 0:
            else if ( orientation == -1 )
            {
                return m_spCutElementHandler->get_intersection_point(Intersect, vrtPosIn, vrtPosOut, prtIndex, alphaOut);
            }
            else
                UG_THROW("in InterfaceHandlerLocalDiffusion::get_intersection_point(): orientationInterface not set!\n");
        }
    
    //////////////////////////////////////////////////////////
    /// InterfaceHandlerLocalBase methods
    //////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////
    /// central method: update_elem():
    ///     it is called during TFVGeom:update() and computes the new cornders
    ///     i.e. 'vCornerCoords' of the cut element. Based on that ALL local
    ///     TFVGeom-computations (compute ip's, normals, gradients, ...) follow
    ///     as for the standard case
        bool update(GridObject* elem, const MathVector<TWorldDim>* vCornerCoords);
    
    /// re-implement the base class method, since the cut element data does not need to
    //   be recomputed => herein: call of 'get_element_modus() instead of 'compute_element_modus()'
 	    bool update_elem(GridObject* elem, const MathVector<TWorldDim>* vCornerCoords)
 	    {
 	    	compute_and_set_prtIndex(elem);
	    	return update(elem, vCornerCoords);
	    }

    // called by FV1FTGeom::update() to fill data 'm_vBF' with the faces on
    //  the interface
        void update_inner_boundary(const MathVector<TWorldDim>* vCornerCoords)
        {
            if ( this->StdFV_assembling() )
                update_inner_boundary_radial_StdFV(vCornerCoords);
            else
                update_inner_boundary_radial(vCornerCoords);
        }
    
        void update_inner_boundary_for2()
        {
            if ( this->StdFV_assembling() )
            {update_inner_boundary_radial_for2_StdFV();}
            else
            {update_inner_boundary_radial_for2();}
        }
    
    // new methods for the update of inner boundary faces:
    //  -> compute the data 'm_vRadialAtIP' and 'm_vRadialAtCo'
    //  -> in case of a moving particle, the radial direction is needed
    //      for the computation of the stresses
    //  -> this data is used in local-to-global-mapper
        void update_inner_boundary_radial(const MathVector<TWorldDim>* vCornerCoords);
        void update_inner_boundary_radial_StdFV(const MathVector<TWorldDim>* vCornerCoords);

    
    // REMARK: _for2()-methods not finally tested!!
        void update_inner_boundary_radial_for2();
        void update_inner_boundary_radial_for2_StdFV();

    
    //////////////////////////////////////////////////////////
    ///  original class methods + helper methods for 'update()'
    //////////////////////////////////////////////////////////
    
    /// collects all corners of the cut element
        int CollectCorners_FlatTop_3d(GridObject* elem);
        int get_cutMode(std::vector<Vertex*> vVertex);
        int CollectCorners_FlatTop_Prism3(GridObject* elem);
        int CollectCorners_FlatTop_Prism4(GridObject* elem);
        int CollectCorners_FlatTop_Pyramid(GridObject* elem);
        int CollectCorners_FlatTop_originalTet(GridObject* elem);
    
    /// !!! not implemented new => taken from the base class:
/*
        void ResortQuadrilateral_for2(std::vector<MathVector<dim> > vQuadriCorners);
     
        int CollectCorners_StdFV_for2(GridObject* elem);
        int CollectCorners_FlatTop_2d_for2(GridObject* elem);
        int CollectTriangle_and_Quadri_for2(GridObject* elem, Vertex* vrtInside);
        int CollectQuadrilateral_besideTri_for2(GridObject* elem);
        int CollectQuadrilateral_for2(GridObject* elem);
     
    	void compute_cut_element_data_for2(GridObject* elem);
*/

    ///////////////////////////////////////////////////////////////
    /// forwarding to 'm_spCutElementHandler'
    ///////////////////////////////////////////////////////////////

        int get_prtIndex()          { return m_spCutElementHandler->get_prtIndex(); }
        int getPrtIndex(size_t dof) { return m_spCutElementHandler->get_prtIndex(dof); }
    
        void compute_and_set_prtIndex(GridObject* elem)
        { m_spCutElementHandler->compute_and_set_prtIndex(elem); }

    ///////////////////////////////////////////////////////////////
    /// access to the particle velocities, stored in the
    ///     'ParticleProvider' class
    
        MathVector<dim> get_transSol(size_t prtIndex, size_t timeSeriesInd)
            { return m_spCutElementHandler->get_transSol(prtIndex, timeSeriesInd); }
        MathVector<dim> get_rotSol(size_t prtIndex, size_t timeSeriesInd)
            { return m_spCutElementHandler->get_rotSol(prtIndex, timeSeriesInd); }
    
        MathMatrix<dim,dim> get_rotationMat(MathVector<dim> radialVector)
            { return m_spCutElementHandler->get_rotationMat(radialVector); }


	///////////////////////////////////////////////////////////////
	/// new methods
	///////////////////////////////////////////////////////////////

    /// for boundary computations
        bool is_boundary_face_for2(const size_t sideID);

    /// acces to 'm_vQuadriCorners_for2'
        const std::vector<MathVector<dim> > quadriCorners() const { return m_vQuadriCorners_for2; }
    
    // switches the order of storage in 'm_vQuadriCorners_for2' and 'this->m_vQuadriOrigID'
    // --> earlier needed during 'CollectQuadrilateral_besideTri_for2()'
    // --> currently not needed anymore
        void switch_order();
    

   	/// called during 'ParticleMapper::modify_LocalData()':
    //  --> resizes the local data (jacobian and defect) due to
    //      the potentially increased number of corners of a cut element
    	void resize_local_indices(LocalVector& locU);
    	void resize_local_indices(LocalVector& locU, size_t numCo);


    ///////////////////////////////////////////////////////////////
    /// new setter methods
    ///////////////////////////////////////////////////////////////

		void set_local_couplings_jac(LocalMatrix rotJ_ind, LocalMatrix rotJ_rot)
			{ m_rotJ_ind = rotJ_ind; m_rotJ_rot = rotJ_rot; }

		void set_local_couplings_def(LocalVector rotD)
			{ m_rotD = rotD;}

    ///////////////////////////////////////////////////////////////
    /// new getter methods
    ///////////////////////////////////////////////////////////////

        SmartPtr<ParticleProvider<dim> > get_particles() { return m_spCutElementHandler->get_particles(); }
        size_t num_particles() const { return m_spCutElementHandler->num_particles();}

	    number get_density(int prtIndex) { return m_spCutElementHandler->get_density(prtIndex); }
	    number get_density_fluid() { return m_fluidDensity; }
	    number get_kinVisc_fluid() { return m_fluidKinVisc; }
        MathVector<dim> get_center(int prtIndex){ return m_spCutElementHandler->get_center(prtIndex);}

        SmartPtr<CutElementHandler_FlatTop<dim> > get_cutElementHandler() { return m_spCutElementHandler; }

	/// access to single entry of 'm_vRadialAtIP'
	    MathVector<dim> radial_at_ip(size_t i) { return m_vRadialAtIP[i]; }
	/// access to single entry of 'm_vRadialAtCo'
	    MathVector<dim> radial_at_co(size_t i) { return m_vRadialAtCo[i]; }


		std::vector<interfaceBF>& get_boundary_faces() { return m_vBF; }
		const LocalIndices& get_local_indices() const  { return m_ind; }

  
    ///////////////////////////////////////////////////////////////
	/// used in local to global mapper:
		number get_rotJ_ind(size_t fct1, size_t dof1, size_t fct2, size_t dof2)
			{ return m_rotJ_ind(fct1, dof1, fct2, dof2); }
		number get_rotJ_rot(size_t fct1, size_t dof1, size_t fct2, size_t dof2)
			{ return m_rotJ_rot(fct1, dof1, fct2, dof2); }
		number get_rotD(size_t fct, size_t dof)
			{ return m_rotD(fct, dof); }


    ///////////////////////////////////////////////////////////////
	/// data for Particle
	///////////////////////////////////////////////////////////////

	   	std::vector<interfaceBF> m_vBF;					// updated during FV1FTGeom::update_inner_boundary_faces()
 
    /// fuid parameter imported by Constructor()
     	number m_fluidDensity;
     	number m_fluidKinVisc;

    /// local data for assembling:
     	LocalMatrix m_rotJ_ind;
     	LocalMatrix m_rotJ_rot;
     	LocalVector m_rotD;

    /// size of local algebra for cut element: 'm_numFct' x 'm_numCo'
     	size_t m_numFct;
    
    /// number of corners of cut element
     	size_t m_numCo;

    /// new local algebra for resized cut element
     	LocalIndices m_ind;

    /// radial data on the interface is needed for the computation of
    //  the rotational component of the stresses
   		std::vector<MathVector<dim> > m_vRadialAtIP; 	// global position of scvf on interface
   		std::vector<MathVector<dim> > m_vRadialAtCo; 	// global position of scv on interface

        std::vector<MathVector<dim> > m_vQuadriCorners_for2;

    /// flag for computation of 'm_vRadialAtIP' during call of 'update_inner_boundary_faces()'
     	const bool m_bRadial_forMassEq_equals_Normal; 				// default = true

    
    /// flag indiating if local data needs to be updated
        bool m_bBndDataNeeded;							// default = false;
                                                        // flag used e.g. in 'FV1FTGeometry::update()'


		SmartPtr<CutElementHandler_FlatTop<dim> > m_spCutElementHandler;
 };


}// end namespace ug

#include "interface_handler_particle_tools.h"
#include "interface_handler_particle_impl.h"

#endif /* INTERFACE_HANDLER_LOCAL_PARTICLE_H_ */
