/*
 * interface_handler_local.h
 *
 *  Created on: 15.01.2015
 *      Author: suze
 */

#ifndef INTERFACE_HANDLER_FLAT_TOP_H_
#define INTERFACE_HANDLER_FLAT_TOP_H_

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


  		InterfaceHandlerLocalParticle(SmartPtr<CutElementHandlerFlatTop<dim> > cutElementHandler,
  		 				number fluidDensity, number fluidKinVisc);


    	virtual ~InterfaceHandlerLocalParticle()	{}

    //////////////////////////////////////////////////////////////////
    /// virtual base class methods which need to be implemented:
    //////////////////////////////////////////////////////////////////
    
    // used for method 'CollectCorners_':
        number get_LSvalue_byPosition(MathVector<dim> vrtPos)
        {
        // ToTo
        //            if ( m_prtIndex == -1 )
        //               UG_THROW("'get_LSvalue_byPosition()': value of m_prtIndex not valid!\n");        
            return m_spCutElementHandler->get_LSvalue_byPosition(vrtPos, m_prtIndex);
        }
        number get_LSvalue(Vertex* vrt, const int prtIndex)
        {
            UG_THROW("Attention in 'get_LSvalue(vrt, prtIndex)': check if second argument = prtIndex is valid!...prtIndex = " << prtIndex << "\n");
            return m_spCutElementHandler->get_LSvalue(vrt, prtIndex);
        }

    
        bool get_intersection_point(MathVector<dim>& Intersect, Vertex* vrtOutsideCirc, Vertex* vrtInsideCirc)
        {
            const MathVector<dim>& vrtPosOut = this->m_aaPos[vrtOutsideCirc];
            const MathVector<dim>& vrtPosIn  = this->m_aaPos[vrtInsideCirc];
        
            if ( this->m_orientationInterface == 1 )
                return m_spCutElementHandler->get_intersection_point(Intersect, vrtPosOut, vrtPosIn, m_prtIndex);
        // inverse order of 'vrtPosOut' and 'vrtPosIn' for call of 'get_intersection_point()'
        // to avoid error for alpha < 0:
            else if ( this->m_orientationInterface == -1 )
            {
                if ( m_prtIndex == -1 ) UG_THROW("'get_intersection_point()': value of m_prtIndex not valid!\n");
                return m_spCutElementHandler->get_intersection_point(Intersect, vrtPosIn, vrtPosOut, m_prtIndex);
            }
            else
                UG_THROW("in InterfaceHandlerLocalDiffusion::get_intersection_point(): m_orientationInterface not set!\n");
        }
    
        bool get_intersection_point(MathVector<dim>& Intersect, Vertex* vrtOutsideCirc, Vertex* vrtInsideCirc, std::vector<number>& alphaOut)
        {
            const MathVector<dim>& vrtPosOut = this->m_aaPos[vrtOutsideCirc];
            const MathVector<dim>& vrtPosIn  = this->m_aaPos[vrtInsideCirc];
        
            if ( this->m_orientationInterface == 1 )
                return m_spCutElementHandler->get_intersection_point(Intersect, vrtPosOut, vrtPosIn, m_prtIndex, alphaOut);
        // inverse order of 'vrtPosOut' and 'vrtPosIn' for call of 'get_intersection_point()'
        // to avoid error for alpha < 0:
            else if ( this->m_orientationInterface == -1 )
            {
                if ( m_prtIndex == -1 ) UG_THROW("'get_intersection_point()': value of m_prtIndex not valid!\n");
                return m_spCutElementHandler->get_intersection_point(Intersect, vrtPosIn, vrtPosOut, m_prtIndex, alphaOut);
            }
            else
                UG_THROW("in InterfaceHandlerLocalDiffusion::get_intersection_point(): m_orientationInterface not set!\n");
        }
    
    //////////////////////////////////////////////////////////
    /// original class methods
    //////////////////////////////////////////////////////////

 	    bool update_elem(GridObject* elem, const MathVector<TWorldDim>* vCornerCoords)
 	    {
 	    	set_prtIndex(elem);
	    	return update(elem, vCornerCoords);
	    }

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
    
    // new methods
        void update_inner_boundary_radial(const MathVector<TWorldDim>* vCornerCoords);
        void update_inner_boundary_radial_for2();
        void update_inner_boundary_radial_for2_StdFV();
        void update_inner_boundary_radial_old();
        void update_inner_boundary_radial_StdFV(const MathVector<TWorldDim>* vCornerCoords);

    //////////////////////////////////////////////////////////
    ///  original class methods + helper methods for 'update()'
    //////////////////////////////////////////////////////////
    
    /// collects all corners of the flat top element
        int CollectCorners_FlatTop_3d(GridObject* elem);
        int get_cutMode(std::vector<Vertex*> vVertex);
        int CollectCorners_FlatTop_Prism3(GridObject* elem);
        int CollectCorners_FlatTop_Prism4(GridObject* elem);
        int CollectCorners_FlatTop_Pyramid(GridObject* elem);
        int CollectCorners_FlatTop_originalTet(GridObject* elem);
    
    /// !!! nicht in Base!
    /*void ResortQuadrilateral_for2(std::vector<MathVector<dim> > vQuadriCorners);
     
        int CollectCorners_StdFV_for2(GridObject* elem);
        int CollectCorners_FlatTop_2d_for2(GridObject* elem);
        int CollectTriangle_and_Quadri_for2(GridObject* elem, Vertex* vrtInside);
        int CollectQuadrilateral_besideTri_for2(GridObject* elem);
        int CollectQuadrilateral_for2(GridObject* elem);
     
    	void compute_flat_top_data_for2(GridObject* elem);
     
     */
    
	//////////////////////////////////////////////////////////////////////////
	/// methods, which need to be implemented originating from the base class
	//////////////////////////////////////////////////////////////////////////

		bool is_FTVertex(Vertex* vrt)
		{
			if ( m_prtIndex == -1 ) {return m_spCutElementHandler->is_FTVertex(m_prtIndex, vrt);}
			else 					return m_spCutElementHandler->is_FTVertex(vrt); }
		bool is_FTVertex(Vertex* vrt, size_t vrtIndex)
		{ return is_FTVertex(vrt); }


		bool is_FTVertex(int& prtIndex, Vertex* vrt)
		{ return m_spCutElementHandler->is_outsideFluid(prtIndex, vrt);	}

		bool is_OutsideVertex(Vertex* vrt, size_t vrtIndex)
		{ 	if ( m_prtIndex == -1 ) return m_spCutElementHandler->is_FTVertex(m_prtIndex, vrt);
			else					return m_spCutElementHandler->is_OutsideVertex(vrt, vrtIndex); }


		bool is_nearInterfaceVertex(Vertex* vrt, size_t vrtIndex)
 		{ return m_spCutElementHandler->is_nearInterfaceVertex(vrt, vrtIndex); }
		bool is_nearInterfaceVertex(Vertex* vrt)
 		{ return is_nearInterfaceVertex(vrt, 0); }



	///////////////////////////////////////////////////////////////
	/// new methods
	///////////////////////////////////////////////////////////////

    /// for boundary computations
        bool is_boundary_face_for2(const size_t sideID);

    /// acces to 'm_vQuadriCorners_for2'
        const std::vector<MathVector<dim> > quadriCorners() const { return m_vQuadriCorners_for2; }
    
    /// needed during 'CollectQuadrilateral_besideTri_for2()':
        void switch_order();
    
    	void set_prtIndex(GridObject* elem);

   	/// called during 'ParticleMapper::modify_LocalData()':
    	void resize_local_indices(LocalVector& locU);
    	void resize_local_indices(LocalVector& locU, size_t numCo);
    
        SmartPtr<CutElementHandlerFlatTop<dim> > get_cutElementHandler() { return m_spCutElementHandler; }

    /// updates 'm_elemModus'; called by FV1FTGeom::update()!!
        bool update(GridObject* elem, const MathVector<TWorldDim>* vCornerCoords);	// = preprocess() of flat_top.h
    
    /// should be called by 'CutElementHandler::update()' directly => no call via this method!
        ElementModus compute_element_modus(GridObject* elem, const int interfaceOrientation)
        { return m_spCutElementHandler->compute_element_modus(elem, interfaceOrientation); }
    
    /// called by 'InterfaceHandlerLocalParticle::update()'
        ElementModus get_element_modus(GridObject* elem)
        { return m_spCutElementHandler->get_element_modus(elem); }


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

	    number get_density(int prtIndex) { return m_spCutElementHandler->get_density(prtIndex); }
	    number get_density_fluid() { return m_fluidDensity; }
	    number get_kinVisc_fluid() { return m_fluidKinVisc; }
        MathVector<dim> get_center(int prtIndex){ return m_spCutElementHandler->get_center(prtIndex);}

	/// access to single entry of 'm_vRadialAtIP'
	    MathVector<dim> radial_at_ip(size_t i) { return m_vRadialAtIP[i]; }
	/// access to single entry of 'm_vRadialAtCo'
	    MathVector<dim> radial_at_co(size_t i) { return m_vRadialAtCo[i]; }

	    SmartPtr<ParticleProvider<dim> > get_particles() { return m_spCutElementHandler->get_particles(); }

		size_t num_particles() const { return m_spCutElementHandler->num_particles();}

		std::vector<interfaceBF>& get_boundary_faces() { return this->m_vBF; }

		const LocalIndices& get_local_indices() const { return m_ind; }

		int get_prtIndex()
		{ 	if ( m_prtIndex == -1 ) UG_THROW("InterfaceHandlerLocalParticle::apper::get_prtIndex(): prtIndex not set!(): prtIndex not set!\n");
			return m_prtIndex;
		}
 		int getPrtIndex(size_t dof)
		{
			int prtIndex = m_spCutElementHandler->get_prtIndex(dof);
			return prtIndex;
		}

	/// used in mapper:
		number get_rotJ_ind(size_t fct1, size_t dof1, size_t fct2, size_t dof2)
			{ return m_rotJ_ind(fct1, dof1, fct2, dof2); }
		number get_rotJ_rot(size_t fct1, size_t dof1, size_t fct2, size_t dof2)
			{ return m_rotJ_rot(fct1, dof1, fct2, dof2); }

		number get_rotD(size_t fct, size_t dof)
			{ return m_rotD(fct, dof); }

	///////////////////////////////////////////////////////////////
	/// get solution values (used by 'immersedbnd_cond.h'

		MathVector<dim> get_transSol(size_t prtIndex, size_t timeSeriesInd)
			{ return m_spCutElementHandler->get_transSol(prtIndex, timeSeriesInd); }
		MathVector<dim> get_rotSol(size_t prtIndex, size_t timeSeriesInd)
			{ return m_spCutElementHandler->get_rotSol(prtIndex, timeSeriesInd); }

		MathMatrix<dim,dim> get_rotationMat(MathVector<dim> radialVector)
			{ return m_spCutElementHandler->get_rotationMat(radialVector); }


    ///////////////////////////////////////////////////////////////
	/// data for Particle
	///////////////////////////////////////////////////////////////

	   	std::vector<interfaceBF> m_vBF;					// updated during FV1FTGeom::update_inner_boundary_faces()

	/// actual particle index: set during 'm_spCutElementHandler->compute_element_modus()'!!!
		int m_prtIndex;

        std::map<Vertex*, int> m_vPrtIndices;

    /// fuid parameter imported by Constructor()
     	number m_fluidDensity;
     	number m_fluidKinVisc;

    /// local data for assembling:
     	LocalMatrix m_rotJ_ind;
     	LocalMatrix m_rotJ_rot;
     	LocalVector m_rotD;

    /// size of local algebra for flat top element: 'm_numFct' x 'm_numCo'
     	size_t m_numFct;
    
    /// number of corners of flat top element
     	size_t m_numCo;

    /// new local algebra for resized flat top element
     	LocalIndices m_ind;

   		std::vector<MathVector<dim> > m_vRadialAtIP; 	// global position of scvf on interface
   		std::vector<MathVector<dim> > m_vRadialAtCo; 	// global position of scv on interface

    /// flag for computation of 'm_vRadialAtIP' during call of 'update_inner_boundary_faces()'
     	const bool m_bRadial_forMassEq_equals_Normal; 				// default = true

    
    /// flag indiating if local data needs to be updated
        bool m_bBndDataNeeded;							// default = false;
                                                        // flag used e.g. in 'FV1FTGeometry::update()'

	    std::vector<MathVector<dim> > m_vQuadriCorners_for2;

		SmartPtr<CutElementHandlerFlatTop<dim> > m_spCutElementHandler;
 };


}// end namespace ug

#include "interface_handler_particle_tools.h"
#include "interface_handler_particle_impl.h"

#endif /* INTERFACE_HANDLER_FLAT_TOP_H_ */
