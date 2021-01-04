/*
 * interface_handler_diffusion.h
 *
 *  Created on: 15.01.2015
 *      Author: suze
 */

#ifndef INTERFACE_HANDLER_LOCAL_DIFFUSION_H_
#define INTERFACE_HANDLER_LOCAL_DIFFUSION_H_


#include "lib_grid/multi_grid.h"

#include "../interface_handler_base.h"


namespace ug{
    
template <int TDim, int TWorldDim, class TInterfaceHandler>
class DimFV1CutGeometry;
    
template <int TWorldDim>
class InterfaceHandlerLocalDiffusion : public InterfaceHandlerLocalBase<TWorldDim>
{
//    typedef InterfaceHandlerLocalBase<TWorldDim> base_class;

	public:
	/// world dimension
		static const int dim = TWorldDim;

	///	max number of geometric objects in a dimension
	// 	(most objects in 1 dim, i.e. number of edges, but +1 for 1D)
        static const int maxMid = DimFV1CutGeometry<dim, dim, InterfaceHandlerLocalDiffusion<dim> >::maxNumSCVF + 1;
 
	///	used traits
		typedef fv1_dim_traits<dim, dim> traits;

	/// used boundary face type
        typedef typename DimFV1CutGeometry<dim, dim, InterfaceHandlerLocalDiffusion<dim> >::BF interfaceBF;

	///	Type of position coordinates
		typedef MathVector<dim> position_type;

	///	Type of Position Attachment
		typedef Attachment<position_type> position_attachment_type;

	///	Type of Accessor to the Position Data Attachment
	 	 typedef Grid::VertexAttachmentAccessor<position_attachment_type> position_accessor_type;


	 	InterfaceHandlerLocalDiffusion(SmartPtr<DiffusionInterfaceProvider<dim> > interfaceProvider,
                                       SmartPtr<CutElementHandler_TwoSided<dim> > cutElementHandler);
    
        virtual ~InterfaceHandlerLocalDiffusion()	{}
    
    //////////////////////////////////////////////////////////////////
    /// virtual methods in 'IInterfaceHandlerLocal'
    //////////////////////////////////////////////////////////////////
 
    // simple forwarding to the associated 'CutElementHandler' class:

        bool is_onInterfaceVertex(Vertex* vrt, size_t vrtIndex)
            {	return m_spCutElementHandler->is_onInterfaceVertex(vrt, vrtIndex);}
    
        bool is_OutsideVertex(Vertex* vrt, size_t vrtIndex)
            {	return m_spCutElementHandler->is_OutsideVertex(vrt, vrtIndex);}
    
        bool is_nearInterfaceVertex(Vertex* vrt, size_t vrtIndex)
            {	return m_spCutElementHandler->is_nearInterfaceVertex(vrt, vrtIndex);}
    
        bool is_nearInterface(Vertex* vrt)
            {	return m_spCutElementHandler->is_nearInterface(vrt);}
    
    //////////////////////////////////////////////////////////////////
    /// virtual methods of 'InterfaceHandlerLocalBase'
    //////////////////////////////////////////////////////////////////
 
    // forwarding to m_spCutElementHandler:
        number get_LSvalue_byPosition(MathVector<dim> vrtPos)
            { return m_spCutElementHandler->get_LSvalue_byPosition(vrtPos); }
    
        bool get_intersection_point(MathVector<dim>& Intersect, Vertex* vrtOutsideCirc,
                                    Vertex* vrtInsideCirc);
        bool get_intersection_point(MathVector<dim>& Intersect, Vertex* vrtOutsideCirc,
                                    Vertex* vrtInsideCirc, std::vector<number>& alphaOut);

    
    ///	recomputes the 'vCornerCoords' of the cut element and derives 'm_roid'
        void compute_cut_element_data(GridObject* elem);

        void compute_and_set_prtIndex(GridObject* elem)
            { m_spCutElementHandler->compute_and_set_prtIndex(elem); }
    
    /// central method called by 'compute_cut_element_data()' to compute 'vCornerCoords'
        int CollectCorners_FlatTop_2d(GridObject* elem);
    
    /// helper method within 'CollectCorners_FlatTop_2d' to bring the newly computed
    //      'vCornerCoords' in correct order (counter clockwise)
        void Resort_RealCornerID();
    
        void Collect_Data_Nitsche(GridObject* elem);
    

    //////////////////////////////////////////////////////////////////////////////
    /// initialize interface boundary conditions via lua-call
    //////////////////////////////////////////////////////////////////////////////
    
        void set_source_data_lua(const number interfaceSourceValue)
            { m_interfaceSource = interfaceSourceValue; m_luaSource_isSet = true; }
        void set_jump_data_lua(const number interfaceJumpValue)
            { m_interfaceJump = interfaceJumpValue; m_luaJump_isSet = true; }
        void set_jump_grad_data_lua(const MathVector<2>& interfaceJumpGradValue)
            { m_interfaceJumpGrad[0] = interfaceJumpGradValue[0];
              m_interfaceJumpGrad[1] = interfaceJumpGradValue[1]; m_luaJumpGrad_isSet = true; }
        void set_diffusion_coeff_data_lua(const MathVector<2>& diffusionCoeffs)
            { m_diffusionCoeff[0] = diffusionCoeffs[0]; m_diffusionCoeff[1] = diffusionCoeffs[1];
              m_luaDiffusion_isSet = true; }
    
    
    //////////////////////////////////////////////////////////
    /// getter methods
    //////////////////////////////////////////////////////////


    /// access to m_vBF (needed during 'FV1CutGeom::update_inner_boundary_faces()')
         std::vector<interfaceBF>& get_boundary_faces() { return m_vBF; }
    /// called during 'update()'-method of Fv1CutGeom
        void clear_boundary_faces() { m_vBF.clear(); }
 
        const LocalIndices& get_local_indices() const { return m_ind; }

        MathVector<dim> get_center(int prtIndex){ return m_spCutElementHandler->get_center(prtIndex);}

    // methods for writing and getting GLOBAL data of interface nodes
        size_t get_or_insert_vertex(const MathVector<dim>& vrtPos);
        size_t get_or_insert_vertex(Vertex* vrt);
        bool find_vrtPos(const MathVector<dim>& vrtPos);

    	size_t get_index_for_global_index_Nitsche(const size_t i);
    	size_t get_global_index_Nitsche(const size_t i) { return m_verticesGlobalIndex[i]; }
    	const size_t get_num_NitscheDoFs() { return m_MapInserted_Nitsche.size(); }

    // used for Nitsche
		void set_Nitsche(bool bNitsche){ m_bNitsche = bNitsche;}
		bool get_Nitsche(){ return m_bNitsche;}
    
	// ---> used during 'diffusion_interface/diffusion_interface.h:initialize_interface_Nitsche()':
	    size_t get_or_insert_indexPair_Nitsche(const size_t global_index);

    
	////////////////////////////////////////////////////////////////////////////////////
    /// methods called by the Elem Disc 'ConvectionDiffusionFV1_cutElem'
    ////////////////////////////////////////////////////////////////////////////////////

    /////////////////////////////////////////////////////////////////////////////////////////////
    // (A) 'integral'-methods: for the computation of the l2-error within
    // 'ConvectionDiffusionFV1_cutElem::add_l2error_A_elem()'
    
        void L2Error_add(const number value) { m_L2Error += value; }
        void L2Error_init(){ m_L2Error = 0.0;}
        number get_L2Error(){ return m_L2Error;}
    
    /////////////////////////////////////////////////////////////////////////////////////////////
    /// (B) further methods
    
 		void resize_local_data(LocalVector locU);

        const bool get_boolian_for_diffusion();

        void set_DoF_tag(const bool bFactor2_for_DoFIndex, const ReferenceObjectID roid)
        {   if ( roid == ROID_TRIANGLE )           {m_shift_DoFIndex_tri  = bFactor2_for_DoFIndex; return;}
            else if ( roid == ROID_QUADRILATERAL ) {m_shift_DoFIndex_quad = bFactor2_for_DoFIndex; return;}
            else {UG_THROW("in InterfaceHandlerLocalDiffusion::set_DoF_tag_tri(): invalid roid!\n");}
        }

  		void set_bScaleDoFs(bool bScaleDoF) { m_scaleDoFs = bScaleDoF; }
        bool get_bScaleDoFs() { return m_scaleDoFs; }

        const bool get_bNearInterface() const {return m_bNearInterface;}
        void set_bNearInterface(bool bNearInterface) { m_bNearInterface = bNearInterface;}
 
    /////////////////////////////////////////////////////////////////////////////////////////////
    // (C) access methods for assemgling of the local defect and jacobian for
    // the elem Disc 'ConvectionDiffusionFV1_cutElem':
    //  --> during call of 'add_jac_A_elem()' and 'add_def_A_elem()'
    
        void reset_defect_on_interface(LocalVector& locD, const size_t size);
        void reset_jacobian_on_interface(LocalMatrix& locJ, const size_t size);

        void set_jacobian(const LocalMatrix locJ, const ReferenceObjectID roid)
        {   if ( roid == ROID_TRIANGLE )           m_locJ_tri = locJ;
            else if ( roid == ROID_QUADRILATERAL ) m_locJ_quad = locJ;
            else {UG_THROW("in InterfaceHandlerLocalDiffusion::set_jacobian(): invalid roid!\n");}
        }
    
        void set_defect(const LocalVector locD, const ReferenceObjectID roid)
        {   if ( roid == ROID_TRIANGLE )           m_locD_tri = locD;
            else if ( roid == ROID_QUADRILATERAL ) m_locD_quad = locD;
            else {UG_THROW("in InterfaceHandlerLocalDiffusion::set_defect(): invalid roid!\n");}
        }
    
    
        const size_t get_numVerticesPos() const                   { return m_verticesPos.size(); }
        const MathVector<dim> get_VerticesPos(size_t index) const { return m_verticesPos[index]; }
        const double get_sol(size_t index) const                  { return m_verticesValue[index]; }
        const double get_sol(size_t index, size_t fct) const      { return m_verticesValue[index]; }

    // also called by class 'DiffusionInterfaceMapper'
        LocalMatrix& get_local_jacobian(const ReferenceObjectID roid)
        {   if ( roid == ROID_TRIANGLE ) return m_locJ_tri;
            else if ( roid == ROID_QUADRILATERAL ) return m_locJ_quad;
            else {UG_THROW("in InterfaceHandlerLocalDiffusion::get_local_jacobian(): invalid roid!\n");}
        }
    
        LocalVector& get_local_defect(const ReferenceObjectID roid)
        {   if ( roid == ROID_TRIANGLE ) return m_locD_tri;
            else if ( roid == ROID_QUADRILATERAL ) return m_locD_quad;
            else {UG_THROW("in InterfaceHandlerLocalDiffusion::get_local_defect(): invalid roid!\n");}
        }

        LocalVector& get_local_solution(const ReferenceObjectID roid)
        {   if ( roid == ROID_TRIANGLE ) return m_locU_tri;
            else if ( roid == ROID_QUADRILATERAL ) return m_locU_quad;
            else {UG_THROW("in InterfaceHandlerLocalDiffusion::get_local_solution(): invalid roid!\n");}
        }
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////
    // (D) access methods for local assembling of the boundary conditions on the
    // immersed interface into the defect and jacobian for the elem Disc
    // 'ConvectionDiffusionFV1_cutElem'
    //      --> during call of 'add_jac_A_elem()' and 'add_def_A_elem()'

    // diffusion can only be set as values, with constant diff coefficient at each side of the interface.
    // no _impl() methods provided for analytic lua-functions
        number get_diffusion();
        number get_diffusion(const bool bElementIsOutside);
    
        double get_jump(const MathVector<dim> position);
        double get_jump_impl(const MathVector<dim> position);
    
        double get_jump_grad(const MathVector<dim> position);
        double get_jump_grad_impl(const MathVector<dim> position);
    
        double get_source(const MathVector<dim> position);
        double get_source_impl(const MathVector<dim> position);

    // setter methods (called by ConvectionDiffusionFV1_cutElem::get_local_data() )
		void set_local_sol(LocalVector& solU, const size_t size, const LocalVector& lvec,
                           const int orientation);

        void set_jump_values     (LocalVector& jumpOut, LocalIndices ind, const size_t size);
        void set_jump_grad_values(LocalVector& jumpGradOut, LocalIndices ind, const size_t size);
        void set_source(const std::vector<double> sourceIm, LocalVector& sourceOut, LocalIndices ind,
                        const size_t size, const bool bElementIsCut);
    
        bool lies_onInterface_tri(const size_t newID);
        bool lies_onInterface_quad(const size_t newID);
        bool lies_onInterface_size(const size_t newID, size_t size);

    
   //////////////////////////////////////////////////////////////////////////////
   /// methods called by class 'DiffusionInterfaceMapper'
   //////////////////////////////////////////////////////////////////////////////

		const size_t get_index_shift_tri() const { return m_shift_DoFIndex_tri; }
		const size_t get_index_shift_quad() const { return m_shift_DoFIndex_quad; }

    
    // LOCAL access to 'VertexModus' via 'm_vvVertexMode' data in 'CutElementHandler' class
        bool check_vertex_modus(VertexModus vrtModus, size_t vrtIndex, const int interfaceOrientation)
            {	return m_spCutElementHandler->check_vertex_modus(vrtModus, vrtIndex, interfaceOrientation);}

    // writes solution of the interface nodes (i.e. of the new DoFs) from
    // the global vector 'vec' into this->m_verticesValue-array:
        void set_interface_values(const std::vector<double > verticesValues);
    
   
    //////////////////////////////////////////////////////////////////////////////
    /// access to entries of the GLOBAL index of interface nodes 'm_vRealCornerID'
    //////////////////////////////////////////////////////////////////////////////
    
    // get the 'real index':
    // case1: node lies on interface => returns the entry location within
    //          'InterfaceHandlerDiffusion::m_MapNearVertices'
    // case2: node lies on an original mesh node: => returns the usual,
    //          local index of vertex
    
    /// called by own class 'InterfaceHandlerLocalDiffusion'
        const size_t real_index(size_t i) const { return m_vRealCornerID[i]; }
        const size_t real_index_size(size_t i, const size_t size) const
            {   if      ( size == 3 ) return m_vRealCornerID_tri[i];
                else if ( size == 4 ) return m_vRealCornerID_quad[i];
                else    UG_THROW("in 'real_index_size()': wrong size (should be 3 or 4!): " << size << "n");
            }
    
    // called by class 'DiffusionInterfaceMapper':
        const size_t real_index_tri(size_t i) const { return m_vRealCornerID_tri[i]; }
        const size_t real_index_quad(size_t i) const { return m_vRealCornerID_quad[i]; }

    
    //////////////////////////////////////////////////////////
    /// output
    //////////////////////////////////////////////////////////
    
        void print_Nitsche_Data();
        void print_CutElementData();
    
    //////////////////////////////////////////////////////////
    /// acces methods for Nitsche
    //////////////////////////////////////////////////////////
    
        number vAlpha(size_t i, size_t j)           { return m_vAlpha[i][j]; }
        MathVector<dim> vIntersectionPnts(size_t i) { return m_vIntersectionPnts[i]; }
        MathMatrix<dim+1,dim+1> vShapeValues()      { return m_vShapeValues; }
        MathVector<dim> NormalToFace()              { return m_NormalToFace; }
        number Gamma()                              { return m_Gamma; }
        number Area()                               { return m_Area; }
        number AreaOrig()                           { return m_AreaOrig; }
        number AreaScale()                          { return m_AreaScale; }
        number IntegralGamma(size_t i)              { return m_vIntegralGamma[i]; }


   //////////////////////////////////////////////////////////////////////////////
   /// class member
   //////////////////////////////////////////////////////////////////////////////

    // m_vBF stores the boundary faces of the immersed boundary;
    //      --> updated during FV1CutGeom::update_inner_boundary_faces()
	    std::vector<interfaceBF> m_vBF;

    // storage of new vertices at the immersed interface
	    std::map<MathVector<dim>, size_t> m_MapInserted;    // <position of interface node, index in map>
	    std::vector<MathVector<dim> > m_verticesPos;
	    std::vector<double> m_verticesValue;

	// RealCornerID: for access do solution via DoFRef
    // 2 cases: (1) corner lies on interface or (2) corner lies on original grid node
    // case1: contains the entry counter of the node within
    //          'InterfaceHandlerDiffusion::m_MapNearVertices' (usually > 4)
    // case2:contains the usual, 'local' index of vertex
	    std::vector<size_t> m_vRealCornerID;
	    std::vector<size_t> m_vRealCornerID_tri;
	    std::vector<size_t> m_vRealCornerID_quad;

    // m_vInterfaceID_tri/_quad := LOCAL index of interface nodes (< 4)
	    std::vector<size_t> m_vInterfaceID_tri;
	    std::vector<size_t> m_vInterfaceID_quad;

	/// size of local algebra for cut element: 'm_numFct' x 'm_numCo'
		size_t m_numFct;
	/// number of corners of cut element
		size_t m_numCo;

	/// new local algebra for resized cut element
		LocalIndices m_ind;

	// local data for assembling on element being cut into triangle and quadrilateral:
		LocalMatrix m_locJ_tri;
		LocalMatrix m_locJ_quad;
		LocalVector m_locD_tri;
		LocalVector m_locD_quad;
		LocalVector m_locU_tri;
		LocalVector m_locU_quad;

        bool m_scaleDoFs;           // if m_scaleDoFs = true: 2 new DoFs will be placed for each interface node (for jump in value)
        number m_L2Error;
        bool m_bNearInterface;
    
	// scale factor for access to DoFIndex on triangle or quadri as cut element
	// --> for call during 'add_local_def/jac_to_global_interface()':
		bool m_shift_DoFIndex_tri;
		bool m_shift_DoFIndex_quad;
    
    ///////////////////////////////////////////////////////////////
    /// data for assembling interface boundary condition
    
    // boolians for boundary condition data:
    // if true: constant values are given via lua-call
    // if false: position-dependent values are given via inline-function
        bool m_luaSource_isSet;
        bool m_luaJump_isSet;
        bool m_luaJumpGrad_isSet;
        bool m_luaDiffusion_isSet;
    
    // if boolians above are true: constant values are given via lua-call and stored here
        number m_interfaceSource;
        number m_interfaceJump;
        MathVector<2> m_interfaceJumpGrad;
        MathVector<2> m_diffusionCoeff;

    ///////////////////////////////////////////////////////////////
    /// data for Nitsche -> 'Collect_Data_Nitsche()':

		bool m_bNitsche;
	    std::map<size_t, size_t> m_MapInserted_Nitsche;
		std::vector<std::vector<number> > m_vAlpha;
		std::vector<MathVector<dim> > m_vIntersectionPnts;
		MathMatrix<dim+1,dim+1> m_vShapeValues;     // dim+1 = number of vertices for simplicial mesh!
		MathVector<dim> m_NormalToFace;
		number m_Gamma;
		MathVector<dim> m_insidePnt;
		number m_Area;
		number m_AreaOrig;
		number m_AreaScale;
		std::vector<number> m_vIntegralGamma;

    // filled during 'immersed_interface/interface_handler_local_base_tools.h: get_or_insert_indexPair_Nitsche()'
    // AND get_or_insert_indexPair_Nitsche() called during 'diffusion_interface/diffusion_interface.h:initialize_interface_Nitsche()'
    // and used for local-to-global mapper
		std::vector<size_t> m_verticesGlobalIndex;

    ///////////////////////////////////////////////////////////////
    /// essential members for further acces
    
	/// contains radius, center and orientation of all given particles
		SmartPtr<DiffusionInterfaceProvider<dim> > m_spInterfaceProvider;
    
    /// computes and contains 'ElementModus', 'VertexModus' and handles access
		SmartPtr<CutElementHandler_TwoSided<dim> > m_spCutElementHandler;

};

    
}// end namespace ug


#include "interface_handler_diffusion_tools.h"
#include "interface_handler_diffusion_impl.h"

#endif /* INTERFACE_HANDLER_LOCAL_DIFFUSION_H_ */
