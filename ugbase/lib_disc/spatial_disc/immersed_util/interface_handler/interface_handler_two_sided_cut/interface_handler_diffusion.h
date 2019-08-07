/*
 * interface_handler_local.h
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
                                       SmartPtr<CutElementHandlerImmersed<dim> > cutElementHandler);
        virtual ~InterfaceHandlerLocalDiffusion()	{}


    	//////////////////////////////////////////////////////////
    	/// getter methods
    	//////////////////////////////////////////////////////////


    /// access to single entry of 'm_vRealCornerOfigID' (in loc to glob for diffusion)
    	const size_t real_index(size_t i) const { return m_vRealCornerID[i]; }
    	const size_t real_index_size(size_t i, const size_t size) const
    	{ if ( size == 3 ) return m_vRealCornerID_tri[i];
    	  if ( size == 4 ) return m_vRealCornerID_quad[i];
    	  UG_THROW("in 'real_index_size()': wrong size (should be 3 or 4!): " << size << "n");
    	}
    	const size_t real_index_tri(size_t i) const { return m_vRealCornerID_tri[i]; }
    	const size_t real_index_quad(size_t i) const { return m_vRealCornerID_quad[i]; }

    	size_t get_index_for_global_index_Nitsche(const size_t i);
    	size_t get_global_index_Nitsche(const size_t i) { return m_verticesGlobalIndex[i]; }
    	const size_t get_num_NitscheDoFs() { return m_MapInserted_Nitsche.size(); }

		void set_Nitsche(bool bNitsche){ m_bNitsche = bNitsche;}
		bool get_Nitsche(){ return m_bNitsche;}

    
        number get_LSvalue_byPosition(MathVector<dim> vrtPos){};
        number get_LSvalue(Vertex* vrt, const int prtIndex){};

    // ToDo: needed method?
	// used for method 'InterfaceHandlerLocal::get_intersection_point()':
		bool get_intersection_point(MathVector<dim>& Intersect, Vertex* vrtOutsideCirc, Vertex* vrtInsideCirc)
		{
			const MathVector<dim>& vrtPosOut = this->m_aaPos[vrtOutsideCirc];
			const MathVector<dim>& vrtPosIn  = this->m_aaPos[vrtInsideCirc];

			if ( this->m_orientationInterface == 1 )
				return this->m_spInterfaceProvider->get_intersection_point(Intersect, vrtPosOut, vrtPosIn, 0);
		// inverse order of 'vrtPosOut' and 'vrtPosIn' for call of 'get_intersection_point()'
		// to avoid error for alpha < 0:
			else if ( this->m_orientationInterface == -1 )
				return this->m_spInterfaceProvider->get_intersection_point(Intersect, vrtPosIn, vrtPosOut, 0);
			else
				UG_THROW("in InterfaceHandlerLocalDiffusion::get_intersection_point(): m_orientationInterface not set!\n");
		}

		bool get_intersection_point(MathVector<dim>& Intersect, Vertex* vrtOutsideCirc, Vertex* vrtInsideCirc, std::vector<number>& alphaOut)
		{
			const MathVector<dim>& vrtPosOut = this->m_aaPos[vrtOutsideCirc];
			const MathVector<dim>& vrtPosIn  = this->m_aaPos[vrtInsideCirc];

			if ( this->m_orientationInterface == 1 )
				return this->m_spInterfaceProvider->get_intersection_point(Intersect, vrtPosOut, vrtPosIn, 0, alphaOut);
		// inverse order of 'vrtPosOut' and 'vrtPosIn' for call of 'get_intersection_point()'
		// to avoid error for alpha < 0:
			else if ( this->m_orientationInterface == -1 )
				return this->m_spInterfaceProvider->get_intersection_point(Intersect, vrtPosIn, vrtPosOut, 0, alphaOut);
			else
				UG_THROW("in InterfaceHandlerLocalDiffusion::get_intersection_point(): m_orientationInterface not set!\n");
		}

   	///////////////////////////////////////////////////////////////
   	/// base class methods: adjusted!
   	///////////////////////////////////////////////////////////////

	    int CollectCorners_FlatTop_2d(GridObject* elem);
	    void Collect_Data_Nitsche(GridObject* elem);

	  ///	updates all std::vectors-data and especially derives 'm_roid'
	    void compute_flat_top_data(GridObject* elem);

	///////////////////////////////////////////////////////////////
	/// new methods
	///////////////////////////////////////////////////////////////

	// calls base class method, BUT: also clears m_vBF! --> m_vBF NOT a member of base class!
    	bool update_elem(GridObject* elem, const MathVector<TWorldDim>* vCornerCoords);

    /// called by 'InterfaceHandlerLocalParticle::update()'
    	ElementModus compute_element_modus(GridObject* elem, const int interfaceOrientation)
     	{ return m_spCutElementHandler->compute_element_modus(elem, interfaceOrientation); }

    /// called by 'InterfaceHandlerLocalParticle::update()'
    	ElementModus get_element_modus(GridObject* elem)
    	{ return m_spCutElementHandler->get_element_modus(elem); }

		bool is_FTVertex(Vertex* vrt, size_t vrtIndex)
		{	return m_spCutElementHandler->is_FTVertex(vrt, vrtIndex);}

	 	bool is_OutsideVertex(Vertex* vrt, size_t vrtIndex)
		{	return m_spCutElementHandler->is_OutsideVertex(vrt, vrtIndex);}

	     bool is_nearInterfaceVertex(Vertex* vrt, size_t vrtIndex)
		{	return m_spCutElementHandler->is_nearInterfaceVertex(vrt, vrtIndex);}

        bool is_nearInterface(Vertex* vrt)
        {	return m_spCutElementHandler->is_nearInterface(vrt);}
    
   	///////////////////////////////////////////////////////////////
    
        bool check_vertex_modus(VertexModus vrtModus, size_t vrtIndex)
        {	return m_spCutElementHandler->check_vertex_modus(vrtModus, vrtIndex);}
        bool check_vertex_modus(VertexModus vrtModus, size_t vrtIndex, const int interfaceOrientation)
        {	return m_spCutElementHandler->check_vertex_modus(vrtModus, vrtIndex, interfaceOrientation);}
    

   	///////////////////////////////////////////////////////////////


	// used for diffusion
	    size_t get_or_insert_vertex(const MathVector<dim>& vrtPos);
        size_t get_or_insert_vertex(Vertex* vrt);
        bool find_vrtPos(const MathVector<dim>& vrtPos);

	// used for Nitsche
	// ---> used during 'diffusion_interface/diffusion_interface.h:initialize_interface_Nitsche()':
	    size_t get_or_insert_indexPair_Nitsche(const size_t global_index);

	    void Resort_RealCornerID();

	    bool lies_onInterface_tri(const size_t newID);
	    bool lies_onInterface_quad(const size_t newID);
	    bool lies_onInterface_size(const size_t newID, size_t size);

    /// access to m_vBF (needed during 'FV1CutGeom::update_inner_boundary_faces()')
	    const interfaceBF* boundary_faces() const { return &m_vBF[0]; }
	    std::vector<interfaceBF>& get_boundary_faces() { return m_vBF; }

    	const LocalIndices& get_local_indices() const { return m_ind; }

		void set_jac_bool(bool bJac) { m_jac_tag = bJac; }
		bool get_jac_bool() { return m_jac_tag; }

		void set_bScaleDoFs(bool bScaleDoF) { m_scaleDoFs = bScaleDoF; }
		bool get_bScaleDoFs() { return m_scaleDoFs; }

 		const size_t get_numVerticesPos() const { return m_verticesPos.size(); }
 		const MathVector<dim> get_VerticesPos(size_t index) const { return m_verticesPos[index]; }
 		const double get_sol(size_t index) const { return m_verticesValue[index]; }
 		const double get_sol(size_t index, size_t fct) const { return m_verticesValue[index]; }

		void add_to_integral(const number value) { m_integral += value; }
		void init_integral(){ m_integral = 0.0;}
		number get_integral(){ return m_integral;}

 		void print_Nitsche_Data();
 		void print_InterfaceIDdata();

 		number vAlpha(size_t i, size_t j) { return m_vAlpha[i][j]; }
 		MathVector<dim> vIntersectionPnts(size_t i) { return m_vIntersectionPnts[i]; }
 		MathMatrix<dim+1,dim+1> vShapeValues() { return m_vShapeValues; }
 		MathVector<dim> NormalToFace() { return m_NormalToFace; }
 		number Gamma() { return m_Gamma; }
 		number Area() { return m_Area; }
 		number AreaOrig() { return m_AreaOrig; }
 		number AreaScale() { return m_AreaScale; }
 		number IntegralGamma(size_t i) { return m_vIntegralGamma[i]; }

	//////////////////////////////////////////////////////////////////////////////
    /// methods called during 'convection_diffusion_fv1.cpp' -> 'geo.___()' :
    //////////////////////////////////////////////////////////////////////////////

 		void resize_local_data(LocalVector locU);

		void set_DoF_tag_tri(const bool bFactor2_for_DoFIndex)
		{ m_shift_DoFIndex_tri = bFactor2_for_DoFIndex; return; }
		void set_DoF_tag_quad(const bool bFactor2_for_DoFIndex)
		{ m_shift_DoFIndex_quad = bFactor2_for_DoFIndex; return; }

        bool get_vertex_modus_tri(size_t dof, const int orientation)
        {  return m_spCutElementHandler->get_vertex_modus_orient(dof, 1); }
        bool get_vertex_modus_quad(size_t dof, const int orientation)
        {  return m_spCutElementHandler->get_vertex_modus_orient(dof, -1); }
    

		void reset_defect_on_interface(LocalVector& locD, const size_t size);
		void reset_jacobian_on_interface(LocalMatrix& locJ, const size_t size);

		void set_local_sol(LocalVector& solU, const size_t size, const LocalVector& lvec, const int orientation);

    // diffusion can only be set as values, with constant diff coefficient at each side of the interface.
    // no _impl() methods provided for analytic lua-functions
        number get_diffusion();
        number get_diffusion(const bool bElementIsOutside);
    
    // interface data (source, jump, jump_grad) can be defined as values or lua functions
		LocalVector set_jump_values(LocalIndices ind, const size_t size);
		LocalVector set_jump_grad_values(LocalIndices ind, const size_t size);
 		LocalVector set_source(const std::vector<double> sourceIm, LocalIndices ind, const size_t size, const bool bElementIsCut);
    // 'set_source' instance used for navier stokes elem disc
        LocalVector set_source(LocalIndices ind, const size_t size, const bool bElementIsCut)
        {LocalVector dummy; return dummy;}


		void set_jacobian_tri(const LocalMatrix locJ) { m_locJ_tri = locJ; }
		void set_jacobian_quad(const LocalMatrix locJ){ m_locJ_quad = locJ; }

		void set_defect_tri(const LocalVector locD) { m_locD_tri = locD; }
		void set_defect_quad(const LocalVector locD){ m_locD_quad = locD; }

		LocalVector& get_local_solution_tri()  { return m_locU_tri; }
		LocalVector& get_local_solution_quad() { return m_locU_quad; }

	// also called in mapper
		LocalMatrix& get_local_jacobian_tri()  { return m_locJ_tri; }
		LocalMatrix& get_local_jacobian_quad() { return m_locJ_quad; }

		LocalVector& get_local_defect_tri()  { return m_locD_tri; }
		LocalVector& get_local_defect_quad() { return m_locD_quad; }

        const bool get_bNearInterface() const {return m_bNearInterface;}
        void set_bNearInterface(bool bNearInterface) { m_bNearInterface = bNearInterface;}
    
        const bool get_boolian_for_diffusion();
    
   //////////////////////////////////////////////////////////////////////////////
   /// methods called during mapper
   //////////////////////////////////////////////////////////////////////////////

		const size_t get_index_shift_tri() const { return m_shift_DoFIndex_tri; }
		const size_t get_index_shift_quad() const { return m_shift_DoFIndex_quad; }


	//////////////////////////////////////////////////////////////////////////////
    /// methods called ONLY in '_impl.cpp':
    //////////////////////////////////////////////////////////////////////////////

        void set_source_lua(const number interfaceSourceValue) { m_interfaceSource = interfaceSourceValue; }
        void set_jump_lua(const number interfaceJumpValue) { m_interfaceJump = interfaceJumpValue; }
        void set_jump_grad_lua(const MathVector<2>& interfaceJumpGradValue)
            { m_interfaceJumpGrad[0] = interfaceJumpGradValue[0]; m_interfaceJumpGrad[1] = interfaceJumpGradValue[1];}
        void set_diffusion_coeff_lua(const MathVector<2>& diffusionCoeffs)
            { m_diffusionCoeff[0] = diffusionCoeffs[0]; m_diffusionCoeff[1] = diffusionCoeffs[1];}

        bool check_interface_data(const bool bBndFct);
		
        double get_jump(const MathVector<dim> position);
        double get_jump_impl(const MathVector<dim> position);


		double get_jump_grad(const MathVector<dim> position);
        double get_jump_grad_impl(const MathVector<dim> position);

		double get_source(const MathVector<dim> position);
        double get_source_impl(const MathVector<dim> position);

	// writes solution of global vector vec into this->m_verticesValue-array:
		void write_solution(const std::vector<double > verticesValues);

        MathVector<dim> get_center(int prtIndex){ return m_spCutElementHandler->get_center(prtIndex);}
    

   //////////////////////////////////////////////////////////////////////////////
   /// members:
   //////////////////////////////////////////////////////////////////////////////

	    std::vector<interfaceBF> m_vBF;			// updated during FV1CutGeom::update_inner_boundary_faces()

	    std::map<MathVector<dim>, size_t> m_MapInserted;

	    std::vector<MathVector<dim> > m_verticesPos;
	    std::vector<double> m_verticesValue;

	// the REalCornerID can also be an index > 3 and > 4, since it is then the index
	// of the InterfaceVertex; see 'm_verticesPos' and 'm_verticesValue'
	    std::vector<size_t> m_vRealCornerID;
	    std::vector<size_t> m_vRealCornerID_tri;
	    std::vector<size_t> m_vRealCornerID_quad;

	    std::vector<size_t> m_vInterfaceID_tri;
	    std::vector<size_t> m_vInterfaceID_quad;

		bool m_scaleDoFs;
		number m_integral;
		bool m_jac_tag;


	/// size of local algebra for flat top element: 'm_numFct' x 'm_numCo'
		size_t m_numFct;
	/// number of corners of flat top element
		size_t m_numCo;

	/// new local algebra for resized flat top element
		LocalIndices m_ind;

	// local data for assembling:
		LocalMatrix m_locJ_tri;
		LocalMatrix m_locJ_quad;
		LocalVector m_locD_tri;
		LocalVector m_locD_quad;
		LocalVector m_locU_tri;
		LocalVector m_locU_quad;

        bool m_bNearInterface;
    
	// scale factor for access to DoFIndex on triangle or quadri as cut element
	// --> for call during 'add_local_def/jac_to_global_interface()':
		bool m_shift_DoFIndex_tri;
		bool m_shift_DoFIndex_quad;

     ///////////////////////////////////////////////////////////////
	 /// data for Nitsche -> 'Collect_Data_Nitsche()':
	 ///////////////////////////////////////////////////////////////

		bool m_bNitsche;
	    std::map<size_t, size_t> m_MapInserted_Nitsche;
		std::vector<std::vector<number> > m_vAlpha;
		std::vector<MathVector<dim> > m_vIntersectionPnts;
		MathMatrix<dim+1,dim+1> m_vShapeValues;  // dim+1 = number of vertices for simplicial mesh!
		MathVector<dim> m_NormalToFace;
		number m_Gamma;
		MathVector<dim> m_insidePnt;
		number m_Area;
		number m_AreaOrig;
		number m_AreaScale;
		std::vector<number> m_vIntegralGamma;

 	  // filled during 'moving_interface/interface_handler_local_base_tools.h: get_or_insert_indexPair_Nitsche()'
	  // AND get_or_insert_indexPair_Nitsche() called during 'diffusion_interface/diffusion_interface.h:initialize_interface_Nitsche()'
	  // and used for local-to-global mapper
		std::vector<size_t> m_verticesGlobalIndex;


    // boolians for boundary data
        bool m_bBndFct;
        number m_interfaceSource;
        number m_interfaceJump;
        MathVector<2> m_interfaceJumpGrad;
        MathVector<2> m_diffusionCoeff;

    
	/// contains radius, center and density of all given particles
		SmartPtr<DiffusionInterfaceProvider<dim> > m_spInterfaceProvider;
		SmartPtr<CutElementHandlerImmersed<dim> > m_spCutElementHandler;

};

    
}// end namespace ug


#include "interface_handler_diffusion_tools.h"
#include "interface_handler_diffusion_impl.h"

#endif /* INTERFACE_HANDLER_LOCAL_DIFFUSION_H_ */
