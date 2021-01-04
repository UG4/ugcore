/*
 * interface_handler_local.h
 *
 *  Created on: 15.01.2015
 *      Author: suze
 */

#ifndef INTERFACE_HANDLER_BASE_H_
#define INTERFACE_HANDLER_BASE_H_

#include "lib_grid/multi_grid.h"
#include "../cut_element_handler/cut_element_handler.h"
//#include "../cut_element_handler/cut_element_handler_flat_top.h"
//#include "../cut_element_handler/cut_element_handler_two_sided.h"
#include "../interface_provider/interface_provider_particle.h"

namespace ug{

template <int TDim, int TWorldDim, class TInterfaceHandler>
class DimFV1CutGeometry;
    
    
////////////////////////////////////////////////////////////////////////
//

class IInterfaceHandlerLocal
{

	public:
   	IInterfaceHandlerLocal() : m_elemModus(INSIDE_DOM){};

   	virtual ~IInterfaceHandlerLocal(){}

    /// access to element modus
 	ElementModus elementModus(){ return m_elemModus; }

    ////////////////////////////////////////////////////////////////////////
 	/// virtual class contains all methods which need to be implemented by
 	/// the class, derived from 'InterfaceHandlerLocalBase'
    ////////////////////////////////////////////////////////////////////////

    /// computes the element modus and writes it to 'm_elemModus'
 	virtual ElementModus get_element_modus(GridObject* elem)
 	{ 	UG_THROW("in 'IInterfaceHandlerLocal': Missing virtual method: get_element_modus()."); }

 	virtual ElementModus compute_element_modus(GridObject* elem)
 	{ 	UG_THROW("in 'IInterfaceHandlerLocal': Missing virtual method: compute_element_modus()."); }

 	virtual bool is_onInterfaceVertex(Vertex* vrt, size_t vrtIndex)
 	{ 	UG_THROW("in 'IInterfaceHandlerLocal': Missing virtual method: is_onInterfaceVertex()."); }

 	virtual bool is_OutsideVertex(Vertex* vrt, size_t vrtIndex)
 	{   UG_THROW("in 'IInterfaceHandlerLocal': Missing virtual method: is_OutsideVertex()."); }

 	virtual bool is_nearInterfaceVertex(Vertex* vrt, size_t vrtIndex)
 	{   UG_THROW("in 'IInterfaceHandlerLocal': Missing virtual method: is_nearInterfaceVertex()."); }


 protected:
 	/// is actual element inside or outside domain or cut by the interface?
 	ElementModus m_elemModus;

};

template <int TWorldDim>
class InterfaceHandlerLocalBase : public IInterfaceHandlerLocal
{

public:
/// world dimension
	static const int dim = TWorldDim;

///	max number of geometric objects in a dimension
// 	(most objects in 1 dim, i.e. number of edges, but +1 for 1D)
	static const int maxMid = 2*dim + 1;

///	used traits
	typedef fv1_dim_traits<dim, dim> traits;

///	Type of position coordinates
	typedef MathVector<dim> position_type;

///	Type of Position Attachment
	typedef Attachment<position_type> position_attachment_type;

///	Type of Accessor to the Position Data Attachment
 	 typedef Grid::VertexAttachmentAccessor<position_attachment_type> position_accessor_type;


	 InterfaceHandlerLocalBase(SmartPtr<CutElementHandler_TwoSided<dim> > cutElementHandler);
	 InterfaceHandlerLocalBase(SmartPtr<CutElementHandler_FlatTop<dim> > cutElementHandler);

    virtual ~InterfaceHandlerLocalBase(){}

    ///////////////////////////////////////////////////////////////////////////
    /// virtual base class methods:
    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////
    /// central method: update_elem():
    ///     it is called during TFVGeom:update() and computes the new cornders
    ///     i.e. 'vCornerCoords' of the cut element. Based on that ALL local
    ///     TFVGeom-computations (compute ip's, normals, gradients, ...) follow
    ///     as for the standard case
    
    // central method called by TFVGeom:update()
        virtual bool update_elem(GridObject* elem, const MathVector<TWorldDim>* vCornerCoords);
    
    ///	recomputes the 'vCornerCoords' of the cut element and derives 'm_roid'
        virtual void compute_cut_element_data(GridObject* elem);
    
    /// sets data needed for usual computations in case of INSIDE_DOM: 'm_vCornerCoords' and 'm_roid'
        virtual void set_element_data(GridObject* elem, const MathVector<TWorldDim>* vCornerCoords,
                                       ReferenceObjectID roid);
    /// computes the 'ElementModus' based on the corners of 'elem' and also sets some BoolMarker
    //  and cut element data ('VertexModus') in 'CutElementHander' storage
        virtual ElementModus compute_element_modus(GridObject* elem)
            { return m_spCutElementHandler->compute_element_modus(elem); }
    
        virtual ElementModus get_element_modus(GridObject* elem)
            { return m_spCutElementHandler->get_element_modus(elem); }

    /// updates inner boundary data (called during FV1CutGeom::update())
        virtual void update_inner_boundary(const MathVector<TWorldDim>* vCornerCoords){};
        virtual void update_inner_boundary_for2(){};

    
    
    ////////////////////////////////////////////////////////////////
    ///  helper methods for 'compute_cut_element_data()': _tools.h
    ///  ALL methods perform:
    ///     (1) computation of the new corners of a cut element
    ///             by computing the intersection with an edge
    ///     (2) adding the positions to 'vCornerCoords' data of
    ///         the finite element/volume geometry data storage
    ////////////////////////////////////////////////////////////////
    
    /// collects new corners of the cut element
        int CollectCorners_StdFV(GridObject* elem);
        int CollectCorners_FlatTop_2d(GridObject* elem);
        int CollectCorners_FlatTop_3d(GridObject* elem);
        int get_cutMode(std::vector<Vertex*> vVertex);
        int CollectCorners_FlatTop_Prism3(GridObject* elem);
        int CollectCorners_FlatTop_Prism4(GridObject* elem);
        int CollectCorners_FlatTop_Pyramid(GridObject* elem);
        int CollectCorners_FlatTop_originalTet(GridObject* elem);
    
    /// helper methods called during 'CollectCorners_FlatTop_2d' to bring the newly computed
    //      'vCornerCoords' in correct order (counter clockwise)
        void ResortQuadrilateral(std::vector<std::pair<MathVector<dim>, size_t > > vInsideCorners,
                                 std::vector<std::pair<MathVector<dim>, size_t > > vOutsideCorners,
                                 MathVector<dim> normalDir);
        bool isCCW(std::vector<MathVector<dim> > vCornerCoords, MathVector<dim> normal);

    // newly computed intersectino points get collected during 'CollectCorners_...()' only,
    //  if they were not collected before
        bool isIncluded(std::vector<MathVector<dim> > vCheckList, MathVector<dim> checkPoint);
    
    //////////////////////////////////////////////////////////////////////
    ///  central methods to compute the intersectin with an edge
    //////////////////////////////////////////////////////////////////////
  
    // computes the 'level-set' value of a point w.r.t. the interface, in order to
    //  derive its 'VertexModus', i.e. INSIDE, OUTSIDE, ON_INTERFACE
        virtual number get_LSvalue_byPosition(MathVector<dim> vrtPos) = 0;
    
        number get_LSvalue(Vertex* vrt)
            { return get_LSvalue_byPosition(m_aaPos[vrt]); }
    
        virtual bool get_intersection_point(MathVector<dim>& Intersect, Vertex* vrtOutsideCirc,
                                        Vertex* vrtInsideCirc) = 0;
        virtual bool get_intersection_point(MathVector<dim>& Intersect, Vertex* vrtOutsideCirc,
                                        Vertex* vrtInsideCirc, std::vector<number>& alphaOut)= 0;

    
    //////////////////////////////////////////////////////////
    ///  further helper methods in _tools.h
    //////////////////////////////////////////////////////////
    
    /// for boundary computations
        bool lies_onInterface(const size_t newID);
    
        bool remapped_fromInterface(const size_t origID);
        bool remapped_fromInterface(const size_t origID, size_t& get_interfaceID);
        bool is_boundary_face(const size_t sideID);
    
        size_t get_vertex_index(Vertex* vrt, GridObject* elem);

    //////////////////////////////////////////////////////////
    /// access methods
    //////////////////////////////////////////////////////////

	/// access to roid
	    ReferenceObjectID roid(){ return m_roid; }

   /// access to 'm_vCornerCoords'
	    const std::vector<MathVector<dim> > corners() const { return m_vCornerCoords; }
	    const MathVector<dim>* pCorners() const { return &m_vCornerCoords[0]; }
    
     /// access to single entry of 'm_vCornerCoords'
	    const MathVector<dim> corner(size_t i) { return m_vCornerCoords[i]; }
    
	/// access to 'm_vInterfaceID'
	    std::vector<size_t> interface_id_all() { return m_vInterfaceID; }
    
	/// access to single entry of 'm_vInterfaceID'
	    size_t interface_id(size_t i) { return m_vInterfaceID[i]; }
    
	/// access to single entry of 'm_vOriginalCornerID'
	    const size_t corner_orig(size_t i) const { return m_vOriginalCornerID[i]; }

	/// access to single entry of 'm_vNOInterfaceID'
	    size_t NOinterface_id(size_t i) { return m_vNOInterfaceID[i]; }

	    size_t numCo() { return m_vCornerCoords.size(); }

    
    //////////////////////////////////////////////////////////
    /// setter methods
    //////////////////////////////////////////////////////////
    
    /// writes member 'm_roid'
        void set_roid_2d();
        void set_roid_3d();
    
    /// called during modify_LocalData()->set_QuadriSol() for case CUT_BY_2_INTERFACE
        void set_QuadriSol(LocalVector locD, LocalVector locU)
        { m_quadriLocD = locD; m_quadriLocU = locU;}
    
	    void set_StdFV_assembling(bool bValue) { m_bUseStdFVAssembling = bValue; }
	    bool StdFV_assembling() { return m_bUseStdFVAssembling; }
  
        void set_print_cutElemData(bool bValue) { m_bPrintCutElemDataToFile = bValue; }
        bool print_cutElment_data() { return m_bPrintCutElemDataToFile; }

    // called by FV1Geom during geo.update() of ElemDisc:
        void set_orientation(const int orientation) { m_spCutElementHandler->set_orientation(orientation); }
        int get_orientation()                       { return m_spCutElementHandler->get_orientation(); }

    
    //////////////////////////////////////////////////////////
    /// output
    //////////////////////////////////////////////////////////

        void print_CutElementData();

    
	///////////////////////////////////////////////////////////////
    /// base members
	///////////////////////////////////////////////////////////////

	//private:
	    SmartPtr<MultiGrid> m_spMG;
	    position_attachment_type m_aPos;		///<Position Attachment
	    position_accessor_type m_aaPos;			///<Accessor


	/// new reference object id != grid objectreference_object_id()!!
	    ReferenceObjectID m_roid;

	/// all corners of cut element (old + new = intersection of interface with element edgees)
	    std::vector<MathVector<dim> > m_vCornerCoords;

	/// local indices of interface corners
    /// (in case of a remapping onto the original grid nodes, it is the local
    ///  index of the node accross the interface, lying on the same edge as
    ///     the interface node
	    std::vector<size_t> m_vInterfaceID;

    /// local indices of all not-interface corners
	    std::vector<size_t> m_vNOInterfaceID;

	/// maps the new corner indices to the original indices:
	/// 	-> vOriginalCornerID[i] = j => the j-th corner of the original
 	///			element is associated with the i-th corner of the new element
	    std::vector<size_t> m_vOriginalCornerID;


	/// 4 corners forming the quadrilateral of an element cut by 2 particles
	/// ---> REMARK: made sure by CollectCorner-methods, that corner 1,2 belong to prtIndex 0
	/// 				and corner 3,4 belong to prtIndex 1
	    std::vector<size_t> m_vQuadriOrigID;
	    LocalVector m_quadriLocD;
	    LocalVector m_quadriLocU;

    /// for a circular interface and orientation = 1, the
    ///     inner of the circle belongs to the outside domain
	    int m_orientationInterface;                         // default = 1


    /// flag for call of CollectCorners_StdFV()
	    bool m_bUseStdFVAssembling;							// default = false

    /// flag for print of cut-element data into file
        bool m_bPrintCutElemDataToFile;                     // default = false
    
	/// contains radius, center and density of all given particles
	    SmartPtr<CutElementHandlerBase<dim> > m_spCutElementHandler;

};



}// end namespace ug

#include "interface_handler_base_tools.h"
#include "interface_handler_base_impl.h"

#endif /* INTERFACE_HANDLER_BASE_H_ */
