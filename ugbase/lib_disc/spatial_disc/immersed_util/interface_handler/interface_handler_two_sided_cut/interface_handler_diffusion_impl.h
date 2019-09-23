/*
 * interface_handler_local_impl.h
 *
 *  Created on: 19.01.2015
 *      Author: suze
 */

#ifndef INTERFACE_HANDLER_LOCAL_DIFFUSION_IMPL_H_
#define INTERFACE_HANDLER_LOCAL_DIFFUSION_IMPL_H_


#include "interface_handler_diffusion_user_data.h"


namespace ug{
    
    
///////////////////////////////////////////////////////////////
/// methods for class 'IInterfaceHandlerLocal'
///////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////
/// methods for class 'InterfaceHandlerLocalDiffusion'
///////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////
//	Constructor
///////////////////////////////////////////////////////////////
template <int TWorldDim>
InterfaceHandlerLocalDiffusion<TWorldDim>::
InterfaceHandlerLocalDiffusion(SmartPtr<DiffusionInterfaceProvider<dim> > interfaceProvider,
                               SmartPtr<CutElementHandler_TwoSided<dim> > cutElementHandler) :
    InterfaceHandlerLocalBase<dim>(cutElementHandler),
    m_numFct(0),
    m_numCo(0),
    m_scaleDoFs(false),
    m_L2Error(0.0),
    m_bNearInterface(false),
    m_shift_DoFIndex_tri(false),
    m_shift_DoFIndex_quad(false),
    m_luaSource_isSet(false),
    m_luaJump_isSet(false),
    m_luaJumpGrad_isSet(false),
    m_luaDiffusion_isSet(false),
    m_interfaceSource(0.0),
    m_interfaceJump(0.0),
    m_interfaceJumpGrad(0.0, 0.0),
    m_diffusionCoeff(0.0, 0.0),
    m_bNitsche(false),
    m_Gamma(0.0),
    m_insidePnt(0.0),
    m_Area(0.0),
    m_AreaOrig(0.0),
    m_AreaScale(0.0),
    m_spInterfaceProvider(interfaceProvider),
    m_spCutElementHandler(cutElementHandler)
{
    m_vBF.clear();
    m_MapInserted.clear();
    m_verticesPos.clear();
    m_verticesValue.clear();
    m_vRealCornerID.clear();
    m_vRealCornerID_tri.clear();
    m_vRealCornerID_quad.clear();
    m_vInterfaceID_tri.clear();
    m_vInterfaceID_quad.clear();
    
    m_MapInserted_Nitsche.clear();
    m_vAlpha.clear();
    m_vIntersectionPnts.clear();
    m_verticesGlobalIndex.clear();
}
//////////////////////////////////////////////////////////////////
/// virtual methods in 'InterfaceHandlerLocalBase'
//////////////////////////////////////////////////////////////////
 
template<int TWorldDim>
bool InterfaceHandlerLocalDiffusion<TWorldDim>::
get_intersection_point(MathVector<dim>& Intersect, Vertex* vrtOutsideCirc, Vertex* vrtInsideCirc)
{
    const int orientation = this->get_orientation();
    const int prtIndex = m_spCutElementHandler->get_prtIndex();

    if ( prtIndex == -1 ) UG_THROW("'get_intersection_point()': value of prtIndex not valid!\n");
  
    const MathVector<dim>& vrtPosOut = this->m_aaPos[vrtOutsideCirc];
    const MathVector<dim>& vrtPosIn  = this->m_aaPos[vrtInsideCirc];
        
    if ( orientation == 1 )
        return this->m_spInterfaceProvider->get_intersection_point(Intersect, vrtPosOut, vrtPosIn, prtIndex);
// inverse order of 'vrtPosOut' and 'vrtPosIn' for call of 'get_intersection_point()'
// to avoid error for alpha < 0:
    else if ( orientation == -1 )
        return this->m_spInterfaceProvider->get_intersection_point(Intersect, vrtPosIn, vrtPosOut, prtIndex);
    else
        UG_THROW("in InterfaceHandlerLocalDiffusion::get_intersection_point(): m_orientationInterface not set!\n");
}
  
template<int TWorldDim>
bool InterfaceHandlerLocalDiffusion<TWorldDim>::
get_intersection_point(MathVector<dim>& Intersect, Vertex* vrtOutsideCirc, Vertex* vrtInsideCirc, std::vector<number>& alphaOut)
{
    const int orientation = this->get_orientation();
    const int prtIndex = m_spCutElementHandler->get_prtIndex();
    
    if ( prtIndex == -1 ) UG_THROW("'get_intersection_point()': value of prtIndex not valid!\n");
    
    const MathVector<dim>& vrtPosOut = this->m_aaPos[vrtOutsideCirc];
    const MathVector<dim>& vrtPosIn  = this->m_aaPos[vrtInsideCirc];
        
    if ( orientation == 1 )
        return this->m_spInterfaceProvider->get_intersection_point(Intersect, vrtPosOut, vrtPosIn, prtIndex, alphaOut);
// inverse order of 'vrtPosOut' and 'vrtPosIn' for call of 'get_intersection_point()'
// to avoid error for alpha < 0:
    else if ( orientation == -1 )
        return this->m_spInterfaceProvider->get_intersection_point(Intersect, vrtPosIn, vrtPosOut, prtIndex, alphaOut);
    else
        UG_THROW("in InterfaceHandlerLocalDiffusion::get_intersection_point(): m_orientationInterface not set!\n");
}
    
///////////////////////////////////////////////////////////////
//	new methods for Diffusion
///////////////////////////////////////////////////////////////


template<int TWorldDim>
void InterfaceHandlerLocalDiffusion<TWorldDim>::
resize_local_data(LocalVector locU)
{
	LocalIndices ind = locU.get_indices();

// resize for cut triangle
	ind.resize_dof(0, 3);
	m_locU_tri.resize(ind);
	m_locD_tri.resize(ind);
	m_locJ_tri.resize(ind);

// if node is near interface: triangle is cut into 2 triangles:
	size_t local_dimension = 4;
	if ( get_bNearInterface() )
		local_dimension = 3;

// resize for cut quadrilateral
	ind.resize_dof(0, local_dimension);
	m_locU_quad.resize(ind);
	m_locD_quad.resize(ind);
	m_locJ_quad.resize(ind);
    
	return;
}

template<int TWorldDim>
void InterfaceHandlerLocalDiffusion<TWorldDim>::
set_interface_values(const std::vector<double > verticesValues)
{
	this->m_verticesValue.clear();

	for (size_t i = 0; i < verticesValues.size(); ++i)
 		this->m_verticesValue.push_back(verticesValues[i]);


// through error:
	if ( this->m_verticesValue.size() != verticesValues.size() )
	{
		UG_LOG("m_verticesValue.size(): " << this->m_verticesValue.size() << "\n");
		UG_LOG("verticesValues.size(): " << verticesValues.size() << "\n");
		UG_THROW("in InterfaceHandlerLocalDiffusion::set_interface_values: wrong size of m_verticesValue!\n");
	}

}

template<int TWorldDim>
void InterfaceHandlerLocalDiffusion<TWorldDim>::
reset_defect_on_interface(LocalVector& locD, const size_t size)
{
	if ( size > locD.num_all_dof(0) )
    {
        UG_LOG("in 'reset_defect_on_interface()': size = " << size << ", locD.size = " << locD.num_all_dof(0) << "\n");
		UG_THROW("in 'reset_defect_on_interface()': size = " << size << ", locD.size = " << locD.num_all_dof(0) << " => claimed size is NOT equal to size of solution vector!\n");
    }
// loop and set solution in 'solU_tri':
	for (size_t dof = 0; dof < locD.num_all_dof(0); ++dof)
	{
	// if dof_real is index of m_vertex: get solution from class:
		if ( this->lies_onInterface_size(dof, size) )
			locD.value(0, dof) = 0.0;
	}
}

template<int TWorldDim>
void InterfaceHandlerLocalDiffusion<TWorldDim>::
reset_jacobian_on_interface(LocalMatrix& locJ, const size_t size)
{

	if ( size > locJ.num_all_row_dof(0) )
    {
        UG_LOG("in 'reset_jacobian_on_interface()': size = " << size << ", locJ.size = " << locJ.num_all_row_dof(0) << "\n");
		UG_THROW("in 'reset_jacobian_on_interface()': size = " << size << ", locJ.size = " << locJ.num_all_row_dof(0) << " => claimed size is NOT equal to size of solution vector!\n");
    }
// loop and set solution in 'solU_tri':
	for (size_t dof1 = 0; dof1 < locJ.num_all_row_dof(0); ++dof1)
	{
		if ( this->lies_onInterface_size(dof1, size) )
		{
		// erase all col-values of chosen row dof1:
			for (size_t dof2 = 0; dof2 < locJ.num_all_col_dof(0); ++dof2)
				locJ.value(0, dof1, 0, dof2) = 0.0;
		}

	}

}

template<int TWorldDim>
const bool InterfaceHandlerLocalDiffusion<TWorldDim>::
get_boolian_for_diffusion()
{
    if( this->elementModus() == OUTSIDE_DOM)
    {
        return true;
    }
    else if( this->elementModus() == INSIDE_DOM)
    {
        return false;
    }
    else if( this->elementModus() == CUT_BY_INTERFACE)
    {
        return false;
    }
    else
    {
        UG_THROW("in InterfaceHandlerLocalDiffusion:get_boolian_for_diffusion: no valid boolian found!\n");
        return false;
    }
        
    return false;
}

///////////////////////////////////////////////////////////////
/// hard coded boundary conditions
///////////////////////////////////////////////////////////////
    
template<int dim>
inline double get_jump_value_ex6(const MathVector<dim> position, const int orientation)
{
    if ( orientation == 1)
        return 0.0;
    
    return exp(position[0])*cos(position[1]);
}
    
template<int dim>
inline double get_jump_grad_value_ex6(const MathVector<dim> position, const int orientation)
{
    if ( orientation == -1)
        return 0.0;
    
    double x = position[0];
    double y = position[1];
    
    return 2.0*exp(x)*(y*sin(y)-x*cos(y));
}
   
    
template <int dim>
inline double get_source_kappa(const MathVector<dim> position, const MathVector<dim> center, const int orientation)
{
        
    if ( orientation == 1)
        return  16.0*16.0;
    
    if ( orientation != -1)
        UG_THROW("wrong orientation!\n");
    
    double dist_x = position[0] - center[0];
    double dist_y = position[1] - center[1];
    double dist = sqrt(dist_x*dist_x+dist_y*dist_y);
            
    return 200*16*dist*dist;
        
    
}
    
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

// _impl()-methods used by user to specify the function
// depending on position
    
template<int TWorldDim>
double InterfaceHandlerLocalDiffusion<TWorldDim>::
get_jump_impl(const MathVector<dim> position)
{
// provide class member data
    const int orientation = this->get_orientation();
    
// call inline function:
    return get_jump_value_ex6<dim>(position, orientation);
    
}
    
template<int TWorldDim>
double InterfaceHandlerLocalDiffusion<TWorldDim>::
get_jump_grad_impl(const MathVector<dim> position)
{
// provide class member data
    const int orientation = this->get_orientation();
        
// call inline function:
    return get_jump_grad_value_ex6<dim>(position, orientation);
 }
    
template<int TWorldDim>
double InterfaceHandlerLocalDiffusion<TWorldDim>::
get_source_impl(const MathVector<dim> position)
{
// provide class member data
    const MathVector<dim> center = get_center(0);
    const int orientation = this->get_orientation();
    
// call inline function:
    return get_source_kappa<dim>(position, center, orientation);
}
 
    
// methods calling the flag to decide whether constant lua user data is
// appplied of position-dependent inline c++ funtion _impl() (see above)
    
template<int TWorldDim>
double InterfaceHandlerLocalDiffusion<TWorldDim>::
get_jump(const MathVector<dim> position)
{
    number jump_value = 0.0;

    if (m_luaJump_isSet)
    {
        jump_value = m_interfaceJump;
    }
    else
    // if the boundary conditions are NOT given by lua call, an own 'inline' implementation is needed:
        jump_value = get_jump_impl(position);

	return jump_value;
}


template<int TWorldDim>
double InterfaceHandlerLocalDiffusion<TWorldDim>::
get_jump_grad(const MathVector<dim> position)
{
    number jump_grad_value = 0.0;
    
    if (m_luaJumpGrad_isSet)
    {
        if ( this->get_orientation() == 1 )
            jump_grad_value = m_interfaceJumpGrad[0];
        else
            jump_grad_value = m_interfaceJumpGrad[1];
    }
    else
    {
    // if the boundary conditions are NOT given by lua call, an own 'inline' implementation is needed:
        jump_grad_value = get_jump_grad_impl(position);
    }
    
    
    return jump_grad_value;
}


template<int TWorldDim>
double InterfaceHandlerLocalDiffusion<TWorldDim>::
get_source(const MathVector<dim> position)
{
    number source_value = 0.0;
    
    if ( m_luaSource_isSet)
        source_value = m_interfaceSource;
    else
    {
    // if the boundary conditions are NOT given by lua call, an own 'inline' implementation is needed:
        source_value = get_source_impl(position);
    }
    
    return source_value;
}

// the diffusion only has the option for constant values
    
template<int TWorldDim>
number InterfaceHandlerLocalDiffusion<TWorldDim>::
get_diffusion()
{
    const int orientation = this->get_orientation();
        
// orientation == 1 corresponds to OUTSIDE of the circle (if the interface is a circle)
    if ( orientation == 1 )
        return m_diffusionCoeff[0];
        
    return m_diffusionCoeff[1];
        
}
    
template<int TWorldDim>
number InterfaceHandlerLocalDiffusion<TWorldDim>::
get_diffusion(const bool bElementIsOutside)
{
// bElementIsOutside corresponds to INSIDE of the circle ( if the interface is a circle)
    if ( bElementIsOutside )
        return m_diffusionCoeff[1];
        
    return m_diffusionCoeff[0];
        
}
    
///////////////////////////////////////////////////////////////
/// setter functions called during elem disc
/// 'ConvectionDiffusionFV1_cutElem' to set bndCond
///////////////////////////////////////////////////////////////

template<int TWorldDim>
void InterfaceHandlerLocalDiffusion<TWorldDim>::
set_source(const std::vector<double> sourceIm, LocalVector& sourceOut, LocalIndices ind,
           const size_t size, const bool bElementIsCut)
{
    ind.resize_dof(0, size);
    sourceOut.resize(ind);
        
// loop and write solution to 'sourceOut':
    for (size_t dof = 0; dof < sourceOut.num_all_dof(0); ++dof)
    {
        if ( !bElementIsCut )
        {
            sourceOut.value(0, dof) = sourceIm[dof];
        }
    // if dof_real is index of m_vertex: get values from InterfaceHandler:
        else if ( this->lies_onInterface_size(dof, size) )
        {
            size_t dof_real = this->real_index_size(dof, size);
            sourceOut.value(0, dof) = get_source(this->get_VerticesPos(dof_real));
        }
        else
        {
            size_t dof_orig = this->corner_orig(dof);
            sourceOut.value(0, dof) = sourceIm[dof_orig];
        }
        
    }
    
    return;
}
    

template<int TWorldDim>
void InterfaceHandlerLocalDiffusion<TWorldDim>::
set_jump_values(LocalVector& jumpOut, LocalIndices ind, const size_t size)
{
    ind.resize_dof(0, size);
    jumpOut.resize(ind);
    
// loop and write values to 'jumpOut':
    for (size_t dof = 0; dof < jumpOut.num_all_dof(0); ++dof)
    {
        
    // if dof_real is index of m_vertex: get solution from class:
        if ( this->lies_onInterface_size(dof, size) )
        {
            size_t dof_real = this->real_index_size(dof, size);
            jumpOut.value(0, dof) = get_jump(this->get_VerticesPos(dof_real));
        }
        else
        {
            jumpOut.value(0, dof) = 0.0;
        }
    }
        
    return;
}



template<int TWorldDim>
void InterfaceHandlerLocalDiffusion<TWorldDim>::
set_jump_grad_values(LocalVector& jumpGradOut, LocalIndices ind, const size_t size)
{
    ind.resize_dof(0, size);
    jumpGradOut.resize(ind);
        
    // loop and set solution in 'solU_tri':
    for (size_t dof = 0; dof < jumpGradOut.num_all_dof(0); ++dof)
    {
            
    // if dof_real is index of m_vertex: get solution from class:
        if ( this->lies_onInterface_size(dof, size) )
        {
            size_t dof_real = this->real_index_size(dof, size);
            jumpGradOut.value(0, dof) = get_jump_grad(this->get_VerticesPos(dof_real));
         }
        else
        {
            jumpGradOut.value(0, dof) = 0.0;
        }
    }
        
    return;
        
}
    
template<int TWorldDim>
void InterfaceHandlerLocalDiffusion<TWorldDim>::
set_local_sol(LocalVector& solU, const size_t size, const LocalVector& lvec, const int orientation)
{
	if ( size > solU.num_all_dof(0) )
    {
        UG_LOG("in 'set_local_sol()': size = " << size << ", solU.size = " << solU.num_all_dof(0) << "\n");
		UG_THROW("in 'set_local_sol()': size = " << size << ", solU.size = " << solU.num_all_dof(0) << " => claimed size is NOT equal to size of solution vector!\n");
    }
// loop and set solution in 'solU':
	for (size_t dof = 0; dof < solU.num_all_dof(0); ++dof)
	{
		size_t dof_real = this->real_index_size(dof, size);

	// if dof_real is index of m_vertex: get solution from the 'm_verticesValue'-array:
		if ( this->lies_onInterface_size(dof, size) )
		{
			solU.value(0, dof) = this->get_sol(dof_real);
		//////////////////////////////////////////////
		// add jump value in case of element, lying on \Omega^-:
		//////////////////////////////////////////////
	//		solU.value(0, dof) += get_jump_value_ex3(this->get_VerticesPos(dof_real));
		}
		else
		{
			solU.value(0, dof) = lvec.value(0,dof_real);
		}
	}
}



///////////////////////////////////////////////////////////////
//	new methods from base class
///////////////////////////////////////////////////////////////

// see mod_elem_flat_top() of flat_top.h
template <int TWorldDim>
void InterfaceHandlerLocalDiffusion<TWorldDim>::
compute_cut_element_data(GridObject* elem)
{
// get new element corners and according element type
	if ( this->StdFV_assembling() ) // Version an Stelle von 'm_bUsualAss = true'
	{
		this->CollectCorners_StdFV(elem);

		if ( dim == 2 ) this->set_roid_2d();
		if ( dim == 3 ) this->set_roid_3d();
	}
	else
	{
		if ( dim == 2 )
		{
			if ( get_Nitsche() )
				{ Collect_Data_Nitsche(elem); this->set_roid_2d(); }
			else
				{ CollectCorners_FlatTop_2d(elem); this->set_roid_2d(); }
 		}

		if ( dim == 3 ) { this->CollectCorners_FlatTop_3d(elem); this->set_roid_3d(); }
	}

 // some checks:
	if ( this->m_vCornerCoords.size() == 0 )
		UG_THROW("m_vCornerCoords.size() = " << this->m_vCornerCoords.size() << "not possible!\n");

	if ( dim == 2 )
		if ( elem->reference_object_id() != ROID_TRIANGLE )
			UG_THROW("Discretisation only coded for triangular elements!\n");

	if ( dim == 3 )
		if ( elem->reference_object_id() != ROID_TETRAHEDRON )
			UG_THROW("Discretisation only coded for tetrahedral elements!\n");

// output computed data to file
    if ( this->print_cutElment_data() )
        print_CutElementData();



}

} // end namespace ug



#endif /* INTERFACE_HANDLER_LOCAL_DIFFUSION_IMPL_H_ */
