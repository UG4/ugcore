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
                               SmartPtr<CutElementHandlerImmersed<dim> > cutElementHandler) :
    InterfaceHandlerLocalBase<dim>(cutElementHandler),
    m_spInterfaceProvider(interfaceProvider),
    m_spCutElementHandler(cutElementHandler),
    m_bNitsche(false),
    m_numFct(0),
    m_numCo(0),
    m_shift_DoFIndex_tri(false),
    m_shift_DoFIndex_quad(false),
    m_interfaceSource(1100.13),
    m_interfaceJump(1100.13),
    m_interfaceJumpGrad(1100.13,1100.13),
    m_diffusionCoeff(1100.13, 1100.13)
{
    m_verticesPos.clear();
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
write_solution(const std::vector<double > verticesValues)
{
	this->m_verticesValue.clear();

	for (size_t i = 0; i < verticesValues.size(); ++i)
 		this->m_verticesValue.push_back(verticesValues[i]);


// through error:
	if ( this->m_verticesValue.size() != verticesValues.size() )
	{
		UG_LOG("m_verticesValue.size(): " << this->m_verticesValue.size() << "\n");
		UG_LOG("verticesValues.size(): " << verticesValues.size() << "\n");
		UG_THROW("in InterfaceHandlerLocalDiffusion::write_solution: wrong size of m_verticesValue!\n");
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
   

    
    template <int dim>
    inline double get_jump_value_ex5(const MathVector<dim> position, const MathVector<dim> center, const int orientation)
    {
        return 0.0;
    }
    
    template <int dim>
    inline double get_jump_grad_value_ex3(const MathVector<dim> position, const MathVector<dim> center, const int orientation)
    {
        return 0.0;
        
        if ( orientation == -1)
            return 0.0;
        
        double absValue = position[0]*position[0] + position[1]* position[1];
        double sum = position[0] + position[1];
        
        double returnValue = 8*(2*absValue - sum)*exp(-absValue);
        
        return returnValue;
    }
    
    template <int dim>
    inline double get_source_kappa(const MathVector<dim> position, const MathVector<dim> center, const int orientation)
    {
        
        if ( orientation == 1)
        {
            return  16.0*16.0;
        }
        else
        {
            if ( orientation != -1)
                UG_THROW("wrong orientation!\n");
            
            double dist_x = position[0] - center[0];
            double dist_y = position[1] - center[1];
            double dist = sqrt(dist_x*dist_x+dist_y*dist_y);
            
            return 200*16*dist*dist;
            
        }
        
    }
    

    
template<int TWorldDim>
double InterfaceHandlerLocalDiffusion<TWorldDim>::
get_jump_impl(const MathVector<dim> position)
{
// provide class member data
    const MathVector<dim> center = get_center(0);
    const int orientation = this->get_orientation();
    
// call inline function:
    return get_jump_value_ex5<dim>(position, center, orientation);
    
}
    
template<int TWorldDim>
double InterfaceHandlerLocalDiffusion<TWorldDim>::
get_jump_grad_impl(const MathVector<dim> position)
{
// provide class member data
    const MathVector<dim> center = get_center(0);
    const int orientation = this->get_orientation();
        
// call inline function:
    return get_jump_grad_value_ex3<dim>(position, center, orientation);
        
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
 
    
    
    
template<int TWorldDim>
double InterfaceHandlerLocalDiffusion<TWorldDim>::
get_jump(const MathVector<dim> position)
{
    number jump_value = 0.0;

    if (m_bBndFct)
    {
    // if the boundary conditions are given ...
        jump_value = get_jump_impl(position);
    }
    else
        jump_value = m_interfaceJump;

	return jump_value;
}


template<int TWorldDim>
double InterfaceHandlerLocalDiffusion<TWorldDim>::
get_jump_grad(const MathVector<dim> position)
{
    number jump_grad_value = 0.0;
    
    if (m_bBndFct)
    {
    // if the boundary conditions are given ...
        jump_grad_value = get_jump_grad_impl(position);
    }
    else
    {
        if ( this->get_orientation() == 1 )
            jump_grad_value = m_interfaceJumpGrad[0];
        else
            jump_grad_value = m_interfaceJumpGrad[1];
    }
    
    return jump_grad_value;
}


template<int TWorldDim>
double InterfaceHandlerLocalDiffusion<TWorldDim>::
get_source(const MathVector<dim> position)
{
    number source_value = 0.0;
    
    if (m_bBndFct)
    {
    // if the boundary conditions are given ...
        source_value = get_source_impl(position);
    }
    else
        source_value = m_interfaceSource;
    
    return source_value;
}

    
///////////////////////////////////////////////////////////////
/// setter functions called during elem disc to set bndCond
///////////////////////////////////////////////////////////////

template<int TWorldDim>
number InterfaceHandlerLocalDiffusion<TWorldDim>::
get_diffusion()
{
    const int orientation = this->get_orientation();

   // orientation == 1 corresponds to outside circle, if interface is a circle
    if ( orientation == 1 ) return m_diffusionCoeff[0];
    else return m_diffusionCoeff[1];

}

template<int TWorldDim>
number InterfaceHandlerLocalDiffusion<TWorldDim>::
get_diffusion(const bool bElementIsOutside)
{        
// m_diffusionCoeff[1] corresponds to inside circle, bElementIsOutside corresponds to outside circle
    if ( bElementIsOutside ) return m_diffusionCoeff[1];
    else return m_diffusionCoeff[0];
        
}
    
template<int TWorldDim>
LocalVector InterfaceHandlerLocalDiffusion<TWorldDim>::
set_source(const std::vector<double> sourceIm, LocalIndices ind, const size_t size, const bool bElementIsCut)
{
	LocalVector source;
	ind.resize_dof(0, size);
	source.resize(ind);

// loop and set solution in 'solU_tri':
	for (size_t dof = 0; dof < source.num_all_dof(0); ++dof)
	{
/*		if ( dof > 2 )
			source.value(0, dof) = sourceIm[2]; // done during 'add_def_elem_local()'
		else
			source.value(0, dof) = sourceIm[dof]; // done during 'add_def_elem_local()'
*/

	// if dof_real is index of m_vertex: get solution from class:
		if ( !bElementIsCut )
		{
  			source.value(0, dof) = sourceIm[dof]; // done during 'add_def_elem_local()'
 		}
		else if ( this->lies_onInterface_size(dof, size) )
		{
			size_t dof_real = this->real_index_size(dof, size);
 			source.value(0, dof) = get_source(this->get_VerticesPos(dof_real));
// 			source.value(0, dof) = get_source_Fedkiw_ex5(this->get_VerticesPos(dof_real));
		}
		else
		{
			size_t dof_orig = this->corner_orig(dof);
		// careful: this->corner() returns corner coords of new corners => use dof-Index, NOT dof_orig-Index?!
 //			source.value(0, dof) = get_source_kappa(this->corner(dof));

			source.value(0, dof) = sourceIm[dof_orig]; // done during 'add_def_elem_local()'
 		}

	}

	return source;
}


template<int TWorldDim>
LocalVector InterfaceHandlerLocalDiffusion<TWorldDim>::
set_jump_values(LocalIndices ind, const size_t size)
{
	LocalVector jump;
	ind.resize_dof(0, size);
	jump.resize(ind);

// loop and set solution in 'solU_tri':
	for (size_t dof = 0; dof < jump.num_all_dof(0); ++dof)
	{

	// if dof_real is index of m_vertex: get solution from class:
		if ( this->lies_onInterface_size(dof, size) )
		{
			size_t dof_real = this->real_index_size(dof, size);
 			jump.value(0, dof) = get_jump(this->get_VerticesPos(dof_real));
		}
		else
		{
			jump.value(0, dof) = 0.0;
		}
	}

	return jump;
}


template<int TWorldDim>
LocalVector InterfaceHandlerLocalDiffusion<TWorldDim>::
set_jump_grad_values(LocalIndices ind, const size_t size)
{
	LocalVector jump_grad;
	ind.resize_dof(0, size);
	jump_grad.resize(ind);

// loop and set solution in 'solU_tri':
	for (size_t dof = 0; dof < jump_grad.num_all_dof(0); ++dof)
	{

	// if dof_real is index of m_vertex: get solution from class:
		if ( this->lies_onInterface_size(dof, size) )
		{
			size_t dof_real = this->real_index_size(dof, size);
 			jump_grad.value(0, dof) = get_jump_grad(this->get_VerticesPos(dof_real));
//			jump_grad.value(0, dof) = get_jump_grad_value_ex5(this->get_VerticesPos(dof_real));
//            jump_grad.value(0, dof) = get_jump_grad_value_kappa_Frei(this->get_VerticesPos(dof_real));
		}
		else
		{
			jump_grad.value(0, dof) = 0.0;
		}
	}

	return jump_grad;

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
// loop and set solution in 'solU_tri':
	for (size_t dof = 0; dof < solU.num_all_dof(0); ++dof)
	{
        
		size_t dof_real = this->real_index_size(dof, size);

	// if dof_real is index of m_vertex: get solution from class:
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
	//		UG_LOG("tri: *value = " << lvec.value(0,dof_real) << "\n");
			solU.value(0, dof) = lvec.value(0,dof_real);
		}
	}
}
    
template <int TWorldDim>
bool InterfaceHandlerLocalDiffusion<TWorldDim>::
check_interface_data(const bool bBndFct)
{
    // check value 1100.13 is set as default during constructor:
    if ( m_interfaceSource == 1100.13 || m_interfaceJump == 1100.13 || m_interfaceJumpGrad[0] == 1100.13 || m_interfaceJumpGrad[1] == 1100.13 || m_diffusionCoeff[0] == 1100.13 || m_diffusionCoeff[1] == 1100.13)
        UG_THROW("interface data not completely specified! call setter functions to define diffusion coefficients, source, jump and jump of the gradient at the interface.\n");
        
    // if the boundary data is given by C++ functions, the values need to be set to zero:
    if ( bBndFct )
    {
        if ( fabs(m_interfaceSource) > 0.000001 || fabs(m_interfaceJump) > 0.000001 || fabs(m_interfaceJumpGrad[0]) > 0.000001|| fabs(m_interfaceJumpGrad[0]) > 0.000001 )
            UG_THROW("if the boundary data is given by C++ functions, the values need to be set to zero\n");
    }
        
    m_bBndFct = bBndFct;
}
    
/*
template <int TWorldDim>
int InterfaceHandlerLocalDiffusion<TWorldDim>::
CollectCorners_FlatTop_2d(GridObject* elem)
{
	//////////////////////////////////////////////
	// 1) fill vector with fluid corners:
	//////////////////////////////////////////////

	this->m_vCornerCoords.clear();
	this->m_vInterfaceID.clear();
	this->m_vOriginalCornerID.clear();

// buffer vectors for (cornerCoords, cornerIndex)
	std::vector<std::pair<MathVector<dim>, size_t > > vOutsideCorners;
	std::vector<std::pair<MathVector<dim>, size_t > > vInsideCorners;
	std::vector<std::pair<MathVector<dim>, size_t > > vNearIntCorners;

//	collect all vertices of the element
	std::vector<Vertex*> vVertex;
	CollectVertices(vVertex, *this->m_spMG, elem);

 	bool isFTVertex = false;
 	for(size_t i = 0; i < vVertex.size(); ++i)
	{
	// remember boolian for check, weather there existe at least one vertex, which is FT!
		isFTVertex = this->is_FTVertex(vVertex[i], i);
		if ( isFTVertex )
			break;
	}

	if ( !isFTVertex )
		UG_THROW("Error in 'CollectCorners_FlatTop_2d': no vertex is FTVertex: should be true for at least 1 vertex!\n");

	//	collect all edges of the element
	std::vector<Edge*> vEdges;
	CollectEdgesSorted(vEdges, *this->m_spMG, elem);

	// loop vertices
	//////////////////////////////////////////////
	// REMARK:
	// order is the same as in 'vCornerCoords', therefore we can be sure, that the
	// order of the new 'vCornerIBCoords' will be consistent with the grid standard
	//////////////////////////////////////////////

	bool bNearInterface = false;
	for(size_t i = 0; i < vVertex.size(); ++i)
	{
		Vertex* vrtRoot = vVertex[i];

		//////////////////////////////////////////////
		// case 1:
		// vertex insideDomain
 		if ( !this->is_FTVertex(vrtRoot, i) )
 		{
			if ( this->is_nearInterfaceVertex(vrtRoot, i) )
				UG_THROW("NearInterface BUT !is_FT => neuerdings Fehler!!....\n");

			this->m_vCornerCoords.push_back(this->m_aaPos[vrtRoot]);
			this->m_vOriginalCornerID.push_back(i);

			vInsideCorners.push_back(std::make_pair(this->m_aaPos[vrtRoot], i));
		}
		//////////////////////////////////////////////
  		// case 2:
		// vertex = FT + ON interface
		// 		=> KEINE Berechnung von 'intersectionPoint' notwendig! -> pushen und alten index pushen

		// REMARK: is_nearInterfaceVerx = false per default, if m_vThresholdOnLevel = 0.0
		else if ( this->is_nearInterfaceVertex(vrtRoot, i) )
		{
 			bNearInterface = true;
 			this->m_vCornerCoords.push_back(this->m_aaPos[vrtRoot]);
 			this->m_vOriginalCornerID.push_back(i);
 			this->m_vInterfaceID.push_back(this->m_vCornerCoords.size()-1);  // attention: push AFTER 'm_vCornerCoords.push_back()'!!

			vOutsideCorners.push_back(std::make_pair(this->m_aaPos[vrtRoot], i));
			vNearIntCorners.push_back(std::make_pair(this->m_aaPos[vrtRoot], i));

		}
		//////////////////////////////////////////////
  		// case 3:
		// vertex 'outsideFluid'
		// 		=> NEUE Position berechen+pushen und alten index pushen
		else
		{
 			//////////////////////////////////////////////////////////////////////////////////////////
			// loop alle edges, die interface schneiden und damit einen neuen intersectionPnt
			// beitragen zum damit assoziierten alten index
			for(size_t e = 0; e < vEdges.size(); ++e)
			{
				Edge* edge = vEdges[e];
				std::vector<Vertex*> vVertexEdge;
				CollectVertices(vVertexEdge, *this->m_spMG, edge);
				if ( vVertexEdge.size() != 2 )
					UG_THROW("error in collecting vertices associated to an edge!....EXIT!...\n");

				Vertex* vrt1 = vVertexEdge[0];
				Vertex* vrt2 = vVertexEdge[1];
				size_t vrtInd1 = get_vertex_index(vrt1, elem);
				size_t vrtInd2 = get_vertex_index(vrt2, elem);

				MathVector<dim> intersectionPnt;

	 		///////////////////////////////////////////////////////////////////
 			// lies vrtRoot on a cutted edge?
		 	///////////////////////////////////////////////////////////////////
			// case1: vrtRoot is intersectionPnt with insideCorner = near_interface_corner => remove!
				if ( this->is_nearInterfaceVertex(vrt2, vrtInd2) || this->is_nearInterfaceVertex(vrt1, vrtInd1) )
				{ bNearInterface = true; continue; }
			 // case2: vert2 = outsideParticle && vrt1 = insideParticle:
				else if ( vrtRoot == vrt1 && !this->is_FTVertex(vrt2, vrtInd2) ){
					this->get_intersection_point(intersectionPnt, vrt2, vrt1);
 				}
			// case3: vrt1 = outsideParticle && vrt2 = insideParticle:
				else if ( vrtRoot == vrt2 && !this->is_FTVertex(vrt1, vrtInd1) )
					this->get_intersection_point(intersectionPnt, vrt1, vrt2);
				else
 					continue;

			// check for correct inersectionPnt
				if ( fabs(this->get_LSvalue_byPosition(intersectionPnt)) > 1e-6  )
					UG_THROW("in 'CollectIBCorners2d()': Error in computation of 'intersectionPnt':\n "
							" intersectionPnt = " << intersectionPnt << "\n distance from interace = " << fabs(get_LSvalue_byPosition(intersectionPnt)) << "\n");

	 		///////////////////////////////////////////////////////////////////
	 		// only push_back, if not included yet!
			// 	-> can be ONLY the case, if the intersectionPoint is a node
	 			if ( ! this->isIncluded(this->m_vCornerCoords, intersectionPnt) )
	 			{

	 				this->m_vCornerCoords.push_back(intersectionPnt);
	 				this->m_vOriginalCornerID.push_back(i);
	 				this->m_vInterfaceID.push_back(this->m_vCornerCoords.size()-1);  // attention: push AFTER 'm_vCornerCoords.push_back()'!!

	 				vOutsideCorners.push_back(std::make_pair(intersectionPnt, i));
   	 			}


 			} // end edge-loop

		} // end else-case

 	} // end vrt-loop

////////////////////////////////////////////////////////////////////////////////////////////
// Postprecessing for quadrilaterals ( <=>  vOutsideCorners == 2 )
// (vInsideCorners.size() == 2) && (bNearInterface)	 => ALL nodes insideFluid, BUT one ON surface
//		=> no Quadrilateral, but Triangle!!
////////////////////////////////////////////////////////////////////////////////////////////
	MathVector<dim> normalDir(0.0);
	if ( (this->m_vCornerCoords.size() == 4) && (!bNearInterface) && (dim == 2) )
		this->ResortQuadrilateral(vInsideCorners, vOutsideCorners, normalDir);
	else if ( bNearInterface )
	{
	// Quadrilateral -> Triangle
		if ( vInsideCorners.size() == 1 ) // case 1
		{
			// do nothing, since re-sorting not necessary...???
		}
	// skip whole element, since only FT points are included
		else if ( vInsideCorners.size() == 0 )
			UG_THROW("in 'CollectCorners_FlatTop_2d()': vInsideCorners.size() "
					"= " << vInsideCorners.size() << "not possible!\n");
	}


	return this->m_vCornerCoords.size();

}


// called by geo.update()!!
template <int TWorldDim>
bool InterfaceHandlerLocalDiffusion<TWorldDim>::
update_elem(GridObject* elem, const MathVector<TWorldDim>* vCornerCoords, int interfaceOrientation)
{
	bool do_update_local = false;
	this->m_vBF.clear();

// computing flat top modus
	this->m_elemModus = compute_element_modus(elem, interfaceOrientation);

	switch(this->m_elemModus)
	{
 		case INSIDE_DOM:	   if ( dim == 2 ) this->set_flat_top_data(elem, vCornerCoords, ROID_TRIANGLE);
 							   if ( dim == 3 ) this->set_flat_top_data(elem, vCornerCoords, ROID_TETRAHEDRON);
							   break;	// usual assembling
		case OUTSIDE_DOM: 	   if ( dim == 2 ) this->set_flat_top_data(elem, vCornerCoords, ROID_TRIANGLE);
							   if ( dim == 3 ) this->set_flat_top_data(elem, vCornerCoords, ROID_TETRAHEDRON);
							   break;	// usual assembling
		case CUT_BY_INTERFACE: this->compute_flat_top_data(elem);
								//if ( m_roid == ROID_PYRAMID )UG_THROW("PYRAMID\n");
							   do_update_local = true;
						  	   break;  // flat top assembling
		default:
			throw(UGError("Error in InterfaceHandlerLocalDiffusion::update(): switch(m_elemModus)!"));
	}

 	return do_update_local;

}
*/

///////////////////////////////////////////////////////////////
//	new methods from base class
///////////////////////////////////////////////////////////////

// see mod_elem_flat_top() of flat_top.h
template <int TWorldDim>
void InterfaceHandlerLocalDiffusion<TWorldDim>::
compute_flat_top_data(GridObject* elem)
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


	this->m_elemModus = this->get_element_modus(elem); // computed via 'compute_element_modus()' during 'update_marker()'

	if ( 0 ){
	if( this->m_elemModus == CUT_BY_INTERFACE)
	{
		UG_LOG("_________________ compute_flat_top_data()_________________\n");

		for ( size_t i = 0; i < this->m_vCornerCoords.size(); ++i )
			UG_LOG("m_vCornerCoords = " <<this-> m_vCornerCoords[i] << "\n");
		for ( size_t i = 0; i < this->m_vOriginalCornerID.size(); ++i )
			UG_LOG("Original: id = " << this->m_vOriginalCornerID[i] << "\n");
		UG_LOG("\n");
		for ( size_t i = 0; i < this->m_vInterfaceID.size(); ++i )
			UG_LOG("Interface: id = " << this->m_vInterfaceID[i] << "\n");
		UG_LOG("\n");
		if ( this->m_roid == ROID_TRIANGLE )
		{
			for ( size_t i = 0; i < this->m_vRealCornerID_tri.size(); ++i )
				UG_LOG("m_vRealCornerID_tri: id = " << m_vRealCornerID_tri[i] << "\n");
			UG_LOG("\n");
 			for ( size_t i = 0; i < m_vInterfaceID_tri.size(); ++i )
 				UG_LOG("Interface tri: id = " << m_vInterfaceID_tri[i] << "\n");
 			UG_LOG("\n");
		}
		if ( this->m_roid == ROID_QUADRILATERAL )
		{
			for ( size_t i = 0; i < m_vRealCornerID_quad.size(); ++i )
				UG_LOG("m_vRealCornerID_quad: id = " << m_vRealCornerID_quad[i] << "\n");
			UG_LOG("\n");
			for ( size_t i = 0; i < m_vInterfaceID_quad.size(); ++i )
				UG_LOG("Interface quad: id = " << m_vInterfaceID_quad[i] << "\n");
			UG_LOG("\n");
		}
	}
	}

}

template <int TWorldDim>
bool InterfaceHandlerLocalDiffusion<TWorldDim>::
update_elem(GridObject* elem, const MathVector<TWorldDim>* vCornerCoords)
{
	bool do_update_local = false;
	m_vBF.clear();

// computing flat top modus
	this->m_elemModus = this->compute_element_modus(elem, this->m_orientationInterface);

	switch(this->m_elemModus)
	{
 		case INSIDE_DOM:	   if ( dim == 2 ) this->set_flat_top_data(elem, vCornerCoords, ROID_TRIANGLE);
 							   if ( dim == 3 ) this->set_flat_top_data(elem, vCornerCoords, ROID_TETRAHEDRON);
							   break;	// usual assembling
		case OUTSIDE_DOM: 	   if ( dim == 2 ) this->set_flat_top_data(elem, vCornerCoords, ROID_TRIANGLE);
							   if ( dim == 3 ) this->set_flat_top_data(elem, vCornerCoords, ROID_TETRAHEDRON);
							   break;	// usual assembling
		case CUT_BY_INTERFACE: compute_flat_top_data(elem);
								//if ( m_roid == ROID_PYRAMID )UG_THROW("PYRAMID\n");
							   do_update_local = true;
						  	   break;  // flat top assembling
		default:
			throw(UGError("Error in InterfaceHandlerLocalBase::update(): switch(m_elemModus)!"));
	}

 	return do_update_local;

}


} // end namespace ug



#endif /* INTERFACE_HANDLER_LOCAL_DIFFUSION_IMPL_H_ */
