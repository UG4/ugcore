/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
 * Author: Susanne Hoellbacher
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#ifndef __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__PARTICLE_TRANSFER_IMPL__
#define __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__PARTICLE_TRANSFER_IMPL__

#include "particle_transfer.h"
#include "lib_disc/reference_element/reference_mapping_provider.h"
#include "lib_disc/local_finite_element/local_finite_element_provider.h"
#include "lib_disc/function_spaces/grid_function_util.h"
#include "lib_grid/algorithms/debug_util.h"								// ElementDebugInfo

////////////////////////////////////////////////////////////////////////////////
// Almost the same implementation as 'ParticleTransfer' (see std_transfer.h)
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// 	Changes compared to the 'ParticleTransfer' class are the call of
// 	'modify_vector()' in the method 'prolongate()' and 'restrict()'.
//
// 	The physical solution values are set in each inactive nodes within
//	particles in order to transfer the valid solutions between the levels.
////////////////////////////////////////////////////////////////////////////////


namespace ug{
    
    template <typename TAlgebra>
    void ParticleAssembleInjectionForP1Lagrange1(typename TAlgebra::matrix_type& mat,
                                                 const DoFDistribution& coarseDD,
                                                 const DoFDistribution& fineDD)
    {
        
        PROFILE_FUNC_GROUP("gmg");
        //  Allow only lagrange P1 functions
        for(size_t fct = 0; fct < fineDD.num_fct(); ++fct)
            if(fineDD.local_finite_element_id(fct).type() != LFEID::LAGRANGE ||
               fineDD.local_finite_element_id(fct).order() != 1)
                UG_THROW("AssembleInjectionForP1Lagrange: "
                         "Interpolation only implemented for Lagrange P1 functions.");
        
        // 	get MultiGrid
        const MultiGrid& grid = *coarseDD.multi_grid();
        
        // 	get number of dofs on different levels
        //	const size_t numFineDoFs = fineDD.num_indices();
        //	const size_t numCoarseDoFs = coarseDD.num_indices();
        
        // 	resize matrix
        //mat.resize_and_clear(numCoarseDoFs, numFineDoFs);
        
        std::vector<size_t> coarseInd, fineInd;
        
        // 	RegularVertex iterators
        typedef DoFDistribution::traits<RegularVertex>::const_iterator const_iterator;
        const_iterator iter, iterBegin, iterEnd;
        
        iterBegin = fineDD.template begin<RegularVertex>();
        iterEnd = fineDD.template end<RegularVertex>();
        
        // 	loop nodes of fine subset
        for(iter = iterBegin; iter != iterEnd; ++iter)
        {
            // 	get father
            GridObject* geomObj = grid.get_parent(*iter);
            Vertex* vert = dynamic_cast<Vertex*>(geomObj);
            
            //	Check if father is RegularVertex
            if(vert != NULL)
            {
                // get global indices
                coarseDD.inner_algebra_indices(vert, coarseInd);
            }
            else continue;
            
            // 	get global indices
            fineDD.inner_algebra_indices(*iter, fineInd);
            
            for(size_t i = 0; i < coarseInd.size(); ++i)
                mat(coarseInd[i], fineInd[i]) = 1.0;
        }
    }
    
    
    template <typename TAlgebra>
    void ParticleAssembleInjectionForP1Lagrange2(typename TAlgebra::matrix_type& P,
                                                 const DoFDistribution& coarseDD,
                                                 const DoFDistribution& fineDD)
    {
        UG_THROW("OK2...\n");
        PROFILE_FUNC_GROUP("gmg");
        // 	allow only lagrange P1 functions
        for(size_t fct = 0; fct < fineDD.num_fct(); ++fct)
            if(fineDD.lfeid(fct).type() != LFEID::LAGRANGE ||
               fineDD.lfeid(fct).order() != 1)
                UG_THROW("AssembleStdProlongationForP1Lagrange:"
                         "Interpolation only implemented for Lagrange P1 functions.");
        
        //  resize matrix
        P.resize_and_clear(fineDD.num_indices(), coarseDD.num_indices());
        
        //  iterators
        const MultiGrid& mg = *coarseDD.multi_grid();
        typedef DoFDistribution::traits<Vertex>::const_iterator const_iterator;
        const_iterator iter, iterBegin, iterEnd;
        
        //  loop subsets on fine level
        std::vector<size_t> vParentIndex, vChildIndex;
        std::vector<DoFIndex> vParentDoF, vChildDoF;
        for(int si = 0; si < fineDD.num_subsets(); ++si)
        {
            iterBegin = fineDD.template begin<Vertex>(si);
            iterEnd = fineDD.template end<Vertex>(si);
            
            //  loop vertices for fine level subset
            for(iter = iterBegin; iter != iterEnd; ++iter)
            {
                //	get element
                Vertex* child = *iter;
                
                //  get father
                GridObject* parent = mg.get_parent(child);
                
                //	check if child contained in coarseDD. This should always be false
                //	for a GridLevel::LEVEL, but might be the case for GridLevel::SURFACE
                //	and an adaptive grid-part used by both dds. In such a case we can
                //	simply set identity.
                if(coarseDD.is_contained(child)){
                    //	get indices
                    coarseDD.inner_algebra_indices(child, vParentIndex);
                    fineDD.inner_algebra_indices(child, vChildIndex);
                    UG_ASSERT(vParentIndex.size() == vChildIndex.size(), "Size mismatch");
                    
                    //	set identity
                    for(size_t i = 0; i < vParentIndex.size(); ++i)
                        P(vChildIndex[i], vParentIndex[i]) = 1.0;
                    
                    //	this child is perfectly handled
                    continue;
                }
                else{
                    //	check if parent exists (this should always be the case, except in
                    //	the case that 'child' is a v-slave)
                    if(!parent) continue;
                    
                    if(!coarseDD.is_contained(parent)){
                        UG_THROW("ParticleTransfer: A parent element is not contained in "
                                 " coarse-dd nor the child element in the coarse-dd. "
                                 "This should not happen.")
                    }
                }
                
                //	type of father
                const ReferenceObjectID roid = parent->reference_object_id();
                
                //	loop all components
                for(size_t fct = 0; fct < fineDD.num_fct(); fct++)
                {
                    //	check that fct defined on subset
                    if(!fineDD.is_def_in_subset(fct, si)) continue;
                    
                    //  get global indices
                    fineDD.inner_dof_indices(child, fct, vChildDoF);
                    
                    //	detect type of father
                    switch(roid)
                    {
                        case ROID_VERTEX:
                        {
                            Vertex* vrt = dynamic_cast<Vertex*>(parent);
                            coarseDD.inner_dof_indices(vrt, fct, vParentDoF);
                            DoFRef(P, vChildDoF[0], vParentDoF[0]) = 1.0;
                        }
                            break;
                        default: UG_THROW("AssembleStdProlongationForP1Lagrange: Element Father"
                                          "is of unsupported type "<<roid);
                    }
                }
            }
        }
    }
    
    
    
    template <typename TDomain, typename TAlgebra>
    void ParticleTransfer<TDomain, TAlgebra>::
    set_solution_to_zero(vector_type& sol, const int levIndex, const int prtIndex, ConstSmartPtr<DoFDistribution> dd)
    {
        
        typedef typename std::vector<grid_base_object*>::iterator ListIter;
        std::vector<grid_base_object*> ElemListLog = m_spParticleHandlerGlobal->m_vvvElemListCut[levIndex][prtIndex];
        
        for(ListIter listIter = ElemListLog.begin(); listIter != ElemListLog.end(); ++listIter)
        {
            //	collect all vertices of the element
            std::vector<Vertex*> vVertex;
            CollectVertices(vVertex, *m_spParticleHandlerGlobal->m_spMG, *listIter);
            
            //	loop vertices
            for(size_t v = 0; v < vVertex.size(); ++v)
            {
                if ( m_spParticleHandlerGlobal->m_spFlatTopVrtMarker->is_marked(vVertex[v]) )
                {
                    // loop velocity DoFs
                    for (size_t fct = 0; fct < dim; ++fct)
                    {
                        //	create multi index
                        std::vector<DoFIndex>  vInd;
                        //	get multi indices
                        if(dd->inner_dof_indices(vVertex[v], fct, vInd) != 1)
                            UG_THROW("Only one index expected.");
                        
                        //	set velocity solution in FlatTopVrt to zero
                        DoFRef(sol, vInd[0]) = 0.0;
                    }
                }
                
            } // end vrt-loop
        } // end cutElem-loop
        
    }
    
    template <typename TDomain, typename TAlgebra>
    void ParticleTransfer<TDomain, TAlgebra>::
    adjust_prolongation(vector_type& corrFine, GridLevel fineLvl,
                        const vector_type& corrCoarse, GridLevel coarseLvl,
                        ConstSmartPtr<ApproximationSpace<TDomain> > spApproxSpace)
    {
        const int coarseIndex = m_spParticleHandlerGlobal->get_Index(coarseLvl);
        const int fineIndex = m_spParticleHandlerGlobal->get_Index(fineLvl);
        
        size_t numPrt = m_spParticleHandlerGlobal->num_particles();
        for(size_t p = 0; p < numPrt; ++p)
        {
            
#ifdef UG_PARALLEL
            std::vector< grid_base_object* > ElemList = m_spParticleHandlerGlobal->m_vvvElemListCut[fineIndex][p];
            if( ElemList.size() == 0 )
                continue;
#endif
            ///////////////////////////////////////////////////////////
            // 	1) projection of the defect of extrapolated ghost vrt:
            ////////////////////////////////////////////////////////////
            
            ConstSmartPtr<DoFDistribution> ddFine = spApproxSpace->dof_distribution(fineLvl);
            set_solution_to_zero(corrFine, fineIndex, p, ddFine);
            
            ///////////////////////////////////////////////////////////
            // 	2) injection for extraDoFs:
            ////////////////////////////////////////////////////////////
            std::vector<DoFIndex> transIndCoarse = m_spParticleHandlerGlobal->get_transInd(coarseIndex, p);
            std::vector<DoFIndex> rotIndCoarse = m_spParticleHandlerGlobal->get_rotInd(coarseIndex, p);
            
            std::vector<DoFIndex> transIndFine = m_spParticleHandlerGlobal->get_transInd(fineIndex, p);
            std::vector<DoFIndex> rotIndFine = m_spParticleHandlerGlobal->get_rotInd(fineIndex, p);
            
            for ( int d = 0; d < dim; ++d )
            {
                DoFRef(corrFine, transIndFine[d]) 	= DoFRef(corrCoarse, transIndCoarse[d]);
                DoFRef(corrFine, rotIndFine[d]) 	= DoFRef(corrCoarse, rotIndCoarse[d]);
            }
        }
        
        
    }
    
    
    template <typename TDomain, typename TAlgebra>
    void ParticleTransfer<TDomain, TAlgebra>::
    adjust_restriction(vector_type& dCoarse, GridLevel coarseLvl,
                       const vector_type& dFine, GridLevel fineLvl,
                       ConstSmartPtr<ApproximationSpace<TDomain> > spApproxSpace)
    {
        // UG_LOG("nothing...\n");
        
        ConstSmartPtr<DoFDistribution> ddCoarse = spApproxSpace->dof_distribution(coarseLvl);
        ConstSmartPtr<DoFDistribution> ddFine = spApproxSpace->dof_distribution(coarseLvl);
        
        const int coarseIndex = m_spParticleHandlerGlobal->get_Index(coarseLvl,ddCoarse);
        const int fineIndex = m_spParticleHandlerGlobal->get_Index(fineLvl,ddFine);
        
        size_t numPrt = m_spParticleHandlerGlobal->num_particles();
        //	loop over all particles and initialize defects:
        for(size_t p = 0; p < numPrt; ++p)
        {
            
#ifdef UG_PARALLEL
            std::vector< grid_base_object* > ElemList = m_spParticleHandlerGlobal->m_vvvElemListCut[fineIndex][p];
            if( ElemList.size() == 0 )
                continue;
#endif
            
            ///////////////////////////////////////////////////////////
            // 		1) set solution outside fluid to zero:
            ////////////////////////////////////////////////////////////
            
            set_solution_to_zero(dCoarse, coarseIndex, p, ddCoarse);
            
            ///////////////////////////////////////////////////////////
            // 	2) injection for extraDoFs:
            ////////////////////////////////////////////////////////////
            std::vector<DoFIndex> transIndCoarse = m_spParticleHandlerGlobal->get_transInd(coarseIndex, p);
            std::vector<DoFIndex> rotIndCoarse = m_spParticleHandlerGlobal->get_rotInd(coarseIndex, p);
            
            std::vector<DoFIndex> transIndFine = m_spParticleHandlerGlobal->get_transInd(fineIndex, p);
            std::vector<DoFIndex> rotIndFine = m_spParticleHandlerGlobal->get_rotInd(fineIndex, p);
            
            for ( int d = 0; d < dim; ++d )
            {
                DoFRef(dCoarse, transIndCoarse[d])	= DoFRef(dFine, transIndFine[d]);
                DoFRef(dCoarse, rotIndCoarse[d])	= DoFRef(dFine, rotIndFine[d]);
            }
            
        }
        
    }
    
    
    template <typename TDomain, typename TAlgebra>
    size_t ParticleTransfer<TDomain, TAlgebra>::
    get_vertex_index(Vertex* vrt, GridObject* elem)
    {
        std::vector<Vertex*> vVertex;
        CollectVertices(vVertex, *m_spParticleHandlerGlobal->m_spMG, elem);
        
        for(size_t i = 0; i < vVertex.size(); ++i)
            if ( vrt == vVertex[i])
                return i;
        
        UG_THROW("in CutElementHandlerImmersed::get_vertex_index: no index found!\n");
    }
    
    
    template <typename TDomain, typename TAlgebra>
    size_t ParticleTransfer<TDomain, TAlgebra>::
    get_vertex_mode(Vertex* coarseVrt0, Vertex* coarseVrt1, Vertex* childVrt )
    {
        //////////////////////////////////////////////////////////
        // alpha = (0.5, 0.5) => mode = 0
        // alpha = (0.0, 0.0) => mode = 1
        // alpha = (1.0, 0.0) => mode = 2
        // alpha = new weight => mode = 3
        //////////////////////////////////////////////////////////
        
        bool bCoarseVrt0_isInside = false;
        bool bCoarseVrt1_isInside = false;
        bool bChildVrt_isInside = false;
        
        
        if ( !m_spParticleHandlerGlobal->m_spOutsideMarker->is_marked(coarseVrt0) && !m_spParticleHandlerGlobal->m_spFlatTopVrtMarker->is_marked(coarseVrt0))
            bCoarseVrt0_isInside = true;
        if ( !m_spParticleHandlerGlobal->m_spOutsideMarker->is_marked(coarseVrt1) && !m_spParticleHandlerGlobal->m_spFlatTopVrtMarker->is_marked(coarseVrt1))
            bCoarseVrt1_isInside = true;
        if ( !m_spParticleHandlerGlobal->m_spOutsideMarker->is_marked(childVrt) && !m_spParticleHandlerGlobal->m_spFlatTopVrtMarker->is_marked(childVrt))
            bChildVrt_isInside = true;
        
        int prtIndex = -1;
        bool isOutside_vrt0 	= m_spParticleHandlerGlobal->is_outsideFluid(prtIndex, coarseVrt0);
        bool isOutside_vrt1 	= m_spParticleHandlerGlobal->is_outsideFluid(prtIndex, coarseVrt1);
        
        if ( m_spParticleHandlerGlobal->m_spOutsideMarker->is_marked(coarseVrt0) && !isOutside_vrt0 )
            UG_THROW("inconsistent 0!\n");
        if ( !m_spParticleHandlerGlobal->m_spOutsideMarker->is_marked(coarseVrt0) && isOutside_vrt0 )
            UG_THROW("inconsistent 1!\n");
        if ( m_spParticleHandlerGlobal->m_spOutsideMarker->is_marked(coarseVrt1) && !isOutside_vrt1 )
            UG_THROW("inconsistent 2!\n");
        if ( !m_spParticleHandlerGlobal->m_spOutsideMarker->is_marked(coarseVrt1) && isOutside_vrt1 )
            UG_THROW("inconsistent 3!\n");
        
        // inside edge => alpha = (0.5, 0.5)
        if ( bCoarseVrt0_isInside && bCoarseVrt1_isInside)
            return 0;
        // flat top edge => alpha = (0.5, 0.5) OR (1.0, 0.0) OR alphaNew!
        // both flat top => alpha = (0.5,0.5) for pressure, for vel: reset later during adjust_prolongation:
        else if ( m_spParticleHandlerGlobal->m_spFlatTopVrtMarker->is_marked(coarseVrt0) && m_spParticleHandlerGlobal->m_spFlatTopVrtMarker->is_marked(coarseVrt1) )
            return 0;
        // outside edge => alpha = (0.0, 0.0)
        else if ( m_spParticleHandlerGlobal->m_spOutsideMarker->is_marked(coarseVrt0) && m_spParticleHandlerGlobal->m_spOutsideMarker->is_marked(coarseVrt1) )
            return 1;
        else if ( !m_spParticleHandlerGlobal->m_spFlatTopVrtMarker->is_marked(coarseVrt0) && !m_spParticleHandlerGlobal->m_spFlatTopVrtMarker->is_marked(coarseVrt1) )
        {UG_THROW("error in 'particle_transfer_impl:get_vertex_mode': one vertex must be flat top!\n");}
        
        // check: exactly one vertex must be inside:
        else if ( !bCoarseVrt0_isInside && !bCoarseVrt1_isInside )
        {UG_THROW("error in 'particle_transfer_impl:get_vertex_mode': one vertex must be inside!\n");}
        // if one vertex is inside, BUT the FlatTopVertex is also ON interface: (0.5, 0.5) !
        else if ( m_spParticleHandlerGlobal->m_spNearInterfaceVrtMarker->is_marked(coarseVrt0) || m_spParticleHandlerGlobal->m_spNearInterfaceVrtMarker->is_marked(coarseVrt1) )
        {
            /*	std::vector<DoFIndex> vParentDoF;
             coarseDD.inner_dof_indices(coarseVrt0, 0, vParentDoF);
             if ( m_spParticleHandlerGlobal->m_spNearInterfaceVrtMarker.is_marked(coarseVrt0) )
             UG_THROW("oho 0: vParentDoF = " << vParentDoF[0] << "\n");
             coarseDD.inner_dof_indices(coarseVrt1, 0, vParentDoF);
             if ( m_spParticleHandlerGlobal->m_spNearInterfaceVrtMarker.is_marked(coarseVrt1) )
             UG_THROW("oho 1: vParentDoF = " << vParentDoF[0] << "\n");*/
            if ( m_spParticleHandlerGlobal->m_spNearInterfaceVrtMarker->is_marked(coarseVrt0) )
                UG_LOG("0 pos: " << m_spParticleHandlerGlobal->m_aaPos[coarseVrt0] << "\n");
            
            if ( m_spParticleHandlerGlobal->m_spNearInterfaceVrtMarker->is_marked(coarseVrt1) )
                UG_LOG("1: pos: " << m_spParticleHandlerGlobal->m_aaPos[coarseVrt1] << "\n");
            return 0;
        }
        // childVertex lies inside => compute new weigthing for prolongation!
        else if ( bChildVrt_isInside )
            return 3;
        // childVertex lies outside AND on flat top edge => constant prologation for pressure and reset during
        // adjust_prolongation for velocity:
        else
            return 2;
        
        UG_THROW("error in 'particle_transfer_impl:get_vertex_mode': no case found!\n");
        
    }

    
//////////////////////////////////////////////////////////
/// Std-Transfer methods
//////////////////////////////////////////////////////////

template <typename TDomain, typename TAlgebra>
void ParticleTransfer<TDomain, TAlgebra>::
assemble_prolongation_p1(matrix_type& P,
                         const DoFDistribution& fineDD,
                         const DoFDistribution& coarseDD)
{
	PROFILE_FUNC_GROUP("gmg");
// 	allow only lagrange P1 functions
	for(size_t fct = 0; fct < fineDD.num_fct(); ++fct)
		if(fineDD.lfeid(fct).type() != LFEID::LAGRANGE ||
			fineDD.lfeid(fct).order() != 1)
			UG_THROW("AssembleStdProlongationForP1Lagrange:"
				"Interpolation only implemented for Lagrange P1 functions.");

//  resize matrix
	P.resize_and_clear(fineDD.num_indices(), coarseDD.num_indices());

    UG_LOG("begin assemble_prolongation_p1:\n");

//  iterators
	const MultiGrid& mg = *coarseDD.multi_grid();
	typedef DoFDistribution::traits<Vertex>::const_iterator const_iterator;
	const_iterator iter, iterBegin, iterEnd;

//  loop subsets on fine level
	std::vector<size_t> vParentIndex, vChildIndex;
	std::vector<DoFIndex> vParentDoF, vChildDoF;
	for(int si = 0; si < fineDD.num_subsets(); ++si)
	{
		iterBegin = fineDD.template begin<Vertex>(si);
		iterEnd = fineDD.template end<Vertex>(si);

	//  loop vertices for fine level subset
		for(iter = iterBegin; iter != iterEnd; ++iter)
		{
		//	get element
			Vertex* child = *iter;

		//  get father
			GridObject* parent = mg.get_parent(child);

		//	check if child contained in coarseDD. This should always be false
		//	for a GridLevel::LEVEL, but might be the case for GridLevel::SURFACE
		//	and an adaptive grid-part used by both dds. In such a case we can
		//	simply set identity.
			if(coarseDD.is_contained(child)){
			//	get indices
				coarseDD.inner_algebra_indices(child, vParentIndex);
				fineDD.inner_algebra_indices(child, vChildIndex);
				UG_ASSERT(vParentIndex.size() == vChildIndex.size(), "Size mismatch");

			//	set identity
				for(size_t i = 0; i < vParentIndex.size(); ++i)
					P(vChildIndex[i], vParentIndex[i]) = 1.0;

			//	this child is perfectly handled
				continue;
			}
			else{
			//	check if parent exists (this should always be the case, except in
			//	the case that 'child' is a v-slave)
				if(!parent) continue;

				if(!coarseDD.is_contained(parent)){
					UG_THROW("ParticleTransfer: Parent element \n"
							<< ElementDebugInfo(mg, parent) <<
							"is not contained in coarse-dd nor the child element\n"
							<< ElementDebugInfo(mg, child) <<
							" in the coarse-dd. This should not happen.")
				}
			}

		//	type of father
			const ReferenceObjectID roid = parent->reference_object_id();

		//	loop all components
			for(size_t fct = 0; fct < fineDD.num_fct(); fct++)
			{
			//	check that fct defined on subset
				if(!fineDD.is_def_in_subset(fct, si)) continue;

			//  get global indices
				fineDD.inner_dof_indices(child, fct, vChildDoF);

			//	detect type of father
				switch(roid)
				{
					case ROID_VERTEX:
					{
						Vertex* vrt = dynamic_cast<Vertex*>(parent);
						coarseDD.inner_dof_indices(vrt, fct, vParentDoF);
						DoFRef(P, vChildDoF[0], vParentDoF[0]) = 1.0;
					}
					break;
					case ROID_EDGE:
                    {
                        Edge* edge = dynamic_cast<Edge*>(parent);
                        std::vector<number> alpha;
                        alpha.resize(2,0.5);
                        
                        bool newWeights1 = false;
                        bool newWeights2 = false;
                        MathVector<dim> intersectionPnt;
                        int prtIndex = -1;
                        
                        bool isOutside_vrt0 	= m_spParticleHandlerGlobal->is_outsideFluid(prtIndex, edge->vertex(0));
                        bool isOutside_vrt1 	= m_spParticleHandlerGlobal->is_outsideFluid(prtIndex, edge->vertex(1));
                        size_t vrtInd0 = get_vertex_index(edge->vertex(0), parent);
                        size_t vrtInd1 = get_vertex_index(edge->vertex(1), parent);
                        
                        /*
                         // output stuff:
                         std::vector<DoFIndex> vParentDoF1;
                         coarseDD.inner_dof_indices(edge->vertex(0), 0, vParentDoF1);
                         UG_LOG("oho 0: vParentDoF1 = " << vParentDoF1[0] << "\n");
                         std::vector<DoFIndex> vParentDoF2;
                         coarseDD.inner_dof_indices(edge->vertex(1), 0, vParentDoF2);
                         UG_LOG("oho 1: vParentDoF2 = " << vParentDoF2[0] << "\n");
                         */
                        size_t vertexMode = get_vertex_mode(edge->vertex(0), edge->vertex(1), child);
                        
                        //bool isOutside_childDoF = m_spParticleHandlerGlobal->is_outsideFluid(prtIndex, child);
                        
                        
                        if ( vertexMode == 1 )
                        {
                            alpha[1] = alpha[0] = 0.0;
                        }
                        else if ( vertexMode == 2 )
                        {
                            if ( isOutside_vrt0 && isOutside_vrt1 )
                                UG_THROW("hmmmm\n");
                            if ( isOutside_vrt0 )
                            {alpha[1] = 0.0; alpha[0] = 1.0;}
                            else
                            {alpha[1] = 1.0; alpha[0] = 0.0;}
                        }
                        // only new weighting, if childDoF insideFluid => get boolian isInside_childDoF
                        else if ( vertexMode == 3 ) //!isOutside_childDoF )
                        {
                            UG_LOG("vertexMode == 3\n");
                            
                            // case1: vrt0 = insideFluid && vrt1 = outsideFluid:
                            if ( !isOutside_vrt0 && isOutside_vrt1 ){
                                
                                coarseDD.inner_dof_indices(edge->vertex(1), fct, vParentDoF);
                                if ( m_spParticleHandlerGlobal->is_nearInterfaceVertex(edge->vertex(1), vrtInd1) )
                                    UG_THROW("1 oho, vParentDoF = " << vParentDoF[0] << "\n");
                                
                                
                                m_spParticleHandlerGlobal->m_spInterfaceProvider->get_intersection_point(intersectionPnt, m_spParticleHandlerGlobal->m_aaPos[edge->vertex(0)], m_spParticleHandlerGlobal->m_aaPos[edge->vertex(1)], prtIndex, alpha);
                                newWeights1 = true;
                                coarseDD.inner_dof_indices(edge->vertex(0), fct, vParentDoF);
                                UG_LOG("1: vertex[0] = " << vParentDoF[0] << "\n");
                                UG_LOG("1: alpha = " << alpha[0] << ", " << alpha[1] << "\n");
                                UG_LOG("VORHER: alpha[i] = " << alpha[0] << ", " << alpha[1] << "\n");
                                
                                alpha[0] = 0.5/alpha[0];
                                alpha[1] = 1.0 - alpha[0];
                                //	UG_THROW("NACHHER: alpha[i] = " << alpha[0] << ", " << alpha[1] << "\n");
                            }
                            // case1: vrt0 = outsideFluid && vrt1 = insideFluid:
                            else if ( isOutside_vrt0 && !isOutside_vrt1 ){
                                
                                coarseDD.inner_dof_indices(edge->vertex(0), fct, vParentDoF);
                                if ( m_spParticleHandlerGlobal->is_nearInterfaceVertex(edge->vertex(0), vrtInd0) )
                                    UG_THROW("2 oho, vParentDoF = " << vParentDoF[0] << "\n");
                                
                                m_spParticleHandlerGlobal->m_spInterfaceProvider->get_intersection_point(intersectionPnt, m_spParticleHandlerGlobal->m_aaPos[edge->vertex(1)], m_spParticleHandlerGlobal->m_aaPos[edge->vertex(0)], prtIndex, alpha);
                                newWeights2 = true;
                                coarseDD.inner_dof_indices(edge->vertex(1), fct, vParentDoF);
                                UG_LOG("2: vertex[0] = " << vParentDoF[0] << "\n");
                                UG_LOG("2: alpha[i] = " << alpha[0] << ", " << alpha[1] << "\n");
                                UG_LOG("VORHER: alpha[i] = " << alpha[0] << ", " << alpha[1] << "\n");
                                
                                // switch also ordering of original alpha!!!
                                alpha[1] = 0.5/alpha[0];
                                alpha[0] = 1.0 - alpha[1];
                                //	UG_THROW("NACHHER: alpha[i] = " << alpha[0] << ", " << alpha[1] << "\n");
                            }
                            else
                                UG_THROW("vertexMode == 3: hmmm\n");
                            
                        } // end if (vertexMode == 3 )
                        else if ( vertexMode != 0 )
                            UG_THROW("in particle_transfer_impl.h: no valid vertexMode computed: vertexMode = " << vertexMode << "\n");
                        
                        
                        for(int i = 0; i < 2; ++i)
                        {
                            Edge* edge = dynamic_cast<Edge*>(parent);
                            coarseDD.inner_dof_indices(edge->vertex(i), fct, vParentDoF);
                            
                            // use prtIndices for prolongation if 'isOutside_parentDoF':
                            bool isOutside_parentDoF = false;
                            
                            if ( newWeights2 || newWeights1 )
                                isOutside_parentDoF = m_spParticleHandlerGlobal->is_outsideFluid(prtIndex, edge->vertex(i));
                            
                            if ( isOutside_parentDoF && fct != dim )
                            {
                                const int coarseIndex = m_spParticleHandlerGlobal->get_Index(coarseDD.grid_level());
                                
                                std::vector<DoFIndex> transIndCoarse = m_spParticleHandlerGlobal->get_transInd(coarseIndex, prtIndex);
                                
                                std::vector<DoFIndex> rotIndCoarse = m_spParticleHandlerGlobal->get_rotInd(coarseIndex, prtIndex);
                                
                                /*					UG_LOG("fct = " << fct << "\n");
                                 UG_LOG("transIndCoarse[0] = " << transIndCoarse[0] << "\n");
                                 UG_LOG("transIndCoarse[1] = " << transIndCoarse[1] << "\n");
                                 UG_LOG("transIndCoarse[fct] = " << transIndCoarse[fct] << "\n");
                                 */
                                
                                // use of 'transIndCoarse[fct]', since 'normally' it is:
                                //		fct == vParentDoF[0][1]
                                DoFRef(P, vChildDoF[0], transIndCoarse[fct]) = alpha[i];
                                DoFRef(P, vChildDoF[0], rotIndCoarse[fct]) = alpha[i];
                                // iff NOT 'normal' case: EXIT!!
                                if ( vParentDoF[0][1] != fct )
                                    UG_THROW("no natural index mapping => use of incex 'fct' not valid for next operation!\n");
                                
                            }
                            else
                                DoFRef(P, vChildDoF[0], vParentDoF[0]) = alpha[i];
                            if ( 0 ) //newWeights2 || newWeights1 )
                            {
                                UG_LOG("intersectionPnt = " << intersectionPnt << "\n");
                                UG_LOG("i = " << i << "vParentDoF[0] = " << vParentDoF[0] << "\n");
                                UG_LOG("i = " << i << "vChildDoF[0] = " << vChildDoF[0] << "\n");
                                UG_LOG("alpha[i] = " << alpha[i] << "\n");
                                
                                UG_LOG("DoFRef(P, vChildDoF[0], vParentDoF[0]) = " << DoFRef(P, vChildDoF[0], vParentDoF[0]) << "\n");
                                
                            }
                            
                        }
                    }
					break;
					case ROID_QUADRILATERAL:
					for(int i = 0; i < 4; ++i)
					{
						Face* face = dynamic_cast<Face*>(parent);
						coarseDD.inner_dof_indices(face->vertex(i), fct, vParentDoF);
						DoFRef(P, vChildDoF[0], vParentDoF[0]) = 0.25;
					}
					break;
					case ROID_HEXAHEDRON:
					for(int i = 0; i < 8; ++i)
					{
						Volume* hexaeder = dynamic_cast<Volume*>(parent);
						coarseDD.inner_dof_indices(hexaeder->vertex(i), fct, vParentDoF);
						DoFRef(P, vChildDoF[0], vParentDoF[0]) = 0.125;
					}
					break;
					default: UG_THROW("AssembleStdProlongationForP1Lagrange: Element father"
									 " is of unsupported type "<< roid << " for "
									 << ElementDebugInfo(mg, child) << ".");
				}
			}
		}
	}
}
/*
template <typename TDomain>
void ProjectGlobalPositionToElem(std::vector<MathVector<TDomain::dim> >& vGlobPos,
                                 GridObject* parent, const TDomain& domain)
{
	const int parentDim = parent->base_object_id();

	// vertex and full dim parent must match
	if(parentDim == 0 || parentDim == TDomain::dim)
		return;

//	get the vertices
	std::vector<MathVector<TDomain::dim> > vCornerCoord;
	switch(parentDim)
	{
		case EDGE:
		{
			CollectCornerCoordinates(vCornerCoord, *static_cast<Edge*>(parent), domain, true);
			MathVector<TDomain::dim> dir;
			VecSubtract(dir, vCornerCoord[1], vCornerCoord[0]);
			for(size_t p = 0; p < vGlobPos.size(); ++p){
				ProjectPointToRay(vGlobPos[p], vGlobPos[p], vCornerCoord[0], dir);
			}
		}
		break;
		case FACE:
		{
			CollectCornerCoordinates(vCornerCoord, *static_cast<Face*>(parent), domain, true);
			MathVector<TDomain::dim> normal;
			MathVector<TDomain::dim> a, b;
			VecSubtract(a, vCornerCoord[1], vCornerCoord[0]);
			VecSubtract(b, vCornerCoord[2], vCornerCoord[0]);
			VecCross(normal, a,b);

			for(size_t p = 0; p < vGlobPos.size(); ++p){
				ProjectPointToPlane(vGlobPos[p], vGlobPos[p], vCornerCoord[0], normal);
			}
		}
		break;
		default: UG_THROW( "Base Object type not found.");
	}
}
*/


template <typename TDomain, typename TAlgebra>
template <typename TChild>
void ParticleTransfer<TDomain, TAlgebra>::
assemble_prolongation(matrix_type& P,
                      const DoFDistribution& fineDD,
                      const DoFDistribution& coarseDD,
                      ConstSmartPtr<TDomain> spDomain)
{
	PROFILE_FUNC_GROUP("gmg");

//  iterators
	MultiGrid& mg = *const_cast<MultiGrid*>(coarseDD.multi_grid().get());
	typedef typename DoFDistribution::traits<TChild>::const_iterator const_iterator;
	const_iterator iter, iterBegin, iterEnd;

//  loop subsets on coarse level
	std::vector<DoFIndex> vParentDoF, vChildDoF;
	std::vector<size_t> vParentIndex, vChildIndex;
	for(int si = 0; si < fineDD.num_subsets(); ++si)
	{
		iterBegin = fineDD.template begin<TChild>(si);
		iterEnd = fineDD.template end<TChild>(si);

	//	check, which cmps to consider on this subset
		std::vector<LFEID> vLFEID;
		std::vector<size_t> vFct;
		for(size_t fct = 0; fct < fineDD.num_fct(); ++fct){
			if(fineDD.max_fct_dofs(fct, TChild::dim, si) == 0) continue;
			vFct.push_back(fct);
			vLFEID.push_back(fineDD.lfeid(fct));
		}
		if(vFct.empty()) continue;

	//  loop elems on coarse level for subset
		for(iter = iterBegin; iter != iterEnd; ++iter)
		{
		//	get child
			TChild* child = *iter;

		//	get parent
			GridObject* parent = mg.get_parent(child);

		//	check if child contained in coarseDD. This should always be false
		//	for a GridLevel::LEVEL, but might be the case for GridLevel::SURFACE
		//	and an adaptive grid-part used by both dds. In such a case we can
		//	simply set identity.
			if(coarseDD.is_contained(child)){
			//	get indices
				coarseDD.inner_algebra_indices(child, vParentIndex);
				fineDD.inner_algebra_indices(child, vChildIndex);
				UG_ASSERT(vParentIndex.size() == vChildIndex.size(), "Size mismatch");

			//	set identity
				for(size_t i = 0; i < vParentIndex.size(); ++i)
					P(vChildIndex[i], vParentIndex[i]) = 1.0;

			//	this child is perfectly handled
				continue;
			}
			else{

			//	check if parent exists (this should always be the case, except in
			//	the case that 'child' is a v-slave)
				if(!parent) continue;

				if(!coarseDD.is_contained(parent)){
					UG_THROW("ParticleTransfer: A parent element is not contained in "
							" coarse-dd nor the child element in the coarse-dd. "
							"This should not happen.")
				}
			}

		//	loop all components
			for(size_t f = 0; f < vFct.size(); f++)
			{
			//	get comp and lfeid
				const size_t fct = vFct[f];
				const LFEID& lfeID = vLFEID[f];

			//  get global indices
				fineDD.inner_dof_indices(child, fct, vChildDoF);

			//	switch space type
				switch(lfeID.type())
				{
					case LFEID::PIECEWISE_CONSTANT:
					{
						coarseDD.dof_indices(parent, fct, vParentDoF);
						UG_ASSERT(vChildDoF.size() == 1, "Must be one.");
						UG_ASSERT(vParentDoF.size() == 1, "Must be one.");

						DoFRef(P, vChildDoF[0], vParentDoF[0]) =  1.0;
					}
					break;

					case LFEID::CROUZEIX_RAVIART:
					{
					//	get dimension of parent
						const int parentDim = parent->base_object_id();
						std::vector<GridObject*> vParent;

					//	check if to interpolate from neighbor elems
						if(parentDim == lfeID.dim()){
							// case: Side inner to parent. --> Parent fine.
							vParent.push_back(parent);
						} else if(parentDim == lfeID.dim() - 1){
							// case: parent is Side. --> Get neighbor elems
							typedef typename TChild::sideof TElem;
							std::vector<TElem*> vElem;
							coarseDD.collect_associated(vElem, parent);
							for(size_t p = 0; p < vElem.size(); ++p)
								vParent.push_back(vElem[p]);

						} else {
							UG_THROW("ParticleTransfer: For CR parent must be full-dim "
									"elem or a side (dim-1). But has dim: "<<parentDim);
						}


					//	global positions of fine dofs
						std::vector<MathVector<TDomain::dim> > vDoFPos;
						InnerDoFPosition(vDoFPos, child, *spDomain, lfeID);

					//	loop contributions from parents
						for(size_t i = 0; i < vParent.size(); ++i)
						{
						//	get coarse indices
							coarseDD.dof_indices(vParent[i], fct, vParentDoF);

						//	get shapes at global positions
							std::vector<std::vector<number> > vvShape;
							ShapesAtGlobalPosition(vvShape, vDoFPos, vParent[i], *spDomain, lfeID);

						//	add restriction
							for(size_t ip = 0; ip < vvShape.size(); ++ip)
								for(size_t sh = 0; sh < vvShape[ip].size(); ++sh)
									DoFRef(P, vChildDoF[ip], vParentDoF[sh]) +=
											(1./vParent.size()) * vvShape[ip][sh];
						}
					}
					break;

					case LFEID::LAGRANGE:
					{
					//	get coarse indices
						coarseDD.dof_indices(parent, fct, vParentDoF);

					//	global positions of child dofs
						std::vector<MathVector<TDomain::dim> > vDoFPos;
						InnerDoFPosition(vDoFPos, child, *spDomain, lfeID);

					//	project
					//	ProjectGlobalPositionToElem(vDoFPos, parent, *spDomain);

					//	get shapes at global positions
						std::vector<std::vector<number> > vvShape;
						ShapesAtGlobalPosition(vvShape, vDoFPos, parent, *spDomain, lfeID);

					//	set restriction
						for(size_t ip = 0; ip < vvShape.size(); ++ip)
							for(size_t sh = 0; sh < vvShape[ip].size(); ++sh)
								DoFRef(P, vChildDoF[ip], vParentDoF[sh]) = vvShape[ip][sh];
					}
					break;

					default:
						UG_THROW("ParticleTransfer: Local-Finite-Element: "<<lfeID<<
						         " is not supported by this Transfer.")

				} // end LFEID-switch
			} // end fct - cmps
		} // end fine - elements
	} // end subset
}

template <typename TDomain, typename TAlgebra>
void ParticleTransfer<TDomain, TAlgebra>::
assemble_prolongation(matrix_type& P,
                      const DoFDistribution& fineDD,
                      const DoFDistribution& coarseDD,
                      ConstSmartPtr<TDomain> spDomain)
{
	//  resize matrix
	P.resize_and_clear(fineDD.num_indices(), coarseDD.num_indices());

	// loop all base types carrying indices on fine elems
	if(fineDD.max_dofs(VERTEX)) assemble_prolongation<Vertex>(P, fineDD, coarseDD, spDomain);
	if(fineDD.max_dofs(EDGE)) assemble_prolongation<Edge>(P, fineDD, coarseDD, spDomain);
	if(fineDD.max_dofs(FACE)) assemble_prolongation<Face>(P, fineDD, coarseDD, spDomain);
	if(fineDD.max_dofs(VOLUME)) assemble_prolongation<Volume>(P, fineDD, coarseDD, spDomain);
}


template <typename TDomain, typename TAlgebra>
template <typename TChild>
void ParticleTransfer<TDomain, TAlgebra>::
assemble_restriction(matrix_type& R,
                     const DoFDistribution& coarseDD,
                     const DoFDistribution& fineDD,
                     ConstSmartPtr<TDomain> spDomain)
{
	PROFILE_FUNC_GROUP("gmg");

//  iterators
	MultiGrid& mg = *const_cast<MultiGrid*>(coarseDD.multi_grid().get());
	typedef typename DoFDistribution::traits<TChild>::const_iterator const_iterator;
	const_iterator iter, iterBegin, iterEnd;

//  loop subsets on coarse level
	std::vector<DoFIndex> vParentDoF, vChildDoF;
	std::vector<size_t> vParentIndex, vChildIndex;
	for(int si = 0; si < fineDD.num_subsets(); ++si)
	{
		iterBegin = fineDD.template begin<TChild>(si);
		iterEnd = fineDD.template end<TChild>(si);

	//	check, which cmps to consider on this subset
		std::vector<LFEID> vLFEID;
		std::vector<size_t> vFct;
		for(size_t fct = 0; fct < fineDD.num_fct(); ++fct){
			if(fineDD.max_fct_dofs(fct, TChild::dim, si) == 0) continue;
			vFct.push_back(fct);
			vLFEID.push_back(fineDD.lfeid(fct));
		}
		if(vFct.empty()) continue;

	//  loop elems on coarse level for subset
		for(iter = iterBegin; iter != iterEnd; ++iter)
		{
		//	get child
			TChild* child = *iter;

		//	get parent
			GridObject* parent = mg.get_parent(child);

		//	check if child contained in coarseDD. This should always be false
		//	for a GridLevel::LEVEL, but might be the case for GridLevel::SURFACE
		//	and an adaptive grid-part used by both dds. In such a case we can
		//	simply set identity.
			if(coarseDD.is_contained(child)){
			//	get indices
				coarseDD.inner_algebra_indices(child, vParentIndex);
				fineDD.inner_algebra_indices(child, vChildIndex);
				UG_ASSERT(vParentIndex.size() == vChildIndex.size(), "Size mismatch");

			//	set identity
				for(size_t i = 0; i < vParentIndex.size(); ++i)
					R(vParentIndex[i], vChildIndex[i]) = 1.0;

			//	this child is perfectly handled
				continue;
			}
			else{

			//	check if parent exists (this should always be the case, except in
			//	the case that 'child' is a v-slave)
				if(!parent) continue;

				if(!coarseDD.is_contained(parent)){
					UG_THROW("ParticleTransfer: A parent element is not contained in "
							" coarse-dd nor the child element in the coarse-dd. "
							"This should not happen.")
				}
			}

		//	loop all components
			for(size_t f = 0; f < vFct.size(); f++)
			{
			//	get comp and lfeid
				const size_t fct = vFct[f];
				const LFEID& lfeID = vLFEID[f];

			//  get global indices
				fineDD.inner_dof_indices(child, fct, vChildDoF);

			//	switch space type
				switch(lfeID.type())
				{
					case LFEID::PIECEWISE_CONSTANT:
					{
						coarseDD.dof_indices(parent, fct, vParentDoF);
						UG_ASSERT(vChildDoF.size() == 1, "Must be one.");
						UG_ASSERT(vParentDoF.size() == 1, "Must be one.");

						DoFRef(R, vParentDoF[0], vChildDoF[0]) =  1.0;
					}
					break;

					case LFEID::CROUZEIX_RAVIART:
					{
					//	get dimension of parent
						const int parentDim = parent->base_object_id();
						std::vector<GridObject*> vParent;

					//	check if to interpolate from neighbor elems
						if(parentDim == lfeID.dim()){
							// case: Side inner to parent. --> Parent fine.
							vParent.push_back(parent);
						} else if(parentDim == lfeID.dim() - 1){
							// case: parent is Side. --> Get neighbor elems
							typedef typename TChild::sideof TElem;
							std::vector<TElem*> vElem;
							coarseDD.collect_associated(vElem, parent);
							for(size_t p = 0; p < vElem.size(); ++p){
							//	NOTE: This is not the transposed of the prolongation
							//		  in adaptive case, since we only restrict to
							//		  covered parts.
								if(mg.num_children<TElem>(vElem[p]) > 0)
									vParent.push_back(vElem[p]);
							}

						} else {
							UG_THROW("ParticleTransfer: For CR parent must be full-dim "
									"elem or a side (dim-1). But has dim: "<<parentDim);
						}


					//	global positions of fine dofs
						std::vector<MathVector<TDomain::dim> > vDoFPos;
						InnerDoFPosition(vDoFPos, child, *spDomain, lfeID);

					//	loop contributions from parents
						for(size_t i = 0; i < vParent.size(); ++i)
						{
						//	get coarse indices
							coarseDD.dof_indices(vParent[i], fct, vParentDoF);

						//	get shapes at global positions
							std::vector<std::vector<number> > vvShape;
							ShapesAtGlobalPosition(vvShape, vDoFPos, vParent[i], *spDomain, lfeID);

						//	add restriction
							for(size_t ip = 0; ip < vvShape.size(); ++ip)
								for(size_t sh = 0; sh < vvShape[ip].size(); ++sh)
									DoFRef(R, vParentDoF[sh], vChildDoF[ip]) +=
											(1./vParent.size()) * vvShape[ip][sh];
						}
					}
					break;

					case LFEID::LAGRANGE:
					{
					//	get coarse indices
						coarseDD.dof_indices(parent, fct, vParentDoF);

					//	global positions of child dofs
						std::vector<MathVector<TDomain::dim> > vDoFPos;
						InnerDoFPosition(vDoFPos, child, *spDomain, lfeID);

					//	get shapes at global positions
						std::vector<std::vector<number> > vvShape;
						ShapesAtGlobalPosition(vvShape, vDoFPos, parent, *spDomain, lfeID);

					//	set restriction
						for(size_t ip = 0; ip < vvShape.size(); ++ip)
                        {
                            
                            for(size_t sh = 0; sh < vvShape[ip].size(); ++sh)
                            {
                                // vertex on coarse grid which is FT can possibly have no child => injection fails!
                                // 	=> project average of neighbouring pressure - which is naturally zero in irrelevant vertices:-)
                                
                                //	global positions of child dofs
                                std::vector<MathVector<TDomain::dim> > vDoFPosParent;
                                InnerDoFPosition(vDoFPosParent, parent, *spDomain, lfeID);
                                if ( vDoFPosParent.size() == 2 )
                                    UG_THROW("vDoFPosParent.size() == 2, but = " << vDoFPosParent.size() << "\n");
                                
                                switch(parent->base_object_id())
                                {
                                    case VERTEX:
                                    {
                                        const char* filename = "parent_vertex";
                                        std::string name(filename);
                                        char ext[50]; sprintf(ext, ".txt");
                                        name.append(ext);
                                        FILE* outputFile = fopen(name.c_str(), "a");
                                        //	create multi index
                                        std::vector<DoFIndex>  vInd;  // ToDo: hier std::vector??
                                        //	get multi indices
                                        Vertex* vrt = static_cast<Vertex*>(parent);
                                        if(coarseDD.dof_indices(vrt, fct, vInd) != 1)
                                            UG_THROW("Only one index expected.");
                                        if ( m_spParticleHandlerGlobal->is_FTVertex(vrt) )
                                        {
                                            fprintf(outputFile,"ParentPos: %e \t %e (index = %lu)\n", vDoFPosParent[0][0], vDoFPosParent[0][1], vInd[0][0]);
                                            DoFRef(R, vParentDoF[sh], vChildDoF[ip]) = vvShape[ip][sh];
                                        }
                                        fclose(outputFile);
                                        
                                    }
                                        break;
                                    case EDGE:
                                    {
                                        
                                        const char* filename = "parent_edge";
                                        std::string name(filename);
                                        char ext[50]; sprintf(ext, ".txt");
                                        name.append(ext);
                                        FILE* outputFile = fopen(name.c_str(), "a");
                                        
                                        Edge* edge = dynamic_cast<Edge*>(parent);
                                        
                                        for ( size_t i = 0; i < 2; ++i )
                                        {
                                            //	create multi index
                                            std::vector<DoFIndex>  vInd;  // ToDo: hier std::vector??
                                            //	get multi indices
                                            if(coarseDD.dof_indices(edge->vertex(i), fct, vInd) != 1)
                                                UG_THROW("Only one index expected.");
                                            if ( m_spParticleHandlerGlobal->is_FTVertex(edge->vertex(i)) && vParentDoF[sh][0] == vInd[0][0] )
                                            {
                                                fprintf(outputFile,"vParentDoF[%lu], vInd[0]: %lu \t %lu \n", i, vParentDoF[sh][0], vInd[0][0]);
                                                DoFRef(R, vParentDoF[sh], vChildDoF[ip]) = vvShape[ip][sh];
                                            }
                                        }
                                        fclose(outputFile);
                                        
                                        
                                    }
                                        break;
                                    case FACE: UG_THROW("Parent type wrong: Is FACE, but should be always VERTEX for triangular grids!"); break;
                                    case VOLUME: UG_THROW("Parent type wrong: Is VOLUME, but should be always VERTEX for triangular grids!"); break;
                                    default: UG_THROW("Base Object type not found.");
                                }
                                //	if ( vParentDoF[sh] )
                                
                            }
                        }
					}
					break;

					default:
						UG_THROW("ParticleTransfer: Local-Finite-Element: "<<lfeID<<
						         " is not supported by this Transfer.")

				} // end LFEID-switch
			} // end fct - cmps
		} // end fine - elements
	} // end subset
}

template <typename TDomain, typename TAlgebra>
void ParticleTransfer<TDomain, TAlgebra>::
assemble_restriction(matrix_type& R,
                     const DoFDistribution& coarseDD,
                     const DoFDistribution& fineDD,
                     ConstSmartPtr<TDomain> spDomain)
{
	//  resize matrix
	R.resize_and_clear(coarseDD.num_indices(), fineDD.num_indices());

	// loop all base types carrying indices on fine elems
	if(fineDD.max_dofs(VERTEX)) assemble_restriction<Vertex>(R, coarseDD, fineDD, spDomain);
	if(fineDD.max_dofs(EDGE)) assemble_restriction<Edge>(R, coarseDD, fineDD, spDomain);
	if(fineDD.max_dofs(FACE)) assemble_restriction<Face>(R, coarseDD, fineDD, spDomain);
	if(fineDD.max_dofs(VOLUME)) assemble_restriction<Volume>(R, coarseDD, fineDD, spDomain);
}


template <typename TDomain, typename TAlgebra>
SmartPtr<typename TAlgebra::matrix_type>
ParticleTransfer<TDomain, TAlgebra>::
prolongation(const GridLevel& fineGL, const GridLevel& coarseGL,
             ConstSmartPtr<ApproximationSpace<TDomain> > spApproxSpace)
{
    UG_LOG("begin prolongation:\n");

	if(fineGL.level() - coarseGL.level() != 1)
		UG_THROW("ParticleTransfer: Can only project between successive level, "
				"but fine = "<<fineGL<<", coarse = "<<coarseGL);

	if(fineGL.type() != coarseGL.type())
		UG_THROW("ParticleTransfer: Can only project between dof distributions of "
				"same type, but fine = "<<fineGL<<", coarse = "<<coarseGL);

	// remove old revisions
	remove_outdated(m_mProlongation, spApproxSpace->revision());

	// key of this restriction
	TransferKey key(coarseGL, fineGL, spApproxSpace->revision());

	// check if must be created
	if(m_mProlongation.find(key) == m_mProlongation.end())
	{
		SmartPtr<matrix_type> P =
				m_mProlongation[key] = SmartPtr<matrix_type>(new matrix_type);

		ConstSmartPtr<DoFDistribution> spCoarseDD = spApproxSpace->dof_distribution(coarseGL);
		ConstSmartPtr<DoFDistribution> spFineDD = spApproxSpace->dof_distribution(fineGL);

		bool P1LagrangeOnly = false;
		if(m_p1LagrangeOptimizationEnabled){
			P1LagrangeOnly = true;
			for(size_t fct = 0; fct < spApproxSpace->num_fct(); ++fct)
				if(spApproxSpace->lfeid(fct).type() != LFEID::LAGRANGE ||
					spApproxSpace->lfeid(fct).order() != 1)
					P1LagrangeOnly = false;
		}

		if(P1LagrangeOnly){
			assemble_prolongation_p1(*P, *spFineDD, *spCoarseDD);
		} else{
			assemble_prolongation(*P, *spFineDD, *spCoarseDD, spApproxSpace->domain());
		}

		for (int type = 1; type < CT_ALL; type = type << 1)
		{
			for (size_t i = 0; i < m_vConstraint.size(); ++i)
			{
				if (m_vConstraint[i]->type() & type)
					m_vConstraint[i]->adjust_prolongation(*P, spFineDD, spCoarseDD, type);
			}
		}

		#ifdef UG_PARALLEL
		P->set_storage_type(PST_CONSISTENT);
		#endif

		write_debug(*P, "P", fineGL, coarseGL);
	}

	return m_mProlongation[key];
}

template <typename TDomain, typename TAlgebra>
SmartPtr<typename TAlgebra::matrix_type>
ParticleTransfer<TDomain, TAlgebra>::
restriction(const GridLevel& coarseGL, const GridLevel& fineGL,
            ConstSmartPtr<ApproximationSpace<TDomain> > spApproxSpace)
{
	if(fineGL.level() - coarseGL.level() != 1)
		UG_THROW("ParticleTransfer: Can only project between successive level, "
				"but fine = "<<fineGL<<", coarse = "<<coarseGL);

	if(fineGL.type() != coarseGL.type())
		UG_THROW("ParticleTransfer: Can only project between dof distributions of "
				"same type, but fine = "<<fineGL<<", coarse = "<<coarseGL);

	// remove old revisions
	remove_outdated(m_mRestriction, spApproxSpace->revision());

	// key of this restriction
	TransferKey key(coarseGL, fineGL, spApproxSpace->revision());

	// check if must be created
	if(m_mRestriction.find(key) == m_mRestriction.end())
	{
		SmartPtr<matrix_type> R =
				m_mRestriction[key] = SmartPtr<matrix_type>(new matrix_type);

		ConstSmartPtr<DoFDistribution> spCoarseDD = spApproxSpace->dof_distribution(coarseGL);
		ConstSmartPtr<DoFDistribution> spFineDD = spApproxSpace->dof_distribution(fineGL);

		if(m_bUseTransposed)
			R->set_as_transpose_of(*prolongation(fineGL, coarseGL, spApproxSpace));
		else
        {
			assemble_restriction(*R, *spCoarseDD, *spFineDD, spApproxSpace->domain());
			ParticleAssembleInjectionForP1Lagrange1<TAlgebra>(*R,*spCoarseDD, *spFineDD);
        }

		#ifdef UG_PARALLEL
		R->set_storage_type(PST_CONSISTENT);
		#endif

		for (int type = 1; type < CT_ALL; type = type << 1)
		{
			for (size_t i = 0; i < m_vConstraint.size(); ++i)
			{
				if (m_vConstraint[i]->type() & type)
					m_vConstraint[i]->adjust_restriction(*R, spCoarseDD, spFineDD, type);
			}
		}

		write_debug(*R, "R", coarseGL, fineGL);
	}

	return m_mRestriction[key];
}

template <typename TDomain, typename TAlgebra>
void ParticleTransfer<TDomain, TAlgebra>::
prolongate(GF& uFine, const GF& uCoarse)
{
	PROFILE_FUNC_GROUP("gmg");

	if(!bCached)
		UG_THROW("ParticleTransfer: currently only cached implemented.");

	const GridLevel& coarseGL = uCoarse.grid_level();
	const GridLevel& fineGL = uFine.grid_level();
	ConstSmartPtr<ApproximationSpace<TDomain> > spApproxSpace = uFine.approx_space();
	if(uCoarse.approx_space() != spApproxSpace)
		UG_THROW("ParticleTransfer: cannot prolongate between grid functions from "
				"different approximation spaces.");

	try{
		//prolongation(fineGL, coarseGL, spApproxSpace)->apply(uFine, uCoarse);
#ifdef UG_PARALLEL
		MatMultDirect(uFine, m_dampProl, *prolongation(fineGL, coarseGL, spApproxSpace), uCoarse);
#else
		prolongation(fineGL, coarseGL, spApproxSpace)->axpy(uFine, 0.0, uFine, m_dampProl, uCoarse);
#endif

	// 	adjust using constraints
		for (int type = 1; type < CT_ALL; type = type << 1)
		{
			for (size_t i = 0; i < m_vConstraint.size(); ++i)
			{
				if (m_vConstraint[i]->type() & type)
					m_vConstraint[i]->adjust_prolongation(uFine, fineGL, uCoarse, coarseGL, type);
			}
		}
        
        adjust_prolongation(uFine, fineGL, uCoarse, coarseGL, spApproxSpace);

	}
	UG_CATCH_THROW("ParticleTransfer:prolongation: Failed for fine = "<<fineGL<<" and "
	               " coarse = "<<coarseGL);

// 	check CR functions
#ifdef UG_PARALLEL
	bool bCROnly = true;
	for(size_t fct = 0; fct < spApproxSpace->num_fct(); ++fct)
		if(spApproxSpace->lfeid(fct).type() != LFEID::CROUZEIX_RAVIART &&
				spApproxSpace->lfeid(fct).type() != LFEID::PIECEWISE_CONSTANT)
			bCROnly = false;

	if(bCROnly){
		ScaleLayoutValues(&uFine, uFine.layouts()->master(), 0.5);
		ScaleLayoutValues(&uFine, uFine.layouts()->slave(), 0.5);
		AdditiveToConsistent(&uFine, uFine.layouts()->master(), uFine.layouts()->slave(),
		                     &uFine.layouts()->comm());
	}
#endif
}

template <typename TDomain, typename TAlgebra>
void ParticleTransfer<TDomain, TAlgebra>::
do_restrict(GF& uCoarse, const GF& uFine)
{
    UG_LOG("___ do_restrict ...\n");

	PROFILE_FUNC_GROUP("gmg");

	if(!bCached)
		UG_THROW("ParticleTransfer: currently only cached implemented.");

	const GridLevel& coarseGL = uCoarse.grid_level();
	const GridLevel& fineGL = uFine.grid_level();
	ConstSmartPtr<ApproximationSpace<TDomain> > spApproxSpace = uFine.approx_space();
	if(uCoarse.approx_space() != spApproxSpace)
		UG_THROW("ParticleTransfer: cannot prolongate between grid functions from "
				"different approximation spaces.");
	try{

		restriction(coarseGL, fineGL, spApproxSpace)->
				apply_ignore_zero_rows(uCoarse, m_dampRes, uFine);

        
//		m_spParticleHandlerGlobal->plotCoarse(uCoarse);
//		m_spParticleHandlerGlobal->plotFine(uFine);

	// 	adjust using constraints
		for (int type = 1; type < CT_ALL; type = type << 1)
		{
			for (size_t i = 0; i < m_vConstraint.size(); ++i)
			{
				if (m_vConstraint[i]->type() & type)
					m_vConstraint[i]->adjust_restriction(uCoarse, coarseGL, uFine, fineGL, type);
			}
		}

        adjust_restriction(uCoarse, coarseGL, uFine, fineGL, spApproxSpace);

	} UG_CATCH_THROW("ParticleTransfer:do_restrict: Failed for fine = "<<fineGL<<" and "
	                 " coarse = "<<coarseGL);
}

template <typename TDomain, typename TAlgebra>
SmartPtr<ITransferOperator<TDomain, TAlgebra> >
ParticleTransfer<TDomain, TAlgebra>::clone()
{
	SmartPtr<ParticleTransfer> op(new ParticleTransfer);
	for(size_t i = 0; i < m_vConstraint.size(); ++i)
		op->add_constraint(m_vConstraint[i]);
	op->set_restriction_damping(m_dampRes);
	op->set_prolongation_damping(m_dampProl);
	op->set_debug(m_spDebugWriter);
	op->enable_p1_lagrange_optimization(p1_lagrange_optimization_enabled());
	op->set_use_transposed(m_bUseTransposed);
    op->set_global_handler(m_spParticleHandlerGlobal);
    
	return op;
}

template <typename TDomain, typename TAlgebra>
void ParticleTransfer<TDomain, TAlgebra>::
write_debug(const matrix_type& mat, std::string name,
            const GridLevel& glTo, const GridLevel& glFrom)
{
	PROFILE_FUNC_GROUP("debug");
//	if no debug writer set, we're done
	if(m_spDebugWriter.invalid()) return;

//	cast dbg writer
	SmartPtr<GridFunctionDebugWriter<TDomain, TAlgebra> > dbgWriter =
			m_spDebugWriter.template cast_dynamic<GridFunctionDebugWriter<TDomain, TAlgebra> >();

//	check success
	if(dbgWriter.invalid()) return;

//	add iter count to name
	name.append("_").append(ToString(glTo.level()));
	name.append("_").append(ToString(glFrom.level()));
	name.append(".mat");

//	write
	GridLevel gridLev = dbgWriter->grid_level();
	dbgWriter->set_grid_levels(glFrom, glTo);
	dbgWriter->write_matrix(mat, name.c_str());
	dbgWriter->set_grid_level(gridLev);
}

} // end namespace ug

#endif /* __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__PARTICLE_TRANSFER_IMPL__ */
