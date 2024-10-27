/*
 * ArteExpandFracs3D.cpp
 *
 *  Created on: 06.10.2024
 *      Author: Markus M. Knodel
 *
 *  * expand fractures using the Arte algorithm, 3D case
 *
 * Author: Markus Knodel, inspired by Arte from Fuchs and Sebastian Reiters code for fracture expansion without Arte
 *
 * implementing a class that gives the basic tools for Arte in 3D
 * might be templated at a later stage to fulfill 2D and 3D purposes, if suitable
 * ( so far 2D case one entire function, not so perfect, but running....)
 *
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

#include <lib_grid/algorithms/extrusion/ArteExpandFracs3D.h>

namespace ug
{

ArteExpandFracs3D::ArteExpandFracs3D(
		Grid & grid, SubsetHandler & sh,
	    std::vector<FractureInfo> const & fracInfos,
		bool useTrianglesInDiamonds, bool establishDiamonds )
	: m_grid(grid),
	  m_sh(sh),
	  m_fracInfos(fracInfos),
	  m_useTrianglesInDiamonds(useTrianglesInDiamonds),
	  m_establishDiamonds(establishDiamonds),
	  m_aaPos(Grid::VertexAttachmentAccessor<APosition>()),
//	  m_facDescr(FaceDescriptor()),
//	  m_volDescr(VolumeDescriptor()),
	  m_fracInfosBySubset(std::vector<FractureInfo>()),
	  //m_sel(Selector()),
	  m_aAdjMarkerVFP(AttVerFracProp()),
	  m_aaMarkVrtVFP( Grid::VertexAttachmentAccessor<AttVerFracProp>()),
//	  m_aAdjMarkerVFP(AttVerFracProp()),
	  m_aaMarkEdgeVFP(Grid::EdgeAttachmentAccessor<AttVerFracProp>()),
	  m_aAdjMarkerB(ABool()),
	  m_aaMarkFaceB(Grid::FaceAttachmentAccessor<ABool>()),
	  m_originalFractureFaces(std::vector<Face*>())
{
	// Notloesung, da copy constructor privat
	m_sel = Selector();
}


ArteExpandFracs3D::~ArteExpandFracs3D()
{
	//  Auto-generated destructor stub
}

bool ArteExpandFracs3D::run()
{
	if( ! initialize() )
		return false;

	if( ! setSelector() )
		return false;

	if( ! attachMarkersFirst() )
		return false;

	if( ! countAndSelectFracBaseNums() )
		return false;

	if( ! assignOrigFracInfos() )
		return false;

	if( ! attachMarkersSecond() )
		return false;


	UG_LOG("under construction " << std::endl);

	if( ! detachMarkers() )
		return false;

	return {};
}

bool ArteExpandFracs3D::initialize()
{
	UG_LOG("initialize " << std::endl);

	//	access position attachment
	if(!m_grid.has_vertex_attachment(aPosition) )
	{
		UG_LOG("Error in ExpandFractures Arte 3D: Missing position attachment");
		return false;
	}

	m_aaPos = Grid::VertexAttachmentAccessor<APosition>(m_grid, aPosition);

	//	make sure that the required options are enabled.

	if( ! m_grid.option_is_enabled(VOLOPT_AUTOGENERATE_FACES) )
	{
		UG_LOG("WARNING in Arte 3D init : grid option VOLOPT_AUTOGENERATE_FACES autoenabled.\n");
		m_grid.enable_options(VOLOPT_AUTOGENERATE_FACES);
	}

	if( ! m_grid.option_is_enabled(FACEOPT_AUTOGENERATE_EDGES) )
	{
		UG_LOG("WARNING in Arte 3D init: grid option FACEOPT_AUTOGENERATE_EDGES autoenabled.\n");
		m_grid.enable_options(FACEOPT_AUTOGENERATE_EDGES);
	}

	//	vectors that allow to access fracture properties by subset index
	m_fracInfosBySubset = std::vector<FractureInfo>( m_sh.num_subsets(), FractureInfo(-1, -1, 0) );

	for( size_t i = 0; i < m_fracInfos.size(); ++i)
	{
		if( m_fracInfos[i].subsetIndex >= m_sh.num_subsets())
		{
			throw(UGError("Bad subsetIndex in given fracInfos."));
		}

		m_fracInfosBySubset[ m_fracInfos[i].subsetIndex] = m_fracInfos[i];

	}


	return true;
}


bool ArteExpandFracs3D::setSelector()
{
	//	Collect surrounding volumes, faces and edges of all fractures in a selector
	//	and select fracture faces, edges and vertices too.

	m_sel = Selector(m_grid);

	m_sel.enable_autoselection(false);
	m_sel.enable_selection_inheritance(true);	//required for DistributeExpansionMarks3D. disabled later on.
	m_sel.enable_strict_inheritance(false);

	return true;
}


bool ArteExpandFracs3D::attachMarkersFirst()
{
	// attachment pair boundary is fracture, number fractures crossing

	m_aAdjMarkerVFP = AttVerFracProp();

	VertexFractureProperties<IndexType> vfp0( false, 0 );
	// default value: no boundary fracture, no fractures crossing

	m_grid.attach_to_vertices_dv( m_aAdjMarkerVFP, vfp0 );
	m_aaMarkVrtVFP = Grid::VertexAttachmentAccessor<AttVerFracProp> ( m_grid, m_aAdjMarkerVFP );

	//m_aAdjMarkerVFP = AttVerFracProp();
	m_grid.attach_to_edges_dv( m_aAdjMarkerVFP, vfp0 );

	m_aaMarkEdgeVFP = Grid::EdgeAttachmentAccessor<AttVerFracProp>( m_grid, m_aAdjMarkerVFP );

	m_aAdjMarkerB = ABool(); // used to know if an face is frac face

	m_grid.attach_to_faces_dv( m_aAdjMarkerB, false );
	m_aaMarkFaceB = Grid::FaceAttachmentAccessor<ABool>( m_grid, m_aAdjMarkerB );




	//  TODO FIXME
	//  das fehlt hier , Analogon 2D Fall!!!!!!!!!!!!! der geht hier eigentlich weiter
	// die Vertizes, Faces und Edges, die mit einer Kluft zu tun haben
	//	using VertFracTrip = VertexFractureTriple<Edge*, Face*, vector3>;
	//	using VecVertFracTrip = std::vector<VertFracTrip>;
	//	VecVertFracTrip vertexNoInfo;

	return true;
}


bool ArteExpandFracs3D::detachMarkers()
{
	m_grid.detach_from_vertices( m_aAdjMarkerVFP );
	m_grid.detach_from_edges( m_aAdjMarkerVFP );
	m_grid.detach_from_faces( m_aAdjMarkerB );

	return true;
}

bool ArteExpandFracs3D::countAndSelectFracBaseNums()
{
	for(size_t i_fi = 0; i_fi < m_fracInfos.size(); ++i_fi )
	{
		int fracInd = m_fracInfos[i_fi].subsetIndex;

		for( FaceIterator iter = m_sh.begin<Face>(fracInd); iter != m_sh.end<Face>(fracInd); ++iter )
		{
			Face* fac = *iter;

			m_sel.select(fac);

			m_aaMarkFaceB[fac] = true;

			for(size_t i = 0; i < fac->num_vertices(); ++i)
			{
				Vertex* vrt = fac->vertex(i);
				m_sel.select(vrt);
				m_aaMarkVrtVFP[vrt]++;

				if( IsBoundaryVertex3D(m_grid, vrt) )
				{
					m_aaMarkVrtVFP[vrt].setIsBndFracVertex();
				}
			}

			std::vector<Edge*> tmpEdges; // used for temporary results.

			CollectEdges( tmpEdges, m_grid, fac );

			for( auto const & edg : tmpEdges )
			{
				m_aaMarkEdgeVFP[edg]++;

				if( IsBoundaryEdge3D( m_grid, edg ) )
					m_aaMarkEdgeVFP[edg].setIsBndFracVertex();

				m_sel.select(edg);
			}
		}
	}

#if 0

	// TODO FIXME das ist was komisches, was von Prof. Reiter da ist, es werden edges gesplittet, für was?

	edges.clear();
	for(EdgeIterator iter = sel.begin<Edge>();
		iter != sel.end<Edge>(); ++iter)
	{
		Edge* e = *iter;
		if(aaMarkVRT[e->vertex(0)] != 2 &&
			aaMarkVRT[e->vertex(1)] != 2 &&
			aaMarkEDGE[e] > 1)
		{
			edges.push_back(e);
		}
	}

	for(size_t i = 0; i < edges.size(); ++i){
		vector3 center = CalculateCenter(edges[i], aaPos);
		RegularVertex* v =	SplitEdge<RegularVertex>(grid, edges[i], false);
		aaPos[v] = center;
		aaMarkVRT[v] = 2;
		sel.select(v);
	//	assign adjacency values for associated selected edges (2 to each)
		for(Grid::AssociatedEdgeIterator iter = grid.associated_edges_begin(v);
			iter != grid.associated_edges_end(v); ++iter)
		{
			if(sel.is_selected(*iter))
				aaMarkEDGE[*iter] = 2;
		}
	}


#endif

	for( VertexIterator iter = m_sel.begin<Vertex>(); iter != m_sel.end<Vertex>(); ++iter)
	{
		Vertex* vrt = *iter;

		bool wahl = true;

		bool isBnd = m_aaMarkVrtVFP[ vrt ].getIsBndFracVertex();
		auto numCrosFrac = m_aaMarkVrtVFP[ vrt ].getNumberFracEdgesInVertex();

		// TODO FIXME das ist richtig für den 2D Fall, aber passt das auch im 3D Fall????
		if( ! isBnd && numCrosFrac == 1 )
		{
			wahl = false;
		}

		// was, wenn numCrossFrac == 0 ist?
		// wieso werden die boundary vrt ausgeschlossen, oder sollen die nicht ausgeschlossen werden?
		// schon im 2D Fall unklar, hier noch wirrer!!! TODO FIXME

		if( wahl )
//		if( m_aaMarkVrtVFP[vrt].getNumberFracEdgesInVertex() > 1 ) // TODO FIXME stimmt das so?
		{
			//	select all associated edges, faces and volumes
			m_sel.select( m_grid.associated_edges_begin(vrt), m_grid.associated_edges_end(vrt) );
			m_sel.select( m_grid.associated_faces_begin(vrt), m_grid.associated_faces_end(vrt) );
			m_sel.select( m_grid.associated_volumes_begin(vrt), m_grid.associated_volumes_end(vrt) );
		}
	}



	// TODO FIXME hier geht es vom 3D Fall Sebastianartig weiter

	return true;
}

bool ArteExpandFracs3D::assignOrigFracInfos()
{
	m_originalFractureFaces.clear();

	for( FaceIterator iter = m_sel.begin<Face>(); iter != m_sel.end<Face>(); ++iter)
	{
		if( m_aaMarkFaceB[*iter] == true )
			m_originalFractureFaces.push_back(*iter);
	}

	m_fracInfosBySubset = std::vector<FractureInfo>( m_sh.num_subsets(), FractureInfo(-1, -1, 0) );

	for(size_t i = 0; i < m_fracInfos.size(); ++i)
	{
		if( m_fracInfos[i].subsetIndex >= m_sh.num_subsets())
		{
			throw(UGError("Bad subsetIndex in given fracInfos."));
		}

		m_fracInfosBySubset[ m_fracInfos[i].subsetIndex ] = m_fracInfos[i];
	}

//	disable selection inheritance to avoid infinite recursion.
	m_sel.enable_selection_inheritance(false);

	return true;
}

bool ArteExpandFracs3D::attachMarkersSecond()
{

	return {};
}




} /* namespace ug */
