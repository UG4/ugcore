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

#include <boost/function.hpp>

#include "lib_grid/algorithms/geom_obj_util/geom_obj_util.h"
#include "lib_grid/callbacks/callbacks.h"
#include "lib_grid/grid/grid_util.h"

#include <stack>
#include <utility>
#include <vector>
#include <type_traits>
#include <limits>
#include <atomic>
#include <cstddef>
#include <bitset>
#include <string>
#include <cmath>

#include "support.h"
#include "support3D.h"


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
	  m_aAdjMarkerVFP(AttVertFracProp()),
	  m_aaMarkVrtVFP( Grid::VertexAttachmentAccessor<AttVertFracProp>()),
//	  m_aAdjMarkerVFP(AttVertFracProp()),
	  m_aaMarkEdgeVFP(Grid::EdgeAttachmentAccessor<AttVertFracProp>()),
	  m_aAdjMarkerB(ABool()),
	  m_aaMarkFaceB(Grid::FaceAttachmentAccessor<ABool>()),
	  m_originalFractureFaces(std::vector<Face*>()),
//	  m_attVrtVec(AttVrtVec()),
//	  m_aaVrtVecVol( Grid::VolumeAttachmentAccessor<AttVrtVec>() ),
	  m_aAdjInfoEdges(AttVecEdge()),
	  m_aAdjInfoFaces(AttVecFace()),
	  m_aAdjInfoVols(AttVecVol()),
	  m_aaVrtInfoAssoEdges( Grid::VertexAttachmentAccessor<AttVecEdge>()),
	  m_aaVrtInfoAssoFaces( Grid::VertexAttachmentAccessor<AttVecFace>()),
	  m_aaVrtInfoAssoVols( Grid::VertexAttachmentAccessor<AttVecVol>()),
	  m_aAdjInfoAVVFT( AttVecVertFracTrip() ),
	  m_aaVrtInfoFraTri(Grid::VertexAttachmentAccessor<AttVecVertFracTrip>()),
//	  m_vrtxFractrQuadrplVec(VrtxFractrQuadrplArte3DVec())
	  m_attVrtVec(AttVrtVec()),
	  m_aaVrtVecVol( Grid::VolumeAttachmentAccessor<AttVrtVec>() ),
	  m_vecCrossVrtInf(std::vector<CrossVertInf>())
{
	// Notloesung, nicht in die erste Initialisierung vor geschweifter Klammer, da copy constructor privat
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

	UG_LOG("initialisiert" << std::endl);

	if( ! setSelector() )
		return false;

	UG_LOG("selektiert" << std::endl);

	if( ! attachMarkers() )
		return false;

	UG_LOG("attached" << std::endl);

	if( ! countAndSelectFracBaseNums() )
		return false;

	UG_LOG("gezaehlt" << std::endl);


	if( ! assignOrigFracInfos() )
		return false;

	UG_LOG("assigniert" << std::endl);


	if( ! establishNewVrtBase() )
		return false;

	UG_LOG("etabliert" << std::endl);


	if( ! generateVertexInfos() )
		return false;

	UG_LOG("generiert" << std::endl);

	return true;

	if( ! createConditionForNewVrtcs() )
		return false;

	UG_LOG("kreiert" << std::endl);

	if( ! loop2EstablishNewVertices() )
		return false;

	UG_LOG("loopiert" << std::endl);

	UG_LOG("under construction " << std::endl);

	if( ! detachMarkers() )
		return false;

	UG_LOG("detachiert" << std::endl);

	return true;
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

//	m_sel = Selector(m_grid);

	m_sel.assign_grid(m_grid);

	m_sel.enable_autoselection(false);
	m_sel.enable_selection_inheritance(true);	//required for select and mark, disabled later
	m_sel.enable_strict_inheritance(false);

//	bool strictInherit = m_sel.strict_inheritance_enabled();
//
//	UG_LOG("strikte Inheritenz ist " << strictInherit << " und false ist " << false << std::endl);

	return true;
}


bool ArteExpandFracs3D::attachMarkers()
{
	// first part

	// attachment pair boundary is fracture, number fractures crossing

	m_aAdjMarkerVFP = AttVertFracProp();

//	support::VertexFracturePropertiesVol<IndexType> vfp0; // false, 0 );
	VertxFracPropts vfp0; // false, 0 );
	// default value: no boundary fracture, no fractures crossing

	m_grid.attach_to_vertices_dv( m_aAdjMarkerVFP, vfp0 );
	m_aaMarkVrtVFP = Grid::VertexAttachmentAccessor<AttVertFracProp> ( m_grid, m_aAdjMarkerVFP );

	//m_aAdjMarkerVFP = AttVertFracProp();
	m_grid.attach_to_edges_dv( m_aAdjMarkerVFP, vfp0 );

	m_aaMarkEdgeVFP = Grid::EdgeAttachmentAccessor<AttVertFracProp>( m_grid, m_aAdjMarkerVFP );

	m_aAdjMarkerB = ABool(); // used to know if an face is frac face

	m_grid.attach_to_faces_dv( m_aAdjMarkerB, false );
	m_aaMarkFaceB = Grid::FaceAttachmentAccessor<ABool>( m_grid, m_aAdjMarkerB );

	// second part

//	m_grid.attach_to_volumes(m_attVrtVec);
//	m_aaVrtVecVol = Grid::VolumeAttachmentAccessor<AttVrtVec>( m_grid, m_attVrtVec);

	std::vector<Edge*> noEdge;
	std::vector<Face*> noFace;
	std::vector<Volume*> noVol;

	m_aAdjInfoEdges = AttVecEdge();
	m_aAdjInfoFaces = AttVecFace();
	m_aAdjInfoVols = AttVecVol();

	m_grid.attach_to_vertices_dv( m_aAdjInfoEdges, noEdge );
	m_aaVrtInfoAssoEdges = Grid::VertexAttachmentAccessor<AttVecEdge>( m_grid, m_aAdjInfoEdges );

	m_grid.attach_to_vertices_dv( m_aAdjInfoFaces, noFace );
	m_aaVrtInfoAssoFaces = Grid::VertexAttachmentAccessor<AttVecFace>( m_grid, m_aAdjInfoFaces );

	m_grid.attach_to_vertices_dv( m_aAdjInfoVols, noVol );
	m_aaVrtInfoAssoVols = Grid::VertexAttachmentAccessor<AttVecVol>( m_grid, m_aAdjInfoVols );


	//  TODO FIXME
	//  das fehlt hier , Analogon 2D Fall!!!!!!!!!!!!! der geht hier eigentlich weiter
	// die Vertizes, Faces und Edges, die mit einer Kluft zu tun haben
	//	using VertFracTrip = VertexFractureTriple<Edge*, Face*, vector3>;
	//	using VecVertFracTrip = std::vector<VertFracTrip>;
	//	VecVertFracTrip vertexNoInfo;

	// AttVecVertFracTrip m_aAdjInfoAVVFT;

	VecVertFracTrip vertexNoInfo;

	m_aAdjInfoAVVFT = AttVecVertFracTrip();

	m_grid.attach_to_vertices_dv( m_aAdjInfoAVVFT, vertexNoInfo );

	m_aaVrtInfoFraTri = Grid::VertexAttachmentAccessor<AttVecVertFracTrip>(m_grid,  m_aAdjInfoAVVFT );


	//	associate a vector of vertices for each volume adjacent to the frac.
	//	An entry later will contain the new vertex, if the
	//	corresponding vertex is an inner fracture vertex, and nullptr if not.

	m_attVrtVec = AttVrtVec();

	m_grid.attach_to_volumes(m_attVrtVec);

	m_aaVrtVecVol = Grid::VolumeAttachmentAccessor<AttVrtVec>(m_grid, m_attVrtVec);


	return true;
}


bool ArteExpandFracs3D::detachMarkers()
{
	m_grid.detach_from_vertices( m_aAdjMarkerVFP );
	m_grid.detach_from_edges( m_aAdjMarkerVFP );
	m_grid.detach_from_faces( m_aAdjMarkerB );

	m_grid.detach_from_vertices( m_aAdjInfoEdges );
	m_grid.detach_from_vertices( m_aAdjInfoFaces );
	m_grid.detach_from_vertices( m_aAdjInfoVols );

	m_grid.detach_from_vertices( m_aAdjInfoAVVFT  );

	m_grid.detach_from_volumes( m_attVrtVec );


	return true;
}

bool ArteExpandFracs3D::countAndSelectFracBaseNums()
{
	UG_LOG("countandselect" << std::endl);

	for(size_t i_fi = 0; i_fi < m_fracInfos.size(); ++i_fi )
	{
		int fracIndSudo = m_fracInfos[i_fi].subsetIndex;

		UG_LOG("sudo ind " << fracIndSudo << std::endl);

		for( FaceIterator iter = m_sh.begin<Face>(fracIndSudo); iter != m_sh.end<Face>(fracIndSudo); ++iter )
		{
			Face* fac = *iter;

//			UG_LOG("Gesicht " << m_aaPos[fac] << std::endl);

			vector3 facCenter = CalculateCenter( fac, m_aaPos );
			UG_LOG("fac center " << facCenter << std::endl);

			for( IndexType i = 0; i < fac->num_vertices(); i++ )
				UG_LOG("Vertex " << i << " -> " << m_aaPos[fac->vertex(i)] << std::endl);

			UG_LOG("alle Vertizes" << std::endl);

			m_sel.select(fac);

			UG_LOG("selektiert msel " << fac << std::endl);

			m_aaMarkFaceB[fac] = true;

			UG_LOG("mark bool " << fac << std::endl);

//			return true;

			std::vector<Edge*> facEdges;

			UG_LOG("kollektiere ecken" << std::endl);

			CollectEdges( facEdges, m_grid, fac );

			IndexType d_anzahlEcken = facEdges.size();

			if( d_anzahlEcken == 0 )
			{
				UG_LOG("keine Ecken " << std::endl);
				return true;
			}

			UG_LOG("Anzahl Ecken " << d_anzahlEcken << std::endl);

			for( auto const & edg : facEdges )
			{
				m_aaMarkEdgeVFP[edg]++;

				if( IsBoundaryEdge3D( m_grid, edg ) )
					m_aaMarkEdgeVFP[edg].setIsBndFracVertex();

				m_sel.select(edg);
			}

			UG_LOG("Ecken gesammelt " << std::endl);

			for(size_t i = 0; i < fac->num_vertices(); ++i)
			{
				Vertex* vrt = fac->vertex(i);
				m_sel.select(vrt);
				// TODO FIXME hier anpassen, herausfinden ob fracture geschlossen
				// oder inneres Ende, und Anzahl der umgebenden fractures bestimmen!!!

				auto & vrtxFracPrps = m_aaMarkVrtVFP[ vrt ];

//				m_aaMarkVrtVFP[vrt]++;
				vrtxFracPrps++;
				// vielleicht auch in nachfolgendem Loop über alle selektierten Vertizes,
				// wo dann die attached faces abgelaufen und dabei die Subdomain Nummer ermittelt wird
				// kann eventuell sogar im Hauptloop gemacht werden, dann könnte man praktisch
				// das alte vertex fracture properties von 2D weiter verwenden, noch zu klären

//				m_aaMarkVrtVFP[vrt].addFractSudo(fracIndSudo);
				vrtxFracPrps.addFractSudo(fracIndSudo);

				if( IsBoundaryVertex3D(m_grid, vrt) )
				{
//					m_aaMarkVrtVFP[vrt].setIsBndFracVertex();
					vrtxFracPrps.setIsBndFracVertex();
				}

				// die Ecken heraus filtern, die mit diesem Vertex assoziert sind
				std::vector<Edge*> attEdg;

				for( auto const & ae: facEdges )
				{
					if( EdgeContains(ae,vrt) )
						attEdg.push_back(ae);
				}

				if( attEdg.size() == 2 )
				{
					EdgePair edgPr( attEdg[0], attEdg[1] );

					AttachedFaceEdgeSudo afes( fac, edgPr, fracIndSudo );

					vrtxFracPrps.addAttachedElem(afes);
				}
				else
				{
					UG_THROW("number of attached edges wrong " << std::endl);
				}

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

	UG_LOG("neuer Beginn" << std::endl);
//	return true;

	for( VertexIterator iter = m_sel.begin<Vertex>(); iter != m_sel.end<Vertex>(); ++iter)
	{
		Vertex* vrt = *iter;

		bool wahl = true;

		auto & vrtxFracPrps = m_aaMarkVrtVFP[ vrt ];

//		bool isBnd = m_aaMarkVrtVFP[ vrt ].getIsBndFracVertex();
		bool isBnd = vrtxFracPrps.getIsBndFracVertex();
//		auto numCrosFrac = m_aaMarkVrtVFP[ vrt ].getNumberFracEdgesInVertex();

		VertxFracPropts::VrtxFracStatus vfpStatus =  vrtxFracPrps.getVrtxFracStatus();

		if( vfpStatus == VertxFracPropts::noFracSuDoAtt )
			UG_THROW("vertex selected and no frac " << std::endl);

		// TODO FIXME das ist richtig für den 2D Fall, aber passt das auch im 3D Fall???? nein, da SudoFrac Zahl nötig
//		if( ! isBnd && numCrosFrac == 1 )
//		if( ! isBnd && vfpStatus == VertxFracPropts::oneFracSuDoAtt && )
		// TODO FIXME was, wenn ein Teil geschlossen ist der fractures, und ein anderer nicht???
		//static_assert(false);


		// bestimmen, ob die vertizes der fracture vertizes von ihrer subdomain komplett umzingelt werden
		// muss vor  hier geklärt werden!!!

//		VecPairSudoBool sudoIsSourrounded;



		bool allClosed = false;

		allClosed = isVrtxSurroundedByFracFaces( vrt, vrtxFracPrps ); //, sudoIsSourrounded );
//			// case boundary: figure out if the vertex is surrounded by frac faces of which two end in
//			// boundary edges, similar the case when the boundary face has itself two
//			// boundary edges, where the vertex is connected to both of them, then it is easy

//		if( ! allClosed )
//			return false;

//		return true;

//		if( ! isBnd )
//		{
//			UG_LOG("test if closed" << std::endl);
//
////			m_sh.assign_subset(vrt,m_sh.num_subsets());
//
////			UG_LOG("test if closed assign" << std::endl);
//
//
//		}
//		else
//		{
//
//
//
////			allClosed =
//		}
//		// TODO FIXME auch bei einer boundary muss das gemacht werden!

		UG_LOG("getestet if closed " << m_aaPos[vrt] << std::endl);


//		return true;

		if( ! allClosed == vrtxFracPrps.getInfoAllFracturesSameClosedState<false>() )
			UG_THROW("da ist was schief gegangen " << std::endl);

		// TODO FIXME ist das so richtig? kann sein, dass das zu vereinfacht ist!!!!!
		// sudo is suourrounded muss übertragen werden TODO FIXME
//		if( ! isBnd && vrtxFracPrps.getInfoAllFracturesSameClosedState<false>() )
		// bei mehreren subdoms eventuell komplizierter, kommt aber hoffentlich am Rand nicht vor......
		if( vrtxFracPrps.getInfoAllFracturesSameClosedState<false>() )
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

			std::vector<Edge*> assoEdg;
			std::vector<Face*> assoFac;
			std::vector<Volume*> assoVol;

			for( std::vector<Edge *>::iterator iterEdg = m_grid.associated_edges_begin(vrt);
											   iterEdg != m_grid.associated_edges_end(vrt);
											   iterEdg++ )
			{
				assoEdg.push_back(*iterEdg);
			}

			for( std::vector<Face *>::iterator iterFac = m_grid.associated_faces_begin(vrt);
											   iterFac != m_grid.associated_faces_end(vrt);
											   iterFac++ )
			{
				assoFac.push_back(*iterFac);
			}

			for( std::vector<Volume *>::iterator iterVol = m_grid.associated_volumes_begin(vrt);
											   	 iterVol != m_grid.associated_volumes_end(vrt);
											   	 iterVol++ )
			{
				assoVol.push_back(*iterVol);
			}

			m_aaVrtInfoAssoEdges[vrt] = assoEdg;
			m_aaVrtInfoAssoFaces[vrt] = assoFac;
			m_aaVrtInfoAssoVols[vrt] = assoVol;


		}
	}

	UG_LOG("vertex Infos Runde eins fertig " << std::endl);


	return true;
}

// herausfinden für Sudo der frac, ob bezüglich dieser sudo die faces geschlossen sind, oder ob ein Fracture End vorliegt
bool ArteExpandFracs3D::isVrtxSurroundedByFracFaces( Vertex * const & vrt, VertxFracPropts & vrtxFracPrps )
//, VecPairSudoBool & sudoSurrounded )
{
	// TODO FIXME wenn an Boundary, dann auch auf closed open unterscheiden - sowohl wenn nur edge an
	// boundary, aber auch wenn ein ganzes face dort, noch unklar, was das bedeutet

	// ganz ähnlich wie im 2D Fall, Loopen, im Kreis drehen, kucken, ob wir vom Anfang ans Ende kommen,
	// und ob das Ende der edges wieder der Anfang der edges ist, da wir uns auf faces drehen

	// case boundary: figure out if the vertex is surrounded by frac faces of which two end in
	// boundary edges, similar the case when the boundary face has itself two
	// boundary edges, where the vertex is connected to both of them, then it is easy


	VecAttachedFaceEdgeSudo vafes = vrtxFracPrps.getAllAttachedElems();

	std::vector<IndexType> sudoList = vrtxFracPrps.getSudoList();

//	for( auto const & sudo : sudoList )
//	{
//		std::vector<Face*> tmpFaces;
//
//		CollectFaces( tmpFace, m_grid, vrt );
//
//	}

	// first compare sudo list, if equal

	std::vector<IndexType> sudosTestList;

	std::vector<bool> sudoFound( sudoList.size(), false );

	UG_LOG("sudo list size " << sudoList.size() << std::endl );

	UG_LOG("vafes list size VA " << vafes.size() << std::endl );


	for( auto const & af : vafes )
	{
		bool found = false;

		IndexType sudo = af.getSudo();

		UG_LOG("sudo af " << sudo << std::endl);

		for( IndexType i = 0; i < sudoList.size(); i++ )
		{
			UG_LOG("sudo list dusso is " << sudoList[i] << std::endl);

			if( sudo == sudoList[i] )
			{
				sudoFound[i] = true;
				found = true;
			}
		}


		if( ! found )
			UG_THROW("sudo nicht gefunden muss aber da sein " << std::endl);
	}

	UG_LOG("alles gefunden " << std::endl);


	for( auto const & sf: sudoFound )
	{
		if( sf == false )
		{
			UG_LOG("Falsch" << std::endl);
			UG_THROW("sudo not found but must be there " << std::endl);
		}
	}

	UG_LOG("alles gefunden Test " << std::endl);


	// sort faces with respect to subdomain - macht das wirklich Sinn, wie umständlich das gemacht wird jetzt?

	// if we arrive here, all sudos found, the entire circle closing can start

	VecPairSudoBool sudoSurrounded;

	for( auto const & sudo : sudoList )
	{
		VecAttachedFaceEdgeSudo vecAttFacSudo;

		for( auto const & attFac : vafes )
		{
			if( sudo == attFac.getSudo() )
			{
				UG_LOG("die sudo gefunden " << sudo << std::endl);
				vecAttFacSudo.push_back(attFac);
			}
		}

		VecAttachedFaceEdgeSudo vecAttFacSudoSort;



		bool isClosed = false;

		bool isBndVrtx = IsBoundaryVertex3D(m_grid,vrt);

		if( ! isBndVrtx )
		{
			UG_LOG("No boundary vertex test closed " << m_aaPos[vrt] << std::endl);

			if( vecAttFacSudo.size() == 1 )
			{
				// no need for further investigations for inner faces
				isClosed = false;
				vecAttFacSudoSort = vecAttFacSudo;
			}
			else
			{
				isClosed = sortElemCircleIsClosed( vecAttFacSudo, vecAttFacSudoSort );
			}

		}
		else // different treatment boundary vertex
		{
			// figure out start face with a boundary edge, and figure out an eventual additional boundary edge
			// if only one face, then check if two attached boundary edges, then true, else false

			if( vecAttFacSudo.size() == 1 )
			{
				AttachedFaceEdgeSudo & singleEntry = vecAttFacSudo[0];

				EdgePair faceEdgs = singleEntry.getLowElm();

				if( IsBoundaryEdge3D(m_grid, faceEdgs.first) && IsBoundaryEdge3D(m_grid, faceEdgs.second ) )
					isClosed = true;

				// very simple
				vecAttFacSudoSort = vecAttFacSudo;

			}
			else
			{
				// figure out a begin face with a boundary edge and figure out another face with a boundary edge

				Edge * beginEdge = nullptr;
				Edge * endEdge = nullptr;

				IndexType startFaceIndx = 0;
				IndexType endFaceIndx = 0;

				for( auto const & afs : vecAttFacSudo )
				{
					Face * fac = afs.getManifElm();

					EdgePair edgs = afs.getLowElm();

					Edge * edgOne = edgs.first;
					Edge * edgTwo = edgs.second;

					if( beginEdge == nullptr )
					{
						if( IsBoundaryEdge3D(m_grid, edgOne) )
						{
							beginEdge = edgTwo;
						}
						else if( IsBoundaryEdge3D(m_grid, edgTwo) )
						{
							beginEdge = edgOne;
						}
						else
						{
							startFaceIndx++;
						}
					}
					else
					{
						if( IsBoundaryEdge3D(m_grid, edgOne) )
						{
							endEdge = edgOne;
						}
						else if( IsBoundaryEdge3D(m_grid, edgTwo) )
						{
							endEdge = edgTwo;
						}
					}

					if( endEdge != nullptr )
						break;

					endFaceIndx++;
				}

				if( beginEdge == nullptr && endFaceIndx == vecAttFacSudo.size() )
//					|| beginEdge == nullptr || endEdge == nullptr )
				{
					UG_LOG("keine boundary edges" << std::endl);

					startFaceIndx = -1;

					if( endEdge != nullptr )
						UG_THROW("Ende nicht null, Anfang null " << std::endl);
				}

				UG_LOG("Boundary vertex test closed " << m_aaPos[vrt] << std::endl);

//				int d_num = 0;
//				for( auto const & afs : vecAttFacSudo )
//				{
//					d_num++;
//					Face * fac = afs.getManifElm();
//					m_sh.assign_subset(fac, m_sh.num_subsets());
//				}
//				UG_LOG("number of surr faces " << d_num << std::endl );

//				m_sh.assign_subset(beginEdge,m_sh.num_subsets());
//				m_sh.assign_subset(endEdge,m_sh.num_subsets());
//				m_sh.assign_subset(vecAttFacSudo[startFaceIndx].getManifElm(),m_sh.num_subsets());

				isClosed = sortElemCircleIsClosed( vecAttFacSudo, vecAttFacSudoSort, startFaceIndx, beginEdge, endEdge );

			}
		}

		UG_LOG("-------------------------------" << std::endl);
		UG_LOG("is closed " << isClosed << " at " << m_aaPos[vrt] << std::endl);
		UG_LOG("-------------------------------" << std::endl);

//		if( isClosed )
//		{
//			m_sh.assign_subset(vrt,3);
//		}
//		else
//		{
//			m_sh.assign_subset(vrt,4);
//		}

		if( vecAttFacSudo.size() != vecAttFacSudoSort.size() )
		{
//			return false;

			UG_THROW("Die Sortierung ist komplett schief gegangen " << std::endl);
		}

		// DEBUG Zeug, später entfernen!!!!!!
//		for( auto const & afss : vecAttFacSudoSort )
//		{
//			Face * fac = afss.getManifElm();
//
//			m_sh.assign_subset(fac, m_sh.num_subsets());
//		}

		PairSudoBool ic( sudo, isClosed );

		sudoSurrounded.push_back( ic );

	}

	vrtxFracPrps.setInfoAllFractureSudosIfClosed(sudoSurrounded);

	bool allClosed = true;

	for( auto const & ic : sudoSurrounded )
	{
		if( ! ic.second )
			allClosed = false;
	}

	return allClosed;
}

bool ArteExpandFracs3D::sortElemCircleIsClosed( VecAttachedFaceEdgeSudo const & vecAttFac,
												VecAttachedFaceEdgeSudo & vecSortedFac,
												int startFacIndexUser,
//												int endFacIndexUser,
//												IndexType startEdgeIndexUser,
//												IndexType endEdgeIndexUser
//												Face * const & startFacUser,
//												Face * const & endFacUser,
												Edge * const & startEdgUser,
												Edge * const & endEdgUser
												)
{
	UG_LOG("Schliesstest" << std::endl);

	IndexType originalSize = vecAttFac.size();

	if( originalSize == 0 )
	{
		UG_THROW("zu klein zum sortieren " << std::endl);
		return false;
	}
	else if ( originalSize == 1 )
	{
		UG_THROW("should not happen size 1, should have been mentioned before " << std::endl);
	}

	UG_LOG("Kopieren zwecks sortieren " << std::endl);

	for( auto const & af : vecAttFac )
	{
		UG_LOG("die sudos innen sind " << af.getSudo() << std::endl);
	}

	VecAttachedFaceEdgeSudo copyVecAttFac = vecAttFac;

	for( auto const & af : copyVecAttFac )
	{
		UG_LOG("die sudos kopiert sind " << af.getSudo() << std::endl);
	}

	IndexType beginIndx = 0;

	UG_LOG("begin Index " << beginIndx << std::endl);

	bool simpleConnectionTest = false;

	Edge * startEdgeForced = nullptr;

	if( startFacIndexUser >= 0 )
	{
		UG_LOG("Veränderung " << startFacIndexUser << std::endl);
		beginIndx = startFacIndexUser;
		simpleConnectionTest = true; // user hopefully did it
	}
	else
	{
		// we need to ensure to start at a fracture which is not in between, in case that circle not closed
		// so ensure that all fracture faces have a fracture face at both edges as neighbour,
		// in principle this is already sufficient to answer the question which we want to know
		// all the rest here is useless playing in principle

		bool broken = false;

		for( IndexType i = 0; i < vecAttFac.size(); i++ )
		{
			IndexType firstSideConnected = 0;
			IndexType secondSideConnected = 0;

			AttachedFaceEdgeSudo afBase = vecAttFac[i];

			Face * faceBase = afBase.getManifElm();
			EdgePair edgPairBase = afBase.getLowElm();

			Edge * edgeBaseOne = edgPairBase.first;
			Edge * edgeBaseTwo = edgPairBase.second;

			for( IndexType j = 0; j < vecAttFac.size(); j++ )
			{
				if( i != j )
				{
					AttachedFaceEdgeSudo afCompr = vecAttFac[j];

					Face * faceCompr = afCompr.getManifElm();
					EdgePair edgPairCompr = afCompr.getLowElm();

					Edge * edgeComprOne = edgPairCompr.first;
					Edge * edgeComprTwo = edgPairCompr.second;

					if( edgeComprOne == edgeBaseOne || edgeComprTwo == edgeBaseOne )
						firstSideConnected++;

					if( edgeComprOne == edgeBaseTwo || edgeComprTwo == edgeBaseTwo )
						secondSideConnected++;
				}

			}

			if( vecAttFac.size() > 2 && ( firstSideConnected > 1 || secondSideConnected > 1 ) )
			{
				UG_THROW("zu oft verbunden " << std::endl );
			}
			else if( vecAttFac.size() == 2 && ( firstSideConnected > 2 || secondSideConnected > 2 ) )
			{
				UG_THROW("zu oft verbunden " << std::endl );
			}
			else if( firstSideConnected == 0 )
			{
				// face is open into one direction, already clear that not closed!!

				beginIndx = i;
				simpleConnectionTest = false;
				startEdgeForced = edgeBaseTwo;
				UG_LOG("forcieren 1 " << std::endl);
				broken = true;
				break;
			}
			else if( secondSideConnected == 0 )
			{
				// face is open into one direction, already clear that not closed!!

				beginIndx = i;
				simpleConnectionTest = false;
				startEdgeForced = edgeBaseOne;
				UG_LOG("forcieren 2 " << std::endl);
				broken = true;
				break;
			}
			else if( firstSideConnected == 1 && secondSideConnected == 1 )
			{
				simpleConnectionTest = true; // as long as a look
			}
			else
			{
				UG_THROW("komischer Verbindungsfall " << std::endl);
			}

			if( broken )
				break;

		}
	}

	UG_LOG("begin Index X " << beginIndx << std::endl);


//	AttachedFaceEdgeSudo initialAFES = *(copyAttFac.begin());
	AttachedFaceEdgeSudo initialAFES = copyVecAttFac[beginIndx];

	IndexType sudo = initialAFES.getSudo();

	UG_LOG("sudo " << beginIndx << " ist " << sudo << std::endl);

	Face * beginFacLoop = initialAFES.getManifElm();

	// TODO FIXME wird das irgendwo verwendet? wieso nicht?
//	Face * endFacLoop = nullptr;

//	if( endFacIndexUser != -1 )
//	{
//		endFacLoop = copyVecAttFac[endFacIndexUser].getManifElm();
//	}

	EdgePair beginEdges = initialAFES.getLowElm();

	Edge * beginEdgeLoop = beginEdges.second;
	Edge * targetedEndEdgeLoop = beginEdges.first; // should be closed! should end at same edge as it begins!

//
//	return true;

//	if( startEdgeIndexUser != -1 )
//	{
//		beginEdgeLoop =
//	}

//	if( startFacUser != nullptr )
//	{
//		beginFacLoop = startFacUser;
//
//	}
//

//	if( ! FaceContains(beginFacLoop,beginEdgeLoop->vertex(0)) )
//		UG_THROW("Sortierung Gesicht hat nicht die gewünschte Ecke " << std::endl);
//
//	if( endFacUser != nullptr )
//	{
//		endFacLoop = endFacUser;
//	}
//

	if( startEdgUser != nullptr )
	{
		beginEdgeLoop = startEdgUser;

		// check if part of the begin face!

		if( ! FaceContains( beginFacLoop, beginEdgeLoop ) )
			UG_THROW("Anfangsgesicht hat keine Anfangsecke " << std::endl);

	}

	if( startEdgeForced != nullptr )
	{
		beginEdgeLoop = startEdgeForced;

		// check if part of the begin face!

		if( ! FaceContains( beginFacLoop, beginEdgeLoop ) )
			UG_THROW("Anfangsgesicht hat keine Anfangsecke forced " << std::endl);

//		m_sh.assign_subset(startEdgeForced, m_sh.num_subsets());
		UG_LOG("forciert " << std::endl);

//		return false;

	}


	if( endEdgUser != nullptr )
	{
		// check if part of end face at end of loop somehow
		targetedEndEdgeLoop = endEdgUser;

		if( startEdgeForced != nullptr )
			UG_THROW("das muss schief gehen, Chaos " << std::endl);

//		if( endFacLoop != nullptr )
//		{
//			if( ! FaceContains( endFacLoop, targetedEndEdgeLoop ) )
//				UG_THROW("Endgesicht hat keine Endecke " << std::endl);
//		}
	}
//

	// DEBUG
//	m_sh.assign_subset(beginFacLoop, m_sh.num_subsets());
//	m_sh.assign_subset(beginEdgeLoop, m_sh.num_subsets());
//	m_sh.assign_subset(targetedEndEdgeLoop, m_sh.num_subsets());


	// Du musst sortieren. Du musst einen Zeitplan machen. Das kann man lernen. Du kannst das selber machen.

	IndexType countedCrossedFaces = 1;

	vecSortedFac.push_back( initialAFES );

	copyVecAttFac.erase( copyVecAttFac.begin() + beginIndx );

//	Face * face2Append = beginFacLoop;
	Edge * startEdge2Append = beginEdgeLoop;

//	m_sh.assign_subset(startEdge2Append,m_sh.num_subsets());

	UG_LOG("while loop anfangen " << std::endl);

	IndexType d_whi = 0;

	while( copyVecAttFac.size() != 0 )
	{

		UG_LOG("in while loop " << d_whi << std::endl);
		d_whi++;

		IndexType foundCommEdg = 0;

		Edge * nextEdge = nullptr;

//		for( auto const & caf : copyVecAttFac )
		for( VecAttachedFaceEdgeSudo::iterator itAttFES  = copyVecAttFac.begin();
											   itAttFES != copyVecAttFac.end();
											   itAttFES++
		)
		{
			AttachedFaceEdgeSudo caf = *itAttFES;

			Face * d_Fac = caf.getManifElm();

//			m_sh.assign_subset(d_Fac,m_sh.num_subsets());

			EdgePair edgPr = caf.getLowElm();

			Edge * edgOne = edgPr.first;
			Edge * edgTwo = edgPr.second;

//			m_sh.assign_subset(edgOne,m_sh.num_subsets());
//			m_sh.assign_subset(edgTwo,m_sh.num_subsets());

//			return true;

			IndexType hasEdge = 0;

			Edge * overNextEdge = nullptr;

			if( edgOne == startEdge2Append )
			{
				nextEdge = edgOne;
				overNextEdge = edgTwo;
				hasEdge++;
			}

			if( edgTwo == startEdge2Append )
			{
				nextEdge = edgTwo;
				overNextEdge = edgOne;
				hasEdge++;
			}

			if( hasEdge > 1 )
				UG_THROW("zu viele Ecken und Kanten " << std::endl);

			if( hasEdge == 1 )
			{
				Face * fac2App = caf.getManifElm();
//				m_sh.assign_subset(fac2App,m_sh.num_subsets());
				EdgePair edgesNextFace( edgOne, edgTwo );
				AttachedFaceEdgeSudo nextAttFES( fac2App, edgesNextFace, sudo );

				vecSortedFac.push_back(nextAttFES);

				copyVecAttFac.erase(itAttFES);
				foundCommEdg++;
				startEdge2Append = overNextEdge;

				break;
			}


		}

		if( foundCommEdg > 1 )
			UG_THROW("Kein Anschluss unter dieser Nummer " << std::endl);

		if( foundCommEdg == 0 )
		{
			UG_LOG("Kreislauf nicht geschlossen " << std::endl);

			if( nextEdge != nullptr )
				UG_THROW("nicht konsistent, null und null " << std::endl);

			break;
		}

		if( nextEdge == nullptr )
		{
			if( foundCommEdg != 0 )
				UG_THROW("nicht konsistent, null und null v2 " << foundCommEdg << " -> " << nextEdge << std::endl);

	//			return false;
		}

	}

	if( originalSize != vecSortedFac.size() )
	{
		UG_THROW("Sortierung hat nicht funktioniert " << std::endl);
		return false;
	}

	if( startEdge2Append != targetedEndEdgeLoop )
	{
		if( simpleConnectionTest )
			UG_THROW("obwohl offen oder vorgegeben trotzdem Ziel nicht erreicht?" << std::endl);

		UG_LOG("Ende nicht erreicht, Kreis nicht geschlossen, innerer Rand vermutlich" << std::endl);
		return false;
	}


	UG_LOG("Kreislauf Faces 3D Test ob Knoten umrandet geschlossen " << std::endl);

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

bool ArteExpandFracs3D::establishNewVrtBase()
{
	//	iterate over all surrounding volumes to enable shifted vertices, this loop taken from SR but shortened

	for( VolumeIterator iterSurrVol = m_sel.volumes_begin(); iterSurrVol != m_sel.volumes_end(); ++ iterSurrVol )
	{
		Volume* sv = *iterSurrVol;

		std::vector<Vertex*> & newVrts = m_aaVrtVecVol[sv];
		newVrts.resize(sv->num_vertices());

		for(size_t iVrt = 0; iVrt < sv->num_vertices(); ++ iVrt )
		{
			newVrts[iVrt] = nullptr;
		}

		// erstmal so tun, als ob keine neuen Vertizes erzeugt werden an den alten Vertizes

	}

	return true;
}

// Analogon zu VertrexFractureInfo in 2D, wo jeder Vertex eine Liste bekommt, wo alle die ihm angehängten
// Ecken, Faces und Volumen gespeichert werden; dazu die Normalen, und vielleicht noch weitere Infos
bool ArteExpandFracs3D::generateVertexInfos()
{
	// TODO FIXME das wird benötigt

	// sowas von der Art als attachement bei den attachments, und dann mit Leben füllen für jeden Vertex
	// in dieser Funktion;
	// vielleicht braucht es auch Edge Infos, oder nur Edge Infos?
//	VecVertFracTrip vertexNoInfo;
//	using AttVecVertFracTrip = Attachment<VecVertFracTrip>;
//	AttVecVertFracTrip aAdjInfoAVVFT;
//	grid.attach_to_vertices_dv( aAdjInfoAVVFT, vertexNoInfo );
//	Grid::VertexAttachmentAccessor<AttVecVertFracTrip> aaVrtInfoFraTri(grid,  aAdjInfoAVVFT );

	// Lebendigmachung in:
	// 	for( auto & fsfpmv : fracSubdom_facePerpendMinVal ) .....
	// von dort lernen!!!!!

	// notwendig: face, normal, volume, edge

	// TODO FIXME das wollen wir nicht, sondern das alte Vertex Fracture Triple
//	m_vrtxFractrQuadrplVec = VrtxFractrQuadrplArte3DVec();
	// TODO FIXME diese Einträge erzeugen

	// TODO FIXME kann vielleicht vereinfacht werden durch einen Loop über alle Vertizes,
	// und das Abfragen der dort schon gespeicherten vertex property Geschichten, sonst
	// wird manches doppelt gemoppelt
	
	for( auto const & fracInf : m_fracInfos )
	{
		IndexType fracSudo = fracInf.subsetIndex;


//		for(EdgeIterator iterEdg = m_sh.begin<Edge>(fracIndSudo); iterEdg != m_sh.end<Edge>(fracIndSudo); iterEdg++ )
//
		for( FaceIterator iterFac = m_sh.begin<Face>(fracSudo); iterFac != m_sh.end<Face>(fracSudo); iterFac++ )
		{

//			VrtxFractrQuadrplArte3D vrtxFractrQuadrpl;
			
			Face* fac = *iterFac;
			
			auto sudoFacInnerLoop = m_sh.get_subset_index(fac);
			
			if( sudoFacInnerLoop != fracSudo )
				UG_THROW("Subdomain Index Fehler 3D " << std::endl);

//			std::vector<Vertex*> verticesFac;
//
//			for( IndexType i = 0; i < fac->num_vertices(); i++ )
//			{
//				verticesFac.push_back( fac->vertex(i) );
//			}


			std::vector<Volume*> assoVols;

			// wo kam denn der Käse her?
//			if( ! m_grid.option_is_enabled(FACEOPT_STORE_ASSOCIATED_VOLUMES) )
//				UG_THROW("How to collect asso vols?" << std::endl);


			if(! m_grid.option_is_enabled(VOLOPT_AUTOGENERATE_FACES) )
			{
				UG_LOG("WARNING grid option VOLOPT_AUTOGENERATE_FACES autoenabled.\n");
				m_grid.enable_options(VOLOPT_AUTOGENERATE_FACES);
			}


//			for( Grid::AssociatedVolumeIterator volIt  = m_grid.associated_volumes_begin(fac);
//												volIt != m_grid.associated_volumes_end(fac);
//												volIt++
//			)
//			{
//				assoVols.push_back(*volIt);
//			}

			CollectVolumes(assoVols, m_grid, fac );
//
//			for( auto const & av : assoVols )
//			{
//				m_sh.assign_subset(av, m_sh.num_subsets());
//			}
//
//			return true;

//			using VolNormPair = std::pair< Volume*, vector3 >;
//
//			using VecVolNormPair = std::vector<VolNormPair>;
//
//			VecVolNormPair vecNormalsAwayVol;

//			UG_LOG("Center Face " << CalculateCenter(fac,m_aaPos) << std::endl);

			for( auto const & kuhVol : assoVols )
			{
				bool facFound = false;
				IndexType numFd = 0;

				for( IndexType iSide = 0; iSide < kuhVol->num_sides(); iSide++ )
				{
					Face * kuhFac = m_grid.get_side(kuhVol, iSide);

//					UG_LOG("Center Kuh Face " << CalculateCenter(kuhFac,m_aaPos) << std::endl);

					// TODO FIXME eigentliches Ziel ist es, den face descriptor des Volumens zu finden,
					// das mit dem face übereinstimmt, alle anderen Seiten des Volumens sind egal
					// Funktion suchen, die ausgehend von einem Face, das ein Volumen begrenzt,
					// den zum Volumen passenden FaceDescriptor findet, also auch richtige Orientierung
					// der Vertices beinhaltet
					// FRAGE TODO FIXME ist ein face descriptor von der Orientierung zum Volumen abhängig
					// oder hängt der nur vom Face ab, das eine vorgegebene Oriertierung hat?

					bool checkCoincide = checkIfFacesVerticesCoincide( kuhFac, fac );

					if( checkCoincide )
					{
						facFound = true;
						numFd++;

						if( kuhFac != fac )
							UG_LOG("Kuh Fac ist nicht fac " << std::endl);

						FaceDescriptor facDescr;

						// TODO FIXME testen, ob der Face Descriptor von der Orientierung abhängt
						// also testen, ob sich der face descriptor ändert, wenn das Volumen
						// auf der einen und auf der anderen Seite des faces hängt
						// deswegen auch die ganze Prozedur mit den kuhFacs, die hoffentlich
						// je nach Volumen anders orientiert sind als das eigentliche Face,
						// aber dieselben Vertices haben, also geometrisch gleich sind, aber anders orientiert!!!!

						// TODO FIXME andere Hergehensweise vielleicht:
						// von m_aaVrtInfoAssoVols ausgehen, darüber loopen, oder die in einen Vektor stecken,
						// wo die Vertices dabei sind, dann kann man sich vielelicht ein paar Klimmzüge sparen,
						// vielleicht aber auch nicht.....

						kuhVol->face_desc( iSide, facDescr );

						vector3 normal;

						CalculateNormal( normal, & facDescr, m_aaPos );

						vector3 facCenter = CalculateCenter( kuhFac, m_aaPos );
						vector3 kuhCenter = CalculateCenter( fac, m_aaPos );
						vector3 kuhVolCenter = CalculateCenter( kuhVol, m_aaPos);

//						UG_LOG("Normale zum face descriptor " << normal << " , " << facCenter << std::endl);
//						UG_LOG("Normale zum Kuhh descriptor " << normal << " , " << kuhCenter << std::endl);
//						UG_LOG("Zentrum des Vol")

//						UG_LOG("fac " << fac << std::endl );
//						UG_LOG("kuh " << kuhFac << std::endl );

						UG_LOG("Normale ist " << normal << " fac " << facCenter
								<< " vol " << kuhVolCenter << std::endl);


//						VolNormPair normalsAwayVol( kuhVol, normal );
//
//						vecNormalsAwayVol.push_back( normalsAwayVol );

						std::vector<Edge*> facEdgs;

						CollectEdges( facEdgs, m_grid, fac);

						for( IndexType iF = 0; iF < fac->num_vertices(); iF++ )
						{
							Vertex * vrt = fac->vertex(iF);

			//				verticesFac.push_back(vrt);

							std::vector<Edge*> edgOfVrtx;

							for( auto const & edg : facEdgs )
							{
								if( EdgeContains(edg, vrt) )
								{
									edgOfVrtx.push_back(edg);
								}
							}

							if( edgOfVrtx.size() == 2 )
							{
								EdgePair commonEdges(edgOfVrtx[0], edgOfVrtx[1]); //  fill values
								// edges commun between face and volume, with the vertex included as well, i.e. two possibilities

								VertFracTrip infoVerticesThisFace( fac, fracSudo, kuhVol, normal, commonEdges );

								// TODO FIXME hier irgendwie graphische Ausgabe von irgendwas

								m_aaVrtInfoFraTri[vrt].push_back( infoVerticesThisFace );

								// DEBUG OUTPUT; REMOVE LATER
//								m_sh.assign_subset(kuhVol,m_sh.num_subsets());
							}
							else
							{
								UG_THROW("Mein Face das hat keine Ecken, keine Ecken hat mein Face" << std::endl);
							}
						}

					}
				}

				if( ! facFound || numFd != 1 )
				{
					UG_THROW("Kein Kuh Volumen gefunden" << std::endl);
					return false;
				}

			}

		}


	}

	return true;
}

bool ArteExpandFracs3D::checkIfFacesVerticesCoincide( Face * const & facOne, Face * const & facTwo )
{

	if( facOne->size() != facTwo->size() )
		return false;

	std::vector<Vertex* > facOneVrtcs, facTwoVrtcs;

	collectFaceVertices( facOneVrtcs, facOne );
	collectFaceVertices( facTwoVrtcs, facTwo );

	for( auto const & vrtOne : facOneVrtcs )
	{
		bool found = false;

		IndexType numFd = 0;

		for( auto const & vrtTwo : facTwoVrtcs )
		{
			if( vrtOne == vrtTwo )
			{
				found = true;
				numFd++;
			}
		}

		if( ! found || numFd != 1 )
			return false;
	}

	return true;
}

bool ArteExpandFracs3D::collectFaceVertices( std::vector<Vertex*> & facVrt, Face * const & fac )
{
	if( fac == nullptr )
		return false;

	facVrt.clear();

	for( IndexType iF = 0; iF < fac->num_vertices(); iF++ )
	{
		Vertex * vrt = fac->vertex(iF);

		facVrt.push_back( vrt );
	}

	return true;
}



// major function of new grid generation, in Keil Style, but functional grid, only the diamonds have to be
// established in additional functionalities independent of this function
bool ArteExpandFracs3D::loop2EstablishNewVertices()
{
	m_vecCrossVrtInf = std::vector<CrossVertInf>();


	// TODO FIXME sowas von der Art wird nötig sein als Vorbereitung für die Diamanten,
	// Infos darin speichern, vielleicht auch noch notwendig, die Kanten zu speichern oder die faces,
	// zu klären im Laufe der Implementation
	// std::vector<CrossVertInf > vecCrossVrtInf;

	// zentraler Loop

	// TODO FIXME vorher noch den attVrtVec und den aaVrtVecFac analog implementieren, das fehlt noch!!!
	// ebenso seine Befüllung, braucht noch eine Funktion dazwischen, die attachments selber in die
	// attach Funktion natürlich
	


	for( VertexIterator iterV = m_sel.begin<Vertex>(); iterV != m_sel.end<Vertex>(); ++ iterV )
	{
		Vertex * oldVrt = *iterV;

		// Position dieses Vertex
		vector3 posOldVrt = m_aaPos[oldVrt];

		// TODO FIXME diese Funktion mit Leben und Analytischer Geometrie 13. Klasse füllen

		VecVertFracTrip & vecVertFracTrip = m_aaVrtInfoFraTri[oldVrt];

		std::vector<Edge*> & allAssoEdges = m_aaVrtInfoAssoEdges[oldVrt];
		std::vector<Face*> & allAssoFaces = m_aaVrtInfoAssoFaces[oldVrt];
		std::vector<Volume*> & allAssoVolumes = m_aaVrtInfoAssoVols[oldVrt];

		UG_LOG("vertex at " << posOldVrt << std::endl);

		auto & vrtxFracPrps = m_aaMarkVrtVFP[ oldVrt ];

		bool vrtxIsBndVrt = vrtxFracPrps.getIsBndFracVertex();

		UG_LOG("is bndry " << vrtxIsBndVrt << std::endl);

		VertxFracPropts::VrtxFracStatus statusThisVrtx = vrtxFracPrps.getVrtxFracStatus();

		if( ! vrtxIsBndVrt )
		{
			if( vrtxFracPrps.getInfoAllFracturesSameClosedState<false>() )
			{
				// gar nix tun, alle offen, innerer Vertex, darf man hier ankommen? NEIN TODO FIXME
				UG_THROW("hier sollten wir nicht angekommen sein " << std::endl);
			}
			// TODO FIXME: was, wenn ein Zwischending, Mischung?

			if( statusThisVrtx == VertxFracPropts::VrtxFracStatus::noFracSuDoAtt )
			{
				UG_THROW("gar keine Frac darf hier nicht ankommen " << std::endl );
			}
			else if( statusThisVrtx == VertxFracPropts::VrtxFracStatus::oneFracSuDoAtt )
			{
				// TODO FIXME erster Fall, eine Fracture, innen, geschlossen, kann eigentlich nur hier ankommen
				UG_LOG("aktuelles Ziel " << std::endl);
			}
			else if( statusThisVrtx == VertxFracPropts::VrtxFracStatus::twoFracSuDoAtt )
			{

			}
			else if( statusThisVrtx == VertxFracPropts::VrtxFracStatus::threeFracSuDoAtt )
			{

			}
			else
			{
				UG_THROW("was für ein Knoten Status???? " << std::endl);
			}
		}
		else // boundary vertex
		{

		}

		// TODO FIXME wichtig: surrounded status closed / open TODO FIXME, sowie Anzahl der schneidenden Fracs
		
	}


	return false;
}

bool ArteExpandFracs3D::createConditionForNewVrtcs()
{

	//	iterate over all surrounding volumes to enable volume changes, this loop taken from SR but shortened
	for(VolumeIterator iterSurrVol = m_sel.volumes_begin(); iterSurrVol != m_sel.volumes_end(); iterSurrVol++ )
	{
		Volume * sv = *iterSurrVol;

		std::vector<Vertex*>& newVrts = m_aaVrtVecVol[sv];
		newVrts.resize(sv->num_vertices());

		for(size_t iVrt = 0; iVrt < sv->num_vertices(); iVrt++ )
		{
			newVrts[iVrt] = nullptr;
		}
			// erstmal so tun, als ob keine neuen Vertizes erzeugt werden an den alten Vertizes
	}


	return true;
}


} /* namespace ug */
