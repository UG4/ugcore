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
	  m_facDescr(FaceDescriptor()),
	  m_volDescr(VolumeDescriptor()),
	  m_tmpEdges(std::vector<Edge*>()),
	  m_tmpFaces(std::vector<Face*>()),
	  m_tmpVols(std::vector<Volume*>())
{

}


ArteExpandFracs3D::~ArteExpandFracs3D()
{
	//  Auto-generated destructor stub
}

bool ArteExpandFracs3D::run()
{
	if( ! initialize() )
		return false;

	UG_LOG("under construction " << std::endl);
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



	return {};
}


} /* namespace ug */
