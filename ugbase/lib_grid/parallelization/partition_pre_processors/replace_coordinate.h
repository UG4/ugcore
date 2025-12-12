/*
 * Copyright (c) 2017:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
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

#ifndef __H__UG_replace_coordinate
#define __H__UG_replace_coordinate

namespace ug {

///	Temporarily replaces the specified coordinate in the given position attachment
/**	\warning 	After 'partitioning_starts' was called, the coordinates of the
 * 				underlying mesh may be changed. Make sure to call
 *				'partitioning_done' when the partitioner is done. Also make sure
 *				to call 'partitioning_done' if errors occurred during partitioning.
 */
template <int dim>
class PPP_ReplaceCoordinate : public IPartitionPreProcessor
{
public:
	using vector_t = MathVector<dim>;
	using apos_t = Attachment<vector_t>;
	using aapos_t = Grid::VertexAttachmentAccessor<apos_t>;

	PPP_ReplaceCoordinate (apos_t aPos, ANumber aNewCoord, int newCoordIndex) :
		m_aPos (aPos),
		m_aNewCoord (aNewCoord),
		m_newCoordIndex (newCoordIndex)
	{}
	~PPP_ReplaceCoordinate() override = default;

	void partitioning_starts (	MultiGrid* mg,
		                          IPartitioner* partitioner) override {
		if (!mg->has_vertex_attachment(m_aPos)){
			UG_LOG("WARNING: Target attachment not found in PartitionPreProcessor "
			       "PPP_ReplaceCoordinate");
			return;
		}

		if (!mg->has_vertex_attachment(m_aNewCoord)){
			UG_LOG("WARNING: Attachment 'aNewCoord' not found in PartitionPreProcessor "
			       "PPP_ReplaceCoordinate");
			return;
		}

		UG_COND_THROW (!mg->has_vertex_attachment(m_aNewCoord),
		               "VertexAttachment missing: aNewCoord");

	//	attach a temporary field where we store the original coordinates
		mg->attach_to_vertices (m_aOrigCoord);

	//	store the current original coordinate of m_aPos in m_aOrigCoord
	//	and replace with new coordinate

		typename aapos_t::ContainerType& positions = *mg->get_attachment_data_container<Vertex> (m_aPos);
		ANumber::ContainerType& origCoords = *mg->get_attachment_data_container<Vertex> (m_aOrigCoord);
		ANumber::ContainerType& newCoords = *mg->get_attachment_data_container<Vertex> (m_aNewCoord);

		UG_COND_THROW (	(positions.size() != origCoords.size())
		                || (positions.size() != newCoords.size()),
		                "Attachment containers have incompatible sizes");

		for(size_t i = 0; i < positions.size(); ++i){
			origCoords[i] = positions[i][m_newCoordIndex];
			positions[i][m_newCoordIndex] = newCoords[i];
		}
	}


	void partitioning_done (	MultiGrid* mg,
		                        IPartitioner* partitioner) override {
		if (!mg->has_vertex_attachment(m_aPos)){
			UG_LOG("WARNING: Target attachment not found in PartitionPreProcessor "
			       "PPP_ReplaceCoordinate");
			return;
		}

		if (!mg->has_vertex_attachment(m_aOrigCoord)){
			UG_LOG("WARNING: Attachment 'aOrigCoord' not found in PartitionPreProcessor "
			       "PPP_ReplaceCoordinate");
			return;
		}

		typename aapos_t::ContainerType& positions
			= *mg->get_attachment_data_container<Vertex> (m_aPos);

		typename ANumber::ContainerType& origCoords
			= *mg->get_attachment_data_container<Vertex> (m_aOrigCoord);

		UG_COND_THROW (	positions.size() != origCoords.size(),
		                "Attachment containers have incompatible sizes");

		for(size_t i = 0; i < positions.size(); ++i){
			positions[i][m_newCoordIndex] = origCoords[i];
		}

		mg->detach_from_vertices (m_aOrigCoord);
	}

private:
	apos_t		m_aPos;
	ANumber		m_aNewCoord;
	ANumber		m_aOrigCoord;
	int			m_newCoordIndex;
};

}//	end of namespace

#endif