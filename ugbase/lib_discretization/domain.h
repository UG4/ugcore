/*
 * domain.h
 *
 *  Created on: 17.12.2009
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__DOMAIN__
#define __H__LIBDISCRETIZATION__DOMAIN__

#include "lib_grid/lg_base.h"

namespace ug{

//	predeclarations
class DistributedGridManager;

/**
 *
 * A Domain collects and exports relevant informations about the
 * physical domain, that will be discretized. It will be used as
 * a template Parameter in several classes to destinguish at compile-time
 * between needed types and parameters.
 *
 * An Implementation of the Domain interface has to fulfill the following requirements:
 *
 * const static int dim = ...
 * typedef ... position_type
 * typedef ... position_attachment_type
 * typedef ... position_accessor_type
 *
 * [ may be extended in future ]
 *
 */



template <int d, typename TGrid, typename TSubsetHandler>
class Domain {
	public:
		// world dimension
		static const int dim = d;

		// type of position coordinates (e.g. MathVector<dim>)
		typedef MathVector<dim> position_type;

		// grid type
		typedef TGrid grid_type;

		// subset handler type
		typedef TSubsetHandler subset_handler_type;

		// type of position attachement (since spacial positions of vertices are attached to the topological grid)
		typedef Attachment<position_type> position_attachment_type;

		// type of accessor to the position data attachement
		typedef Grid::VertexAttachmentAccessor<position_attachment_type> position_accessor_type;

		typedef DistributedGridManager distributed_grid_manager_type;

	public:
		Domain(uint options = GRIDOPT_STANDARD_INTERCONNECTION) :
			m_grid(options), m_sh(m_grid), m_distGridMgr(NULL)
			{
				m_aPos = GetDefaultPositionAttachment<position_attachment_type>();

				// let position accessor access Vertex Coordinates
				if(!m_grid.template has_attachment<VertexBase>(m_aPos))
					m_grid.template attach_to<VertexBase>(m_aPos);

				m_aaPos.access(m_grid, m_aPos);

#ifdef UG_PARALLEL
				m_distGridMgr = new DistributedGridManager(m_grid);
#endif
			}

		~Domain()
		{
			if(m_distGridMgr)
				delete m_distGridMgr;
		}

		inline TGrid& get_grid() {return m_grid;};
		inline const TGrid& get_grid() const {return m_grid;};

		inline TSubsetHandler& get_subset_handler() {return m_sh;};
		inline const TSubsetHandler& get_subset_handler() const {return m_sh;};

		inline int get_dim() const {return d;};

		inline position_attachment_type& get_position_attachment() {return m_aPos;};
		inline const position_attachment_type& get_position_attachment() const {return m_aPos;};

		inline position_accessor_type& get_position_accessor() {return m_aaPos;};
		inline const position_accessor_type& get_position_accessor() const {return m_aaPos;};

		inline DistributedGridManager* get_distributed_grid_manager()	{return m_distGridMgr;}

	protected:
		TGrid m_grid;

		TSubsetHandler m_sh;

		position_attachment_type m_aPos;
		position_accessor_type m_aaPos;

	//	for parallelization only
		DistributedGridManager*	m_distGridMgr;
};

} // end namespace ug


#endif /* __H__LIBDISCRETIZATION__DOMAIN__ */
