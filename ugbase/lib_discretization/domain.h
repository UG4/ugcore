/*
 * domain.h
 *
 *  Created on: 17.12.2009
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__DOMAIN__
#define __H__LIBDISCRETIZATION__DOMAIN__

#include "lib_grid/lib_grid.h"

namespace ug{

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

	public:
		Domain(TGrid& grid, TSubsetHandler& sh, position_attachment_type& aPos) :
			m_grid(grid), m_sh(sh), m_aPos(aPos), m_aaPos(grid, aPos)
			{};

		inline TGrid& get_grid();
		inline TSubsetHandler& get_subset_handler();
		inline uint get_dim();
		inline position_attachment_type& get_position_attachment();
		inline position_accessor_type& get_position_accessor();

	protected:
		TGrid& m_grid;

		TSubsetHandler& m_sh;

		position_attachment_type& m_aPos;
		position_accessor_type m_aaPos;
};

template <int d, typename TGrid, typename TSubsetHandler>
inline
TGrid&
Domain<d, TGrid, TSubsetHandler>::
get_grid()
{
	return m_grid;
}

template <int d, typename TGrid, typename TSubsetHandler>
inline
TSubsetHandler&
Domain<d, TGrid, TSubsetHandler>::
get_subset_handler()
{
	return m_sh;
}

template <int d, typename TGrid, typename TSubsetHandler>
inline
uint
Domain<d, TGrid, TSubsetHandler>::
get_dim()
{
	return d;
}

template <int d, typename TGrid, typename TSubsetHandler>
inline
typename Domain<d, TGrid, TSubsetHandler>::position_attachment_type&
Domain<d, TGrid, TSubsetHandler>::
get_position_attachment()
{
	return m_aPos;
}

template <int d, typename TGrid, typename TSubsetHandler>
inline
typename Domain<d, TGrid, TSubsetHandler>::position_accessor_type&
Domain<d, TGrid, TSubsetHandler>::
get_position_accessor()
{
	return m_aaPos;
}


}


#endif /* __H__LIBDISCRETIZATION__DOMAIN__ */
