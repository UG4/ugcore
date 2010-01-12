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


enum GridType {
	GT_GRID = 1,
	GT_MULTIGRID,
	NUM_GRID_TYPES
};

template <int d>
class Domain {
	public:
		const static int dim = d;
		typedef MathVector<dim> position_type;
		typedef Attachment<position_type> position_attachment_type;
		typedef Grid::VertexAttachmentAccessor<position_attachment_type> position_accessor_type;

	public:
		Domain(Grid& grid, SubsetHandler& sh, position_attachment_type& aaPos);
		Domain(MultiGrid& grid, MGSubsetHandler& sh, position_attachment_type& aaPos);

		inline Grid* get_grid();
		inline MultiGrid* get_multigrid();
		inline GridType get_grid_type();
		inline SubsetHandler* get_subset_handler();
		inline unsigned int get_dim();
		inline position_attachment_type* get_position_attachment();
		inline position_accessor_type get_position_accessor();

	protected:
		GridType _grid_type;
		Grid* _grid;
		MultiGrid* _mg;

		SubsetHandler* _sh;
		MGSubsetHandler* _mg_sh;

		position_attachment_type* _aPos;
		position_accessor_type _aaPos;
};

template <int d>
Domain<d>::Domain(Grid& grid, SubsetHandler& sh, position_attachment_type& aPos) : _aaPos(grid, aPos)
{
	_grid = &grid;
	_grid_type = GT_GRID;
	_sh = &sh;
	_aPos = &aPos;
}

template <int d>
Domain<d>::Domain(MultiGrid& mg, MGSubsetHandler& sh, position_attachment_type& aPos) : _aaPos(mg, aPos)
{
	_mg = &mg;
	_grid_type = GT_MULTIGRID;
	_mg_sh = &sh;
	_aPos = &aPos;
}

template <int d>
inline
GridType Domain<d>::get_grid_type()
{
	return _grid_type;
}

template <int d>
inline
Grid* Domain<d>::get_grid()
{
	return _grid;
}

template <int d>
inline
MultiGrid* Domain<d>::get_multigrid()
{
	return _mg;
}

template <int d>
inline
SubsetHandler* Domain<d>::get_subset_handler()
{
	return _sh;
}

template <int d>
inline
unsigned int Domain<d>::get_dim()
{
	return d;
}

template <int d>
inline
typename Domain<d>::position_attachment_type* Domain<d>::get_position_attachment()
{
	return _aPos;
}

template <int d>
inline
typename Domain<d>::position_accessor_type Domain<d>::get_position_accessor()
{
	return _aaPos;
}


}


#endif /* __H__LIBDISCRETIZATION__DOMAIN__ */
