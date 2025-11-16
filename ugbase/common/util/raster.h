/*
 * Copyright (c) 2016:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG_raster
#define __H__UG_raster

#include "common/math/ugmath_types.h"

namespace ug{

///	Generic raster for arbitrary dimensions.
/** The raster provides access to raster-data in arbitrary dimensions.
 * Common use cases are e.g. 1d/2d/3d image data or density distributions.
 *
 * A raster consists of regularily distributed nodes, which form a structured
 * d-dimensional grid. The space in between those nodes can be interpreted as cells.
 * Through 'min_corner' and 'extension' the spacial properties of the raster are
 * defined. A point in d-dimensional space can thus be mapped to a cell in the raster.
 * This is, e.g., used for interpolation algorithms.
 *
 * Values are stored in nodes. Direct access to the data in such nodes is provided
 * as well as interpolation algorithms that provide values at arbitrary points
 * in d-dimensional space.
 *
 * To ease usage, e.g., in scripting languages, the raster provides a 'selected_node'
 * and a 'cursor'. The first is an index to a currently selected node, the latter
 * is a point in space at which interpolation can be performed.
 * The usage of both is of course optional. The value of the
 * selected_node can be accessed through 'selected_node_value'. The value at
 * the cursor can be interpolated through 'interpolate_at_cursor'. Please note,
 * that the more flexible methods 'node_value(MultiIndex)' and
 * 'interpolate(Coordinate)' are recommended for most use cases.
 *
 * \param T		The template parameter 'T' specifies the underlying data-type.
 *				It has to support the following operations, if 'interpolate'
 *				shall be used:
 *					- T& operator += (const T&)
 *					- T& operator -= (const T&)
 *					- T& operator *= (number)
 *				It furthermore has to support std::numeric_limits<T>::max()
 */

template <class T, int TDIM>
class Raster{
	public:
		class MultiIndex {
			public:
				MultiIndex();
				MultiIndex(size_t i);
				int	dim () const;
				void 		set (size_t i);
				size_t&		operator[] (int d);
				size_t		operator[] (int d) const;

				friend std::ostream& operator << (std::ostream& o, const MultiIndex& mi)
				{
					o << "(";
					for(size_t d = 0; d < TDIM; ++d){
						o << mi[d];
						if(d + 1 < TDIM)
							o << ", ";
					}
					o << ")";
					return o;
				}

			private:
				size_t	m_ind[TDIM];
		};

		class Coordinate {
			public:
				Coordinate();
				Coordinate(number c);
				Coordinate(const MathVector<TDIM, number>& v);

				int	dim () const;
				void 		set (number c);
				number&		operator[] (int d);
				number		operator[] (int d) const;
				Coordinate& operator+= (const Coordinate& c);
				Coordinate& operator-= (const Coordinate& c);
				Coordinate& operator*= (number s);
				
				friend std::ostream& operator << (std::ostream& o, const Coordinate& coord)
				{
					o << "(";
					for(size_t d = 0; d < TDIM; ++d){
						o << coord[d];
						if(d + 1 < TDIM)
							o << ", ";
					}
					o << ")";
					return o;
				}

			private:
				number	m_coord[TDIM];
		};

	///	Creates an empty raster that has to be initialized before use.
	/**	Use the methods 'load_from_asc' or instead 'set_num_nodes', 'create',
	 *	'set_extension', and 'set_min_corner' to initialize the raster.*/
		Raster ();

	///	Creates the new raster by copying the contents of the given raster.
		Raster (const Raster& raster);

	///	Creates a new raster with the specified number of nodes.
	/** Internally calls 'set_num_nodes' and 'create'. The ruster can thus be used
	 *	right after construction. The extension is chosen to be numNodes[d]-1 in each
	 *	dimension, i.e., the distance between neighbored nodes is '1' by default.
	 *	The min-corner is at the origin. You may change those through
	 *	'set_extension' and 'set_min_corner'.*/
		Raster (const MultiIndex& numNodes);
		
	///	Creates a new raster with the specified number of nodes and the specified extension.
	/** Internally calls 'set_num_nodes', 'create', 'set_extension', and 'set_min_corner.
	 *	'minCorner' is optional and defaults to the origin.*/
		Raster (const MultiIndex& numNodes,
				const Coordinate& extension,
				const Coordinate& minCorner = Coordinate(0));

		~Raster ();

		Raster& operator= (const Raster& raster);

		void load_from_asc (const char* filename);
		void save_to_asc (const char* filename) const;

		int dim () const;

	///	sets the number of nodes that shall be used by the raster.
	/**	At least one node is required per dimension. After having set all
	 *	node-numbers, one should call 'create' to actually create the grid.*/
		void set_num_nodes (int dim, size_t num);
		void set_num_nodes (const MultiIndex& mi);

	///	returns the total number of nodes in the raster
		size_t num_nodes_total () const;
	///	returns the number of nodes in the specified dimension
		size_t num_nodes (int dim) const;
	///	returns the number of nodes for each dimension in a multi-index.
		const MultiIndex& num_nodes () const;

	///	creates the raster. Call this method after 'set_num_nodes' has been called for each dimension.
		void create ();

	///	returns the value at the given multi-index (read/write)
		T& node_value (const MultiIndex& mi);

	///	returns the value at the given multi-index (read only)
		T node_value (const MultiIndex& mi) const;

	///	sets the min corner of the raster. Used for interpolation at cursor.
	/** \{ */
		void set_min_corner (int dim, number coord);
		void set_min_corner (const Coordinate& coord);
	/** \} */

	///	returns the min-corner of the raster
		const Coordinate& min_corner () const;

	///	returns the coordinate of the min-corner of the raster for the given dimension
		number min_corner (int dim) const;

	///	sets the extension of the raster. Used for interpolation at cursor.
	/** \{ */
		void set_extension (int dim, number ext);
		void set_extension (const Coordinate& coord);
	/** \} */

	///	returns the extension of the raster
		const Coordinate& extension () const;

	///	returns the extension of the raster for the given dimension
		number extension (int dim) const;

	///	sets the value that shall be considered as 'no-data-value'
		void set_no_data_value(T val);

	///	returns the value that shall be considered 'no-data-value'
		T no_data_value() const;

	///	interpolates the value with the given order at the given coordinate
		T interpolate (const Coordinate& coord, int order) const;

	///	blurs (smoothens) the values by repeatedly averaging between direct neighbors
		void blur(T alpha, size_t iterations);


	/// Creates and runs the specified kernel on all nodes and returns its result
	/** The class TKernel has to feature a default constructor, a type definition 'result_t',
	 * and a method 'result_t result() const'.
	 *
	 * Like all kernels, it furthermore has to specify a method
	 * \code
	 * void operator () (Raster<T, TDIM>& raster,
	 *			   		 const typename Raster<T, TDIM>::MultiIndex& cur);
	 * \endcode
	 */
		template <class TKernel>
		typename TKernel::result_t
		run_on_all();

	/// Runs the specified kernel on all nodes
	/** TKernel can either be a function or a class with an 'operator()'.
	 * The signature should be as follows:
	 * \code
	 * void Func (Raster<T, TDIM>& raster,
	 *			  const typename Raster<T, TDIM>::MultiIndex& cur);
	 * \endcode
	 * or
	 * \code
	 * void operator () (Raster<T, TDIM>& raster,
	 *			   		 const typename Raster<T, TDIM>::MultiIndex& cur);
	 * \endcode
	 */
		template <class TKernel>
		void run_on_all(TKernel& kernel);


	/// Creates and runs the specified kernel on all direct neighbors of a node and returns its result
	/** The class TKernel has to feature a default constructor, a type definition 'result_t',
	 * and a method 'result_t result() const'.
	 *
	 * Like all kernels, it furthermore has to specify a method
	 * \code
	 * void operator () (Raster<T, TDIM>& raster,
	 *			   		 const typename Raster<T, TDIM>::MultiIndex& cur);
	 * \endcode
	 */
		template <class TKernel>
		typename TKernel::result_t
		run_on_nbrs(const MultiIndex& center);

	/// Runs the specified kernel on all direct neighbors of a node
	/** TKernel can either be a function or a class with an 'operator()'.
	 * The signature should be as follows:
	 * \code
	 * void Func (Raster<T, TDIM>& raster,
	 *			  const typename Raster<T, TDIM>::MultiIndex& cur);
	 * \endcode
	 * or
	 * \code
	 * void operator () (Raster<T, TDIM>& raster,
	 *			   		 const typename Raster<T, TDIM>::MultiIndex& cur);
	 * \endcode
	 */
		template <class TKernel>
		void run_on_nbrs(const MultiIndex& center, TKernel& kernel);

	////////////////////////////////////////////////////////////////////////////
	//	Convenience methods for easy scripting below.

	///	Select a node. 'selected_node_value' provides read/write access to the selected node.
	/** \{ */
		void select_node (int dim, size_t index);
		void select_node (const MultiIndex& mi);
	/** \} */

	///	sets the value of the selected node
		void set_selected_node_value (T val);

	///	returns the value of the selected node
		T selected_node_value () const;

	///	Set the coordinate of the cursor. The cursor can be used to interpolate values.
	/** \{ */
		void set_cursor (int dim, number coord);
		void set_cursor (const Coordinate& coord);
	/** \} */

	///	interpolates the value with the given order at the cursor position
		T interpolate_at_cursor (int order) const;

	private:
		template <class TKernel>
		void run_on_all(const MultiIndex& start, TKernel& kernel, int curDim);

		template <class TKernel>
		void run_on_nbrs(const MultiIndex& center, TKernel& kernel, int curDim);

		size_t data_index (
					const MultiIndex& mi,
					int curDim = TDIM - 1,
					size_t curVal = 0) const;

		void update_num_nodes_total();
		void update_cell_extension ();
		void update_cell_extension (int dim);

		T interpolate_linear (
				const MultiIndex& minNodeInd,
				Coordinate& localCoord,
				int curDim = TDIM) const;

		T*			m_data;
		MultiIndex	m_numNodes;
		MultiIndex	m_selNode;
		Coordinate	m_minCorner;
		Coordinate	m_extension;
		Coordinate	m_cellExtension;
		Coordinate	m_cursor;
		size_t		m_numNodesTotal;
		T			m_noDataValue;
};

}//	end of namespace

////////////////////////////////
// include implementation
#include "raster_impl.hpp"

#endif	//__H__UG_raster
