/*
 * Copyright (c) 2011-2022:  G-CSC, Goethe University Frankfurt
 * Author: Lukas Larisch
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


#include "lib_disc/domain.h"
#include "lib_disc/function_spaces/grid_function.h"

#include "lib_disc/io/vtkoutput.h"

namespace ug{

template <typename TDomain, typename TAlgebra>
class GridFunctionOrdering
{
	typedef GridFunction<TDomain, TAlgebra> TGridFunction;
public:
	GridFunctionOrdering(SmartPtr<TGridFunction> spGridFct, const char* name)
	{
		u = spGridFct->clone_without_values();
		auto v = getVector(u);
		//v[1] = 0;
	}


	template <typename TElem>
	void count_sizes(Grid &grid, const TGridFunction &u, int si,
		    int& numVert, int& numElem, int& numConn)
	{
	//	get reference element
		typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;

	//	number of corners of element
		static const int numCo = ref_elem_type::numCorners;

	//	get iterators
		typedef typename IteratorProvider<TGridFunction>::template traits<TElem>::const_iterator const_iterator;
		const_iterator iterBegin = IteratorProvider<TGridFunction>::template begin<TElem>(u, si);
		const_iterator iterEnd = IteratorProvider<TGridFunction>::template end<TElem>(u, si);

	//	loop elements
		for( ; iterBegin != iterEnd; ++iterBegin)
		{
		//	get the element
			TElem *elem = *iterBegin;

		//	count number of elements and number of connections;
		//	handle octahedrons separately by splitting into a top and bottom pyramid
			if(ref_elem_type::REFERENCE_OBJECT_ID != ROID_OCTAHEDRON)
			{
				++numElem;
				numConn += numCo;
			}
			else
			{
			// 	counting top and bottom pyramid
				numElem += 2;
				numConn += 10;
			}

		//	loop vertices of the element
			for(int i = 0; i < numCo; ++i)
			{
			//	get vertex of the element
				Vertex* v = GetVertex(elem,i);

			//	if this vertex has already been counted, skip it
				if(grid.is_marked(v)) continue;

			// count vertex and mark it
				++numVert;
				grid.mark(v);
			}
		}
	}


	template <typename TElem>
	void number_points_elementwise(Grid::VertexAttachmentAccessor<Attachment<int> > &aaVrtIndex,
		                 Grid &grid, const TGridFunction &u, int si, int& n)
	{
	//	get reference element
		typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;
	//	get iterators
		typedef typename IteratorProvider<TGridFunction>::template traits<TElem>::const_iterator const_iterator;
		const_iterator iterBegin = IteratorProvider<TGridFunction>::template begin<TElem>(u, si);
		const_iterator iterEnd = IteratorProvider<TGridFunction>::template end<TElem>(u, si);

	//	loop all elements of the subset
		for( ; iterBegin != iterEnd; ++iterBegin)
		{
		//	get the element
			TElem *elem = *iterBegin;

		//	loop vertices of the element
			for(size_t i = 0; i < (size_t) ref_elem_type::numCorners; ++i)
			{
			//	get vertex of element
				Vertex* v = GetVertex(elem, i);

			//	if vertex has already be handled, skip it
				if(grid.is_marked(v)) continue;

			//	mark the vertex as processed
				grid.mark(v);

			//	number vertex
				aaVrtIndex[v] = n++;
			}
		}
	}

	void create_vtkoutput_ordering(){
//		check functions
		bool bEverywhere = true;
		for(size_t fct = 0; fct < u.num_fct(); ++fct)
		{
		//	check if function is defined everywhere
			if(!u.is_def_everywhere(fct)){
				bEverywhere = false;
				UG_THROW("only serial case implemented!");
			}
		}

//		get the grid associated to the solution
		Grid& grid = *u.domain()->grid();

// 		attach help indices
		typedef ug::Attachment<int> AVrtIndex;
		AVrtIndex aVrtIndex;
		Grid::VertexAttachmentAccessor<AVrtIndex> aaVrtIndex;
		grid.attach_to_vertices(aVrtIndex);
		aaVrtIndex.access(grid, aVrtIndex);

		int dim = DimensionOfSubsets(*u.domain()->subset_handler());

//		counters
		int numVert = 0, numElem = 0, numConn = 0;

		int si = -1;

// 		Count needed sizes for vertices, elements and connections
		try{
//			reset all marks
			grid.begin_marking();

//			switch dimension
			switch(dim)
			{
				case 0: count_sizes<Vertex>(grid, u, si, numVert, numElem, numConn); break;
				case 1: count_sizes<RegularEdge>(grid, u, si, numVert, numElem, numConn);
						count_sizes<ConstrainingEdge>(grid, u, si, numVert, numElem, numConn); break;
				case 2: count_sizes<Triangle>(grid, u, si, numVert, numElem, numConn);
						count_sizes<Quadrilateral>(grid, u, si, numVert, numElem, numConn);
						count_sizes<ConstrainingTriangle>(grid, u, si, numVert, numElem, numConn);
						count_sizes<ConstrainingQuadrilateral>(grid, u, si, numVert, numElem, numConn); break;
				case 3: count_sizes<Tetrahedron>(grid, u, si, numVert, numElem, numConn);
						count_sizes<Pyramid>(grid, u, si, numVert, numElem, numConn);
						count_sizes<Prism>(grid, u, si, numVert, numElem, numConn);
						count_sizes<Octahedron>(grid, u, si, numVert, numElem, numConn);
						count_sizes<Hexahedron>(grid, u, si, numVert, numElem, numConn); break;
				default: UG_THROW(name() << "::create_vtkoutput_ordering: Dimension " << dim << " is not supported.");
			}

//			signal end of marking
			grid.end_marking();
		}
		UG_CATCH_THROW(name() << "::create_vtkoutput_ordering: Can not count piece sizes.");


		MGSubsetHandler& sh = *u.domain()->subset_handler();

//		write grid
		int n = 0;

//		start marking of vertices
		grid.begin_marking();

//		switch dimension
		if(numVert > 0){
			switch(dim){
				case 0: number_points_elementwise<Vertex>(aaVrtIndex, grid, u, si, n); break;
				case 1: number_points_elementwise<RegularEdge>(aaVrtIndex, grid, u, si, n);
					number_points_elementwise<ConstrainingEdge>(aaVrtIndex, grid, u, si, n); break;
				case 2: number_points_elementwise<Triangle>(aaVrtIndex, grid, u, si, n);
					number_points_elementwise<Quadrilateral>(aaVrtIndex, grid, u, si, n);
					number_points_elementwise<ConstrainingTriangle>(aaVrtIndex, grid, u, si, n);
					number_points_elementwise<ConstrainingQuadrilateral>(aaVrtIndex, grid, u, si, n); break;
				case 3: number_points_elementwise<Tetrahedron>(aaVrtIndex, grid, u, si, n);
					number_points_elementwise<Pyramid>(aaVrtIndex, grid, u, si, n);
					number_points_elementwise<Prism>(aaVrtIndex, grid, u, si, n);
					number_points_elementwise<Octahedron>(aaVrtIndex, grid, u, si, n);
					number_points_elementwise<Hexahedron>(aaVrtIndex, grid, u, si, n); break;
				default: UG_THROW(name() << "::create_vtkoutput_ordering: Dimension " << dim << " is not supported.");
			}
		}
//		signal end of marking the grid
		grid.end_marking();
	}

	SmartPtr<TGridFunction> get(){
		return u;
	}

	const char* name() const {return "GridFunctionOrdering";}

private:
	SmartPtr<TGridFunction> u;
};

} //namespace
