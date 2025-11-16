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
#include "lib_disc/function_spaces/grid_function_user_data.h"

#include "lib_disc/io/vtkoutput.h" ////IteratorProvider

namespace ug{

template <typename TDomain, typename TAlgebra>
class GridPointsOrdering
{
	using TGridFunction = GridFunction<TDomain, TAlgebra>;
	using TGridFunctionNumberData = GridFunctionNumberData<GridFunction<TDomain, TAlgebra> >;

	using AVrtIndex = Attachment<int>;

public:
	GridPointsOrdering(SmartPtr<TGridFunction> spGridFct, const char* name)
	{

		m_u = spGridFct->clone_without_values();
		m_name = name;
		m_n = 0;
		m_numVert = 0;
		m_numElem = 0;
		m_numConn = 0;

		create_vtkoutput_ordering();

		std::vector<size_t> index(m_numVert);

		size_t k = 0;

		using const_iterator = typename IteratorProvider<TGridFunction>::template traits<ug::Vertex>::const_iterator;
		const_iterator iterBegin = IteratorProvider<TGridFunction>::template begin<ug::Vertex>(*m_u, -1);
		const_iterator iterEnd = IteratorProvider<TGridFunction>::template end<ug::Vertex>(*m_u, -1);
		for( ; iterBegin != iterEnd; ++iterBegin, ++k)
		{
			ug::Vertex* v = *iterBegin;
			index[k] = m_aaVrtIndex[v];
		}

		k = 0;
		std::vector<DoFIndex> ind(1);

		iterBegin = IteratorProvider<TGridFunction>::template begin<ug::Vertex>(*m_u, -1);
		iterEnd = IteratorProvider<TGridFunction>::template end<ug::Vertex>(*m_u, -1);
		for( ; iterBegin != iterEnd; ++iterBegin, ++k)
		{
			ug::Vertex* v = *iterBegin;

		//	get vector holding all indices on the vertex
			m_u->inner_dof_indices(v, 0, ind);
			DoFRef(*m_u, ind[0]) = index[k];
		}
	}


	template <typename TElem>
	void count_sizes()
	{
//		get the grid associated to the solution
		Grid& grid = *m_u->domain()->grid();

	//	get reference element
		using ref_elem_type = typename reference_element_traits<TElem>::reference_element_type;

	//	number of corners of element
		static constexpr int numCo = ref_elem_type::numCorners;

	//	get iterators
		using const_iterator = typename IteratorProvider<TGridFunction>::template traits<TElem>::const_iterator;
		const_iterator iterBegin = IteratorProvider<TGridFunction>::template begin<TElem>(*m_u, -1);
		const_iterator iterEnd = IteratorProvider<TGridFunction>::template end<TElem>(*m_u, -1);

	//	loop elements
		for( ; iterBegin != iterEnd; ++iterBegin)
		{
		//	get the element
			TElem *elem = *iterBegin;

		//	count number of elements and number of connections;
		//	handle octahedrons separately by splitting into a top and bottom pyramid
			if(ref_elem_type::REFERENCE_OBJECT_ID != ROID_OCTAHEDRON)
			{
				++m_numElem;
				m_numConn += numCo;
			}
			else
			{
			// 	counting top and bottom pyramid
				m_numElem += 2;
				m_numConn += 10;
			}

		//	loop vertices of the element
			for(int i = 0; i < numCo; ++i)
			{
			//	get vertex of the element
				Vertex* v = GetVertex(elem, i);

			//	if this vertex has already been counted, skip it
				if(grid.is_marked(v)) continue;

			// count vertex and mark it
				++m_numVert;
				grid.mark(v);
			}
		}
	}


	template <typename TElem>
	void number_points_elementwise()
	{
//		get the grid associated to the solution
		Grid& grid = *m_u->domain()->grid();

	//	get reference element
		using ref_elem_type = typename reference_element_traits<TElem>::reference_element_type;
	//	get iterators
		using const_iterator = typename IteratorProvider<TGridFunction>::template traits<TElem>::const_iterator;
		const_iterator iterBegin = IteratorProvider<TGridFunction>::template begin<TElem>(*m_u, -1);
		const_iterator iterEnd = IteratorProvider<TGridFunction>::template end<TElem>(*m_u, -1);

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
				m_aaVrtIndex[v] = m_n++;
			}
		}
	}

	void create_vtkoutput_ordering(){
//		check functions
		for(size_t fct = 0; fct < m_u->num_fct(); ++fct)
		{
		//	check if function is defined everywhere
			if(!m_u->is_def_everywhere(fct))
				UG_THROW("only serial case implemented!");
		}

//		get the grid associated to the solution
		Grid& grid = *m_u->domain()->grid();

// 		attach help indices
		AVrtIndex aVrtIndex;
		grid.attach_to_vertices(aVrtIndex);
		m_aaVrtIndex.access(grid, aVrtIndex);

		int dim = DimensionOfSubsets(*m_u->domain()->subset_handler());

// 		Count needed sizes for vertices, elements and connections
		try{
//			reset all marks
			grid.begin_marking();

//			switch dimension
			switch(dim)
			{
				case 0: count_sizes<Vertex>(); break;
				case 1: count_sizes<RegularEdge>();
						count_sizes<ConstrainingEdge>(); break;
				case 2: count_sizes<Triangle>();
						count_sizes<Quadrilateral>();
						count_sizes<ConstrainingTriangle>();
						count_sizes<ConstrainingQuadrilateral>(); break;
				case 3: count_sizes<Tetrahedron>();
						count_sizes<Pyramid>();
						count_sizes<Prism>();
						count_sizes<Octahedron>();
						count_sizes<Hexahedron>(); break;
				default: UG_THROW(name() << "::create_vtkoutput_ordering: Dimension " << dim << " is not supported.");
			}

//			signal end of marking
			grid.end_marking();
		}
		UG_CATCH_THROW(name() << "::create_vtkoutput_ordering: Can not count piece sizes.");

//		start marking of vertices
		grid.begin_marking();

//		switch dimension
		if(m_numVert > 0){
			switch(dim){
				case 0: number_points_elementwise<Vertex>(); break;
				case 1: number_points_elementwise<RegularEdge>();
					number_points_elementwise<ConstrainingEdge>(); break;
				case 2: number_points_elementwise<Triangle>();
					number_points_elementwise<Quadrilateral>();
					number_points_elementwise<ConstrainingTriangle>();
					number_points_elementwise<ConstrainingQuadrilateral>(); break;
				case 3: number_points_elementwise<Tetrahedron>();
					number_points_elementwise<Pyramid>();
					number_points_elementwise<Prism>();
					number_points_elementwise<Octahedron>();
					number_points_elementwise<Hexahedron>(); break;
				default: UG_THROW(name() << "::create_vtkoutput_ordering: Dimension " << dim << " is not supported.");
			}
		}
//		signal end of marking the grid
		grid.end_marking();
	}

	SmartPtr<TGridFunctionNumberData> get(){
		return SmartPtr<TGridFunctionNumberData>(new TGridFunctionNumberData(m_u, m_name));
	}

	const char* name() const {return "GridPointsOrdering";}

private:
	SmartPtr<TGridFunction> m_u;
	const char* m_name;
	size_t m_n;
	int m_numVert;
	int m_numElem;
	int m_numConn;

	Grid::VertexAttachmentAccessor<AVrtIndex> m_aaVrtIndex;
};

} //namespace
