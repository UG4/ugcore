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
#include "lib_disc/reference_element/reference_element.h"

namespace ug{

template <typename TDomain, typename TAlgebra>
class GridFunctionOrdering
{
	using TGridFunction = GridFunction<TDomain, TAlgebra>;
	using TGridFunctionNumberData = GridFunctionNumberData<GridFunction<TDomain, TAlgebra> >;

	using VertexConstIterator = typename TGridFunction::template traits<Vertex>::const_iterator;

public:
	GridFunctionOrdering(SmartPtr<TGridFunction> spGridFct, const char* name)
	{
		m_u = spGridFct->clone_without_values();
		m_name = name;

		std::vector<DoFIndex> ind(1);
		size_t k = 0;
		for(VertexConstIterator iter = m_u->template begin<Vertex>(); iter != m_u->template end<Vertex>(); ++iter)
		{
		//	get vertex
			Vertex* vrt = *iter;

		//	get vector holding all indices on the vertex
			m_u->inner_dof_indices(vrt, 0, ind);
			DoFRef(*m_u, ind[0]) = k++;
		}
	}

	SmartPtr<TGridFunctionNumberData> get(){
		return SmartPtr<TGridFunctionNumberData>(new TGridFunctionNumberData(m_u, m_name));
	}

private:
	SmartPtr<TGridFunction> m_u;
	const char* m_name;
};

} //namespace
