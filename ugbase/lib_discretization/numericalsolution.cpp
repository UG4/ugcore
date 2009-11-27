/*
 * numericalsolution.cpp
 *
 *  Created on: 12.05.2009
 *      Author: andreasvogel
 */

#include "numericalsolution.h"

namespace ug{

		NumericalSolution::NumericalSolution(std::string shortname, std::string longname, std::string description, NumericalSolutionPattern& pattern, Grid& grid)
		{
			m_shortname = shortname;
			m_longname = longname;
			m_description = description;
			m_grid = &grid;

			_pattern = &pattern;

			_GridVector = new Vector;
			_GridVector->create_vector(_pattern->num_doubles());

			return;
		}

		bool NumericalSolution::set_values(ValueFunction fct, int nr_func, SubsetHandler& sh, uint subsetIndex)
		{
			VertexBaseIterator iterBegin, iterEnd, iter;
			iterBegin = sh.begin<VertexBase>(subsetIndex);
			iterEnd = sh.end<VertexBase>(subsetIndex);

			Grid* grid = sh.get_assigned_grid();
			Grid::VertexAttachmentAccessor<APosition> aaPos(*grid, aPosition);

			int index;
			number val;
			MathVector<3> corner;
			std::vector<int> indices;
			std::vector<number> values;

			/* loop over all Vertices */
			int nvertices = 0;
			for(iter = iterBegin; iter != iterEnd; iter++)
			{
				VertexBase *vert = *iter;
				corner = aaPos[vert];

				index = _pattern->get_index(vert, nr_func);
				fct(corner, val);
				values.push_back(val);
				indices.push_back(index);

				nvertices++;
			}

			double* valueArray = new double[indices.size()];
			int* indexArray = new int[indices.size()];

			for(int i=0; i< indices.size(); i++)
			{
				valueArray[i] = values[i];
				indexArray[i] = indices[i];
			}

			if(_GridVector->set_values((int) indices.size(), indexArray, valueArray) != true)
				return false;
			delete valueArray;
			delete indexArray;

			return true;
		}

		Vector* NumericalSolution::GridVector()
		{
			return _GridVector;
		}

		NumericalSolutionPattern* NumericalSolution::get_pattern()
		{
			return _pattern;
		}


		bool NumericalSolution::print()
		{
			std::cout << "shortname: \"" << m_shortname << "\", ";
			std::cout << "longname: \"" << m_longname <<"\", ";
			std::cout << "description: \"" << m_description << "\"" << std::endl;
			return true;
		}

		NumericalSolution::~NumericalSolution()
		{
			delete _GridVector;
		}

		std::string NumericalSolution::get_name(int comp)
		{
			return _pattern->get_name(comp);
		}



}
