/*
 * discretizationscheme.cpp
 *
 *  Created on: 26.06.2009
 *      Author: andreasvogel
 */

#include "discretizationscheme.h"

namespace ug{

//////////////////////
////// Dirichlet
//////////////////////

bool DirichletValues::add_dirichlet_nodes(NumericalSolution& u, int nr_func, DirichletBNDCond* dirichbnd, SubsetHandler& sh, uint subsetIndex)
{
	VertexBaseIterator iterBegin, iterEnd, iter;
	iterBegin = sh.begin<VertexBase>(subsetIndex);
	iterEnd = sh.end<VertexBase>(subsetIndex);

	Grid* grid = sh.get_assigned_grid();
	Grid::VertexAttachmentAccessor<APosition> aaPos(*grid, aPosition);

	int index;
	number val;
	MathVector<3> corner;

	/* loop over all Vertices */
	int nvertices = 0, ndirnodes = 0;
	for(iter = iterBegin; iter != iterEnd; iter++)
	{
		VertexBase *vert = *iter;
		corner = aaPos[vert];
		nvertices++;

		if(IsBoundaryVertex2D(*grid, vert))
		{
			index = (int) u.get_pattern()->get_index(vert, nr_func);
			dirichbnd->BNDValueFunction(corner, val);
			m_vector_values.push_back(val);
			m_vector_indices.push_back(index);
			m_matrixrow_indices.push_back(index);
			ndirnodes++;
		}
	}
	//std::cout << nvertices << " Vertices checked: Nr. Dirichlet nodes found:" << ndirnodes << "(now total: "<< m_vector_indices.size() <<")"<< std::endl;

	return true;
}

bool DirichletValues::set_values(Vector& vec)
{
	double* valueArray = new double[m_vector_indices.size()];
	int* indexArray = new int[m_vector_indices.size()];

	for(int i=0; i<m_vector_indices.size(); i++)
	{
		valueArray[i] = m_vector_values[i];
		indexArray[i] = m_vector_indices[i];
	}

	if(vec.set_values(m_vector_indices.size(), indexArray, valueArray) != true)
		return false;
	delete valueArray;
	delete indexArray;

	std::cout << m_vector_indices.size() << " Boundary nodes set to Dirichlet value" << std::endl;
	return true;
}

bool DirichletValues::set_zero_values(Vector& vec)
{
	double* valueArray = new double[m_vector_indices.size()];
	int* indexArray = new int[m_vector_indices.size()];

	for(int i=0; i<m_vector_indices.size(); i++)
	{
		valueArray[i] = 0.0;
		indexArray[i] = m_vector_indices[i];
	}

	if(vec.set_values(m_vector_indices.size(), indexArray, valueArray) != true)
		return false;
	delete valueArray;
	delete indexArray;

	std::cout << m_vector_indices.size() << " Boundary nodes set to zero value" << std::endl;
	return true;
}


bool DirichletValues::set_rows(Matrix& mat)
{
	int* indexArray = new int[m_matrixrow_indices.size()];

	for(int i=0; i<m_matrixrow_indices.size(); i++)
	{
		indexArray[i] = m_matrixrow_indices[i];
	}

	if(mat.set_dirichletrows(m_matrixrow_indices.size(), indexArray) != true)
		return false;
	delete indexArray;

	std::cout << m_matrixrow_indices.size() << " Matrix rows set to Identity row" << std::endl;
	return true;
}

}
