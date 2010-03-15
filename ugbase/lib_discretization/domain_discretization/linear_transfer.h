/*
 * linear_transfer.h
 *
 *  Created on: 04.12.2009
 *      Author: andreasvogel
 */

#ifndef LINEAR_TRANSFER_H_
#define LINEAR_TRANSFER_H_

#include <iostream>
#include "lib_grid/lib_grid.h"
#include "lib_algebra/lib_algebra.h"
#include "lib_discretization/lib_discretization.h"
#include "common/common.h"

namespace ug {


bool assemble_interpolation(	Matrix& mat, DoFPattern& pattern,
								geometry_traits<VertexBase>::iterator iterBegin,
								geometry_traits<VertexBase>::iterator iterEnd, // iterators for toLevel
								int fromLevel, int toLevel)
{
	if(fromLevel != toLevel - 1)
	{
		std::cout << "ERROR in assemble_interpolation: Can not construct Interpolation from Level " << fromLevel << " to Level " << toLevel << "." << std::endl;
		return false;
	}
	Grid* grid = pattern.get_assigned_subset()->get_assigned_grid();
	MultiGrid* mg = dynamic_cast<MultiGrid*>(grid);
	if(mg == NULL || mg->num_levels() <= 1)
	{
		std::cout << "Pattern is not a MultiGrid pattern or contains only one level. Can not create interpolation" << std::endl;
		return false;
	}

	mat.create_matrix(pattern.num_dofs(toLevel), pattern.num_dofs(fromLevel));

	uint num_func = pattern.num_comp();

	double *values = new double[num_func];
	int *ncols = new int[num_func];
	for(uint i = 0; i < num_func; ++i)
	{
		ncols[i] = 1;
	}
	int *rows = new int[num_func];
	int *cols = new int[num_func];


	for(geometry_traits<VertexBase>::iterator iter = iterBegin; iter != iterEnd; ++iter)
	{
		for(uint nr_func = 0; nr_func < num_func; nr_func++)
		{
			rows[nr_func] = pattern.get_index(*iter, nr_func); // to index
		}

		if(IsBoundaryVertex2D(*grid, *iter)) continue;
		GeometricObject* geomObj = mg->get_parent(*iter);

		VertexBase* vert = dynamic_cast<VertexBase*>(geomObj);
		if(vert != NULL)
		{
			if(IsBoundaryVertex2D(*grid, vert)) continue;
			for(uint nr_func = 0; nr_func < num_func; nr_func++)
			{
				cols[nr_func] = (int) pattern.get_index(vert, nr_func); // from index
				values[nr_func] = 1.0;
			}
			mat.add_values(num_func, ncols, rows, cols, values);
			continue;
		}

		Edge* edge = dynamic_cast<Edge*>(geomObj);
		if(edge != NULL)
		{
			for(int i = 0; i < 2; ++i)
			{
				vert = edge->vertex(i);
				if(IsBoundaryVertex2D(*grid, vert)) continue;
				for(uint nr_func = 0; nr_func < num_func; nr_func++)
				{
					cols[nr_func] = (int) pattern.get_index(vert, nr_func); // from index
					values[nr_func] = 0.5;
				}
				mat.add_values(num_func, ncols, rows, cols, values);
			}
			continue;
		}

		Quadrilateral* quad = dynamic_cast<Quadrilateral*>(geomObj);
		if(quad != NULL)
		{
			for(int i = 0; i < 4; ++i)
			{
				vert = quad->vertex(i);
				if(IsBoundaryVertex2D(*grid, vert)) continue;
				for(uint nr_func = 0; nr_func < num_func; nr_func++)
				{
					cols[nr_func] = (int) pattern.get_index(vert, nr_func); // from index
					values[nr_func] = 0.25;
				}
				mat.add_values(num_func, ncols, rows, cols, values);
			}
			continue;
		}

		std::cout << "ERROR in assemble_interpolation: Element Father not detected." << std::endl;
		return false;
	}

	return true;
}

}

#endif /* LINEAR_TRANSFER_H_ */
