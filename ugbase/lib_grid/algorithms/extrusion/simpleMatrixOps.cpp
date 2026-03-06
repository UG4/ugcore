/*
 * simpleMatrixOps.cpp
 *
 *  Created on: 30.12.2024
 *      Author: Markus Knodel
 */

#include <vector>
#include "simpleMatrixOps.h"

namespace ug
{

namespace simpleMatrOps
{

double determinant_2x2(std::vector<std::vector<double>> const & matrix)
{
    return (matrix[0][0] * matrix[1][1]) - (matrix[0][1] * matrix[1][0]);
}

double determinant_3x3(std::vector<std::vector<double>> const & matrix)
{
	double det = 0;
	det += matrix[0][0] * ((matrix[1][1] * matrix[2][2]) - (matrix[1][2] * matrix[2][1]));
	det -= matrix[0][1] * ((matrix[1][0] * matrix[2][2]) - (matrix[1][2] * matrix[2][0]));
	det += matrix[0][2] * ((matrix[1][0] * matrix[2][1]) - (matrix[1][1] * matrix[2][0]));
	return det;
}

std::vector<double> cramerRule(std::vector<std::vector<double>> const & coefficients, std::vector<double> const & constants)
{

	int n = coefficients.size(); // number of unknowns

	std::vector<double> solutions; // vector to store the solutions

	// Calculate the determinant of the coefficients matrix
	double det_coefficients;

	if (n == 2)
	{
		det_coefficients = determinant_2x2(coefficients);
	}
	else if (n == 3)
	{
		det_coefficients = determinant_3x3(coefficients);
	}
	else
	{
		UG_LOG("Error: Cramer's Rule only for 2x2 and 3x3 matrices" << std::endl);
		return solutions;
	}

	// Calculate the solutions for each unknown
	for (int i = 0; i < n; i++)
	{
		// Create a copy of the coefficients matrix and replace the i-th column with the constants vector
		std::vector<std::vector<double>> temp_matrix = coefficients;

		for (int j = 0; j < n; j++)
		{
			temp_matrix[j][i] = constants[j];
		}

		// Calculate the determinant of the modified matrix
		double det_temp;

		if (n == 2)
		{
			det_temp = determinant_2x2(temp_matrix);
		}
		else
		{
			det_temp = determinant_3x3(temp_matrix);
		}

		// Calculate the solution for the i-th unknown
		double solution = det_temp / det_coefficients;
		solutions.push_back(solution);
	}

	return solutions;
}

}


}



