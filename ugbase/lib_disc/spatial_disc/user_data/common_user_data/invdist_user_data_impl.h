/*
 * Copyright (c) 2015:  G-CSC, Goethe University Frankfurt
 * Author: Dmitry Logashenko
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

/*
 * Implementation of the Inverse Distance Weighting (IDW) interpolation for data sets
 */

namespace ug {

/**
 * Computes the interpolation basing on all the interpolation point.
 */
template <int WDim, typename TPntIterator, typename TData>
void IDWInterpolation<WDim, TPntIterator, TData>::compute
(
	data_type & res, ///< interpolated value
	const MathVector<dim> & pos, ///< geometric position where to interpolate
	t_pnt_iter pnt_beg, ///< the first interpolation point
	t_pnt_iter pnt_end, ///< delimiter of the iterpolation points
	number order, ///< order of the interpolation
	number small_dist ///< distance at which we do not distinguish the points
)
{
	if (pnt_beg == pnt_end)
		UG_THROW ("IDWInterpolation: Cannot interpolate using 0 data points!");
	
	data_type sum = 0;
	number factor = 0;
	for (t_pnt_iter pnt = pnt_beg; pnt != pnt_end; ++pnt)
	{
		number dist;
		if ((dist = VecDistance (pos, pnt->pos)) < small_dist)
		{ /* We are at a data point: */
			res = pnt->value;
			return;
		}
		dist = pow (dist, order);
		data_type value = pnt->value;
		value /= dist;
		sum += value;
		factor += 1 / dist;
	}
	sum /= factor;
	res = sum;
}

/**
 * Computes the interpolation basing on the interpolation points in a given ball.
 */
template <int WDim, typename TPntIterator, typename TData>
void IDWInterpolation<WDim, TPntIterator, TData>::compute
(
	data_type & res, ///< interpolated value
	const MathVector<dim> & pos, ///< geometric position where to interpolate
	number R, ///< radius of the ball (if 0 then the whole space)
	t_pnt_iter pnt_beg, ///< the first interpolation point
	t_pnt_iter pnt_end, ///< delimiter of the iterpolation points
	number order, ///< order of the interpolation
	number small_dist ///< distance at which we do not distinguish the points
)
{
	if (R == 0.0) // a special case: we do not consider the ball
	{
		IDWInterpolation<dim, t_pnt_iter, data_type>::compute
								(res, pos, pnt_beg, pnt_end, order, small_dist);
		return;
	}
	
	if (pnt_beg == pnt_end)
		UG_THROW ("IDWInterpolation: Cannot interpolate using 0 data points!");
	
	data_type sum = 0;
	number factor = 0;
	bool no_data_in_ball = true;
	for (t_pnt_iter pnt = pnt_beg; pnt != pnt_end; ++pnt)
	{
		number dist;
		if ((dist = VecDistance (pos, pnt->pos)) < small_dist)
		{ /* We are at a data point: */
			res = pnt->value;
			return;
		}
		if (dist > R) continue;
		no_data_in_ball = false;
		dist = pow (dist, order);
		data_type value = pnt->value;
		value /= dist;
		sum += value;
		factor += 1 / dist;
	}
	if (no_data_in_ball)
		UG_THROW ("IDWInterpolation: No interpolation points in the ball with R = " << R
					<< " and center at" << pos << ".");
	sum /= factor;
	res = sum;
}

/**
 * Loads interpolation points from a given stream.
 */
template <int WDim, typename TData>
void IDWUserData<WDim, TData>::load_data_from (std::istream & input)
{
	std::string input_line;
	std::istringstream line_stream;
	MathVector<dim, number> pos;
	data_type value;
	
	line_stream.exceptions (std::istream::failbit | std::istream::badbit);
	
	while (! input.eof ())
	{
		if (input.fail ())
			UG_THROW ("IDWUserData: Could not load the interpolation points from the file!");
		std::getline (input, input_line);
		if (input_line.length () == 0)
			continue;
		try
		{
			line_stream.str (input_line);
			read_plain_txt (line_stream, pos);
			line_stream >> value;
		}
		catch (std::istream::failure& e)
		{
			UG_THROW ("IDWUserData: Failed to parse line '" << input_line << "' in the input file.");
		};
		
		this->append (pos, value);
	}
}

/**
 * Loads interpolation points from a given file.
 */
template <int WDim, typename TData>
void IDWUserData<WDim, TData>::load_data_from (const char * file_name)
{
	std::ifstream input (file_name, std::ifstream::in);
	if (input.fail ())
		UG_THROW ("IDWUserData: Cannot open data file '" << file_name << "' for input!");
	this->load_data_from (input);
}

} // end namespace ug

/* End of File */
