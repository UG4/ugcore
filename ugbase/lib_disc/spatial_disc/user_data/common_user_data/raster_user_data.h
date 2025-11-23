/*
 * Copyright (c) 2020:  G-CSC, Goethe University Frankfurt
 * Author: Arne Naegel
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

#ifndef __UG__LIB_DISC__SPATIAL_DISC__RASTER_USER_DATA_H__
#define __UG__LIB_DISC__SPATIAL_DISC__RASTER_USER_DATA_H__

#include "common/common.h"
#include "common/util/raster.h"
#include "common/math/ugmath_types.h"
#include "lib_disc/spatial_disc/user_data/user_data.h"
#include "lib_disc/spatial_disc/user_data/std_glob_pos_data.h"

namespace ug{

///! 
template <int dim>
class RasterUserData : public StdGlobPosData<RasterUserData<dim>, number, dim, void>
{

public:
	using TData = number;
	using base_type = CplUserData<TData, dim>;

	using TRaster = Raster<number, 2>;
	//using TRaster1 = Raster<number, 1>;
	//using TRaster2 = Raster<number, 2>;
	//using TRaster3 = Raster<number, 3>;
	//using base_type = StdGlobPosData<RasterUserData<dim>, number, dim, void>;

	using input_type = SmartPtr<TRaster>;


private:

	SmartPtr<TRaster>  m_spRaster;
	int m_interpOrder;
	number m_rescale;



public:

	RasterUserData(SmartPtr<TRaster> myRaster)
	: m_spRaster(myRaster), m_interpOrder(1), m_rescale(1.0)
	{}

	virtual ~RasterUserData() = default;


	inline void evaluate(number& y, const MathVector<dim>& x, number t, int si) const
	{
		TRaster::Coordinate xcoord;
		xcoord[0] = x[0];
		if (dim>1) { xcoord[1] = x[1];}
		y = m_spRaster->interpolate(xcoord, m_interpOrder)*m_rescale;
	}

	void set_order(int order) {m_interpOrder = order;}
	void set_scale(number alpha) {m_rescale = alpha;}

	bool requires_grid_fct() const override {return false;}
	bool constant() const override {return false;}
	bool continuous() const override {return false;}

};





}//	end of namespace ug



#endif