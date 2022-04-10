/*
 * Copyright (c) 2022:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__LIB_GRID__CUSTOM_PROJECTOR__
#define __H__UG__LIB_GRID__CUSTOM_PROJECTOR__

#include "common/math/misc/math_util.h"
#include "refinement_projector.h"

#include "lib_disc/spatial_disc/user_data/user_data.h"
#include "bindings/lua/lua_user_data.h"

namespace ug{

class CustomProjector : public RefinementProjector{
public:
	CustomProjector(){}

	virtual ~CustomProjector(){}

	typedef MathVector<3> small_vec_t_3;
	typedef MathVector<6> small_vec_t_6;
	typedef MathVector<12> small_vec_t_12;
	typedef UserData<small_vec_t_6, 6> Position_t_6;
	typedef UserData<small_vec_t_12, 12> Position_t_12;
	typedef SmartPtr<Position_t_6> TSpUserData_6;
	typedef SmartPtr<Position_t_12> TSpUserData_12;
	typedef LuaUserData<small_vec_t_6, 6> TLuaUserData_6;
	typedef LuaUserData<small_vec_t_12, 12> TLuaUserData_12;



///	called when a new vertex was created from an old edge.
	virtual number new_vertex(Vertex* vrt, Edge* parent)
	{
		typename Edge::ConstVertexArray vrts = parent->vertices();

		small_vec_t_6 newpos;
		small_vec_t_6 oldpos;
		oldpos[0] = pos(vrts[0])[0];
		oldpos[1] = pos(vrts[0])[1];
		oldpos[2] = pos(vrts[0])[2];
		oldpos[3] = pos(vrts[1])[0];
		oldpos[4] = pos(vrts[1])[1];
		oldpos[5] = pos(vrts[1])[2];
		(*m_spEdgeProjectorFunc)(newpos, oldpos, 0.0f, 0);
		small_vec_t_3 newpos_short;
		newpos_short[0] = newpos[0];
		newpos_short[1] = newpos[1];
		newpos_short[2] = newpos[2];

		set_pos(vrt, newpos_short);

		return 1;
	}

///	called when a new vertex was created from an old face.
	virtual number new_vertex(Vertex* vrt, Face* parent)
	{
		size_t numVrts = parent->num_vertices();
		if(numVrts != 4){
			UG_THROW("custom_projector::new_vertex(Face): Only implemented for quads!");
		}


		typename Edge::ConstVertexArray vrts = parent->vertices();

		small_vec_t_12 newpos;
		small_vec_t_12 oldpos;
		oldpos[0] = pos(vrts[0])[0];
		oldpos[1] = pos(vrts[0])[1];
		oldpos[2] = pos(vrts[0])[2];
		oldpos[3] = pos(vrts[1])[0];
		oldpos[4] = pos(vrts[1])[1];
		oldpos[5] = pos(vrts[1])[2];
		oldpos[6] = pos(vrts[2])[0];
		oldpos[7] = pos(vrts[2])[1];
		oldpos[8] = pos(vrts[2])[2];
		oldpos[9] = pos(vrts[3])[0];
		oldpos[10]= pos(vrts[3])[1];
		oldpos[11]= pos(vrts[3])[2];
		(*m_spQuadProjectorFunc)(newpos, oldpos, 0.0f, 0);

		small_vec_t_3 newpos_short;
		newpos_short[0] = newpos[0];
		newpos_short[1] = newpos[1];
		newpos_short[2] = newpos[2];

		set_pos(vrt, newpos_short);

		return 1;
	}

///	called when a new vertex was created from an old volume.
	virtual number new_vertex(Vertex* vrt, Volume* parent)
	{
		UG_THROW("custom_projector::new_vertex(Volume): Projector works for surfaces only!");
	}

	void set_edge_projector(const char* strEdgeProjectorFunc){
		m_spEdgeProjectorFunc = make_sp(new TLuaUserData_6(strEdgeProjectorFunc));
	}

	void set_quad_projector(const char* strQuadProjectorFunc){
		m_spQuadProjectorFunc = make_sp(new TLuaUserData_12(strQuadProjectorFunc));
	}

private:
	TSpUserData_6 m_spEdgeProjectorFunc;
	TSpUserData_12 m_spQuadProjectorFunc;
};

}

#endif
