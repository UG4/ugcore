/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
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

#include <fstream>
#include <queue>
#include "lib_grid/lg_base.h"
#include "common/profiler/profiler.h"
#include "simple_grid.h"
#include "edge_length_adjustment.h"
#include "lib_grid/refinement/regular_refinement.h"
#include "common/node_tree/node_tree.h"
#include "lib_grid/algorithms/trees/octree.h"
#include "lib_grid/callbacks/callbacks.h"

using namespace std;

namespace ug
{

///	only for debugging purposes!!!
/**	Output value pairs to gnuplot...
 * \{ */
//#define SMOOTH_REMESHING__GPLOT_ENABLED
#ifdef SMOOTH_REMESHING__GPLOT_ENABLED
	using GnuplotData = vector<pair<number, number> >;
	static GnuplotData gplotLengthFac;
	static GnuplotData gplotMinCurvature;
	static GnuplotData gplotAverageCurvature;

	void WriteGnuplotData(const char* filename, const GnuplotData& data)
	{
		ofstream out(filename);
		if(!out)
			return;

		for(size_t i = 0; i < data.size(); ++i)
			out << data[i].first << " " << data[i].second << endl;

		out.close();
	}

	#define GPLOTPOINT(dataName, x, y) dataName.push_back(make_pair<number, number>((x), (y)));
	#define GPLOTSAVE()	{WriteGnuplotData("length_fac.gplot", gplotLengthFac);\
						WriteGnuplotData("min_curvature.gplot", gplotMinCurvature);\
						WriteGnuplotData("average_curvature.gplot", gplotAverageCurvature);}
#else
//	do nothing if SMOOTH_REMESHING__GPLOT_ENABLED is false
	#define GPLOTPOINT(dataName, x, y)
	#define GPLOTSAVE()
#endif
/** \} */


/*
vector3 PNTrianglePos(const vector3& p0, const vector3& p1, const vector3& p2,
					  const vector3& n0, const vector3& n1, const vector3& n2);

vector3 PNTriangleNorm(const vector3& p0, const vector3& p1, const vector3& p2,
					   const vector3& n0, const vector3& n1, const vector3& n2);

vector3 PNCTrianglePos(const vector3& p0, const vector3& p1, const vector3& p2,
						const vector3& n0, const vector3& n1, const vector3& n2,
						const vector3& cn0, const vector3& cn1, const vector3& cn2);

vector3 PNCTriangleNorm(const vector3& p0, const vector3& p1, const vector3& p2,
						const vector3& n0, const vector3& n1, const vector3& n2,
						const vector3& cn0, const vector3& cn1, const vector3& cn2);
*/


class ILocalRemesher{
	public:
		virtual ~ILocalRemesher() = default;
		virtual void smooth_vertex(Vertex* vrt) = 0;
		virtual Vertex* collapse_edge(Edge* edge) = 0;
		virtual Vertex* split_edge(Edge* edge) = 0;
};


class IPatchRemesher : public ILocalRemesher{
	public:
		virtual ~IPatchRemesher();

	///	set the grid which will be remeshed
		void set_grid(Grid& grid, APosition aPos);
		void set_crease_callbacks(Grid::vertex_traits::callback vrtCreaseCallback,
								  Grid::edge_traits::callback edgeCreaseCallback);
		void set_fixed_callbacks(Grid::vertex_traits::callback vrtFixedCallback,
								 Grid::edge_traits::callback edgeFixedCallback);

	///	Adds a new patch consisting of the given elements and associated vertices
	/**	Make sure to specify the source grid before calling this method.*/
		template <typename TElemIterator>
		void add_surface_patch(TElemIterator begin, TElemIterator end);

		virtual void smooth_vertex(Vertex* vrt) = 0;
		virtual Vertex* collapse_edge(Edge* edge) = 0;
		virtual Vertex* split_edge(Edge* edge) = 0;
		virtual Edge* swap_edge(Edge* edge) = 0;

		virtual vector3 vertex_position(Vertex* vrt) = 0;
		virtual vector3 vertex_normal(Vertex* vrt) = 0;

		virtual number approximation_quality(Edge* e) = 0;
		virtual number approximation_quality(Face* f) = 0;
		virtual number element_quality(Face* f) = 0;

	protected:
		virtual void relocate_vertex(Vertex* vrt);

	private:
		class ProjectedPoint{
			int patchID;
			GridObject* elem;
			vector2 barycentricCoords;
		};

		Grid	m_refGrid;
		Grid*	m_remeshGrid;
};



}//	end of namespace
