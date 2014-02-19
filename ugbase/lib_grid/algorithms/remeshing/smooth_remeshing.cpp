//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m01 d22

#include <fstream>
#include <queue>
#include "lib_grid/lg_base.h"
#include "common/profiler/profiler.h"
#include "simple_grid.h"
#include "edge_length_adjustment.h"
#include "lib_grid/algorithms/refinement/regular_refinement.h"
#include "common/node_tree/node_tree.h"
#include "lib_grid/algorithms/trees/octree.h"
#include "lib_grid/algorithms/callback_util.h"

using namespace std;

namespace ug
{

///	only for debugging purposes!!!
/**	Output value pairs to gnuplot...
 * \{ */
//#define SMOOTH_REMESHING__GPLOT_ENABLED
#ifdef SMOOTH_REMESHING__GPLOT_ENABLED
	typedef vector<pair<number, number> > GnuplotData;
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
		virtual ~ILocalRemesher()	{}
		virtual void smooth_vertex(VertexBase* vrt) = 0;
		virtual VertexBase* collapse_edge(EdgeBase* edge) = 0;
		virtual VertexBase* split_edge(EdgeBase* edge) = 0;
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
		template <class TElemIterator>
		void add_surface_patch(TElemIterator begin, TElemIterator end);

		virtual void smooth_vertex(VertexBase* vrt) = 0;
		virtual VertexBase* collapse_edge(EdgeBase* edge) = 0;
		virtual VertexBase* split_edge(EdgeBase* edge) = 0;
		virtual EdgeBase* swap_edge(EdgeBase* edge) = 0;

		virtual vector3 vertex_position(VertexBase* vrt) = 0;
		virtual vector3 vertex_normal(VertexBase* vrt) = 0;

		virtual number approximation_quality(EdgeBase* e) = 0;
		virtual number approximation_quality(Face* f) = 0;
		virtual number element_quality(Face* f) = 0;

	protected:
		virtual void relocate_vertex(VertexBase* vrt);

	private:
		class ProjectedPoint{
			int 				patchID;
			GridObject* 	elem;
			vector2				barycentricCoords;
		};

		Grid	m_refGrid;
		Grid*	m_remeshGrid;
};



}//	end of namespace
