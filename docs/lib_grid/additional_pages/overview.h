//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m09 d14

/**	\page pageLGOverview libGrid - Overview

	- \ref secGrid "Grid"
	- \ref secGeomObjs "Geometric Objects"
	- \ref secObjectCreation "Object Creation"
	- \ref secObjectAccess "Object Access"
	- \ref secAttachments "Attachments"
	- \ref secTools "Tools"
	- \ref secUtil "Util"
	- \ref secAlgorithms "Algorithms"
	- \ref secFileIO "FileIO"

<hr>
\section secGrid Grid
The central class of libGrid is the Grid class.
It handles element creation, neighbourhood managment, data-attachments and observers.
The class MultiGrid is derived from Grid, and adds a hierarchical structure to the grids elements.
	- ug::Grid
	- ug::MultiGrid

<hr>
\section secGeomObjs Geometric Objects
Geometric Objects are the building blocks of a grid.
There are four different basic geometric objects:
	- ug::VertexBase
	- ug::EdgeBase
	- ug::Face
	- ug::Volume
	
Those basic objects are then further specialized:
	- ug::Vertex
	- ug::Edge
	- ug::Triangle, ug::Quadrilateral
	- ug::Tetrahedron, ug::Hexahedron, ug::Prism, ug::Pyramid
	
For hanging-node support (only in 2d in the moment) the following objects are introduced:
	- ug::HangingVertex
	- ug::ConstrainingEdge. ug::ConstrainingVertex


<hr>
\section secObjectCreation Object Creation
Geometric Objects are created through the ug::Grid::create method (or ug::MultiGrid::create).
Since the Grid class only knows about the basic geometric objects,
ug::Grid::create<TGeomObj> is a template method, that takes the type
of the geometric object that is to be created as template argument:

<code>
	using namespace ug;
	...
	Grid g;
//	create vertices
	Vertex* v1 = *g.create<Vertex>();
	Vertex* v2 = *g.create<Vertex>();
	
//	create an edge
	g.create<Edge>(EdgeDescriptor(v1, v2));
</code>

Geometric Objects in a ug::MultiGrid are created in the same way.
The create method takes an optional parameter:
<code>ug::Grid::create(GeometricObject* pParent = NULL);</code>
The parent is used in different ways:
	- ug::SubsetHandler can automatically assign a subset based on the parents subset.
	- ug::Selector can automatically assign the selection-status based on the parents selection status.
	- A ug::MultiGrid inserts elements one level above the parents level.
	- You can use the parent in your own grid-observer specializations (derive from ug::IGridObserver)
	
All those beahviours can be enabled / disabled in the respective classes.


<hr>
\section secObjectAccess Object Access
libGrid uses the technique of iterators for geometric-object access.
A separate iterator-type exists for each object-type:
	- ug::VertexBaseIterator, ug::VertexIterator, ...
	- ug::EdgeBaseIterator, ug::EdgeIterator, ...
	- ug::FaceIterator, ug::TriangleIterator, ug::QuadrilateralIterator, ...
	- ug::TetrahedronIterator, ug::HexahedronIterator, ug::PrismIterator, ug::PyramidIterator
	
You can query a grid for a begin and an end-iterator for each geometric-object type using
ug::Grid::begin and ug::Grid::end. Both method are template methods. The template argument
specifies the type of geometric-object over which you want to iterate.

<code>
//	Let g be a grid that already contains some geometric objects
//	iterate over all vertices
	for(VertexBaseIterator iter = g.begin<VertexBase>();
		iter != g.end<VertexBase(); ++iter)
	{
		VertexBase* v = *iter;
		...
	}

//	iterate over all faces
	for(FaceIterator iter = g.begin<Face>();
		iter != g.end<Face(); ++iter)
	{
		Face* f = *iter;
		...
	}
	
//	iterate over all triangles
	for(TriangleIterator iter = g.begin<Triangle>();
		iter != g.end<Triangle(); ++iter)
	{
		Triangle* t = *iter;
		...
	}
</code>

The same technique can be used to iterate over all Triangles of a subset:
<code>
	using namespace ug;
	...
	Grid g;
	SubsetHandler sh(g);
	...
//	iterate over all triangles in subset 0
	for(TriangleIterator iter = sh.begin<Triangle>(0);
		iter != sh.end<Triangle>(0); ++iter)
	{
		Triangle* t = *iter;
		...
	}
</code>

or a ug::Selector
<code>
	...
	Selector sel(g);
	...
	for(TriangleIterator iter = sel.begin<Triangle>();
		iter != sel.end<Triangle>(); ++iter)
	...
</code>

or a level of a ug::MultiGrid
<code>
	...
	MultiGrid mg;
	...
	for(TriangleIterator iter = mg.begin<Triangle>(0);
		iter != mg.end<Triangle>(0); ++iter)
</code>

or triangles on level l in subset i of a ug::MGSubsetHandler
<code>
	...
	MGSubsetHandler mgsh(mg);
	...
	for(TriangleIterator iter = mgsh.begin<Triangle>(i, l);
		iter != mgsh.end<Triangle>(i, l); ++iter)
	...
</code>

There are even more classes that support this way of geometric-object iteration.
\sa ug::GeometricObjectCollection, \sa ug::MultiLevelGeometricObjectCollection


<hr>
\section secAttachments Attachments
...


<hr>
\section secTools Tools
There are some classes that help tremendously when implementing algorithms:
the Selector and the SubsetHandler.
	- ug::Selector
	- ug::SubsetHandler

<hr>
\section secUtil Util
...

<hr>
\section secAlgorithms Algorithms
...

<hr>
\section secFileIO FileIO
...

*/
