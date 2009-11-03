//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m09 d14

/**	\page pageLGIntroduction libGrid - Introduction
	
	- \ref secAbout "About"
	- \ref secWhatIsIt "What is libGrid?"
	- \ref secFeatures "Main features?"
	- \ref secCompilation "How to compile libGrid?"
	- \ref secUsage "Usage"
	
<hr>
\section secAbout About
libGrid has been written and is still developed by Sebastian Reiter (Sebastian.Reiter@gcsc.uni-frankfurt.de) at the
<i>Goethe Center for Scientific Computing at the University of Frankfurt, Germany </i>(http://www.g-csc.de).
 
Amongst others, the following people contributed to libGrid:<br>
Niklas Antes, Martin Stepniewski, Nicolas Tessore, Andreas Vogel

<hr>
\section secWhatIsIt What is libGrid?
libGrid is a library for the manipulation of 1, 2 and 3 dimensional grids
written in C++.

It is designed in a flexible way that should allow its usage for all
grid-related algorithms.

<hr>
\section secFeatures Main features
libGrid tries to treat grids in a flexible, but at the same time fast and efficient way.
This goal is achieved by separation of topology and data.<br>
While the elements of a grid are stored in linked lists,
the associated data is stored in separate arrays.
This not only improves cache-hits considerably but also allows to add (and remove)
data to the elements of the grid at runtime.
To avoid performance problems, data is not added to elements one by one, but to all
elements of the same type (Vertices / Edges / Faces / Volumes) at once.<br>
Data managment is fully automated. The access to data is simple, typesafe and fast.

A grid in libGrid can consist of vertices, edges, faces and volumes.
Besides others, the following important concrete types are supported:
	- Vertex
	- Edge
	- Triangle, Quadrilateral
	- Tetrahedron, Hexahedron, Prism, Pyramid

New types can be added without much hassle.
Those types can be combined as desired. A grid could for example consist of some
Vertices that are connected by Triangles plus a single edge.
However, since it is sometimes desireable that missing edges (or faces) are created
automatically if a new face or volume is added, grids feature options
through which such behaviours can be turned on and off.
\sa GeometricObjects

Automated access to the neighbourhood of elements is another feature of libGrid.
Neighbours of an element can be collected on the fly, or are optionally stored
for each element. This allows the user to choose between better runtime or
lower memory consumption.
This again can be changed at runtime, which allows algorithms to enable or disable
features as required.

Since libGrid heavily uses template-programming, efficient algorithms can be implemented
with a minimal amount of code. However, while thorough knowledge of template-programming
is beneficial for writing intelligent and efficient code,
only a very basic understanding is required to start programming with libGrid.

LibGrid also features a mechanism to register observers at a grid. Those observers
are then notified if changes to the grids topology are made.
Combined with the possibility to dynamically add data to a grid at runtime, those
observers are a powerful tool to create stable, robust and modular code. 
\sa ug::Selector, \sa ug::SubsetHandler

<hr>
\section secCompilation How to compile libGrid
libGrid can be build using CMAKE (http://www.cmake.org/).
	- download the current source-distribution of cmake
	- build and install it (configure; make; install;)
	- change directory to libGrids root directory (../libGrid)
	- create a new build directory (mkdir build)
	- change directory to the build directory (cd build)
	- call cmake with libGrids path (cmake ../)
	- call make
 
You should now find all generated libs in libGrid/lib and
the executable in libGrid/bin.

<hr>
\section secUsage Usage
...
*/