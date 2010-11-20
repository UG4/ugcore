// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// y10 m01 d25

#ifndef __H__UG3LAYER__GM__
#define __H__UG3LAYER__GM__

#include "lib_grid/lg_base.h"
#include "gm_constants.h"

////////////////////////////////////////////////////////////////////////
//	types
typedef uint UINT;
typedef int INT;

////////////////////////////////////////////////////////////////////////
//	dimension
#ifndef __THREEDIM__
	#ifndef __TWODIM__
		#define __TWODIM__
	#endif
#endif


////////////////////////////////////////////////////////////////////////
// elements
#ifdef __TWODIM__
	typedef ug::Face	element;
#else if __THREEDIM__
	typedef ug::Volume	element;
#endif

typedef ug::Triangle triangle;
typedef ug::Quadrilateral quadrilateral;

typedef ug::Tetrahedron tetrahedron;
typedef ug::Pyramid pyramid;
typedef ug::Prism prism;
typedef ug::Hexahedron hexahedron;

////////////////////////////////////////////////////////////////////////
// geometrical part
/*
typedef struct format                           FORMAT;

typedef union  vertex                           VERTEX;
typedef struct elementlist                      ELEMENTLIST;
typedef struct node                             NODE;
typedef union  element                          ELEMENT;
typedef struct link                             LINK;
typedef struct edge                             EDGE;
typedef union  geom_object                      GEOM_OBJECT;
typedef union  selection_object         SELECTION_OBJECT;
typedef struct grid                             GRID;
typedef struct multigrid                        MULTIGRID;
typedef union object_with_key           KEY_OBJECT;
*/

////////////////////////////////////////////////////////////////////////
//	grid
//	the grid in the ug3layer is not the equivalent of an ug::Grid.
//	instead it represents a level of the multi-grid hierarchy.
struct UG3Grid
{
///	not used by ug4
	UINT control;
///	not used by ug4
	INT attribut;
///	not used by ug4
	INT status;

///	level within the multi-grid structure
    INT level;

///	stores iterators to the elements
	GeometricObjectCollection goc;
};

typedef UG3Grid GRID;

////////////////////////////////////////////////////////////////////////
//	multi-grid
struct UG3MultiGrid
{
///	One grid for each level
	UG3Grid[MAXLEVEL]	grids;

///	depth of the hierarchy
	INT topLevel;

///	the current level
	INT currentLevel;

///	The real multi-grid
	ug::MultiGrid	mg;
};

typedef UG3MultiGrid MULTIGRID;

void SetCurrentGrid(MULTIGRID* pMG);
MULTIGRID* GetCurrentGrid();



////////////////////////////////////////////////////////////////////////
//	macros for grids (ug/gm/gm.h line 3098)
#define GLEVEL(p)		((p)->level)
#define FIRSTELEMENT(p)         (*(p->goc->begin<ELEMENT>()))
#define PFIRSTELEMENT(p)        FIRSTELEMENT(p)
//#define LASTELEMENT(p)          ((p)->lastelement[0])
//#define PLASTELEMENT(p)         LASTELEMENT(p)
#endif
