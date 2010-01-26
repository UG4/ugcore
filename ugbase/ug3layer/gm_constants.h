// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// y10 m01 d25

////////////////////////////////////////////////////////////////////////
//	constants copied from ug3 (ug/gm/gm.h)

#ifndef __H__UG3LAYER__GM_CONSTANTS__
#define __H__UG3LAYER__GM_CONSTANTS__

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	lines 158 - 388
/** @name Some size parameters */
/*@{*/
/** \brief  maximal space dimension                              */
#define DIM_MAX                                 3  
/** \brief  maximal dimension of boundary surface*/
#define DIM_OF_BND_MAX                  2
/** \brief  maximum depth of triangulation               */
#define MAXLEVEL                                32
/** \brief  use 5 bits for object identification */
#define MAXOBJECTS                              32
/** \brief  max number of elements in selection  */
#define MAXSELECTION               100
/*@}*/

/** @name Some size macros for allocation purposes */
/*@{*/
/** \brief max number of sides of an elem       */
#define MAX_SIDES_OF_ELEM               6
/** \brief max number of edges of an element*/
#define MAX_EDGES_OF_ELEM               12
/** \brief max number of corners of an eleme*/
#define MAX_CORNERS_OF_ELEM             8
/** \brief max number of edges of a side    */
#define MAX_EDGES_OF_SIDE               4
/** \brief max number of edges meeting in co*/
#define MAX_EDGES_OF_CORNER             4
/** \brief max number of corners of a side  */
#define MAX_CORNERS_OF_SIDE     4
/** \brief an edge has always two corners.. */
#define MAX_CORNERS_OF_EDGE             2
/** \brief two sides have one edge in common*/
#define MAX_SIDES_OF_EDGE               2
/** \brief max number of sons of an element */  
enum {MAX_SONS = 30};
/** \brief max number of nodes on elem side */
#define MAX_SIDE_NODES                  9
/** \brief max number of son edges of edge  */  
#define MAX_SON_EDGES                   2
/** \brief max #fine sides touching a coarse*/
#define MAX_SIDES_TOUCHING              10

/** \todo Please doc me! */
#define MAX_ELEM_VECTORS                (MAX_CORNERS_OF_ELEM+MAX_EDGES_OF_ELEM+1+MAX_SIDES_OF_ELEM)
/** \brief max number of doubles in a vector or matrix mod 32 */ 
#define MAX_NDOF_MOD_32        256
/** \brief max number of doubles in a vector or matrix */
#define MAX_NDOF 32*MAX_NDOF_MOD_32  
/*@}*/

  
/****************************************************************************/
/*                                                                          */
/* defines for algebra                                                      */
/*                                                                          */
/****************************************************************************/

/** \brief Number of different data types                                    */
#define MAXVOBJECTS                                             4
/** \brief max number of abstract vector types                  */
#define MAXVECTORS                                              4
#if (MAXVECTORS<MAXVOBJECTS)
        #error *** MAXVECTORS must not be smaller than MAXVOBJECTS ***
#endif

/** \brief to indicate type not defined                                 */
#define NOVTYPE                                                 -1
/** \brief max number of geometric domain parts                 */
#define MAXDOMPARTS                                             4

/** \brief transforms type into bitpattern                              */
#define BITWISE_TYPE(t) (1<<(t))

/* derived sizes for algebra */
/** \brief max number of diff. matrix types                 */
#define MAXMATRICES             MAXVECTORS*MAXVECTORS
/** \brief max number of diff. connections              */
#define MAXCONNECTIONS  (MAXMATRICES + MAXVECTORS)

/** \todo Please doc me! */
#define MATRIXTYPE(rt,ct)   ((rt)*MAXVECTORS+(ct))
/** \todo Please doc me! */
#define DIAGMATRIXTYPE(rt)  (MAXMATRICES+rt)

/** \brief Type of geometric entity which a certain vector is attached to */
enum VectorType {NODEVEC,   /**< Vector associated to a node */
                 EDGEVEC,   /**< Vector associated to an edge */
                 ELEMVEC,   /**< Vector associated to an element */
                 SIDEVEC    /**< Vector associated to an element side */
};

/** @name Some constants for abstract vector type names */
/*@{*/
/** \todo Please doc me! */
#define FROM_VTNAME                                             '0'
/** \todo Please doc me! */
#define TO_VTNAME                                               'z'
/** \todo Please doc me! */
#define MAXVTNAMES                                              (1+TO_VTNAME-FROM_VTNAME)
/*@}*/

/** @name Constants for blockvector description (BVD) */
/*@{*/
/** \brief number for "there is no blockvector"; largest number of type BLOCKNUMBER */
#define NO_BLOCKVECTOR ((BLOCKNUMBER)~0)
/** \brief largest admissible blockvector number */
#define MAX_BV_NUMBER (NO_BLOCKVECTOR - 1)
/** \brief largest admissible blockvector level number */
#define MAX_BV_LEVEL UCHAR_MAX
/** \brief Maximum number
   of entries in a BVD; NOTE: the actual available
   number of entries depends on the range of each entry */
#define BVD_MAX_ENTRIES (sizeof(BVD_ENTRY_TYPE)*CHAR_BIT)
/*@}*/ 

/** @name Constants for BLOCKVECTOR */
/*@{*/
/** \brief symbolic value for BVDOWNTYPE */
enum {BVDOWNTYPEVECTOR,
      BVDOWNTYPEBV,
      BVDOWNTYPEDIAG};

/** \brief symbolic value for BVTVTYPE */
enum {BV1DTV,
      BV2DTV};

enum {BVNOORIENTATION, /**< No special orientation for BVORIENTATION */
      BVHORIZONTAL, /**< Vectors form a horizontal line for BVORIENTATION */
      BVVERTICAL /**< Vectors form a vertical line for BVORIENTATION */
};
/*@}*/    
  
/****************************************************************************/
/*                                                                          */
/* various defines                                                          */
/*                                                                          */
/****************************************************************************/

/**      0 = OK as usual */
/** @name result codes of user supplied functions*/ 
/*@{*/
/** \brief coordinate out of range                              */
#define OUT_OF_RANGE                    1
/** \brief configProblem could not init problem */
#define CANNOT_INIT_PROBLEM     1
    
/** \brief Use of GSTATUS (for grids), use power of 2 */
enum {GSTATUS_BDF         = 1,
      GSTATUS_INTERPOLATE = 2,
      GSTATUS_ASSEMBLED   = 4,
      GSTATUS_ORDERED     = 8};
/*@}*/

/** \brief Selection mode */
enum {nodeSelection=1,     /**< Objects selected are nodes */
      elementSelection=2,   /**< Objects selected are elements */
      vectorSelection=3    /**< Objects selected are vectors */
};

/** \brief Possible values for rule in MarkForRefinement */
enum RefinementRule 
    {NO_REFINEMENT = 0,
     COPY = 1,
     RED =  2,
     BLUE = 3,
     COARSE = 4,
#ifdef __TWODIM__
     BISECTION_1 = 5,
     BISECTION_2_Q = 6,
     BISECTION_2_T1 = 7,
     BISECTION_2_T2 = 8,
     BISECTION_3 = 9
#endif
#ifdef __THREEDIM__
     
     TETRA_RED_HEX = 5,
     
     PRISM_BISECT_1_2 = 9,
     PRISM_QUADSECT = 7,
     PRISM_BISECT_HEX0 = 5,
     PRISM_BISECT_HEX1 = 8,
     PRISM_BISECT_HEX2 = 6,
     PRISM_ROTATE_LEFT = 10,
     PRISM_ROTATE_RGHT = 11,
     PRISM_QUADSECT_HEXPRI0 = 14,
     PRISM_RED_HEX = 15,
     
     HEX_BISECT_0_1 = 5,
     HEX_BISECT_0_2 = 6,
     HEX_BISECT_0_3 = 7,
     HEX_TRISECT_0 = 8,
     HEX_TRISECT_5 = 9,
     HEX_QUADSECT_0 = 12,      
     HEX_QUADSECT_1 = 13,
     HEX_QUADSECT_2 = 14,      
     HEX_BISECT_HEXPRI0 = 15,
     HEX_BISECT_HEXPRI1 = 16

#endif
};

/** \brief Values for element class */
enum MarkClass {NO_CLASS,
                    YELLOW_CLASS,
                    GREEN_CLASS,
                    RED_CLASS,
                    SWITCH_CLASS
                    };

    /** \brief Values for node types (relative to the father element of the vertex) */
    enum {CORNER_NODE,
          MID_NODE,
          SIDE_NODE,
          CENTER_NODE,
          LEVEL_0_NODE
    };

/** @name Macros for the multigrid user data space management */
/*@{*/
#define OFFSET_IN_MGUD(id)              (GetMGUDBlockDescriptor(id)->offset)
#define IS_MGUDBLOCK_DEF(id)    (GetMGUDBlockDescriptor(id)!=NULL)
/*@}*/

/* REMARK: TOPNODE no more available since 970411
   because of problems in parallelisation
   to use it in serial version uncomment define
#define TOPNODE(p)              ((p)->iv.topnode)
*/

/** \brief Modes for LexOrderVectorsInGrid */
enum {OV_CARTES,
      OV_POLAR};

#endif
