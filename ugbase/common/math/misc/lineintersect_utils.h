/* Copyright (C) Graham Rhodes, 2001. 
 * All rights reserved worldwide.
 *
 * This software is provided "as is" without express or implied
 * warranties. You may freely copy and compile this source into
 * applications you distribute provided that the copyright text
 * below is included in the resulting source code, for example:
 * "Portions Copyright (C) Graham Rhodes, 2001"
 */
/**************************************************************************************
|
|           File: lineintersect_utils.h
|
|        Purpose: Function prototypes for line segment intersection utility functions
|
|     Book Title: Game Programming Gems II
|
|  Chapter Title: Fast, Robust Intersection of 3D Line Segments
|
|         Author: Graham Rhodes
|
|      Revisions: 05-Apr-2001 - GSR. Original.
|
**************************************************************************************/
#ifndef _lineintersect_utils_h
#define _lineintersect_utils_h

#include "common/types.h"

/// \addtogroup ugbase_math
/// \{

void IntersectLineSegments(const number A1x, const number A1y, const number A1z,
                           const number A2x, const number A2y, const number A2z,
                           const number B1x, const number B1y, const number B1z,
                           const number B2x, const number B2y, const number B2z,
                           bool infinite_lines, number epsilon, number &PointOnSegAx,
                           number &PointOnSegAy, number &PointOnSegAz, number &PointOnSegBx,
                           number &PointOnSegBy, number &PointOnSegBz, number &NearestPointX,
                           number &NearestPointY, number &NearestPointZ, number &NearestVectorX,
                           number &NearestVectorY, number &NearestVectorZ, bool &true_intersection);

void FindNearestPointOnLineSegment(const number A1x, const number A1y, const number A1z,
                                   const number Lx, const number Ly, const number Lz,
                                   const number Bx, const number By, const number Bz,
                                   bool infinite_line, number epsilon_squared, number &NearestPointX,
                                   number &NearestPointY, number &NearestPointZ,
                                   number &parameter);

void FindNearestPointOfParallelLineSegments(number A1x, number A1y, number A1z,
                                            number A2x, number A2y, number A2z,
                                            number Lax, number Lay, number Laz,
                                            number B1x, number B1y, number B1z,
                                            number B2x, number B2y, number B2z,
                                            number Lbx, number Lby, number Lbz,
                                            bool infinite_lines, number epsilon_squared,
                                            number &PointOnSegAx, number &PointOnSegAy, number &PointOnSegAz,
                                            number &PointOnSegBx, number &PointOnSegBy, number &PointOnSegBz);

void AdjustNearestPoints(number A1x, number A1y, number A1z,
                         number Lax, number Lay, number Laz,
                         number B1x, number B1y, number B1z,
                         number Lbx, number Lby, number Lbz,
                         number epsilon_squared, number s, number t,
                         number &PointOnSegAx, number &PointOnSegAy, number &PointOnSegAz,
                         number &PointOnSegBx, number &PointOnSegBy, number &PointOnSegBz);

// end group ugbase_math
/// \}

#endif
