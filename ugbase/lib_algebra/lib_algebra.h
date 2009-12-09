/*
 * lib_algebra.h
 *
 *  Created on: 02.07.2009
 *      Author: andreasvogel
 */

#ifndef __H__LIB_ALGEBRA__
#define __H__LIB_ALGEBRA__

/*
#include "hypre_algebra/hyprematrix.h"
#include "hypre_algebra/hyprevector.h"
#include "hypre_algebra/hyprelinearsolver.h"
*/

#include "arne_algebra/arnematrix.h"
#include "arne_algebra/arnevector.h"
#include "arne_algebra/arnelinearsolver.h"

namespace ug {

/*
typedef HYPREboomerAMG LinearSolver;
typedef HypreMatrix Matrix;
typedef HypreVector Vector;
*/


typedef ArneJacobi LinearSolver;
typedef ArneMatrix Matrix;
typedef ArneVector Vector;


}

#endif /* __H__LIB_ALGEBRA__ */
