/*
 * simpleMatrixOps.h
 *
 *  Created on: 30.12.2024
 *      Author: Markus Knodel
 */

#ifndef UGCORE_UGBASE_LIB_GRID_ALGORITHMS_EXTRUSION_SIMPLEMATRIXOPS_H_
#define UGCORE_UGBASE_LIB_GRID_ALGORITHMS_EXTRUSION_SIMPLEMATRIXOPS_H_

// #include "lib_grid/lg_base.h"
#include <vector>

namespace ug {
namespace simpleMatrOps {

std::vector<double> cramerRule(std::vector<std::vector<double>> const & coefficients, std::vector<double> const & constants);

}
}



#endif