/*
 * histogramm.h
 *
 *  Created on: 20.09.2013
 *      Author: mrupp
 */

#ifndef HISTOGRAMM_H_
#define HISTOGRAMM_H_

#include <vector>
#include <string>

namespace ug{


std::string HistogrammString(std::vector<double> values);


std::string DistributionPercentage(std::vector<double> values);
}
#endif /* HISTOGRAMM_H_ */
