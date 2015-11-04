
#ifndef __H__UG__LIB_DISC__DOF_MANAGER__LEXORDER__
#define __H__UG__LIB_DISC__DOF_MANAGER__LEXORDER__

#include <vector>
#include <utility> // for pair

#include "lib_disc/function_spaces/approximation_space.h"

namespace ug{

template<int dim>
void ComputeLexicographicOrder(std::vector<size_t>& vNewIndex,
                               std::vector<std::pair<MathVector<dim>, size_t> >& vPos);

/// orders the dof distribution using Cuthill-McKee
template <typename TDomain>
void OrderLexForDofDist(SmartPtr<DoFDistribution> dd, ConstSmartPtr<TDomain> domain);

/// orders the all DofDistributions of the ApproximationSpace using lexicographic order
template <typename TDomain>
void OrderLex(ApproximationSpace<TDomain>& approxSpace, const char *order);

} // end namespace ug

#endif /* __H__UG__LIB_DISC__DOF_MANAGER__LEXORDER__ */
