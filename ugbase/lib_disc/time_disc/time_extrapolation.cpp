#include "time_extrapolation.h"
#include "lib_algebra/cpu_algebra_types.h"

namespace ug {
#ifdef UG_CPU_1
template class AitkenNevilleTimex< CPUAlgebra::vector_type >;
#endif
#ifdef UG_CPU_2
template class AitkenNevilleTimex< CPUBlockAlgebra<2>::vector_type  >;
#endif
#ifdef UG_CPU_3
template class AitkenNevilleTimex< CPUBlockAlgebra<3>::vector_type  >;
#endif
#ifdef UG_CPU_4
template class AitkenNevilleTimex< CPUBlockAlgebra<4>::vector_type >;
#endif
}
