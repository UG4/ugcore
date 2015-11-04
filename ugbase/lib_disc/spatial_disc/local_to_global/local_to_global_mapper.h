
#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__LOCAL_TO_GLOBAL_MAPPER__
#define __H__UG__LIB_DISC__SPATIAL_DISC__LOCAL_TO_GLOBAL_MAPPER__

// extern headers
#include <vector>

// intern headers
#include "lib_disc/dof_manager/dof_distribution.h"

namespace ug{

/// interface for definition of special LocalToGlobal mappings
/**
 * \tparam	TAlgebra			type of Algebra
 */
template <typename TAlgebra>
class ILocalToGlobalMapper
{
	public:
	///	Algebra type
		typedef TAlgebra algebra_type;

	///	Type of algebra matrix
		typedef typename algebra_type::matrix_type matrix_type;

	///	Type of algebra vector
		typedef typename algebra_type::vector_type vector_type;

	public:
	///	default Constructor
		ILocalToGlobalMapper() {};

	///	send local entries to global matrix
		virtual void add_local_vec_to_global(vector_type& vec, const LocalVector& lvec,
				ConstSmartPtr<DoFDistribution> dd) = 0;

	///	send local entries to global rhs
		virtual void add_local_mat_to_global(matrix_type& mat, const LocalMatrix& lmat,
				ConstSmartPtr<DoFDistribution> dd) = 0;

	///	modifies local solution vector for adapted defect computation
		virtual void modify_LocalSol(LocalVector& vecMod, const LocalVector& lvec, ConstSmartPtr<DoFDistribution> dd) = 0;

	///	virtual destructor
		virtual ~ILocalToGlobalMapper() {};

};


} // end namespace ug

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__LOCAL_TO_GLOBAL_MAPPER__*/
