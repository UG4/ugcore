
#ifndef __H__UG__LIB_DISC__LOCAL_SHAPE_FUNCTION_SET__LAGRANGE__LAGRANGE_LOCAL_DOF__
#define __H__UG__LIB_DISC__LOCAL_SHAPE_FUNCTION_SET__LAGRANGE__LAGRANGE_LOCAL_DOF__

#include "lib_disc/local_finite_element/local_dof_set.h"

namespace ug{

/// returns number of DoFs on element type for order p
size_t LagrangeNumDoFs(const ReferenceObjectID elem, const size_t p);

///	returns number of DoFs Subelement for an element type and order p
size_t LagrangeNumDoFOnSub(const ReferenceObjectID elem,
                           const ReferenceObjectID sub, const size_t p);


///////////////////////////////////////////////////////////////////////////////
// Lagrange LocalDoFSets
///////////////////////////////////////////////////////////////////////////////

/// Lagrange DoF Set
template <typename TRefElem>
class LagrangeLDS : public LocalDoFSet
{
	public:
	///	constructor
		LagrangeLDS(size_t order = 1);

	///	sets the order
		void set_order(size_t order);

	///	returns the type of reference element
		ReferenceObjectID roid() const {return TRefElem::REFERENCE_OBJECT_ID;}

	///	returns the total number of DoFs on the finite element
		size_t num_dof() const {return m_vLocalDoF.size();};

	///	returns the number of DoFs on a sub-geometric object type
		size_t num_dof(ReferenceObjectID roid) const;

	///	returns the dof storage
		const LocalDoF& local_dof(size_t dof) const {return m_vLocalDoF[dof];}

	///	returns if the local dof position are exact
		bool exact_position_available() const {return true;};

	protected:
		size_t p;		///< order
		std::vector<LocalDoF> m_vLocalDoF; 	///< association to geom obj
};

} //namespace ug

#endif /* __H__UG__LIB_DISC__LOCAL_SHAPE_FUNCTION_SET__LAGRANGE__LAGRANGE_LOCAL_DOF__ */
