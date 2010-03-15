/*
 * local_dof_pattern.h
 *
 *  Created on: 17.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__LOCAL_DOF_PATTERN__
#define __H__LIBDISCRETIZATION__LOCAL_DOF_PATTERN__

#include "../reference_element/reference_elements.h"

namespace ug{

/**
 * Gives number of dofs for a Reference Element and its sub elements
 */
template <typename TRefElem>
class LocalDoFPattern
{
	public:
		// constructor (init with 0 dofs on every ref element)
		LocalDoFPattern();

		// sets the number of dofs
		void set_num_dofs(ReferenceElementType type, int num);

		// returns the total number of dofs on the finite element
		inline int total_num_dofs() const;

		// returns the number of dof on a geometric object
		inline int num_dofs(ReferenceElementType type) const;

	private:
		int m_num_dof[NUM_REFERENCE_ELEMENTS];
		int m_total_num_dofs;
};

/*
 * Gives number of dofs for all Reference Elements used in a grid.
 * Only possible, if distribution is equal on all elements
 */
class ContinuousDoFPattern
{
	public:
		// constructor (init with 0)
		ContinuousDoFPattern();

		// copy constructor
		ContinuousDoFPattern(const ContinuousDoFPattern& item);

		// assignment
		ContinuousDoFPattern& operator=(const ContinuousDoFPattern& item);

		// add
		ContinuousDoFPattern& operator+=(const ContinuousDoFPattern& item);

		// subtract
		ContinuousDoFPattern& operator-=(const ContinuousDoFPattern& item);

		// add a LocalDoFPattern
		template <typename TRefElem>
		bool add_local_dof_pattern(const LocalDoFPattern<TRefElem>& p);

		// returns the number of dofs needed on a Reference Element
		inline int total_num_dofs(ReferenceElementType type) const;

		// returns the number of dofs needed on a Reference Element
		inline int num_dofs(ReferenceElementType type) const;

		// comparison operator
		friend bool operator==(const ContinuousDoFPattern& lhs, const ContinuousDoFPattern& rhs);

	private:
		int m_num_dofs[NUM_REFERENCE_ELEMENTS];
		int m_total_num_dofs[NUM_REFERENCE_ELEMENTS];
		int m_dim;
};

bool operator==(const ContinuousDoFPattern& lhs, const ContinuousDoFPattern& rhs);

}

#include "local_dof_pattern_impl.h"

#endif /* __H__LIBDISCRETIZATION__LOCAL_DOF_PATTERN__ */
