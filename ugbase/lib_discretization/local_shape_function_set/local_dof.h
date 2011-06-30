/*
 * local_dof.h
 *
 *  Created on: 17.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__LOCAL_DOF__
#define __H__LIBDISCRETIZATION__LOCAL_DOF__

#include "../reference_element/reference_element.h"

namespace ug{

/**
 * This class is used to store for a single degree of freedom (DoF) the location
 * within an element. For continuous finite elements the DoFs are usually
 * associated with a sub-geometric object of the element itself (e.g. a vertex).
 * This can be requested from this class.
 */
class LocalDoF
{
	public:
	///	constructor
	/**
	 * Create a pair describing the position of a DoF within the reference element.
	 *
	 * \param[in]	dim		dimension of sub-geometric object
	 * \param[in]	id		number of sub-geometric object (in the numbering
	 * 						used by the reference element)
	 * \param[in]	offset	if several DoFs are associated with the same
	 * 						sub-geometric object the offset specifies the number
	 * 						within all DoFs on that geometric object
	 */
		LocalDoF(int dim, size_t id, size_t offset)
			: m_dim(dim), m_id(id), m_offset(offset)
		{}

	///	returns the dimension of associated geometric object
		inline int dim(size_t sh) const {return m_dim;}

	///	returns the index for the geometric object (w.r.t reference element numbering)
		inline size_t id(size_t sh) const {return m_id;}

	///	returns the offset for the geometric object
		inline size_t offset(size_t sh) const {return m_offset;}

	protected:
	///	dimension of sub-geometric object
		int m_dim;

	///	id of sub-geometric object in counting of reference element
		size_t m_id;

	///	offset if several DoFs associated to the same geometric object
		size_t m_offset;
};

/**
 * This class provides the interface for the storage of degrees of freedom
 * on a finite element.
 *
 * \tparam 	TRefElem	Reference Element type
 */
template <typename TRefElem>
class ILocalDoFPattern
{
	public:
	///	returns the total number of dofs on the finite element
		virtual size_t num_sh() const;

	///	returns the DoFs storage
		virtual const LocalDoF& storage(size_t sh) const;

	///	returns if DoFs are associated with objects of the dimension
		virtual bool storage_use(int dim) const;
};

/**
 * This class provides a wrapper class into the ILocalDoFPattern interface in
 * order to make it available in that context by paying the price of
 * virtual functions
 */
template <typename TImpl>
class ILocalDoFPatternWrapper
	: public ILocalDoFPattern<typename TImpl::reference_element_type>,
	  private TImpl
{
	private:
	///	implementation type
		typedef TImpl ImplType;

	public:
	///	\copydoc ug::ILocalDoFPattern::num_sh()
		virtual size_t num_sh() const
		{
			return ImplType::num_sh();
		}

	///	\copydoc ug::ILocalDoFPattern::storage()
		virtual const LocalDoF& storage(size_t sh) const
		{
			return ImplType::storage(sh);
		}

	///	\copydoc ug::ILocalDoFPattern::storage_use()
		virtual bool storage_use(int dim) const
		{
			return ImplType::storage_use(dim);
		}

};

} // end namespace ug

#endif /* __H__LIBDISCRETIZATION__LOCAL_DOF__ */
