/*
 * local_dof_set.h
 *
 *  Created on: 17.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__LOCAL_FINITE_ELEMENT__LOCAL_DOF_SET__
#define __H__UG__LIB_DISC__LOCAL_FINITE_ELEMENT__LOCAL_DOF_SET__

#include <vector>
#include <map>

#include "local_finite_element_id.h"
#include "lib_disc/reference_element/reference_element_traits.h"
#include "lib_grid/grid/geometric_base_objects.h"

namespace ug{

/// \ingroup lib_disc_local_finite_elements
/// @{

/**
 * This class is used to store for a single degree of freedom (DoF) the location
 * within an element. For continuous finite elements the DoFs are usually
 * associated with a sub-geometric object of the element itself (e.g. a vertex).
 * This can be requested from this class, which stores the dimension of the
 * sub-element the DoF is located on, the id of the sub-element (w.r.t. to the
 * numbering in the reference elements) and an offset > 0 if there are more than
 * one DoFs associated with the same sub-element.
 */
class LocalDoF
{
	public:
	///	default constructor
		LocalDoF() : m_dim(-1), m_id(0), m_offset(0) {}

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

	///	sets the values
		void set(int dim, size_t id, size_t offset)
		{
			m_dim = dim; m_id = id; m_offset = offset;
		}

	///	returns the dimension of associated geometric object
		inline int dim() const {return m_dim;}

	///	returns the index for the geometric object (w.r.t reference element numbering)
		inline size_t id() const {return m_id;}

	///	returns the offset for the geometric object
		inline size_t offset() const {return m_offset;}

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
 */
class LocalDoFSet
{
	public:
	///	returns the reference dimension
		int dim() const;

	///	returns the Reference object id of the corresponding grid object
		virtual ReferenceObjectID roid() const = 0;

	///	type of shape functions
		virtual LFEID type() const = 0;

	///	returns the total number of dofs on the finite element
	/// \{
				size_t num_dof() const {return num_sh();}
		virtual size_t num_sh() const;
	/// \}

	///	returns the number of DoFs on a sub-geometric object type
		virtual size_t num_dof(ReferenceObjectID roid) const = 0;

	///	returns the DoFs storage
		virtual const LocalDoF& local_dof(size_t dof) const = 0;

	///	returns the number of DoFs on a sub-geometric object of dim and id
		size_t num_dof(int d, size_t id) const;

	///	virtual destructor
		virtual ~LocalDoFSet() {};
};

/**
 * Local DoF Set also providing to local position of the dofs (iff available)
 */
template <int TDim>
class DimLocalDoFSet : public LocalDoFSet
{
	public:
	/// constructor
		DimLocalDoFSet();

	///	returns if the local dof position are exact
		virtual bool exact_position_available() const = 0;

	///	local position of DoF i
	/**
	 * This function returns the local position of a DoF if possible.
	 * \param[in] 	i		number of DoF
	 * \param[out]	pos		Position of DoF
	 * \retval		true 	if position exists
	 * \retval		false 	if no meaningful position available
	 */
		virtual bool position(size_t i, MathVector<TDim>& pos) const = 0;
};

/// @}

/// writes to the output stream
std::ostream& operator<<(std::ostream& out,	const LocalDoFSet& v);

/**
 * Intersection of local dof sets
 */
class CommonLocalDoFSet
{
	///	indicate not set value
		enum{NOT_SPECIFIED = -1};

	public:
	///	constructor
		CommonLocalDoFSet() {clear();}

	///	reset all numbers of dofs to not set
		void clear();

	///	add a local dof set to the intersection
		void add(const LocalDoFSet& set);

	///	number of dofs on a reference element type
		size_t num_dof(ReferenceObjectID roid) const {
			UG_ASSERT(m_vNumDoF[roid] >= 0, "No DoF info for "<<roid);
			return m_vNumDoF[roid];
		}

	protected:
		int m_vNumDoF[NUM_REFERENCE_OBJECTS];
};

/// writes to the output stream
std::ostream& operator<<(std::ostream& out,	const CommonLocalDoFSet& v);


////////////////////////////////////////////////////////////////////////////////
//	Provider
////////////////////////////////////////////////////////////////////////////////

// LocalDoFSetProvider
/** class to provide local DoF sets
 *
 *	This class provides references to Local DoF sets.
 *	It is implemented as a Singleton.
 */
class LocalDoFSetProvider {
	private:
	// 	disallow private constructor
		LocalDoFSetProvider(){};

	// disallow copy and assignment (intentionally left unimplemented)
		LocalDoFSetProvider(const LocalDoFSetProvider&);
		LocalDoFSetProvider& operator=(const LocalDoFSetProvider&);

	// 	private destructor
		~LocalDoFSetProvider();

	// 	Singleton provider
		static LocalDoFSetProvider& inst()
		{
			static LocalDoFSetProvider myInst;
			return myInst;
		};

	private:
	//	creation of nedelec element set for a reference element
		template <typename TRefElem>
		static void create_set(const LFEID& id);

	// 	creation of set
		static void create_set(ReferenceObjectID roid, const LFEID& id);

	// 	creation of set
		static void create_set(const LFEID& id);

	//	type of map holding dof sets for a reference object id
		typedef std::map<LFEID, std::vector<ConstSmartPtr<LocalDoFSet> > > RoidMap;

	//	map holding dof sets for a reference object id
		static RoidMap m_mRoidDoFSet;

	//	type of map holding common dof set for roid of same dimension
		typedef std::map<LFEID, CommonLocalDoFSet> CommonMap;

	//	map holding common dof set for roid of same dimension
		static CommonMap m_mCommonDoFSet;

	public:
	/** register a local DoF set for a given reference element type
	 * This function is used to register a Local Shape Function set for an element
	 * type and the corresponding local DoF set id.
	 *
	 * \param[in]		id 		Identifier for local DoF set
	 * \param[in]		set		Local Shape Function Set to register
	 * \return			bool	true iff registration successful
	 */
		static void register_set(const LFEID& id, ConstSmartPtr<LocalDoFSet> set);

	/** unregister a local DoF set for a given reference element type
	 * This function is used to unregister a Local Shape Function set for an element
	 * type and the corresponding local DoF set id from this Provider.
	 *
	 * \param[in]		id 		Identifier for local DoF set
	 * \return			bool	true iff removal successful
	 */
		static bool unregister_set(const LFEID& id);

	///	returns the common dof set for all reference objects of a dimension
		static const CommonLocalDoFSet& get(const LFEID& id, bool bCreate = true);

	///	returns the local DoF set base for an id
		static const LocalDoFSet& get(ReferenceObjectID roid, const LFEID& id, bool bCreate = true);
};


} // end namespace ug

#endif /* __H__UG__LIB_DISC__LOCAL_FINITE_ELEMENT__LOCAL_DOF_SET__ */
