/*
 * reference_element_mapping.h
 *
 *  Created on: 13.05.2009
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__REFERENCE_ELEMENT__REFERENCE_ELEMENT_MAPPING__
#define __H__LIBDISCRETIZATION__REFERENCE_ELEMENT__REFERENCE_ELEMENT_MAPPING__

#include <cassert>
#include <iostream>
#include <sstream>
#include "common/common.h"
#include "common/math/ugmath.h"
#include "lib_grid/lg_base.h"

namespace ug{

/**
 * This class describes the mapping from a reference element into the real
 * (physical) world. The mapping is initialized by the physical positions of
 * the vertices of the real world element. The order of those points must be
 * given as indicated by the corresponding reference element.
 *
 * Let \f$R\f$ be the reference element and \f$T\f$ be the element. Then, the
 * reference mapping is a mapping:
 * \f[
 * 	\phi:	R \mapsto T
 * \f]
 *
 * \tparam	TRefElem		reference element
 * \tparam	TWorldDim		world dimension
 */
template <typename TRefElem, int TWorldDim>
class ReferenceMapping
{
	public:
	///	world dimension (range space dimension)
		static const int worldDim = TWorldDim;

	///	reference dimension (domain space dimension)
		static const int dim = TRefElem::dim;

	///	flag if mapping is linear (i.e. Jacobian does not depend on x)
		static const bool isLinear = false;

	public:
	///	Default Constructor
		ReferenceMapping();

	///	Constructor setting the corners of the element
		ReferenceMapping(const MathVector<worldDim>* vCorner);

	///	refresh mapping for new set of corners
		void update(const MathVector<worldDim>* vCorner);

	///	map local coordinate to global coordinate
		void local_to_global(MathVector<worldDim>& globPos,
		                     const MathVector<dim> locPos) const;

	///	returns transposed of jacobian
		void jacobian_transposed(MathMatrix<dim, worldDim>& JT,
		                         const MathVector<dim> locPos) const;

	///	returns transposed of the inverse of the jacobian
		void jacobian_transposed_inverse(MathMatrix<worldDim, dim>& JTInv,
		                                 const MathVector<dim> locPos) const;

	///	returns the determinate of the jacobian
		number jacobian_det(const MathVector<dim> locPos) const;
};


///	virtual base class for reference mappings
/**
 * This class is the base class for reference mappings in order to make them
 * selectable through a provider (on the price of virtual functions).
 *
 * \tparam	TDim		reference element dimensino
 * \tparam	TWorldDim	(physical) world dimension
 */
template <int TDim, int TWorldDim = TDim>
class DimReferenceMapping
{
	public:
	///	world dimension (range space dimension)
		static const int worldDim = TWorldDim;

	///	reference dimension (domain space dimension)
		static const int dim = TDim;

	public:
	///	returns if mapping is affine
		virtual bool is_linear() const = 0;

	///	refresh mapping for new set of corners
		virtual void update(const MathVector<worldDim>* vCorner) = 0;

	///	map local coordinate to global coordinate
		virtual void local_to_global(MathVector<worldDim>& globPos,
		                             const MathVector<dim> locPos) const = 0;

	///	map n local coordinate to global coordinate
		virtual void local_to_global(MathVector<worldDim>* globPos,
									 const MathVector<dim>* locPos, size_t n) const = 0;

	///	returns transposed of jacobian
		virtual void jacobian_transposed(MathMatrix<dim, worldDim>& JT,
		                                 const MathVector<dim> locPos) const = 0;

	///	returns transposed of jacobian for n local positions
		virtual void jacobian_transposed(MathMatrix<dim, worldDim>* JT,
										 const MathVector<dim>* locPos, size_t n) const = 0;

	///	returns transposed of the inverse of the jacobian
		virtual void jacobian_transposed_inverse(MathMatrix<worldDim, dim>& JTInv,
		                                         const MathVector<dim> locPos) const = 0;

	///	returns transposed of the inverse of the jacobian for n local positions
		virtual void jacobian_transposed_inverse(MathMatrix<worldDim, dim>* JTInv,
												 const MathVector<dim>* locPos, size_t n) const = 0;

	///	returns the determinate of the jacobian
		virtual number jacobian_det(const MathVector<dim> locPos) const = 0;

	///	returns the determinate of the jacobian for n local positions
		virtual void jacobian_det(number* det, const MathVector<dim>* locPos, size_t n) const = 0;
};


/// wrapper of a ReferenceElementMapping into the virtual base class
template <typename TRefMapping>
class DimReferenceMappingWrapper
	: public DimReferenceMapping<TRefMapping::dim, TRefMapping::worldDim>,
	  	  	  TRefMapping
{
	public:
	///	world dimension (range space dimension)
		static const int worldDim = TRefMapping::worldDim;

	///	reference dimension (domain space dimension)
		static const int dim = TRefMapping::dim;

	public:
	///	returns if mapping is affine
		virtual bool is_linear() const {return TRefMapping::isLinear;}

	///	refresh mapping for new set of corners
		virtual void update(const MathVector<worldDim>* vCorner)
		{
			TRefMapping::update(vCorner);
		}

	///	map local coordinate to global coordinate
		virtual void local_to_global(MathVector<worldDim>& globPos,
		                             const MathVector<dim> locPos) const
		{
			TRefMapping::local_to_global(globPos, locPos);
		}

	///	map n local coordinate to global coordinate
		virtual void local_to_global(MathVector<worldDim>* globPos,
									 const MathVector<dim>* locPos, size_t n) const
		{
			for(size_t ip = 0; ip < n; ++ip)
				TRefMapping::local_to_global(globPos[ip], locPos[ip]);
		}

	///	returns transposed of jacobian
		virtual void jacobian_transposed(MathMatrix<dim, worldDim>& JT,
		                                 const MathVector<dim> locPos) const
		{
			TRefMapping::jacobian_transposed(JT, locPos);
		}

	///	returns transposed of jacobian for n local positions
		virtual void jacobian_transposed(MathMatrix<dim, worldDim>* JT,
										 const MathVector<dim>* locPos, size_t n) const
		{
			if(TRefMapping::isLinear)
			{
				UG_ASSERT(n > 0, "No point specified in mapping.");
				TRefMapping::jacobian_transposed(JT[0], locPos[0]);
				for(size_t ip = 1; ip < n; ++ip)
					JT[ip] = JT[0];
			}
			else
			{
				for(size_t ip = 0; ip < n; ++ip)
					TRefMapping::jacobian_transposed(JT[ip], locPos[ip]);
			}
		}

	///	returns transposed of the inverse of the jacobian
		virtual void jacobian_transposed_inverse(MathMatrix<worldDim, dim>& JTInv,
		                                         const MathVector<dim> locPos) const
		{
			TRefMapping::jacobian_transposed_inverse(JTInv, locPos);
		}

	///	returns transposed of the inverse of the jacobian for n local positions
		virtual void jacobian_transposed_inverse(MathMatrix<worldDim, dim>* JTInv,
												 const MathVector<dim>* locPos, size_t n) const
		{
			if(TRefMapping::isLinear)
			{
				UG_ASSERT(n > 0, "No point specified in mapping.");
				TRefMapping::jacobian_transposed_inverse(JTInv[0], locPos[0]);
				for(size_t ip = 1; ip < n; ++ip)
					JTInv[ip] = JTInv[0];
			}
			else
			{
				for(size_t ip = 0; ip < n; ++ip)
					TRefMapping::jacobian_transposed_inverse(JTInv[ip], locPos[ip]);
			}
		}

	///	returns the determinate of the jacobian
		virtual number jacobian_det(const MathVector<dim> locPos) const
		{
			return TRefMapping::jacobian_det(locPos);
		}

	///	returns the determinate of the jacobian for n local positions
		virtual void jacobian_det(number* det, const MathVector<dim>* locPos, size_t n) const
		{
			if(TRefMapping::isLinear)
			{
				UG_ASSERT(n > 0, "No point specified in mapping.");
				det[0] = TRefMapping::jacobian_det(locPos[0]);
				for(size_t ip = 1; ip < n; ++ip)
					det[ip] = det[0];
			}
			else
			{
				for(size_t ip = 0; ip < n; ++ip)
					det[ip] = TRefMapping::jacobian_det(locPos[ip]);
			}
		}
};


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

struct UG_ERROR_ReferenceMappingMissing : public UGFatalError
{
	UG_ERROR_ReferenceMappingMissing(int dim_, int worldDim_, ReferenceObjectID roid_)
		: UGFatalError(""), dim(dim_), worldDim(worldDim_), roid(roid_)
	{
		std::stringstream ss; ss << "ReferenceMapping not found for "<<roid<<
									" from R^"<<dim<<" to R^"<<worldDim;
		UGFatalError::set_msg(ss.str());
	}
	int dim;
	int worldDim;
	ReferenceObjectID roid;
};

/// class to provide reference mappings
/**
 *	This class provides references mappings. It is implemented as a Singleton.
 */
class ReferenceMappingProvider {
	private:
	// 	disallow private constructor
		ReferenceMappingProvider();

	// disallow copy and assignment (intentionally left unimplemented)
		ReferenceMappingProvider(const ReferenceMappingProvider&);
		ReferenceMappingProvider& operator=(const ReferenceMappingProvider&);

	// 	private destructor
		~ReferenceMappingProvider(){};

	// 	Singleton provider
		static ReferenceMappingProvider& inst()
		{
			static ReferenceMappingProvider myInst;
			return myInst;
		};

	//	This is very dirty implementation, since casting to void. But, it is
	//	efficient, easy and typesafe. Maybe it should be changed to something
	//	inherent typesafe (not using casts) some day
	//	holding all mappings (worldDim x dim x roid)
		void* m_vvvMapping[4][4][NUM_REFERENCE_OBJECTS];

	//	casts void to map
		template <int TDim, int TWorldDim>
		DimReferenceMapping<TDim, TWorldDim>* get_mapping(ReferenceObjectID roid)
		{
			UG_STATIC_ASSERT(TDim <= 3, only_implemented_for_ref_dim_smaller_equal_3);
			UG_STATIC_ASSERT(TWorldDim <= 3, only_implemented_for_ref_dim_smaller_equal_3);
			UG_ASSERT(roid < NUM_REFERENCE_OBJECTS, "Roid specified incorrectly.");
			UG_ASSERT(roid >= 0, "Roid specified incorrectly.");
			return reinterpret_cast<DimReferenceMapping<TDim, TWorldDim>*>
						(m_vvvMapping[TDim][TWorldDim][roid]);
		}

	//	casts map to void
		template <int TDim, int TWorldDim>
		void set_mapping(ReferenceObjectID roid, DimReferenceMapping<TDim, TWorldDim>& map)
		{
			m_vvvMapping[TDim][TWorldDim][roid] = reinterpret_cast<void*>(&map);
		}

	public:
	///	returns a reference to a DimReferenceMapping
	/**
	 * This class returns a reference mapping for a ReferenceObjectID. The
	 * reference dimension and the world dimension must be chosen as
	 * template arguments. An exception is throw if such an mapping does not
	 * exist.
	 *
	 * \param[in]	roid		Reference Object ID
	 * \tparam		TDim		reference element dimension
	 * \tparam		TWorldDim	(physical) world dimension
	 */
		template <int TDim, int TWorldDim>
		static DimReferenceMapping<TDim, TWorldDim>& get(ReferenceObjectID roid)
		{
			DimReferenceMapping<TDim, TWorldDim>* pMap = inst().get_mapping<TDim, TWorldDim>(roid);
			if(!pMap) throw(UG_ERROR_ReferenceMappingMissing(TDim, TWorldDim, roid));
			else return *pMap;
		}
};


} // end namespace ug

#endif /* __H__LIBDISCRETIZATION__REFERENCE_ELEMENT__REFERENCE_ELEMENT_MAPPING__ */
