/*
 * trialspace.h
 *
 *  Created on: 12.05.2009
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__TRIALSPACE__
#define __H__LIBDISCRETIZATION__TRIALSPACE__

#include "referenceelement.h"
#include "dofhandler.h"
#include "common/math/ugmath.h"
#include "lib_grid/lib_grid.h"

namespace ug {

enum TrialSpaceType {
	TST_INVALID = 0,
	TST_P1CONFORM
};


template <typename TElem>
class TrialSpace
{
public:
	static const std::size_t RefDim = geometry_traits<TElem>::BASE_OBJECT_TYPE_ID;

public:
	virtual bool evaluateShape(int nrShapeFct, MathVector< RefDim > locPos, number& value) = 0;
	virtual bool evaluateShape(int nrShapeFct, MathVector< RefDim >[], number values[], int n) = 0;
	virtual bool evaluateShapeGrad(int nrShapeFct, MathVector< RefDim > locPos, MathVector< RefDim >& value) = 0;
	virtual bool positionOfDoF(int nrShapeFct, MathVector< RefDim >& value) = 0;
	virtual uint order() = 0;
	virtual ~TrialSpace<TElem>()
	{};
};

template <typename TTrialSpace>
class trialspace_traits
{};

template <typename TElem>
class P1conform : public TrialSpace<TElem>{
protected:
		static const std::size_t RefDim = TrialSpace<TElem>::RefDim;
		static const std::size_t nsh = reference_element_traits<TElem>::NumberCorners;

public:
		virtual bool evaluateShape(int nrShapeFct, MathVector< RefDim > locPos, number& value);
		virtual bool evaluateShape(int nrShapeFct, MathVector< RefDim > locPos[], number values[], int n);
		virtual bool evaluateShapeGrad(int nrShapeFct, MathVector< RefDim > locPos, MathVector< RefDim >& value);
		virtual bool positionOfDoF(int nrShapeFct, MathVector< RefDim >& value);
		virtual uint order()
		{
			return m_order;
		}

		virtual ~P1conform<TElem>()
		{};
	private:
		static const uint m_order = 1;

};

template<>
template<typename TElem>
class trialspace_traits< P1conform<TElem> >
{
	public:

};


// Singleton, holding all Trial Spaces available
template <typename TElem>
class TrialSpaces {

	private:

		static TrialSpaces& inst()
		{
			static TrialSpaces myInst;
			return myInst;
		};

		// private constructor
		TrialSpaces()
		{};

		inline static TrialSpace<TElem>& get_TrialSpace(TrialSpaceType type)
		{
			static P1conform<TElem> P1Conform;

			if(type == TST_P1CONFORM)
				return P1Conform;
		}

	public:
		static TrialSpace<TElem>& TrialSpace(TrialSpaceType type)
		{
			return inst().get_TrialSpace(type);
		}
};


} /* end namespace libDiscretization */



#endif /* __H__LIBDISCRETIZATION__TRIALSPACE__ */
