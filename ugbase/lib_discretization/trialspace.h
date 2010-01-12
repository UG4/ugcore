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
#include <cassert>

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
	virtual bool evaluateShape(int nrShapeFct, const MathVector< RefDim >& locPos, number& value) const = 0;
	virtual bool evaluateShape(int nrShapeFct, MathVector< RefDim >[], number values[], int n) const = 0;
	virtual bool evaluateShapeGrad(int nrShapeFct, const MathVector< RefDim >& locPos, MathVector< RefDim >& value) const = 0;
	virtual bool positionOfDoF(int nrShapeFct, MathVector< RefDim >& value) const = 0;
	virtual uint order() const = 0;
	virtual uint num_dofs() const = 0;
	virtual ~TrialSpace()
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
		virtual bool evaluateShape(int nrShapeFct, const MathVector< RefDim >& locPos, number& value) const;
		virtual bool evaluateShape(int nrShapeFct, MathVector< RefDim > locPos[], number values[], int n) const;
		virtual bool evaluateShapeGrad(int nrShapeFct, const MathVector< RefDim >& locPos, MathVector< RefDim >& value) const;
		virtual bool positionOfDoF(int nrShapeFct, MathVector< RefDim >& value) const;
		virtual uint order() const { return _order;	}
		virtual uint num_dofs() const { return nsh;	}

		virtual ~P1conform()
		{};
	private:
		static const uint _order = 1;

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

		inline static const TrialSpace<TElem>& get_TrialSpace(TrialSpaceType type)
		{
			static P1conform<TElem> P1Conform;

			if(type == TST_P1CONFORM)
				return P1Conform;
			else
				assert(0 && "TrialSpaceType not implememted. Aborting. \n");
		}

	public:
		static const TrialSpace<TElem>& TrialSpace(TrialSpaceType type)
		{
			return inst().get_TrialSpace(type);
		}
};


} /* end namespace libDiscretization */



#endif /* __H__LIBDISCRETIZATION__TRIALSPACE__ */
