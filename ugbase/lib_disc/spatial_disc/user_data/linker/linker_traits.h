
#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__DATA_LINKER_TRAITS__
#define __H__UG__LIB_DISC__SPATIAL_DISC__DATA_LINKER_TRAITS__

#include "common/common.h"

namespace ug{

////////////////////////////////////////////////////////////////////////////////
// Linker Traits to allow generic programming
////////////////////////////////////////////////////////////////////////////////

/// Linker Traits
template <typename TData, typename TDataIn, typename TRet = TData>
struct linker_traits
{
///	computes out += s * in1 (with appropriate '*')
	static void mult_add(TData& out, const TData& in1, const TDataIn& s);
};

template <>
struct linker_traits<number, number>
{
	static void mult_add(number& out, const number& in1, const number& s)
	{
		out += in1 * s;
	}
};

template <std::size_t dim>
struct linker_traits< MathVector<dim>, number >
{
	static void mult_add(MathVector<dim>& out,
	                     const MathVector<dim>& in1,
						 const number& s)
	{
		VecScaleAppend(out, s, in1);
	}
};

template <std::size_t dim>
struct linker_traits< MathVector<dim>, MathMatrix<dim,dim> >
{
	static void mult_add(MathVector<dim>& out,
	                     const MathVector<dim>& in1,
						 const MathMatrix<dim,dim>& s)
	{
		MatVecMultAppend(out, s, in1);
	}
};

template <std::size_t dim>
struct linker_traits< MathVector<dim>,  MathVector<dim>, number >
{
	static void mult_add(number& out,
	                     const MathVector<dim>& in1,
						 const MathVector<dim>& s)
	{
		out = VecDot(s, in1);
	}
};

template <std::size_t dim>
struct linker_traits< MathMatrix<dim,dim>, number >
{
	static void mult_add(MathMatrix<dim,dim>& out,
	                     const MathMatrix<dim,dim>& in1,
						 const number& s)
	{
		MatScaleAppend(out, s, in1);
	}
};

template <std::size_t dim>
struct linker_traits< MathTensor<4,dim>, number >
{
	static void mult_add(MathTensor<4,dim>& out,
	                     const MathTensor<4, dim>& in1,
						 const number& s)
	{
		out = in1;
		UG_THROW("linker_traits for MathTensor4 not implemented");
	}
};

} // end namespace ug

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__DATA_LINKER_TRAITS__ */
