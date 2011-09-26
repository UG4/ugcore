// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// Somewhen back in 2005-2007

#ifndef __UG__LIB_GRID_NVEC_MATH_TEMPLATES__
#define __UG__LIB_GRID_NVEC_MATH_TEMPLATES__

namespace ug{
namespace libgrid_simplealg{

template <class Number>
void NVecAdd(Number* vOut, Number* v1, Number* v2, Number fac2, int n)
{
	for(int i = 0; i < n; ++i)
		vOut[i] = v1[i] + fac2 * v2[i];
}

template <class Number>
void NVecSub(Number* vOut, Number* v1, Number* v2, Number fac2, int n)
{
	for(int i = 0; i < n; ++i)
		vOut[i] = v1[i] - fac2 * v2[i];
}

template <class Number>
Number NVecDot(Number* v1, Number* v2, int n)
{
	Number r = 0;
	for(int i = 0; i < n; ++i)
		r += v1[i] * v2[i];
	return r;
}

template <class Number>
void NVecCopy(Number* vDest, Number* vSrc, int n)
{
	for(int i = 0; i < n; ++i)
		vDest[i] = vSrc[i];
}

template <class Number>
void NVecMult(Number* vOut, Number s, Number* v, int n)
{
	for(int i = 0; i < n; ++i)
		vOut[i] = s * v[i];
}

template <class Number>
Number NVecNorm2(Number* v, int n)
{
	Number r = 0;
	for(int i = 0; i < n; ++i)
		r += v[i] * v[i];
	return (Number)sqrt(r);
}

}// end of namespace
}// end of namespace
#endif
