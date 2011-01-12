/**
 * \file amg_debug.h
 *
 * \author Martin Rupp
 *
 * \date 24.11.10
 *
 * Goethe-Center for Scientific Computing 2010.
 */


#ifndef __H__LIB_ALGEBRA__FAMG_SOLVER__FAMG_TESTVECTORS_H__
#define __H__LIB_ALGEBRA__FAMG_SOLVER__FAMG_TESTVECTORS_H__

namespace ug {


static void get_testvector_constant(stdvector<double> &testvector, stdvector<size_t> &N, double &i_value, size_t i)
{
	for(size_t j=0; j<N.size(); j++)
		testvector[j] = 1.0;
	i_value = 1.0;
}


/*
template<typename value_type, typename pos_type>
void get_testvector_value_special(value_type &value, const pos_type &myPos, const pos_type &referencePos)
{
	BlockRef(value, TcomponentOut) = myPos[TcomponentIn] - referencePos[TcomponentIn];
}


inline void get_testvector_xx(DenseVector<VariableArray1<double> > &testvector, stdvector<size_t> &N2, int i_index)
{
	return get_testvector_special<0, 0>(testvector, N2, i_index);
}
inline void get_testvector_xy(DenseVector<VariableArray1<double> > &testvector, stdvector<size_t> &N2, int i_index)
{
	return get_testvector_special<0, 1>(testvector, N2, i_index);
}
...
2d :
(x 0), (y 0), (0 x), (0 y)
3d :
(x 0 0) (0 y 0) (0 0 z)  ..?

*/

}

#endif // __H__LIB_ALGEBRA__FAMG_SOLVER__FAMG_TESTVECTORS_H__
