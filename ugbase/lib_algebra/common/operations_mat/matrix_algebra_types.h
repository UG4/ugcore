
#ifndef MATRIX_ALGEBRA_TYPES_H_
#define MATRIX_ALGEBRA_TYPES_H_

namespace ug{

/*
suppose you have your algebra with your vector class yourvector and your matrix class yourmatrix.
when your matrix class is rowwise multiplicable, and you want to do so, set
template<>
class matrix_algebra_type_traits<yourclass>
{
	static const matrix_algebra_type type = MATRIX_USE_ROW_FUNCTIONS;
};
and add a function
inline void yourmatrix::mat_mult_add_row(size_t row, yourvector::value_type &dest, double beta1, const yourvector &v) const;
to your class.

If you cannot or dont want to use rowwise multiplication, set
template<>
class matrix_algebra_type_traits<yourclass>
{
	static const matrix_algebra_type type = MATRIX_USE_GLOBAL_FUNCTIONS;
};

and implement (at least)
inline void MatMultAddDirect(yourvector &dest, const number &beta1, const yourmatrix &A1, const yourvector &w1);
inline void MatMultDirect(yourvector &dest, const number &beta1, const yourmatrix &A1, const yourvector &w1,
							number &alpha1, const yourvector &v1);
Other functions like MatMultAdd(dest, beta1, A1, w1, beta2, A2, w2, alpha1, v1) are then implemented through your functions.
(see class mat_operations_class<vector_t, matrix_t, MATRIX_USE_GLOBAL_FUNCTIONS> for all functions).

*/

enum matrix_algebra_type
{
	MATRIX_USE_ROW_FUNCTIONS,
	MATRIX_USE_GLOBAL_FUNCTIONS,
	MATRIX_USE_OPERATORS,
	MATRIX_USE_MEMBER_FUNCTIONS
};

template<typename vector_t, typename matrix_t, int type>
struct mat_operations_class;

template<typename T>
struct matrix_algebra_type_traits
{
	static const int type = MATRIX_USE_OPERATORS;
};


}

#endif /* MATRIX_ALGEBRA_TYPES_H_ */
