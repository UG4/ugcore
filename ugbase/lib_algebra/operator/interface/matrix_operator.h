
#ifndef __H__LIB_ALGEBRA__OPERATOR__INTERFACE__MATRIX_OPERATOR__
#define __H__LIB_ALGEBRA__OPERATOR__INTERFACE__MATRIX_OPERATOR__

#include "linear_operator.h"
#include "lib_algebra/common/operations_mat/matrix_algebra_types.h"

namespace ug{


///////////////////////////////////////////////////////////////////////////////
// Matrix based linear operator
///////////////////////////////////////////////////////////////////////////////

template <typename M, typename X, typename Y = X>
class MatrixOperator :	public virtual ILinearOperator<X,Y>,
						public M
{
	public:
	// 	Domain space
		typedef X domain_function_type;

	// 	Range space
		typedef Y codomain_function_type;

	// 	Matrix type
		typedef M matrix_type;

	public:
	// 	Init Operator J(u)
		virtual void init(const X& u) {}

	// 	Init Operator L
		virtual void init() {}

	// 	Apply Operator f = L*u (e.g. d = J(u)*c in iterative scheme)
		virtual void apply(Y& f, const X& u) {matrix_type::apply(f,u);}

	// 	Apply Operator, i.e. f = f - L*u;
		virtual void apply_sub(Y& f, const X& u) {matrix_type::matmul_minus(f,u);}

	// 	Access to matrix
		virtual M& get_matrix() {return *this;};
};

template<typename M, typename X, typename Y>
struct matrix_algebra_type_traits<MatrixOperator<M, X, Y> >
{
	enum
	{
		type = matrix_algebra_type_traits<M>::type
	};
};

} // end namespace ug
#endif /* __H__LIB_ALGEBRA__OPERATOR__INTERFACE__MATRIX_OPERATOR__ */
