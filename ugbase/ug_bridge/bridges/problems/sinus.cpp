

#include "../../ug_bridge.h"
#include "common/common.h"
#include "lib_discretization/lib_discretization.h"

namespace ug
{
namespace bridge
{


namespace{

template <size_t dim>
inline void DiffTensor(MathMatrix<dim, dim>& D, const MathVector<dim>& x, number time = 0.0)
{
	for(size_t i = 0; i < dim; ++i){
		for(size_t j = 0; j < dim; ++j){
			D[i][j] = 0;
		}
		D[i][i] = 1.0;
	}
}

template <size_t dim>
inline void ConvVel(MathVector<dim>& v, const MathVector<dim>& x, number time = 0.0)
{
	for(size_t i = 0; i < dim; ++i){
		v[i] = 0;
	}
}

template <size_t dim>
inline void Reaction(number& value, const MathVector<dim>& x, number time = 0.0)
{
	value = 0.0;
}

template <size_t dim>
inline void f(number& value, const MathVector<dim>& x, number time = 0.0)
{
	double s = 2*M_PI;
	value = 0;
	for(size_t i = 0; i < dim; ++i)
		value += s*s*sin(s*x[i]);
}

}

template <int dim>
class SinusConvDiffUserFunction : public IConvDiffUserFunction<2>
{
	public:
		typedef void (*Diff_Tensor_fct)(MathMatrix<dim,dim>&, const MathVector<dim>&, number);
		typedef void (*Conv_Vel_fct)(MathVector<dim>&, const MathVector<dim>&, number);
		typedef void (*Reaction_fct)(number&, const MathVector<dim>&, number);
		typedef void (*Rhs_fct)(number&, const MathVector<dim>&, number);

	public:
		virtual Diff_Tensor_fct get_diff_tensor_function() const {return &DiffTensor<dim>;}
		virtual Conv_Vel_fct get_velocity_function() const {return &ConvVel<dim>;}
		virtual Reaction_fct get_reaction_function() const {return &Reaction<dim>;}
		virtual Rhs_fct get_rhs_function() const {return &f<dim>;}
};

namespace{

template <int dim>
bool BNDCond( number& value, const MathVector<dim>& x, number time = 0.0)
{
	double s = 2*M_PI;
	value = 0;
	for(size_t i = 0; i < dim; ++i)
		value += sin(s*x[i]);

	return true;
}

}

template <int dim>
class SinusDirichletBoundaryFunction : public DirichletBoundaryFunction<dim>
{
	public:
	//	Function Type
		typedef bool (*Boundary_fct)(number&, const MathVector<dim>&, number);

	public:
		virtual Boundary_fct get_bnd_function() const
		{
			return &BNDCond<dim>;
		}
};

void RegisterSinusUserFunctions(Registry& reg)
{
//	DirichletBoundaryFunction
	{
		reg.add_class_<SinusDirichletBoundaryFunction<2>, DirichletBoundaryFunction<2> >("SinusDirichletBoundaryFunction2d")
			.add_constructor();
	}

//	ConvDiffUserFunction
	{
		reg.add_class_<IConvDiffUserFunction<2> >("ConvDiffUserFunction2d");

		reg.add_class_<SinusConvDiffUserFunction<2>, IConvDiffUserFunction<2> >("SinusConvDiffUserFunction2d")
			.add_constructor();
	}
}

} // end namespace
} // end namepace
