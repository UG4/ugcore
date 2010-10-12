

#include "../../ug_bridge.h"
#include "common/common.h"
#include "lib_discretization/lib_discretization.h"

namespace ug
{
namespace bridge
{

namespace {

bool c_BNDCond( number& value, const MathVector<2>& x, number time = 0.0)
{
	if(x[1] == 150)
		if(x[0] > 150 && x[0] < 450)
		{
			value = 1.;
			return true;
		}
	if(x[1] == 0)
	{
		value = 0.;
		return true;
	}

	// else neumann
	value = 0.;
	return false;
}


bool p_BNDCond( number& value, const MathVector<2>& x, number time = 0.0)
{
	// check if one of the upper corners, then dirichlet
	if(x[1] == 150 && (x[0] == 0.0 || x[0] == 600))
	{
		value = 9810 * (150 - x[1]);
		return true;
	}

	// else neumann value
	value = 0.0;
	return false;
}

} // end unnnamed namespace

class c_ElderDirichletBoundaryFunction : public DirichletBoundaryFunction<2>
{
	public:
	//	Function Type
		typedef bool (*Boundary_fct)(number&, const MathVector<2>&, number);

	public:
		virtual Boundary_fct get_bnd_function() const
		{
			return &c_BNDCond;
		}
};

class p_ElderDirichletBoundaryFunction : public DirichletBoundaryFunction<2>
{
	public:
	//	Function Type
		typedef bool (*Boundary_fct)(number&, const MathVector<2>&, number);

	public:
		virtual Boundary_fct get_bnd_function() const
		{
			return &p_BNDCond;
		}
};

namespace {

inline void Porosity(number& n)
{
	n = 0.1;
}
inline void Viscosity(number& visco, number c)
{
	visco = 1e-3;
}
inline void Density(number& density, number c)
{
	density = 1e3 + 0.2e3 * c;
}
inline void D_Density(number& Ddensity, number c)
{
	Ddensity = 0.2e3;
}
inline void Mol_Diff(MathMatrix<2,2>& D)
{
	D[0][0] = 3.565e-6;
	D[1][0] = 0.0;
	D[0][1] = 0.0;
	D[1][1] = 3.565e-6;
}
inline void Permeability(MathMatrix<2,2>& K)
{
	number perm = 4.845e-13;
	//perm = 1e-11;

	K[0][0] = perm;
	K[1][0] = 0.0;
	K[0][1] = 0.0;
	K[1][1] = perm;
}

inline void Gravity(MathVector<2>& gravity)
{
	gravity[0] = 0.0;
	gravity[1] = -9.81;
}


class ElderUserFunction : public IDensityDrivenFlowUserFunction<2>
{
	public:
		static const int dim = 2;

	//	Function Types
		typedef void (*Pososity_fct)(number&);
		typedef void (*Viscosity_fct)(number&, number);
		typedef void (*Density_fct)(number&, number);
		typedef void (*D_Density_fct)(number&, number);
		typedef void (*Mol_Diff_Tensor_fct)(MathMatrix<dim,dim>&);
		typedef void (*Permeability_Tensor_fct)(MathMatrix<dim,dim>&);
		typedef void (*Gravity_fct)(MathVector<dim>&);

	public:
		virtual Pososity_fct get_porosity_function() const {return &Porosity;}
		virtual Viscosity_fct get_viscosity_function() const {return &Viscosity;}
		virtual Density_fct get_density_function() const {return &Density;}
		virtual D_Density_fct get_d_density_function() const {return &D_Density;}
		virtual Mol_Diff_Tensor_fct get_mol_diff_tensor_function() const {return &Mol_Diff;}
		virtual Permeability_Tensor_fct get_perm_tensor_function() const {return &Permeability;}
		virtual Gravity_fct get_gravity_function() const {return Gravity;}
};

} // end unnamed namespace

namespace {

bool c_InitValues(const MathVector<2>& x, number& res)
{
	res = 0.0;

	if(x[1] == 150)
		if(x[0] > 150 && x[0] < 450)
			res = 1.0;

	return true;
}

bool p_InitValues(const MathVector<2>& x, number& res)
{
	res = 9810 * (150 - x[1]);
	return true;
}

} // end unnamed namespace


template <typename TGridFunction>
class InterpolateElder
{
	public:
	bool invoke(TGridFunction& u)
	{
		if(!ug::InterpolateFunction<TGridFunction>(&c_InitValues, u, 0)) return false;
		if(!ug::InterpolateFunction<TGridFunction>(&p_InitValues, u, 1)) return false;
		return true;
	}
};


void RegisterElderUserFunctions(Registry& reg)
{

//	DirichletBoundaryFunction
	{
		reg.add_class_<c_ElderDirichletBoundaryFunction, DirichletBoundaryFunction<2> >("cElderDirichletBoundaryFunction2d")
			.add_constructor();
		reg.add_class_<p_ElderDirichletBoundaryFunction, DirichletBoundaryFunction<2> >("pElderDirichletBoundaryFunction2d")
			.add_constructor();
	}

//	DensityDrivenUserFunction
	{
		reg.add_class_<IDensityDrivenFlowUserFunction<2> >("IDensityDrivenFlowUserFunction2d");

		reg.add_class_<ElderUserFunction, IDensityDrivenFlowUserFunction<2> >("ElderUserFunction2d")
			.add_constructor();
	}

//	Interpolation
	{
		typedef Domain<2, MultiGrid, MGSubsetHandler> domain_type;
		typedef P1ConformDoFDistribution dof_distribution_type;
		typedef MartinAlgebra algebra_type;

#ifdef UG_PARALLEL
		typedef ParallelGridFunction<GridFunction<domain_type, dof_distribution_type, algebra_type> > T;
#else
		typedef GridFunction<domain_type, dof_distribution_type, algebra_type> T;
#endif

		reg.add_class_<InterpolateElder<T> >("InterpolateElder")
			.add_constructor()
			.add_method("invoke", &InterpolateElder<T>::invoke);
	}


}

} // end namespace
} // end namepace
