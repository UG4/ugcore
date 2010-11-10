

#include "ug_bridge/ug_bridge.h"
#include "common/common.h"
#include "lib_discretization/lib_discretization.h"

namespace ug
{
namespace bridge
{

namespace {


static number porosity = 0.1;
static number viscosity = 1e-3;

inline void Porosity(number& n)
{
	n = porosity;
}
inline void Viscosity(number& visco, number c)
{
	visco = viscosity;
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

		void set_values(number porosity_, number viscosity_)
		{
			porosity = porosity_;
			viscosity = viscosity_;
		}
};

} // end unnamed namespace

void RegisterElderUserFunctions(Registry& reg, const char* parentGroup)
{
	const char* grp = parentGroup;

//	DensityDrivenUserFunction
	{
		reg.add_class_<ElderUserFunction, IDensityDrivenFlowUserFunction<2> >("ElderUserFunction2d", grp)
			.add_constructor()
			.add_method("set_values|interactive=false", &ElderUserFunction::set_values,
						"", "Porosity||invokeOnChange=true#Viscosity||invokeOnChange=true");
	}
}

} // end namespace
} // end namepace
