

#include "ug_bridge/ug_bridge.h"
#include "common/common.h"
#include "lib_discretization/lib_discretization.h"

namespace ug
{
namespace bridge
{

namespace {


static number viscosity = 1e-3;


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


class ElderUserFunction : public IDensityDrivenFlowUserFunction<2>
{
	public:
		static const int dim = 2;

	//	Function Types
		typedef void (*Viscosity_fct)(number&, number);
		typedef void (*Density_fct)(number&, number);
		typedef void (*D_Density_fct)(number&, number);

	public:
		virtual Viscosity_fct get_viscosity_function() const {return &Viscosity;}
		virtual Density_fct get_density_function() const {return &Density;}
		virtual D_Density_fct get_d_density_function() const {return &D_Density;}

		void set_values(number viscosity_)
		{
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
						"", "Viscosity||invokeOnChange=true");
	}
}

} // end namespace
} // end namepace
