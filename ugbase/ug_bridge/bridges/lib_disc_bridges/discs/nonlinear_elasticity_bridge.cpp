/*
 * nonlinear_elasticity_bridge.cpp
 *
 *  Created on: 08.06.2011
 *      Authors: andreasvogel
 */

// extern headers
#include <iostream>
#include <sstream>

// include bridge
#include "../../../ug_bridge.h"
#include "../../../registry.h"

// lib_algebra includes
#include "lib_algebra/algebra_selector.h"
#include "lib_algebra/algebra_types.h"
#include "lib_algebra/operator/operator_util.h"
#include "lib_algebra/operator/operator_interface.h"
#include "lib_algebra/operator/operator_inverse_interface.h"

// lib_disc includes
#include "lib_discretization/domain.h"
#include "lib_discretization/function_spaces/grid_function.h"
#include "lib_discretization/function_spaces/approximation_space.h"
#include "lib_discretization/dof_manager/p1conform/p1conform.h"

// fe1_nonlinear_elasticity includes
#include "lib_discretization/spatial_discretization/elem_disc/nonlinear_elasticity/fe1_nonlinear_elasticity.h"


namespace ug
{

namespace bridge
{


template <int dim>
class ElasticityTensorUserData
	: public IUserData<MathTensor<4,dim>, dim>
{
	///	Base class type
		typedef IUserData<MathTensor<4,dim>, dim> base_type;

		using base_type::num_series;
		using base_type::num_ip;
		using base_type::ip;
		using base_type::time;
		using base_type::value;

	public:
	//	Functor Type
		typedef typename base_type::functor_type functor_type;

	//	return functor
		virtual functor_type get_functor() const {return *this;}

	public:
	///	Constructor
		ElasticityTensorUserData() {}

	///	virtual destructor
		virtual ~ElasticityTensorUserData()	{}

	///	evaluates the data at a given point and time
		void operator() (MathTensor<4,dim>& C, const MathVector<dim>& x, number time = 0.0)
		{
			//filling the ElasticityTensor
			C.set(0.0);

			C[0][0][0][0] = 5.;
			C[0][0][1][2] = 12.;
		}
	///	prints Elasticity Tensor C at point x and time t
		void test_evaluate(number x, number y, number z, number t)
		{
			MathTensor<4,dim> C;
			MathVector<dim> v;
			v.x = x;
			if(v.size() > 1)
				v[1] = y;
			if(v.size() > 2)
				v[2] = z;

			this->operator()(C, v, t);
			//UG_LOG("Tensor: " << C);

		}
	///	implement as a IPData
		virtual void compute(bool computeDeriv = false)
		{
			for(size_t s = 0; s < num_series(); ++s)
				for(size_t i = 0; i < num_ip(s); ++i)
				{
					this->operator() (	value(s,i),
										ip(s, i),
										time());
				}
		}
};


template <typename TDomain, typename TAlgebra, typename TDoFDistribution>
void RegisterNonlinearElasticityObjects(Registry& reg, const char* parentGroup)
{
//	typedef domain
	typedef TDomain domain_type;
	typedef TDoFDistribution dof_distribution_type;
	typedef TAlgebra algebra_type;
	typedef typename algebra_type::vector_type vector_type;
	typedef typename algebra_type::matrix_type matrix_type;
	static const int dim = domain_type::dim;

	std::stringstream grpSS; grpSS << parentGroup << "/" << dim << "d";
	std::string grp = grpSS.str();

#ifdef UG_PARALLEL
		typedef ParallelGridFunction<GridFunction<domain_type, dof_distribution_type, algebra_type> > function_type;
#else
		typedef GridFunction<domain_type, dof_distribution_type, algebra_type> function_type;
#endif

/////////////////////////////////////////////////////////////////////////////
// Elem Disc
/////////////////////////////////////////////////////////////////////////////


	typedef FE1NonlinearElasticityElemDisc<domain_type, algebra_type> T;
	typedef IDomainElemDisc<domain_type, algebra_type> TBase;
	std::stringstream ss; ss << "FE1NonlinearElasticity" << dim << "d";
	reg.add_class_<T, TBase >(ss.str().c_str(), grp.c_str())
		.add_constructor()
		.add_method("set_elasticity_tensor", &T::set_elasticity_tensor);
		//methods, which are used in script-file
}


template <typename TAlgebra, typename TDoFDistribution>
bool RegisterNonlinearElasticityDisc(Registry& reg, const char* parentGroup)
{
	typedef TDoFDistribution dof_distribution_type;
	typedef TAlgebra algebra_type;
	typedef typename algebra_type::vector_type vector_type;
	typedef typename algebra_type::matrix_type matrix_type;

	try
	{
	//	get group string
		std::string grp = parentGroup; grp.append("/Discretization");

#ifdef UG_DIM_1
	//	Domain dependend part 1D
		{
			typedef Domain<1, MultiGrid, MGSubsetHandler> domain_type;
			RegisterNonlinearElasticityObjects<domain_type, algebra_type, dof_distribution_type>(reg, grp.c_str());
		}
#endif

#ifdef UG_DIM_2
	//	Domain dependend part 2D
		{
			typedef Domain<2, MultiGrid, MGSubsetHandler> domain_type;
			RegisterNonlinearElasticityObjects<domain_type, algebra_type, dof_distribution_type>(reg, grp.c_str());
		}
#endif

#ifdef UG_DIM_3
	//	Domain dependend part 3D
		{
			typedef Domain<3, MultiGrid, MGSubsetHandler> domain_type;
			RegisterNonlinearElasticityObjects<domain_type, algebra_type, dof_distribution_type>(reg, grp.c_str());
		}
#endif
	}
	catch(UG_REGISTRY_ERROR_RegistrationFailed ex)
	{
		UG_LOG("### ERROR in RegisterNonlinearElasticityDisc: "
				"Registration failed (using name " << ex.name << ").\n");
		return false;
	}

	return true;
}


template <int dim>
bool RegisterElasticityTensor(Registry& reg, const char* parentGroup)
{
	std::stringstream grpSS; grpSS << parentGroup << "/" << dim << "d";
	std::string grp = grpSS.str();

	typedef ElasticityTensorUserData<dim> T;
	typedef IUserData<MathTensor<4,dim>, dim> TBase;
	std::stringstream ss; ss << "ElasticityTensorUserData" << dim << "d";
	reg.add_class_<T, TBase >(ss.str().c_str(), grp.c_str())
		.add_constructor()
		.add_method("operator", &T::operator ())
		.add_method("test_evaluate", &T::test_evaluate);
		//methods, which are used in script-file

	return true;
}

bool RegisterElasticityTensor(Registry& reg, const char* parentGroup)
{
	try
	{
	//	get group string
		std::string grp = parentGroup; grp.append("/Discretization");

#ifdef UG_DIM_1
		RegisterElasticityTensor<1>(reg, parentGroup);
#endif
#ifdef UG_DIM_2
		RegisterElasticityTensor<2>(reg, parentGroup);
#endif
#ifdef UG_DIM_3
		RegisterElasticityTensor<3>(reg, parentGroup);
#endif
	}
	catch(UG_REGISTRY_ERROR_RegistrationFailed ex)
	{
		UG_LOG("### ERROR in RegisterElasticityTensor: "
				"Registration failed (using name " << ex.name << ").\n");
		return false;
	}
	return true;
}


bool RegisterDynamicNonlinearElasticityDisc(Registry& reg, int algebra_type, const char* parentGroup)
{
	bool bReturn = true;

	switch(algebra_type)
	{
	case eCPUAlgebra:		 		bReturn &= RegisterNonlinearElasticityDisc<CPUAlgebra, P1ConformDoFDistribution>(reg, parentGroup); break;
//	case eCPUBlockAlgebra2x2: 		bReturn &= RegisterNonlinearElasticityDisc<CPUBlockAlgebra<2>, GroupedP1ConformDoFDistribution>(reg, parentGroup); break;
//	case eCPUBlockAlgebra3x3: 		bReturn &= RegisterNonlinearElasticityDisc<CPUBlockAlgebra<3>, GroupedP1ConformDoFDistribution>(reg, parentGroup); break;
//	case eCPUBlockAlgebra4x4: 		bReturn &= RegisterNonlinearElasticityDisc<CPUBlockAlgebra<4>, GroupedP1ConformDoFDistribution>(reg, parentGroup); break;
//	case eCPUVariableBlockAlgebra: 	bReturn &= RegisterNonlinearElasticityDisc<CPUVariableBlockAlgebra, GroupedP1ConformDoFDistribution>(reg, parentGroup); break;
	default: UG_ASSERT(0, "Unsupported Algebra Type");
				UG_LOG("Unsupported Algebra Type requested.\n");
				return false;
	}

	bReturn &= RegisterElasticityTensor(reg, parentGroup);

	return bReturn;
}

}//	end of namespace ug
}//	end of namespace interface
