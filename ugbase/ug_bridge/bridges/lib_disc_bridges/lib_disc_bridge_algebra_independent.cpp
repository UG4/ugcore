/*
 * lib_disc_bridge_algebra_independent.cpp
 *
 *  Created on: 14.02.2011
 *      Author: andreasvogel
 */

// include brigde
#include "../../ug_bridge.h"

// common includes
#include "common/math/ugmath.h"
#include "common/math/misc/math_util.h"

// Lagrange Function Space
#include "lib_discretization/local_shape_function_set/lagrange/lagrange.h"

// P1ConformFunctionPattern
#include "lib_discretization/dof_manager/p1conform/p1conform.h"

namespace ug
{
namespace bridge
{

template <int dim, int p>
number func1(const MathVector<dim>& x)
{
	number res = 0.;
	for(size_t d = 0; d < dim; ++d)
		res += x[d];// pow(x[d], p) + 0.5;
	return res;
}

template <int dim, int p>
MathVector<dim> Dfunc1(const MathVector<dim>& x)
{
	MathVector<dim> res;
	for(size_t d = 0; d < dim; ++d)
		res[d] = 1.0;//p * pow(x[d], p-1);

	return res;
}

template <typename TRefElem, int p>
bool TestLagrangeSpacesElem(number& maxDiff)
{
	static const double Small = 1e-7;

	LagrangeLSFS<TRefElem, p> space;
	typedef LagrangeLSFS<TRefElem, p> type;

	UG_LOG("Testing Reference Element of Order "<< type::order << "\n");

	UG_LOG("  Dimension:    " << type::dim << "\n");
	UG_LOG("  Num Shapes:   " << type::nsh << "\n");

	if(type::nsh != space.num_sh())
	{
		UG_LOG("ERROR: Mismatch between type::nsh and num_sh()\n");
		return false;
	}

	UG_LOG("  Multi Indices are: ");
	if(space.num_sh() < 10)
	{
		UG_LOG("\n");
		for(size_t i = 0; i < space.num_sh(); ++i)
		{
			typename type::position_type pos;
			space.position(i, pos);
			UG_LOG("      " << i << ":  ");
			UG_LOG(space.multi_index(i));
			UG_LOG(" at " << pos << "\n");
		}
	}
	else
		UG_LOG(" [output skipped]\n");

//	Multi Indices
	UG_LOG("  Checking MultiIndex mappings ... ");
	for(size_t i = 0; i < space.num_sh(); ++i)
	{
		typename type::multi_index_type ind = space.multi_index(i);
		if(space.index(ind) != i)
		{
			UG_LOG("ERROR: MultiIndex map and back is not identity"
					" for index "<< i << " and MultiIndex " << ind << ".\n");
			return false;
		}
	}
	UG_LOG("PASSED.\n");

//	Shapes
	UG_LOG("  Checking shapes at dofs are kronecker delta ... ");
	for(size_t i = 0; i < space.num_sh(); ++i)
	{
		typename type::position_type pos;
		if(!space.position(i, pos))
		{
			UG_LOG("ERROR: Cannot get position of shape " << i << "\n");
		}

		for(size_t j = 0; j < space.num_sh(); ++j)
		{
			if(i == j )
			{
				number diff = fabs(space.shape(j, pos)-1.0);
				maxDiff = std::max(maxDiff, diff);

				if(diff > Small)
				{
					UG_LOG("ERROR: Shape "<<j<<" at support point "<<i<<" is"
							" not 1.0 but "<<std::scientific <<space.shape(j, pos)<<"\n");
					return false;
				}
			}
			else
			{
				number diff = fabs(space.shape(j, pos));
				maxDiff = std::max(maxDiff, diff);
				if(diff > Small)
				{
					UG_LOG("ERROR: Shape "<<j<<" at support point "<<i<<" is"
							" not 0.0 but "<<std::scientific <<space.shape(j, pos)<<"\n");
					return false;
				}
			}
		}
	}
	UG_LOG("PASSED.\n");

//	sum of shapes at arbitrary point
	typename type::position_type pos;

	UG_LOG("  Testing Sum of Shapes at arbitrary points ... ");
	for(size_t k = 0; k < 20; ++k)
	{
		pos[0] = urand(0.0, 1.0);
		if(type::dim >= 2)
		{
			pos[1] = urand(0.0, 1.0 - pos[0]);
		}
		if(type::dim == 3)
		{
			pos[2] = urand(0.0, 1.0 - pos[0] - pos[1]);
		}

		std::vector<number> vShape(type::nsh);

		space.shapes(&vShape[0], pos);
		number sum = 0, sum2 = 0;
		for(size_t sh = 0; sh < space.num_sh(); ++sh)
		{
			sum += space.shape(sh, pos);
			sum2 += vShape[sh];
		}
		number diff = fabs(sum - 1.0);
		maxDiff = std::max(maxDiff, diff);
		if(diff > Small)
		{
			UG_LOG("ERROR: Sum of Shapes not 1.0 but differ " << diff <<".\n");
			return false;
		}

		if(fabs(sum - sum2) > Small)
		{
			UG_LOG("ERROR: Vector sum and sum not equal. \n");
			return false;
		}
	}
	UG_LOG("PASSED.\n");

	UG_LOG("  Testing interpolation of function ... ");

	static const int dim = type::dim;
	// interpolate
	std::vector<number> vInter(type::nsh);
	for(size_t sh = 0; sh < space.num_sh(); ++sh)
	{
		typename type::position_type pos;
		space.position(sh, pos);

		vInter[sh] = func1<dim, p>(pos);
	}

	// eval at arbitrary points
	for(size_t k = 0; k < 20; ++k)
	{
		pos[0] = urand(0.0, 1.0);
		if(type::dim >= 2)
		{
			pos[1] = urand(0.0, 1.0 - pos[0]);
		}
		if(type::dim == 3)
		{
			pos[2] = urand(0.0, 1.0 - pos[0] - pos[1]);
		}

		std::vector<number> vShape(type::nsh);
		std::vector<MathVector<dim> > vGrad(type::nsh);

		space.shapes(&vShape[0], pos);
		space.grads(&vGrad[0], pos);

		number sum = 0, sum2 = 0;
		MathVector<dim> grad, grad2;
		VecSet(grad, 0.0); VecSet(grad2, 0.0);
		for(size_t sh = 0; sh < space.num_sh(); ++sh)
		{
//			UG_LOG(sh << " Adding: int="<<vInter[sh] << ", grad=" << vGrad[sh] << "\n");
//			UG_LOG(sh << " Adding: int="<<vInter[sh] << ", shape=" << vShape[sh] << "\n");
			sum += vInter[sh] * space.shape(sh, pos);
			sum2 += vInter[sh] * vShape[sh];
			VecScaleAppend(grad, vInter[sh], space.grad(sh, pos));
			VecScaleAppend(grad2, vInter[sh], vGrad[sh]);
		}

		number diff = fabs(sum - func1<dim, p>(pos));
		maxDiff = std::max(maxDiff, diff);
		if(diff > Small)
		{
			UG_LOG("ERROR: Sum of Shapes not ok but differ " << diff <<".\n");
			return false;
		}

		if(fabs(sum - sum2) > Small)
		{
			UG_LOG("ERROR: Vector sum and sum not equal. \n");
			return false;
		}

		MathVector<dim> vecDiff;
		number diffGrad;

		VecSubtract(vecDiff, grad2, Dfunc1<dim, p>(pos));
		diffGrad = VecTwoNorm(vecDiff);
		if(diffGrad > Small)
		{
			MathVector<dim> tmpGrad = Dfunc1<dim, p>(pos);
			UG_LOG("ERROR: Sum of Grad (from vector) not ok but differ " << diffGrad <<".\n");
			UG_LOG("       Grad is " << grad2 << ", exact is " << tmpGrad << ".\n");
			return false;
		}

		VecSubtract(vecDiff, grad, Dfunc1<dim, p>(pos));
		diffGrad = VecTwoNorm(vecDiff);
		if(diffGrad > Small)
		{
			MathVector<dim> tmpGrad = Dfunc1<dim, p>(pos);
			UG_LOG("ERROR: Sum of Grad not ok but differ " << diffGrad <<":\n");
			UG_LOG("       Grad is " << grad << ", exact is " << tmpGrad << ".\n");
			return false;
		}

	}
	UG_LOG("PASSED.\n");

	UG_LOG("\n");
	return true;
}

template <int p>
bool TestLagrangeSpacesOrder(number& maxDiff)
{
	bool bRet = true;

	UG_LOG(" -----  EDGE  ----- \n");
	bRet &= TestLagrangeSpacesElem<ReferenceEdge, p>(maxDiff);

	UG_LOG(" -----  QUADRILATERAL  ----- \n");
	bRet &= TestLagrangeSpacesElem<ReferenceQuadrilateral, p>(maxDiff);

	UG_LOG(" -----  HEXAHEDRON  ----- \n");
	bRet &= TestLagrangeSpacesElem<ReferenceHexahedron, p>(maxDiff);

	UG_LOG(" -----  TRIANGLE  ----- \n");
//	bRet &= TestLagrangeSpacesElem<ReferenceTriangle, p>(maxDiff);

	UG_LOG(" -----  TETRAHEDRON  ----- \n");
//	bRet &= TestLagrangeSpacesElem<ReferenceTetrahedron, p>(maxDiff);

	UG_LOG(" -----  PRISM  ----- \n");
//	bRet &= TestLagrangeSpacesElem<ReferencePrism, p>(maxDiff);

	UG_LOG(" -----  PYRAMID  ----- \n");
	bRet &= TestLagrangeSpacesElem<ReferencePyramid, p>(maxDiff);

	UG_LOG(" #### COMMON RESULT for order " << p <<
	       ": Test " << (bRet?"SUCCEDED":"FAILED") << "\n");

	return bRet;
}

bool TestLagrangeSpaces()
{
	bool bRet = true;

	number maxDiff = 0.0;
	bRet &= TestLagrangeSpacesOrder<1>(maxDiff);
/*	bRet &= TestLagrangeSpacesOrder<2>(maxDiff);
	bRet &= TestLagrangeSpacesOrder<3>(maxDiff);
	bRet &= TestLagrangeSpacesOrder<4>(maxDiff);
	bRet &= TestLagrangeSpacesOrder<5>(maxDiff);
	bRet &= TestLagrangeSpacesOrder<6>(maxDiff);
	bRet &= TestLagrangeSpacesOrder<7>(maxDiff);
	bRet &= TestLagrangeSpacesOrder<8>(maxDiff);
*/
	UG_LOG(" #### COMMON RESULT"
	       ": Test " << (bRet?"SUCCEDED":"FAILED") << "\n");

	UG_LOG("Max Difference from exact value: " << maxDiff << "\n");

	return bRet;
}
/*
void SetDebugLevel(const char* t, int level)
{
	std::string tag(t);

	if(tag == "LIB_DISC_MULTIGRID")
	{
		UG_SET_DEBUG_LEVEL(LIB_DISC_MULTIGRID, level);
	}
	else
	{
		UG_LOG("In SetDebugLevel: Tag '" << tag << "'not found.\n");
	}
}*/

bool RegisterStaticLibDiscretizationInterface(Registry& reg, const char* parentGroup)
{
	try
	{
	//	get group string
		std::string grp = parentGroup; grp.append("/Discretization");

	//	FunctionGroup
		{
			reg.add_class_<FunctionGroup>("FunctionGroup", grp.c_str())
				.add_constructor()
				.add_method("clear", &FunctionGroup::clear)
				.add_method("set_function_pattern", &FunctionGroup::set_function_pattern)
				.add_method("add_function", (bool (FunctionGroup::*)(const char*))&FunctionGroup::add);
		}

	//	FunctionPattern
		{
			typedef FunctionPattern T;
			reg.add_class_<T>("FunctionPattern", grp.c_str())
				.add_method("clear", &T::clear)
				.add_method("add_fct_on_subset", (bool (T::*)(const char*, const char*, int, const char*))&T::add_fct_on_subset)
				.add_method("add_fct", (bool (T::*)(const char*, const char*, int))&T::add_fct);
		}

	//  Debug function
		//reg.add_function("SetDebugLevel", &SetDebugLevel, grp.c_str());

	//	TestLagrangeSpaces
		reg.add_function("TestLagrangeSpaces", &TestLagrangeSpaces, grp.c_str());

#ifdef UG_PARALLEL
	//	IDomainDecompositionInfo, StandardDomainDecompositionInfo
		{
		typedef pcl::IDomainDecompositionInfo Tbase;
		reg.add_class_<Tbase>("IDomainDecompositionInfo", grp.c_str());
		typedef pcl::StandardDomainDecompositionInfo T;
		reg.add_class_<T, Tbase>("StandardDomainDecompositionInfo", grp.c_str())
			.add_constructor()
			.add_method("map_proc_id_to_subdomain_id", &T::map_proc_id_to_subdomain_id)
			.add_method("set_num_subdomains",          &T::set_num_subdomains)
			.add_method("get_num_subdomains",          &T::get_num_subdomains);
		}
#endif

	}
	catch(UG_REGISTRY_ERROR_RegistrationFailed ex)
	{
		UG_LOG("### ERROR in RegisterLibDiscretizationInterface: "
				"Registration failed (using name " << ex.name << ").\n");
		return false;
	}

	return true;
}



} // end namespace bridge
} // end namespace ug
