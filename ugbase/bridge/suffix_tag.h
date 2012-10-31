/*
 * suffix_tag.h
 *
 *  Created on: 24.05.2012
 *      Author: andreasvogel
 */

#ifndef __H__UG_BRIDGE__SUFFIX_TAG__
#define __H__UG_BRIDGE__SUFFIX_TAG__


namespace ug{
namespace bridge{

////////////////////////////////////////////////////////////////////////////////
//	Dimension - Suffix and Tag
////////////////////////////////////////////////////////////////////////////////

/// returns the dim-suffix for a domain (e.g. "3d")
template <int dim>
std::string GetDimensionSuffix()
{
//	the dimension suffix
	std::stringstream ss; ss << dim << "d";
	return ss.str();
}

/// returns the dim-tag for a domain (e.g. "dim=3d")
template <int dim>
std::string GetDimensionTag()
{
	return std::string("dim=").append(GetDimensionSuffix<dim>()).append(";");
}

/// returns dim tag at runtime (e.g. "dim=3d")
inline std::string GetDimensionTag(int dim)
{
//	the dimension tag
	std::stringstream ss; ss << "dim=" << dim << "d;";
	return ss.str();
}

////////////////////////////////////////////////////////////////////////////////
//	Domain - Suffix and Tag
////////////////////////////////////////////////////////////////////////////////

/// returns the dim-suffix for a domain (e.g. "3d")
template <typename TDomain>
std::string GetDomainSuffix(){return GetDimensionSuffix<TDomain::dim>();}

/// returns the dim-tag for a domain (e.g. "dim=3d")
template <typename TDomain>
std::string GetDomainTag() {return GetDimensionTag<TDomain::dim>();}

////////////////////////////////////////////////////////////////////////////////
//	Algebra - Suffix and Tag
////////////////////////////////////////////////////////////////////////////////

/// returns the algebra-suffix (e.g. "CPU3", "CPUFlex")
inline std::string GetAlgebraSuffix(const AlgebraType& algType)
{
//	the algebra suffix
	std::stringstream ss;

//	add type
	if(algType.type() == AlgebraType::CPU) ss << "CPU";
	else if(algType.type() == AlgebraType::CRS) ss << "CRS";
	else UG_THROW("Unknown algebra type.");

//	add blocktype
	if(algType.blocksize() == AlgebraType::VariableBlockSize) ss << "Variable";
	else ss << algType.blocksize();

	return ss.str();
}

/// returns the algebra-suffix (e.g. "CPU3", "CPUFlex")
template <typename TAlgebra>
std::string GetAlgebraSuffix()
{
	return GetAlgebraSuffix(TAlgebra::get_type());
}

/// returns the algebra-suffix (e.g. "alg=CPU3", "alg=CPUVariable")
inline std::string GetAlgebraTag(const AlgebraType& algType)
{
//	the algebra suffix
	std::stringstream ss; ss << "alg=";

//	add type
	if(algType.type() == AlgebraType::CPU) ss << "CPU";
	else if(algType.type() == AlgebraType::CRS) ss << "CRS";
	else UG_THROW("Unknown algebra type.");

//	add blocktype
	if(algType.blocksize() == AlgebraType::VariableBlockSize) ss << "Variable;";
	else ss << algType.blocksize() << ";";

	return ss.str();
}


/// returns the algebra-suffix (e.g. "alg=CPU3", "alg=CPUVariable")
template <typename TAlgebra>
std::string GetAlgebraTag()
{
	return GetAlgebraTag(TAlgebra::get_type());
}

////////////////////////////////////////////////////////////////////////////////
//	Dimension + Algebra - Suffix and Tag
////////////////////////////////////////////////////////////////////////////////

/// returns the algebra-dim-suffix for a domain (e.g. "3dCPU1")
template <int dim, typename TAlgebra>
std::string GetDimensionAlgebraSuffix()
{
	std::string dimAlgSuffix = GetDimensionSuffix<dim>();
	dimAlgSuffix.append(GetAlgebraSuffix<TAlgebra>());
	return dimAlgSuffix;
}

/// returns the dim-tag for a domain (e.g. "dim=3d;alg=CPU1;")
template <int dim, typename TAlgebra>
std::string GetDimensionAlgebraTag()
{
	std::string dimAlgTag = GetDimensionTag<dim>();
	dimAlgTag.append(GetAlgebraTag<TAlgebra>());
	return dimAlgTag;
}
/// returns dim tag at runtime (e.g. "dim=3d;alg=CPU1;")
inline std::string GetDimensionAlgebraTag(int dim, const AlgebraType& algType)
{
	std::string dimAlgTag = GetDimensionTag(dim);
	dimAlgTag.append(GetAlgebraTag(algType));
	return dimAlgTag;
}

////////////////////////////////////////////////////////////////////////////////
//	Domain + Algebra - Suffix and Tag
////////////////////////////////////////////////////////////////////////////////

/// returns the dim-suffix for a domain (e.g. "3dCPU1")
template <typename TDomain, typename TAlgebra>
std::string GetDomainAlgebraSuffix(){return GetDimensionAlgebraSuffix<TDomain::dim, TAlgebra>();}

/// returns the dim-tag for a domain (e.g. "dim=3d;alg=CPU1;")
template <typename TDomain, typename TAlgebra>
std::string GetDomainAlgebraTag(){return GetDimensionAlgebraTag<TDomain::dim, TAlgebra>();}

}//	end bridge
}//	end ug

#endif /* __H__UG_BRIDGE__SUFFIX_TAG__ */
