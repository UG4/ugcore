/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#ifndef __H__UG_BRIDGE__SUFFIX_TAG__
#define __H__UG_BRIDGE__SUFFIX_TAG__


namespace ug{
namespace bridge{

/// \addtogroup bridge
/// \{

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
template<typename TAlgebraTypeType>
inline std::string GetAlgebraSuffix(const TAlgebraTypeType& algType)
{
//	the algebra suffix
	std::stringstream ss;

//	add type
	if(algType.type() == TAlgebraTypeType::CPU) ss << "CPU";
	else if(algType.type() == TAlgebraTypeType::GPU) ss << "GPU";
	else UG_THROW("Unknown algebra type.");

//	add blocktype
	if(algType.blocksize() == TAlgebraTypeType::VariableBlockSize) ss << "Variable";
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
template<typename TAlgebraTypeType>
inline std::string GetAlgebraTag(const TAlgebraTypeType& algType)
{
//	the algebra suffix
	std::stringstream ss; ss << "alg=";

//	add type
	if(algType.type() == TAlgebraTypeType::CPU) ss << "CPU";
	else if(algType.type() == TAlgebraTypeType::GPU) ss << "GPU";
	else UG_THROW("Unknown algebra type.");

//	add blocktype
	if(algType.blocksize() == TAlgebraTypeType::VariableBlockSize) ss << "Variable;";
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
template<typename TAlgebraTypeType>
inline std::string GetDimensionAlgebraTag(int dim, const TAlgebraTypeType& algType)
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

// end group bridge
/// \}

}//	end bridge
}//	end ug

#endif /* __H__UG_BRIDGE__SUFFIX_TAG__ */
