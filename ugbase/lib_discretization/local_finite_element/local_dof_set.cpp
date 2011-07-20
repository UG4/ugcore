/*
 * local_dof_set.cpp
 *
 *  Created on: 18.07.2011
 *      Author: andreasvogel
 */

#include "local_dof_set.h"
#include "lagrange/lagrange_local_dof.h"

namespace ug{

LocalDoFSetProvider::LocalDoFSetProvider()
{
	static bool init = false;

	if(!init)
	{
	//	remember initialization
		init = true;

	//	register all element types that allow higher orders
		if(!init_standard_sets<ReferenceVertex>())
			throw(UGFatalError("Cannot register standard Edge local dof sets."));
		if(!init_standard_sets<ReferenceEdge>())
			throw(UGFatalError("Cannot register standard Edge local dof sets."));
		if(!init_standard_sets<ReferenceTriangle>())
			throw(UGFatalError("Cannot register standard Triangle local dof sets."));
		if(!init_standard_sets<ReferenceQuadrilateral>())
			throw(UGFatalError("Cannot register standard Quadrilateral local dof sets."));
		if(!init_standard_sets<ReferenceTetrahedron>())
			throw(UGFatalError("Cannot register standard Tetrahedron local dof sets."));
		if(!init_standard_sets<ReferencePrism>())
			throw(UGFatalError("Cannot register standard Prism local dof sets."));
		if(!init_standard_sets<ReferenceHexahedron>())
			throw(UGFatalError("Cannot register standard Hexahedron local dof sets."));
		if(!init_standard_sets<ReferencePyramid>())
			throw(UGFatalError("Cannot register standard Pyramid local dof sets."));
	}
}

template <typename TRefElem>
bool LocalDoFSetProvider::init_standard_sets()
{
	static ILocalDoFSetWrapper<LagrangeLDS<TRefElem, 1> > sLagrangeP1;
	if(!register_set(LFEID(LFEID::LAGRANGE, 1), sLagrangeP1)) return false;

	static ILocalDoFSetWrapper<LagrangeLDS<TRefElem, 2> > sLagrangeP2;
	if(!register_set(LFEID(LFEID::LAGRANGE, 2), sLagrangeP2)) return false;

	static ILocalDoFSetWrapper<LagrangeLDS<TRefElem, 3> > sLagrangeP3;
	if(!register_set(LFEID(LFEID::LAGRANGE, 3), sLagrangeP3)) return false;

//	done
	return true;
}


const ILocalDoFSet& LocalDoFSetProvider::get(ReferenceObjectID roid, LFEID id)
{
//	init provider and search for identifier
	RoidMap::const_iterator iter = inst().m_mRoidDoFSet.find(id);

//	if not found
	if(iter == m_mRoidDoFSet.end())
	{
		UG_LOG("ERROR in 'LocalDoFSetProvider::get': "
				"Unknown LocalDoFSet for type "<<id<<" requested.\n");
		throw(UGFatalError("LocalDoFSet for Finite Element Space unknown"));
	}

//	get vector
	const std::vector<const ILocalDoFSet*>& vBase = iter->second;

//	check that dof set registered
	if(vBase[roid] == NULL)
	{
		UG_LOG("ERROR in 'LocalDoFSetProvider::get': "
				"Unknown LocalDoFSet for type "<<id<<" requested for"
				" Reference Element type " <<roid<<".\n");
		throw(UGFatalError("Trial Space type unknown"));
	}

//	return dof set
	return *(vBase[roid]);
}

const CommonLocalDoFSet& LocalDoFSetProvider::get(int dim, LFEID id)
{
//	init provider and search for identifier
	CommonMap::const_iterator iter = inst().m_mCommonDoFSet.find(id);

//	if not found
	if(iter == m_mCommonDoFSet.end())
	{
		UG_LOG("ERROR in 'LocalDoFSetProvider::get': "
				"Unknown LocalDoFSet for type "<<id<<" and dim "<<dim<<" requested.\n");
		throw(UGFatalError("LocalDoFSet for Finite Element Space unknown"));
	}

//	get dimension
	const std::vector<CommonLocalDoFSet> vCommonDoFSet = iter->second;

//	return the common set
	return vCommonDoFSet.at(dim);
}


bool LocalDoFSetProvider::register_set(LFEID id, const ILocalDoFSet& set)
{
//	reference object id
	const ReferenceObjectID roid = set.roid();

//	get vector of types
	std::vector<const ILocalDoFSet*>& vBase = m_mRoidDoFSet[id];

//	resize vector
	vBase.resize(NUM_REFERENCE_OBJECTS, NULL);

//	check that no space has been previously registered to this place
	if(vBase[roid])
	{
		UG_LOG("ERROR in 'LocalDoFSetProvider::register_set()': "
				"LocalDoFSet already registered for type: "<<id<<" and "
				" Reference element type "<<roid<<".\n");
		return false;
	}

//	if ok, add
	vBase[roid] = &set;

//	get common dof set for the dimension
	std::vector<CommonLocalDoFSet>& vCommonSet = m_mCommonDoFSet[id];

//	resize
	vCommonSet.resize(4);

//	add this local dof set
	if(!vCommonSet[set.dim()].add(set))
	{
	//	write error message
		UG_LOG("ERROR in 'LocalDoFSetProvider::register_set()': "
				"Cannot build CommonLocalDoFSet for type: "<<id<<" when adding "
				" Reference element type "<<roid<<".\n");

		UG_LOG("CommonLocalDoFSet is:\n" << vCommonSet[set.dim()]);
		UG_LOG("LocalDoFSet is:\n" << set);

	//	remove set
		m_mCommonDoFSet.erase(id);

	//	return error flag
		return false;
	}

//	done
	return true;
}

std::map<LFEID, std::vector<const ILocalDoFSet*> >
LocalDoFSetProvider::m_mRoidDoFSet = std::map<LFEID, std::vector<const ILocalDoFSet*> >();

std::map<LFEID, std::vector<CommonLocalDoFSet> >
LocalDoFSetProvider::m_mCommonDoFSet = std::map<LFEID, std::vector<CommonLocalDoFSet> >();

////////////////////////////////////////////////////////////////////////////////
// 	CommonLocalDoFSet
////////////////////////////////////////////////////////////////////////////////

void CommonLocalDoFSet::clear()
{
//	set all dofs to not specified
	for(int i = 0; i < NUM_REFERENCE_OBJECTS; ++i)
		m_vNumDoF[i] = NOT_SPECIFIED;
}

///	add a local dof set to the intersection
bool CommonLocalDoFSet::add(const ILocalDoFSet& set)
{
	for(int i = 0; i < NUM_REFERENCE_OBJECTS; ++i)
	{
	//	get roid
		const ReferenceObjectID roid = (ReferenceObjectID) i;

	//	skip if reference element dim of set lower than ref element
		if(set.dim() < ReferenceElementDimension(roid))
			continue;

	//	skip if same dimension but different roid
		if(set.dim() == ReferenceElementDimension(roid))
			if(set.roid() != roid)
				continue;

	//	check if space specifies for roid
		if(set.num_dof(roid) == NOT_SPECIFIED) continue;

	//	check if already value set and iff the same
		if(m_vNumDoF[i] != NOT_SPECIFIED)
			if(m_vNumDoF[i] != set.num_dof(roid))
			{
				UG_LOG("ERROR in 'LocalDoFSetIntersection::add': "
						" Values does not match ("<<m_vNumDoF[i]<<" <-> "
						<< set.num_dof(roid)<<").\n");
				return false;
			}

	//	set value if not already set
		m_vNumDoF[i] = set.num_dof(roid);
	}

//	done
	return true;
}

/// writes to the output stream
std::ostream& operator<<(std::ostream& out,	const CommonLocalDoFSet& v)
{
	for(int i = 0; i < NUM_REFERENCE_OBJECTS; ++i)
	{
		ReferenceObjectID roid = (ReferenceObjectID) i;

		out << roid << ": \t\t" << v.num_dof(roid) << "\n";
	}
	return out;
}

/// writes to the output stream
std::ostream& operator<<(std::ostream& out,	const ILocalDoFSet& v)
{
	for(int i = 0; i < NUM_REFERENCE_OBJECTS; ++i)
	{
		ReferenceObjectID roid = (ReferenceObjectID) i;

		out << roid << ": \t\t" << v.num_dof(roid) << "\n";
	}
	return out;
}

} // end namespace ug
