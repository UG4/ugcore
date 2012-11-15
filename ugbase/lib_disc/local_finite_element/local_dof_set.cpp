/*
 * local_dof_set.cpp
 *
 *  Created on: 18.07.2011
 *      Author: andreasvogel
 */

#include "local_dof_set.h"
#include "lagrange/lagrange_local_dof.h"
#include "crouzeix-raviart/crouzeix_raviart_local_dof.h"
#include "piecewise_constant/piecewise_constant_local_dof.h"
#include "mini/mini_local_dof.h"
#include "nedelec/nedelec_local_dof.h"
#include "lib_disc/reference_element/reference_element_util.h"

namespace ug{

LocalDoFSetProvider::
~LocalDoFSetProvider()
{
//	free all created sets
	for(size_t i = 0; i < m_vCreated.size(); ++i)
		delete m_vCreated[i];
};

template <typename TRefElem>
void LocalDoFSetProvider::create_lagrange_set(size_t order)
{
//	create lagrange set
	LagrangeLDS<TRefElem>* setLagrange = new LagrangeLDS<TRefElem>(order);

//	remember created set for delete in destructor
	m_vCreated.push_back(setLagrange);

//	register the set
	try{
		register_set(LFEID(LFEID::LAGRANGE, order), *setLagrange);
	}
	UG_CATCH_THROW("Unable to register LagrangeLDS");
}

void LocalDoFSetProvider::create_lagrange_sets(size_t order)
{
	create_lagrange_set<ReferenceVertex>(order);
	create_lagrange_set<ReferenceEdge>(order);
	create_lagrange_set<ReferenceTriangle>(order);
	create_lagrange_set<ReferenceQuadrilateral>(order);
	create_lagrange_set<ReferenceTetrahedron>(order);
	create_lagrange_set<ReferencePyramid>(order);
	create_lagrange_set<ReferencePrism>(order);
	create_lagrange_set<ReferenceHexahedron>(order);
}

template <typename TRefElem>
void LocalDoFSetProvider::create_crouzeix_raviart_sets()
{
//	create Crouzeix-Raviart set
	CrouzeixRaviartLDS<TRefElem>* setCrouzeixRaviart = new CrouzeixRaviartLDS<TRefElem>();

//	remember created set for delete in destructor
	m_vCreated.push_back(setCrouzeixRaviart);

//	register the set
	try{
		register_set(LFEID(LFEID::CROUZEIX_RAVIART, 1), *setCrouzeixRaviart);
	}
	UG_CATCH_THROW("Unable to register CrouzeixRaviartLDS");
}

void LocalDoFSetProvider::create_crouzeix_raviart_sets()
{
	create_crouzeix_raviart_sets<ReferenceVertex>();
	create_crouzeix_raviart_sets<ReferenceEdge>();
	create_crouzeix_raviart_sets<ReferenceTriangle>();
	create_crouzeix_raviart_sets<ReferenceQuadrilateral>();
	create_crouzeix_raviart_sets<ReferenceTetrahedron>();
	create_crouzeix_raviart_sets<ReferencePyramid>();
	create_crouzeix_raviart_sets<ReferencePrism>();
	create_crouzeix_raviart_sets<ReferenceHexahedron>();
}


template <typename TRefElem>
void LocalDoFSetProvider::create_piecewise_constant_sets()
{
//	create piecewise constant set
	PiecewiseConstantLDS<TRefElem>* setPiecewiseConstant = new PiecewiseConstantLDS<TRefElem>();

//	remember created set for delete in destructor
	m_vCreated.push_back(setPiecewiseConstant);

//	register the set
	try{
		register_set(LFEID(LFEID::PIECEWISE_CONSTANT, 0), *setPiecewiseConstant);
	}
	UG_CATCH_THROW("Unable to register PiecewiseConstantLDS");
}

void LocalDoFSetProvider::create_piecewise_constant_sets()
{
	create_piecewise_constant_sets<ReferenceVertex>();
	create_piecewise_constant_sets<ReferenceEdge>();
	create_piecewise_constant_sets<ReferenceTriangle>();
	create_piecewise_constant_sets<ReferenceQuadrilateral>();
	create_piecewise_constant_sets<ReferenceTetrahedron>();
	create_piecewise_constant_sets<ReferencePyramid>();
	create_piecewise_constant_sets<ReferencePrism>();
	create_piecewise_constant_sets<ReferenceHexahedron>();
}

template <typename TRefElem>
void LocalDoFSetProvider::create_mini_bubble_sets()
{
//	create piecewise constant set
	MiniBubbleLDS<TRefElem>* setMiniBubble = new MiniBubbleLDS<TRefElem>();

//	remember created set for delete in destructor
	m_vCreated.push_back(setMiniBubble);

//	register the set
	try{
		register_set(LFEID(LFEID::MINI, 1), *setMiniBubble);
	}
	UG_CATCH_THROW("Unable to register MiniBubbleLDS");
}

void LocalDoFSetProvider::create_mini_bubble_sets()
{
	create_mini_bubble_sets<ReferenceVertex>();
	create_mini_bubble_sets<ReferenceEdge>();
	create_mini_bubble_sets<ReferenceTriangle>();
	create_mini_bubble_sets<ReferenceQuadrilateral>();

	create_mini_bubble_sets<ReferenceTetrahedron>();
	create_mini_bubble_sets<ReferencePyramid>();
	create_mini_bubble_sets<ReferencePrism>();
	create_mini_bubble_sets<ReferenceHexahedron>();

}

template <typename TRefElem>
void LocalDoFSetProvider::create_nedelec_sets()
{
//	create nedelec element set
	NedelecLDS<TRefElem>* setNedelec = new NedelecLDS<TRefElem>();

//	remember created set for delete in destructor
	m_vCreated.push_back(setNedelec);

//	register the set
	try{
		register_set(LFEID(LFEID::PIECEWISE_CONSTANT, 0), *setNedelec);
	}
	UG_CATCH_THROW("Unable to register NedelecLDS");
}

void LocalDoFSetProvider::create_nedelec_sets()
{
	create_nedelec_sets<ReferenceVertex>();
	create_nedelec_sets<ReferenceEdge>();
	create_nedelec_sets<ReferenceTriangle>();
	create_nedelec_sets<ReferenceQuadrilateral>();
	create_nedelec_sets<ReferenceTetrahedron>();
	create_nedelec_sets<ReferencePyramid>();
	create_nedelec_sets<ReferencePrism>();
	create_nedelec_sets<ReferenceHexahedron>();
}

void LocalDoFSetProvider::create_set(const LFEID& id)
{
	if(id.type() == LFEID::LAGRANGE)
	{
		create_lagrange_sets(id.order());
	}
	if(id.type() == LFEID::CROUZEIX_RAVIART)
	{
		create_crouzeix_raviart_sets();
	}
	if(id.type() == LFEID::PIECEWISE_CONSTANT)
	{
		create_piecewise_constant_sets();
	}
	if(id.type() == LFEID::MINI)
	{
		create_mini_bubble_sets();
	}

}


const ILocalDoFSet& LocalDoFSetProvider::get(ReferenceObjectID roid, LFEID id, bool bCreate)
{
//	init provider and search for identifier
	RoidMap::const_iterator iter = inst().m_mRoidDoFSet.find(id);

//	if not found
	if(iter == m_mRoidDoFSet.end())
	{
		if(bCreate)
		{
			create_set(id);
			return get(roid, id, false);
		}

		UG_LOG("ERROR in 'LocalDoFSetProvider::get': "
				"Unknown LocalDoFSet for type "<<id<<" requested.\n");
		throw(UGError("LocalDoFSet for Finite Element Space unknown"));
	}

//	get vector
	const std::vector<const ILocalDoFSet*>& vBase = iter->second;

//	check that dof set registered
	if(vBase[roid] == NULL)
	{
		if(bCreate)
		{
			create_set(id);
			return get(roid, id, false);
		}

		UG_LOG("ERROR in 'LocalDoFSetProvider::get': "
				"Unknown LocalDoFSet for type "<<id<<" requested for"
				" Reference Element type " <<roid<<".\n");
		throw(UGError("Trial Space type unknown"));
	}

//	return dof set
	return *(vBase[roid]);
}

const CommonLocalDoFSet& LocalDoFSetProvider::get(int dim, LFEID id, bool bCreate)
{
//	init provider and search for identifier
	CommonMap::const_iterator iter = inst().m_mCommonDoFSet.find(id);

//	if not found
	if(iter == m_mCommonDoFSet.end())
	{
		if(bCreate)
		{
			create_set(id);
			return get(dim, id, false);
		}

		UG_LOG("ERROR in 'LocalDoFSetProvider::get': "
				"Unknown LocalDoFSet for type "<<id<<" and dim "<<dim<<" requested.\n");
		throw(UGError("LocalDoFSet for Finite Element Space unknown"));
	}

//	get dimension
	const std::vector<CommonLocalDoFSet>& vCommonDoFSet = iter->second;

//	return the common set
	return vCommonDoFSet.at(dim);
}


void LocalDoFSetProvider::register_set(LFEID id, const ILocalDoFSet& set)
{
//	reference object id
	const ReferenceObjectID roid = set.roid();

//	get vector of types
	std::vector<const ILocalDoFSet*>& vBase = m_mRoidDoFSet[id];

//	resize vector
	vBase.resize(NUM_REFERENCE_OBJECTS, NULL);

//	check that no space has been previously registered to this place
	if(vBase[roid])
		UG_THROW("LocalDoFSetProvider::register_set(): "
				"LocalDoFSet already registered for type: "<<id<<" and "
				" Reference element type "<<roid<<".");

//	if ok, add
	vBase[roid] = &set;

//	get common dof set for the dimension
	std::vector<CommonLocalDoFSet>& vCommonSet = m_mCommonDoFSet[id];

//	resize
	vCommonSet.resize(4);

//	add this local dof set
	try{
		vCommonSet[set.dim()].add(set);
	}
	catch(UGError& err)
	{
	//	remove set
		m_mCommonDoFSet.erase(id);

	//	write error message
		std::stringstream ss;
		ss<<"LocalDoFSetProvider::register_set(): "
				"Cannot build CommonLocalDoFSet for type: "<<id<<" when adding "
				" Reference element type "<<roid<<".\n"<<
				"CommonLocalDoFSet is:\n" << vCommonSet[set.dim()]<<
				"LocalDoFSet is:\n" << set;
		err.push_msg(ss.str(),__FILE__,__LINE__);
		throw(err);
	}
}

std::map<LFEID, std::vector<const ILocalDoFSet*> >
LocalDoFSetProvider::m_mRoidDoFSet = std::map<LFEID, std::vector<const ILocalDoFSet*> >();

std::map<LFEID, std::vector<CommonLocalDoFSet> >
LocalDoFSetProvider::m_mCommonDoFSet = std::map<LFEID, std::vector<CommonLocalDoFSet> >();

std::vector<ILocalDoFSet*>
LocalDoFSetProvider::m_vCreated = std::vector<ILocalDoFSet*>();

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
void CommonLocalDoFSet::add(const ILocalDoFSet& set)
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
				UG_THROW("LocalDoFSetIntersection::add: "
						" Values does not match ("<<m_vNumDoF[i]<<" <-> "
						<< set.num_dof(roid)<<").");
			}

	//	set value if not already set
		m_vNumDoF[i] = set.num_dof(roid);
	}
}

/// writes to the output stream
std::ostream& operator<<(std::ostream& out,	const CommonLocalDoFSet& v)
{
	for(int i = 0; i < NUM_REFERENCE_OBJECTS; ++i)
	{
		ReferenceObjectID roid = (ReferenceObjectID) i;

		out << std::setw(14) << roid << ":   " << v.num_dof(roid) << "\n";
	}
	return out;
}

/// writes to the output stream
std::ostream& operator<<(std::ostream& out,	const ILocalDoFSet& v)
{
	for(int i = 0; i < NUM_REFERENCE_OBJECTS; ++i)
	{
		ReferenceObjectID roid = (ReferenceObjectID) i;

		out << std::setw(14) << roid << ":   " << v.num_dof(roid) << "\n";
	}
	return out;
}

} // end namespace ug
