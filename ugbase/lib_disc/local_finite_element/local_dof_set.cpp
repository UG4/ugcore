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
{};

template <typename TRefElem>
void LocalDoFSetProvider::create_set(const LFEID& id)
{
	switch(id.type()){
		case LFEID::LAGRANGE:
			register_set(id, ConstSmartPtr<LocalDoFSet>(new LagrangeLDS<TRefElem>(id.order())));
			return;
		case LFEID::PIECEWISE_CONSTANT:
			register_set(id, ConstSmartPtr<LocalDoFSet>(new PiecewiseConstantLDS<TRefElem>));
			return;
		case LFEID::CROUZEIX_RAVIART:
			register_set(id, ConstSmartPtr<LocalDoFSet>(new CrouzeixRaviartLDS<TRefElem>));
			return;
		case LFEID::MINI:
			register_set(id, ConstSmartPtr<LocalDoFSet>(new MiniBubbleLDS<TRefElem>));
			return;
		case LFEID::NEDELEC:
			register_set(id, ConstSmartPtr<LocalDoFSet>(new NedelecLDS<TRefElem>));
			return;
		default: return;
	}
}

void LocalDoFSetProvider::
create_set(ReferenceObjectID roid, const LFEID& id)
{
	try{
	//	switch type
		switch(roid)
		{
			case ROID_VERTEX:		create_set<ReferenceVertex>(id); return;
			case ROID_EDGE:			create_set<ReferenceEdge>(id); return;
			case ROID_TRIANGLE:		create_set<ReferenceTriangle>(id); return;
			case ROID_QUADRILATERAL:create_set<ReferenceQuadrilateral>(id); return;
			case ROID_TETRAHEDRON:	create_set<ReferenceTetrahedron>(id); return;
			case ROID_PRISM:		create_set<ReferencePrism>(id); return;
			case ROID_PYRAMID:		create_set<ReferencePyramid>(id); return;
			case ROID_HEXAHEDRON:	create_set<ReferenceHexahedron>(id); return;
			default: return;
		}
	}
	UG_CATCH_THROW("LocalDoFSetProvider: Creation of set "<<id<<
					" for "<<roid<<" failed.");
}
void LocalDoFSetProvider::
create_set(const LFEID& id)
{
	for(int roid = 0; roid < NUM_REFERENCE_OBJECTS; ++roid)
		create_set((ReferenceObjectID)roid, id);
}

const LocalDoFSet& LocalDoFSetProvider::
get(ReferenceObjectID roid, const LFEID& id, bool bCreate)
{
//	init provider and search for identifier
	RoidMap::const_iterator iter = inst().m_mRoidDoFSet.find(id);

//	if not found
	if(iter == m_mRoidDoFSet.end()){
		if(bCreate){
			create_set(id);
			return get(roid, id, false);
		}
		UG_THROW("LocalDoFSetProvider: Cannot create LocalDoFSet for finite "
				"element type "<<id<<" and "<<roid);
	}

//	get vector
	const std::vector<ConstSmartPtr<LocalDoFSet> >& vBase = iter->second;

//	check that dof set registered
	if(vBase[roid].invalid()){
		if(bCreate){
			create_set(id);
			return get(roid, id, false);
		}
		UG_THROW("LocalDoFSetProvider: Cannot create LocalDoFSet for finite "
				"element type "<<id<<" and "<<roid);
	}

//	return dof set
	return *(vBase[roid]);
}

const CommonLocalDoFSet& LocalDoFSetProvider::
get(const LFEID& id, bool bCreate)
{
//	init provider and search for identifier
	CommonMap::const_iterator iter = inst().m_mCommonDoFSet.find(id);

//	if not found
	if(iter == m_mCommonDoFSet.end()){
		if(bCreate)	{
			create_set(id);
			return get(id, false);
		}
		UG_THROW("LocalDoFSetProvider: Cannot create CommonLocalDoFSet for "<<id);
	}

//	return the common set
	return iter->second;
}


void LocalDoFSetProvider::register_set(const LFEID& id, ConstSmartPtr<LocalDoFSet> set)
{
//	reference object id
	const ReferenceObjectID roid = set->roid();

//	get vector of types
	std::vector<ConstSmartPtr<LocalDoFSet> >& vBase = m_mRoidDoFSet[id];

//	resize vector
	vBase.resize(NUM_REFERENCE_OBJECTS, NULL);

//	if ok, add
	vBase[roid] = set;

//	for creation of CommonLocalDoFSet: skip if not the dimension of the space
	if(set->dim() != id.dim()) return;

//	add this local dof set
	try{
		m_mCommonDoFSet[id].add(*set);
	}
	catch(UGError& err)
	{
	//	write error message
		std::stringstream ss;
		ss<<"LocalDoFSetProvider::register_set(): "
				"Cannot build CommonLocalDoFSet for type: "<<id<<" when adding "
				" Reference element type "<<roid<<".\n"<<
				"CommonLocalDoFSet is:\n" << m_mCommonDoFSet[id]<<
				"LocalDoFSet is:\n" << *set;
		err.push_msg(ss.str(),__FILE__,__LINE__);

	//	remove set
		m_mCommonDoFSet.erase(id);

		throw(err);
	}
}

std::map<LFEID, std::vector<ConstSmartPtr<LocalDoFSet> > >
LocalDoFSetProvider::m_mRoidDoFSet =
		std::map<LFEID, std::vector<ConstSmartPtr<LocalDoFSet> > >();

std::map<LFEID, CommonLocalDoFSet>
LocalDoFSetProvider::m_mCommonDoFSet = std::map<LFEID, CommonLocalDoFSet>();

////////////////////////////////////////////////////////////////////////////////
// 	LocalDoFSet
////////////////////////////////////////////////////////////////////////////////

int LocalDoFSet::dim() const {
	return ReferenceElementDimension(roid());
}

size_t LocalDoFSet::num_sh() const
{
	size_t sum = 0;
	for(int i = 0; i < NUM_REFERENCE_OBJECTS; ++i)
		sum += num_dof((ReferenceObjectID)i);
	return sum;
}

size_t LocalDoFSet::num_dof(int d, size_t id) const
{
	static const ReferenceElement& rRefElem = ReferenceElementProvider::get(roid());
	return num_dof(rRefElem.roid(d, id));
}

////////////////////////////////////////////////////////////////////////////////
// 	DimLocalDoFSet
////////////////////////////////////////////////////////////////////////////////

template <int TDim>
DimLocalDoFSet<TDim>::DimLocalDoFSet(){
}

template class DimLocalDoFSet<1>;
template class DimLocalDoFSet<2>;
template class DimLocalDoFSet<3>;

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
void CommonLocalDoFSet::add(const LocalDoFSet& set)
{
	for(int i = 0; i < NUM_REFERENCE_OBJECTS; ++i)
	{
	//	get roid
		const ReferenceObjectID roid = (ReferenceObjectID) i;

	//	do not override values of same dimension
		if(ReferenceElementDimension(roid) == set.dim()
			&& roid != set.roid()) continue;

	//	check if already value set and iff the same
		if(m_vNumDoF[i] != NOT_SPECIFIED)
			if(m_vNumDoF[i] != (int)set.num_dof(roid))
				UG_THROW("LocalDoFSetIntersection::add: Adding DoF-Spezification "
						"for "<<roid<<" as Subelement of Space for "<<set.roid()<<
						": Values does not match ("<<m_vNumDoF[i]<<" <-> "
						<< set.num_dof(roid)<<").");

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
std::ostream& operator<<(std::ostream& out,	const LocalDoFSet& v)
{
	for(int i = 0; i < NUM_REFERENCE_OBJECTS; ++i)
	{
		ReferenceObjectID roid = (ReferenceObjectID) i;

		out << std::setw(14) << roid << ":   " << v.num_dof(roid) << "\n";
	}
	return out;
}

} // end namespace ug
