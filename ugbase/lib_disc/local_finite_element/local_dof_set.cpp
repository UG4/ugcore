/*
 * local_dof_set.cpp
 *
 *  Created on: 18.07.2011
 *      Author: andreasvogel
 */

#include "local_dof_set.h"
#include "lib_disc/reference_element/reference_element_util.h"

namespace ug{

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

bool LocalDoFSet::operator==(const LocalDoFSet& v) const
{
	if(roid() != v.roid()) return false;
	if(dim() != v.dim()) return false;

	if(num_dof() != v.num_dof()) return false;

	for(int r = 0; r < NUM_REFERENCE_OBJECTS; ++r){
		const ReferenceObjectID roid = (ReferenceObjectID)r;
		if(num_dof(roid) != v.num_dof(roid)) return false;
	}

	for(size_t dof = 0; dof < num_dof(); ++dof)
		if(local_dof(dof) != v.local_dof(dof)) return false;

	return true;
}

////////////////////////////////////////////////////////////////////////////////
// 	DimLocalDoFSet
////////////////////////////////////////////////////////////////////////////////

template <int TDim>
DimLocalDoFSet<TDim>::DimLocalDoFSet(){
}

template <int TDim>
bool DimLocalDoFSet<TDim>::operator==(const DimLocalDoFSet<TDim>& v) const
{
	if(*dynamic_cast<const LocalDoFSet*>(this) != *dynamic_cast<const LocalDoFSet*>(&v))
		return false;

	if(exact_position_available() != v.exact_position_available()) return false;

	for(size_t dof = 0; dof < num_dof(); ++dof){
		MathVector<TDim> pos, vpos;
		position(dof, pos);
		v.position(dof, vpos);

		if(pos != vpos) return false;
	}

	return true;
}

template class DimLocalDoFSet<0>;
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

std::ostream& operator<<(std::ostream& out,	const LocalDoF& v)
{
	out <<"("<<v.dim()<<","<<v.id()<<","<<v.offset()<<")";
	return out;
}

std::ostream& operator<<(std::ostream& out,	const CommonLocalDoFSet& v)
{
	for(int i = 0; i < NUM_REFERENCE_OBJECTS; ++i)
	{
		ReferenceObjectID roid = (ReferenceObjectID) i;

		out << std::setw(14) << roid << ":   ";
		if(v.num_dof(roid) == CommonLocalDoFSet::NOT_SPECIFIED)
			out << "(not specified)";
		else
			out << v.num_dof(roid);
		out << "\n";
	}
	return out;
}

std::ostream& operator<<(std::ostream& out,	const LocalDoFSet& v)
{
	for(int i = 0; i < NUM_REFERENCE_OBJECTS; ++i)
	{
		ReferenceObjectID roid = (ReferenceObjectID) i;
		out << std::setw(14) << roid << ":   " << v.num_dof(roid) << "\n";
	}
	return out;
}

template <int dim>
std::ostream& operator<<(std::ostream& out,	const DimLocalDoFSet<dim>& v)
{
	out << *dynamic_cast<const LocalDoFSet*>(&v);

	for(size_t dof = 0; dof < v.num_dof(); ++dof){
		out << "DoF "<<std::setw(3)<<dof<<": ";
		MathVector<dim> pos; v.position(dof, pos);
		out << v.local_dof(dof) << ", " << pos << "\n";
	}
	return out;
}

template std::ostream& operator<<(std::ostream& out, const DimLocalDoFSet<0>& v);
template std::ostream& operator<<(std::ostream& out, const DimLocalDoFSet<1>& v);
template std::ostream& operator<<(std::ostream& out, const DimLocalDoFSet<2>& v);
template std::ostream& operator<<(std::ostream& out, const DimLocalDoFSet<3>& v);

} // end namespace ug
