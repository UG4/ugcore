/*
 * local_finite_element_provider.cpp
 *
 *  Created on: 17.02.2010
 *      Author: andreasvogel
 */

#include "local_finite_element_provider.h"
#include "lib_disc/reference_element/reference_element_util.h"
#include "lib_disc/reference_element/reference_mapping_provider.h"

// include spaces
#include "lagrange/lagrangep1.h"
#include "lagrange/lagrange.h"
#include "crouzeix-raviart/crouzeix_raviart.h"
#include "piecewise_constant/piecewise_constant.h"
#include "mini/mini.h"
#include "nedelec/nedelec.h"


namespace ug{

////////////////////////////////////////////////////////////////////////////////
//	Wrapper for LocalShapeFunctionSets
////////////////////////////////////////////////////////////////////////////////

/// wrapper class implementing the LocalShapeFunctionSet interface
/**
 * This class wrappes a class passed by the template argument into the
 * virtual ILocalShapeFunctionSet interface and makes it thus usable in that
 * context on the price of virtual functions.
 *
 * \tparam 	TImpl		Implementation of a Local Shape Function Set
 */
template <typename TImpl>
class LocalShapeFunctionSetWrapper
	: public LocalShapeFunctionSet<TImpl::dim,
	  	  	  	  	  	  	  	   typename TImpl::shape_type,
	  	  	  	  	  	  	  	   typename TImpl::grad_type>,
	  public TImpl
{
	/// Implementation
		typedef TImpl ImplType;

	public:
	///	reference element dimension
		static const int dim = TImpl::dim;

	///	Shape type
		typedef typename ImplType::shape_type shape_type;

	///	Gradient type
		typedef typename ImplType::grad_type grad_type;

	public:
	///	constructor
		LocalShapeFunctionSetWrapper(){}

	public:
	///	\copydoc ug::LocalDoFSet::type()
		virtual ReferenceObjectID roid() const {return ImplType::roid();}

	///	\copydoc ug::LocalDoFSet::num_sh()
		virtual size_t num_sh() const {return ImplType::num_sh();}

	///	\copydoc ug::LocalDoFSet::num_dof()
		virtual size_t num_dof(ReferenceObjectID roid) const {return ImplType::num_dof(roid);}

	///	\copydoc ug::LocalDoFSet::local_dof()
		virtual const LocalDoF& local_dof(size_t dof) const {return ImplType::local_dof(dof);}

	///	\copydoc ug::LocalDoFSet::exact_position_available()
		virtual bool exact_position_available() const {return ImplType::exact_position_available();}

	///	\copydoc ug::LocalShapeFunctionSet::position()
		virtual bool position(size_t i, MathVector<dim>& pos) const
		{
			return ImplType::position(i, pos);
		}

	public:
	///	\copydoc ug::LocalShapeFunctionSet::continuous()
		virtual bool continuous() const {return ImplType::continuous();}

	///	\copydoc ug::LocalShapeFunctionSet::shape()
		virtual shape_type shape(size_t i, const MathVector<dim>& x) const
		{
			return ImplType::shape(i, x);
		}

	///	\copydoc ug::LocalShapeFunctionSet::shape()
		virtual void shape(shape_type& s, size_t i, const MathVector<dim>& x) const
		{
			typedef BaseLSFS<TImpl, TImpl::dim,
				  	  	  	  	  	  	  	   typename TImpl::shape_type,
				  	  	  	  	  	  	  	   typename TImpl::grad_type> baseType;

				  	  	  	  	  	  	  	   baseType::shape(s, i, x);
		}

	///	\copydoc ug::LocalShapeFunctionSet::shapes()
		virtual void shapes(shape_type* vShape, const MathVector<dim>& x) const
		{
			ImplType::shapes(vShape, x);
		}

	///	\copydoc ug::LocalShapeFunctionSet::shapes()
		virtual void shapes(std::vector<shape_type>& vShape, const MathVector<dim>& x) const
		{
			ImplType::shapes(vShape, x);
		}

	///	\copydoc ug::LocalShapeFunctionSet::shapes()
		virtual void shapes(std::vector<std::vector<shape_type> >& vvShape,
							const std::vector<MathVector<dim> >& vLocPos) const
		{
			ImplType::shapes(vvShape, vLocPos);
		}

	///	\copydoc ug::LocalShapeFunctionSet::grad()
		virtual void grad(grad_type& g, size_t i, const MathVector<dim>& x) const
		{
			ImplType::grad(g, i, x);
		}

	///	\copydoc ug::LocalShapeFunctionSet::grads()
		virtual void grads(grad_type* vGrad, const MathVector<dim>& x) const
		{
			ImplType::grads(vGrad, x);
		}

	///	\copydoc ug::LocalShapeFunctionSet::grads()
		virtual void grads(std::vector<grad_type>& vGrad, const MathVector<dim>& x) const
		{
			ImplType::grads(vGrad, x);
		}

	///	\copydoc ug::LocalShapeFunctionSet::grads()
		virtual void grads(std::vector<std::vector<grad_type> >& vvGrad,
						   const std::vector<MathVector<dim> >& vLocPos) const
		{
			ImplType::grads(vvGrad, vLocPos);
		}
};

////////////////////////////////////////////////////////////////////////////////
// 	SubLocalDoFSet
////////////////////////////////////////////////////////////////////////////////

/**
 * Intersection of local dofs
 */
template <int dim>
class SubLocalDoFSet : public DimLocalDoFSet<dim>
{
	public:
		template <int setDim>
		SubLocalDoFSet(const ReferenceObjectID roid,
					   const DimLocalDoFSet<setDim>& set)
		   : m_roid(roid),
			 m_bExactPos(set.exact_position_available()),
			 m_bInit(false),
			 m_vNumDoF(NUM_REFERENCE_OBJECTS, 0)
		{
			if(ReferenceElementDimension(roid) != dim)
				UG_THROW("SubLocalDoFSet: templated Reference Dimension "<<dim<<
				         " does not match Reference Element dimension of "<<roid);

			const DimReferenceElement<setDim>& setRefElem =
					ReferenceElementProvider::get<setDim>(set.roid());

		// loop all subs of roid type
			for(size_t id = 0; id < setRefElem.num(dim); ++id){
				if(setRefElem.roid(dim, id) != m_roid) continue;

			// 	get mapping to sub
				std::vector<MathVector<setDim> > vCorner;
				for(size_t co = 0; co < setRefElem.num(dim, id, 0); ++co){
					vCorner.push_back(setRefElem.corner(setRefElem.id(dim, id, 0, co)));
				}
				const DimReferenceMapping<dim, setDim>& map =
							ReferenceMappingProvider::get<dim, setDim>(m_roid, vCorner);

				std::vector<size_t> vNumDoF(NUM_REFERENCE_OBJECTS, 0);
				std::vector<LocalDoF> vLocalDoF;
				std::vector<MathVector<dim> > vLocalPos;

			//	loop subselements of sub
				for(int subDim = 0; subDim <= dim; ++subDim){

					for(size_t i = 0; i < setRefElem.num(dim, id, subDim); ++i){
						const size_t subID = setRefElem.id(dim, id, subDim, i);

						const ReferenceObjectID sroid =  setRefElem.roid(subDim, subID);
						vNumDoF[sroid] = set.num_dof(sroid);

					//	get local DoFs on sub
						std::vector<size_t> vLocalDoFID;
						for(size_t dof = 0; dof < set.num_dof(); ++dof){
							const LocalDoF& localDoF = set.local_dof(dof);
							if(localDoF.dim() != subDim) continue;
							if(localDoF.id() != subID) continue;

							vLocalDoFID.push_back(dof);
						}

					//	loop sub DoFs in correct order
						for(size_t offset = 0; offset < vLocalDoFID.size(); ++offset){

							const LocalDoF* pLocalDoF = NULL;
							size_t localDoFID;
							for(size_t dof = 0; dof < vLocalDoFID.size(); ++dof){
								if(set.local_dof(vLocalDoFID[dof]).offset() == offset){
									pLocalDoF = &set.local_dof(vLocalDoFID[dof]);
									localDoFID = vLocalDoFID[dof];
								}
							}

							if(pLocalDoF == NULL)
								UG_THROW("SubLocalDoFSet: Cannot find local dof "
										"with offset "<<offset<<", but must "
										"exist. Check Implementation of LocalDoFSet ");

						//	map local dof
							MathVector<dim> locPos;
							MathVector<setDim> globPos;
							set.position(localDoFID, globPos);

							try{
								map.global_to_local(locPos, globPos);
							}
							catch(UGError& err){
								std::stringstream ss;
								ss<<"SubLocalDoFSet: Cannot find local position for "
								"global Position "<<globPos<<" of the "<<localDoFID<<
								"'th LocalDoF on "<<set.roid()<<
								" for the "<<id<<"'th "<<roid<<" and its "
								"Sub-"<<sroid<<" with id "<<subID<<" with "
								"corners:\n";
								for(size_t i = 0; i < vCorner.size(); ++i)
									ss << " "<<i<<": "<<vCorner[i] << "\n";

								err.push_msg(ss.str(),__FILE__,__LINE__);
								throw(err);
							}

						//	add
							vLocalDoF.push_back(LocalDoF(subDim, i, offset));
							vLocalPos.push_back(locPos);

							this->set(vNumDoF, vLocalDoF, vLocalPos);
						}
					}
				}
			}
		}

	public:
		virtual ReferenceObjectID roid() const {return m_roid;};
		virtual size_t num_dof(ReferenceObjectID roid) const {return m_vNumDoF[roid];}
		virtual const LocalDoF& local_dof(size_t dof) const {return m_vLocalDoF[dof];}
		virtual bool exact_position_available() const {return m_bExactPos;}
		virtual bool position(size_t i, MathVector<dim>& pos) const
		{
			pos = m_vLocalPos[i];
			return exact_position_available();
		}

	protected:
		void set(const std::vector<size_t>& vNumDoF,
				 const std::vector<LocalDoF>& vLocalDoF,
				 const std::vector<MathVector<dim> >& vLocalPos)
		{
			if(m_bInit){
				if(vNumDoF != m_vNumDoF) UG_THROW("NumDoF mismatch");
				if(vLocalDoF != m_vLocalDoF) UG_THROW("vLocalDoF mismatch");
				if(vLocalPos != m_vLocalPos) UG_THROW("vLocalPos mismatch");
			}
			else{
				m_vNumDoF = vNumDoF;
				m_vLocalDoF = vLocalDoF;
				m_vLocalPos = vLocalPos;
			}
		}


	protected:
		const ReferenceObjectID m_roid; ///< Reference ID this DoF Set is for
		const bool m_bExactPos;

		bool m_bInit;
		std::vector<size_t> m_vNumDoF;
		std::vector<LocalDoF> m_vLocalDoF; ///< Local DoFs of this set
		std::vector<MathVector<dim> > m_vLocalPos; ///< Local Positions of DoFs
};

////////////////////////////////////////////////////////////////////////////////
//	Provider for LocalShapeFunctionSets
////////////////////////////////////////////////////////////////////////////////

template <typename TRefElem>
void LocalFiniteElementProvider::create_lagrange_set(const LFEID& id)
{
//	reference object id
	static const ReferenceObjectID roid = TRefElem::REFERENCE_OBJECT_ID;
	static const int dim = TRefElem::dim;

//	only order >= 1 available
	if(id.order() < 1) return;

//	if refdim == dim create the space
	if(dim == id.dim()){
		ConstSmartPtr<LocalShapeFunctionSet<dim> > set;
		switch(id.order()){
			case 1: set = ConstSmartPtr<LocalShapeFunctionSet<dim> >
						  (new LocalShapeFunctionSetWrapper<LagrangeP1<TRefElem> >);
					break;
			case 2: set = ConstSmartPtr<LocalShapeFunctionSet<dim> >
						  (new LocalShapeFunctionSetWrapper<LagrangeLSFS<TRefElem,2> >);
					break;
			break;
			case 3: set = ConstSmartPtr<LocalShapeFunctionSet<dim> >
						  (new LocalShapeFunctionSetWrapper<LagrangeLSFS<TRefElem,3> >);
					break;
			break;
			case 4: set = ConstSmartPtr<LocalShapeFunctionSet<dim> >
						  (new LocalShapeFunctionSetWrapper<LagrangeLSFS<TRefElem,4> >);
					break;
			break;

			default:{
				SmartPtr<LocalShapeFunctionSetWrapper<FlexLagrangeLSFS<TRefElem> > >
					sSetFlexLagrange(new LocalShapeFunctionSetWrapper<FlexLagrangeLSFS<TRefElem> >);
				sSetFlexLagrange->set_order(id.order());
				set = sSetFlexLagrange;
			}
		}
		register_set(id, set);
		return;
	}
//	if refdim < dim, the restriction of to the subelement is the lagrange space
//	of the refdim
	else if (dim < id.dim()){
	//  get (and maybe create) lagrange space for dim == refdim (always exist)
		ConstSmartPtr<LocalShapeFunctionSet<dim> > set =
			getptr<dim>(roid, LFEID(id.type(), dim, id.order()));

	//	register this space as space on subelement
		register_set(id, set);
	}

	// else: for refdim > dim there exist no restrictions
}

void LocalFiniteElementProvider::
create_lagrange_set(ReferenceObjectID roid, const LFEID& id)
{
	switch(roid){
		case ROID_EDGE:			create_lagrange_set<ReferenceEdge>(id); return;
		case ROID_TRIANGLE:		create_lagrange_set<ReferenceTriangle>(id); return;
		case ROID_QUADRILATERAL:create_lagrange_set<ReferenceQuadrilateral>(id); return;
		case ROID_TETRAHEDRON:	create_lagrange_set<ReferenceTetrahedron>(id); return;
		case ROID_PRISM:		create_lagrange_set<ReferencePrism>(id); return;
		case ROID_HEXAHEDRON:	create_lagrange_set<ReferenceHexahedron>(id); return;
		case ROID_PYRAMID:
			// only space available for order 1
			if(id.order() != 1) return;
			{
				register_set(LFEID(LFEID::LAGRANGE, 3, 1),
							 ConstSmartPtr<LocalShapeFunctionSet<3> >
							(new LocalShapeFunctionSetWrapper<LagrangeP1<ReferencePyramid> >));
			}
		default: return;
	}
}

template <typename TRefElem>
void LocalFiniteElementProvider::create_mini_bubble_set(const LFEID& id)
{
//	reference object id
	static const ReferenceObjectID roid = TRefElem::REFERENCE_OBJECT_ID;
	static const int dim = TRefElem::dim;

//	if refdim == dim create the space
	if(dim == id.dim()){
		if(id == LFEID(LFEID::MINI, dim, 1))
			register_set(id,  ConstSmartPtr<LocalShapeFunctionSet<dim> >(
						 new LocalShapeFunctionSetWrapper<MiniBubbleLSFS<TRefElem> >));
		return;
	}
//	if refdim < dim, the restriction of to the subelement is the lagrange space
//	of the refdim
	else if (dim < id.dim()){
	//  get (and maybe create) lagrange space for dim == refdim (always exist)
		ConstSmartPtr<LocalShapeFunctionSet<dim> > set =
			getptr<dim>(roid, LFEID(LFEID::LAGRANGE, dim, 1));

	//	register this space as space on subelement
		register_set(id, set);
	}

	// else: for refdim > dim there exist no restrictions
}

void LocalFiniteElementProvider::
create_mini_bubble_set(ReferenceObjectID roid, const LFEID& id)
{
	switch(roid){
		case ROID_EDGE:			create_mini_bubble_set<ReferenceEdge>(id); return;
		case ROID_TRIANGLE:		create_mini_bubble_set<ReferenceTriangle>(id); return;
		case ROID_QUADRILATERAL:create_mini_bubble_set<ReferenceQuadrilateral>(id); return;
		case ROID_TETRAHEDRON:	create_mini_bubble_set<ReferenceTetrahedron>(id); return;
		// no spaces implemented for other 3d elems
		default: return;
	}
}

template <typename TRefElem>
void LocalFiniteElementProvider::create_nedelec_set(const LFEID& id)
{
	static const int dim = TRefElem::dim;
	if(id == LFEID(LFEID::NEDELEC, dim, 1))
		register_set(id, ConstSmartPtr<LocalShapeFunctionSet<dim, MathVector<dim>, MathMatrix<dim,dim> > >(
					 new LocalShapeFunctionSetWrapper<NedelecLSFS<TRefElem> >));
	return;
}

void LocalFiniteElementProvider::
create_nedelec_set(ReferenceObjectID roid, const LFEID& id)
{
	switch(roid){
		case ROID_TRIANGLE:		create_nedelec_set<ReferenceTriangle>(id); return;
		case ROID_TETRAHEDRON:	create_nedelec_set<ReferenceTetrahedron>(id); return;
		default: return;
	}
}

template <typename TRefElem>
void LocalFiniteElementProvider::create_piecewise_constant_set(const LFEID& id)
{
	static const int dim = TRefElem::dim;
	if(id == LFEID(LFEID::PIECEWISE_CONSTANT, dim, 0))
		register_set(id, ConstSmartPtr<LocalShapeFunctionSet<dim> >(
					 new LocalShapeFunctionSetWrapper<PiecewiseConstantLSFS<TRefElem> >));
	return;
}

void LocalFiniteElementProvider::
create_piecewise_constant_set(ReferenceObjectID roid, const LFEID& id)
{
	switch(roid){
		case ROID_EDGE:			create_piecewise_constant_set<ReferenceEdge>(id); return;
		case ROID_TRIANGLE:		create_piecewise_constant_set<ReferenceTriangle>(id); return;
		case ROID_QUADRILATERAL:create_piecewise_constant_set<ReferenceQuadrilateral>(id); return;
		case ROID_TETRAHEDRON:	create_piecewise_constant_set<ReferenceTetrahedron>(id); return;
		case ROID_PRISM:		create_piecewise_constant_set<ReferencePrism>(id); return;
		case ROID_PYRAMID:		create_piecewise_constant_set<ReferencePyramid>(id); return;
		case ROID_HEXAHEDRON:	create_piecewise_constant_set<ReferenceHexahedron>(id); return;
		default: return;
	}
}


template <typename TRefElem>
void LocalFiniteElementProvider::create_crouxeiz_raviart_set(const LFEID& id)
{
	static const int dim = TRefElem::dim;
	if(id == LFEID(LFEID::CROUZEIX_RAVIART, dim, 1))
		register_set(id, ConstSmartPtr<LocalShapeFunctionSet<dim> >(
					 new LocalShapeFunctionSetWrapper<CrouzeixRaviartLSFS<TRefElem> >));
	return;
}

void LocalFiniteElementProvider::
create_crouxeiz_raviart_set(ReferenceObjectID roid, const LFEID& id)
{
	switch(roid){
		// no space for dim <= 1
		case ROID_TRIANGLE:		create_crouxeiz_raviart_set<ReferenceTriangle>(id); return;
		case ROID_QUADRILATERAL:create_crouxeiz_raviart_set<ReferenceQuadrilateral>(id); return;
		case ROID_TETRAHEDRON:	create_crouxeiz_raviart_set<ReferenceTetrahedron>(id); return;
		case ROID_PRISM:		create_crouxeiz_raviart_set<ReferencePrism>(id); return;
		case ROID_PYRAMID:		create_crouxeiz_raviart_set<ReferencePyramid>(id); return;
		case ROID_HEXAHEDRON:	create_crouxeiz_raviart_set<ReferenceHexahedron>(id); return;
		default: return;
	}
}

void LocalFiniteElementProvider::
create_set(ReferenceObjectID roid, const LFEID& id)
{
	try{
		switch(id.type())
		{
			case LFEID::LAGRANGE: 			create_lagrange_set(roid, id); break;
			case LFEID::MINI:   			create_mini_bubble_set(roid, id); break;
			case LFEID::NEDELEC:			create_nedelec_set(roid, id); break;
			case LFEID::PIECEWISE_CONSTANT:	create_piecewise_constant_set(roid, id); break;
			case LFEID::CROUZEIX_RAVIART:	create_crouxeiz_raviart_set(roid, id); break;
			default: return;
		}
	}
	UG_CATCH_THROW("LocalFiniteElementProvider: Creation of set "<<id<<
					" for "<<roid<<" failed.");
}

void LocalFiniteElementProvider::
create_set(const LFEID& id)
{
	switch(id.type())
	{
	//	this spaces can create also sub-elem spaces, since they are continuous
		case LFEID::LAGRANGE:
		case LFEID::MINI:

			for(int r = 0; r < NUM_REFERENCE_OBJECTS; ++r){
				const ReferenceObjectID roid = (ReferenceObjectID) r;
				const int dim = ReferenceElementDimension(roid);

				if(dim <= id.dim())
					create_set(roid, id);
			}
			break;

	//	this spaces can not create also sub-elem spaces, since they are discontinuous
		case LFEID::PIECEWISE_CONSTANT:
		case LFEID::CROUZEIX_RAVIART:
		case LFEID::NEDELEC:

			for(int r = 0; r < NUM_REFERENCE_OBJECTS; ++r){
				const ReferenceObjectID roid = (ReferenceObjectID) r;
				const int dim = ReferenceElementDimension(roid);

				if(dim == id.dim())
					create_set((ReferenceObjectID)roid, id);
			}
			break;

		default: return;
	}
}

template <int rdim, int dim>
void LocalFiniteElementProvider::
create_sub_dof_set(ReferenceObjectID roid, const LFEID& id)
{
	if(dim != id.dim())
		UG_THROW("Dimension must match here, internal error. ("<<dim<<","<<id.dim()<<")");

	for(int r = 0; r < NUM_REFERENCE_OBJECTS; ++r){
		const ReferenceObjectID elemRoid = (ReferenceObjectID) r;
		const int elemDim = ReferenceElementDimension(elemRoid);

		if(elemDim != dim) continue;

		ConstSmartPtr<DimLocalDoFSet<dim> > set = get_dof_ptr<dim>(elemRoid, id);
		if(set.invalid()) continue;

		try{
			ConstSmartPtr<DimLocalDoFSet<rdim> > subSet =
				ConstSmartPtr<DimLocalDoFSet<rdim> >(new SubLocalDoFSet<rdim>(roid, *set) );

			ConstSmartPtr<DimLocalDoFSet<rdim> > givenSubSet = get_dof_ptr<rdim>(roid, id, false);

			if(givenSubSet.valid()){
				if(*givenSubSet != *subSet)
					UG_THROW("LocalFiniteElementProvider::create_sub_dof_set: "
							"Creating DimLocalDoFSet for "<<roid<<" and "<<id<<
							". Already registered Set does not match with computed,"
							" but this indicates that the LocalDoFSets have not "
							"been implemented correctly. Check implementation."
							"Sets are: Given:"<<*givenSubSet<<"\nvs. New:\n"<<*subSet);
			}
			else {
				register_set<rdim>(id, subSet);
			}
		}
		UG_CATCH_THROW("LocalFiniteElementProvider::create_sub_dof_set: Cannot "
						"create SubDoFSet for "<<roid<<" and space "<<id<<" when "
						"trying to deduce from LocalDoFSet for "<<elemRoid<<".");
	}
}

void LocalFiniteElementProvider::
create_dof_set(ReferenceObjectID roid, const LFEID& id)
{
	const int dim = id.dim();
	const int rdim = ReferenceElementDimension(roid);

	if(rdim == dim) {
		create_set(roid, id);
		return;
	}

	if(rdim > dim)
		UG_THROW("Wrong dimension for SubDoFs: "<<rdim<<", "<<dim);

	switch(dim){
		case 1:
			switch(rdim){
				case 0: create_sub_dof_set<0, 1>(roid, id); return;
				default: return;
			}
		case 2:
			switch(rdim){
				case 0: create_sub_dof_set<0, 2>(roid, id); return;
				case 1: create_sub_dof_set<1, 2>(roid, id); return;
				default: return;
			}
		case 3:
			switch(rdim){
				case 0: create_sub_dof_set<0, 3>(roid, id); return;
				case 1: create_sub_dof_set<1, 3>(roid, id); return;
				case 2: create_sub_dof_set<2, 3>(roid, id); return;
				default: return;
			}
		default: return;
	}
}


LocalFiniteElementProvider::
LocalFiniteElementProvider()
{};

LocalFiniteElementProvider::
~LocalFiniteElementProvider()
{};

const LocalDoFSet& LocalFiniteElementProvider::
get_dofs(ReferenceObjectID roid, const LFEID& id, bool bCreate)
{
//	init provider and search for identifier
	typedef std::map<LFEID, LocalDoFSets> Map;
	Map::const_iterator iter = inst().m_mLocalDoFSets.find(id);

//	if not found
	if(iter == m_mLocalDoFSets.end() || (iter->second)[roid].invalid()){
		if(bCreate){
			create_dof_set(roid, id);
			return get_dofs(roid, id, false);
		}
		UG_THROW("LocalFiniteElementProvider: Cannot create LocalDoFSet for finite "
				"element type "<<id<<" and "<<roid);
	}

//	return dof set
	return *((iter->second)[roid]);
}

const CommonLocalDoFSet& LocalFiniteElementProvider::
get_dofs(const LFEID& id, bool bCreate)
{
//	init provider and search for identifier
	typedef std::map<LFEID, CommonLocalDoFSet> Map;
	Map::const_iterator iter = inst().m_mCommonDoFSet.find(id);

//	if not found
	if(iter == m_mCommonDoFSet.end()){
		if(bCreate)	{
			create_set(id);
			return get_dofs(id, false);
		}
		UG_THROW("LocalFiniteElementProvider: Cannot create CommonLocalDoFSet for "<<id);
	}

//	return the common set
	return iter->second;
}

bool LocalFiniteElementProvider::continuous(const LFEID& id, bool bCreate)
{
	std::map<LFEID, bool>::iterator iter = m_mContSpace.find(id);
	if(iter == m_mContSpace.end())
	{
	//	try to create the set
		if(bCreate){
			create_set(id);
			return continuous(id, false);
		}

		UG_THROW("LocalFiniteElementProvider::continuous: no shape function "
				"set "<<id<<" registered.");
	}

	return (*iter).second;
}

void LocalFiniteElementProvider::register_set(const LFEID& id, ConstSmartPtr<LocalDoFSet> set)
{
//	reference object id
	const ReferenceObjectID roid = set->roid();

//	register
	m_mLocalDoFSets[id][roid] = set;

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
		ss<<"LocalFiniteElementProvider::register_set(): "
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

std::map<LFEID, LocalFiniteElementProvider::LocalDoFSets>
LocalFiniteElementProvider::m_mLocalDoFSets =
std::map<LFEID, LocalFiniteElementProvider::LocalDoFSets>();

std::map<LFEID, CommonLocalDoFSet>
LocalFiniteElementProvider::m_mCommonDoFSet =
std::map<LFEID, CommonLocalDoFSet>();

std::map<LFEID, bool>
LocalFiniteElementProvider::m_mContSpace =
std::map<LFEID, bool>();

} // namespace ug

