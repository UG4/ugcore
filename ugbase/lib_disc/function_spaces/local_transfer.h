/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
 * Authors: Andreas Vogel, Christian Wehner
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

#ifndef __H__UG__LIB_DISC__FUNCTION_SPACS__LOCAL_TRANSFER__
#define __H__UG__LIB_DISC__FUNCTION_SPACS__LOCAL_TRANSFER__

#include "local_transfer_interface.h"
#include "lib_disc/local_finite_element/local_finite_element_provider.h"
#include "lib_disc/domain_util.h"
#include "lib_disc/reference_element/reference_mapping_provider.h"
#include "lib_disc/function_spaces/dof_position_util.h"

namespace ug {


template <typename TDomain>
class PiecewiseConstantElemTransfer
	: public ElemProlongationBase<TDomain, PiecewiseConstantElemTransfer<TDomain> >
	, public ElemRestrictionBase<TDomain, PiecewiseConstantElemTransfer<TDomain> >
{
	public:
		explicit PiecewiseConstantElemTransfer(const LFEID& lfeid) : m_lfeid(lfeid) {}

		~PiecewiseConstantElemTransfer() override = default;

		bool perform_prolongation_on(GridBaseObjectId gbo) override
		{
			if(m_lfeid.dim() == gbo) return true;
			return false;
		}

		template <typename TElem>
		void prolongate(TElem* parent,
		                TransferValueAccessor& vValueChild,
		                TransferValueAccessor& vValueParent)
		{
			const MultiGrid* mg = IElemProlongation<TDomain>::m_spGrid.get();
			const int numChild = mg->num_children<TElem>(parent);

			vValueParent.access_inner(parent);
			for(int c = 0; c < numChild; ++c){
				TElem* child = mg->get_child<TElem>(parent, c);

				vValueChild.access_inner(child);
				vValueChild[0] = vValueParent[0];
			}
		}

		bool perform_restriction_on(GridBaseObjectId gbo) override
		{
			if(m_lfeid.dim() == gbo) return true;
			return false;
		}

		template <typename TElem>
		void do_restrict(TElem* parent,
		                TransferValueAccessor& vValueChild,
		                TransferValueAccessor& vValueParent)
		{
			const MultiGrid* mg = IElemRestriction<TDomain>::m_spGrid.get();
			const int numChild = mg->num_children<TElem>(parent);

			vValueParent.access_inner(parent);
			UG_ASSERT(vValueParent.size()==1, "Exactly one DoF for Piecewise-Constant")
			vValueParent[0] = 0.0;

			for(int c = 0; c < numChild; ++c){
				TElem* child = mg->get_child<TElem>(parent, c);

				vValueChild.access_inner(child);
				UG_ASSERT(vValueChild.size()==1, "Exactly one DoF for Piecewise-Constant")
				vValueParent[0] += vValueChild[0];
			}
			vValueParent[0] *= (1./numChild);
		}

	protected:
		LFEID m_lfeid;
};

template <typename TDomain>
class P1LagrangeElemTransfer
	: public ElemProlongationBase<TDomain, P1LagrangeElemTransfer<TDomain> >
	, public ElemRestrictionBase<TDomain, P1LagrangeElemTransfer<TDomain> >
{
	public:
		P1LagrangeElemTransfer(const LFEID& lfeid) : m_lfeid(lfeid) {}

		~P1LagrangeElemTransfer() override = default;

		bool perform_prolongation_on(GridBaseObjectId gbo) override
		{
			if(m_lfeid.dim() < gbo) return false;
			return true;
		}

		template <typename TElem>
		void prolongate(TElem* parent,
		                TransferValueAccessor& vValueChild,
		                TransferValueAccessor& vValueParent)
		{
			const MultiGrid* mg = IElemProlongation<TDomain>::m_spGrid.get();
			const int numChild = mg->num_children<Vertex>(parent);
			if(numChild == 0) return;
			if(numChild != 1) UG_THROW("Max child Vertex must be 1");

			Vertex* child = mg->get_child<Vertex>(parent, 0);

			vValueChild.access_inner(child);
			vValueParent.access_closure(parent);

			if (vValueChild.size() == 0)
				return;	// current fct not defined on child

			const int numVrt = vValueParent.size();
			UG_ASSERT(numVrt > 0, "No function value found on parent vertices.");

			vValueChild[0] = vValueParent[0];
			for(int i = 1; i < numVrt; ++i)
				vValueChild[0] += vValueParent[i];
			vValueChild[0] *= 1.0/numVrt;
		}

		bool perform_restriction_on(GridBaseObjectId gbo) override
		{
			if(gbo == VERTEX) return true;
			return false;
		}

		// the following line silences -Woverloaded-virtual
		using ElemRestrictionBase<TDomain, P1LagrangeElemTransfer<TDomain> >::do_restrict;
		void do_restrict(GridObject* parent,
		                 TransferValueAccessor& vValueChild,
		                 TransferValueAccessor& vValueParent)
		{
			UG_THROW("Should not be called.")
		}

		void do_restrict(Vertex* parent,
		                 TransferValueAccessor& vValueChild,
		                 TransferValueAccessor& vValueParent)
		{
			const MultiGrid* mg = IElemRestriction<TDomain>::m_spGrid.get();
			const int numChild = mg->num_children<Vertex>(parent);
			if(numChild != 1) UG_THROW("Num child Vertex must be 1");

			Vertex* child = mg->get_child<Vertex>(parent, 0);

			vValueChild.access_inner(child);
			vValueParent.access_closure(parent);

			if (vValueParent.size() == 0)
				return; // current fct not defined on parent

			UG_ASSERT(vValueChild.size() > 0, "No function value found on child vertex.");

			vValueParent[0] = vValueChild[0];
		}

	protected:
		LFEID m_lfeid;
};


template <typename TDomain>
class StdLagrangeElemTransfer
	: public ElemProlongationBase<TDomain,StdLagrangeElemTransfer<TDomain> >
	, public ElemRestrictionBase<TDomain,StdLagrangeElemTransfer<TDomain> >
{
	public:
	///	world dimension
		static constexpr int dim = TDomain::dim;

	public:
		StdLagrangeElemTransfer(const LFEID& lfeid) : m_lfeid(lfeid) {}

		~StdLagrangeElemTransfer() override = default;

		bool perform_prolongation_on(GridBaseObjectId gbo) override
		{
			if(m_lfeid.order() == 1 && gbo != VERTEX) return false;
			if(m_lfeid.dim() < gbo) return false;
			return true;
		}


		template <typename TParent, typename TChild>
		void prolongate(TParent* parent,
		                TransferValueAccessor& vValueChild,
		                TransferValueAccessor& vValueParent,
		                const LocalShapeFunctionSet<TParent::dim>& lsfs,
		                const std::vector<MathVector<dim> >& vCornerParent,
		                const ReferenceObjectID parentRoid)
		{
			const MultiGrid* mg = IElemProlongation<TDomain>::m_spGrid.get();
			static constexpr int pdim = TParent::dim;

		//	 number of children of the base type
			const int numChild = mg->num_children<TChild>(parent);

		//	loop all children
			for(int c = 0; c < numChild; ++c)
			{
			//	get child
				TChild* child = mg->get_child<TChild>(parent, c);

			//	access the values
				vValueChild.access_inner(child);

			//	global positions of fine dofs
				std::vector<MathVector<dim> > vDoFPos;
				InnerDoFPosition(vDoFPos, child, *IElemProlongation<TDomain>::m_spDomain, m_lfeid);

			//	get Reference Mapping
				DimReferenceMapping<pdim, dim>& map =
					ReferenceMappingProvider::get<pdim, dim>(parentRoid, vCornerParent);

			//	get local position of DoF
				std::vector<MathVector<pdim> > vLocPos(vDoFPos.size(), 0.0);
				map.global_to_local(vLocPos, vDoFPos);

			//	evaluate coarse shape fct at fine local point
				std::vector<std::vector<number> > vvShape;
				lsfs.shapes(vvShape, vLocPos);

			//	sum up interpolated values
				for(size_t csh = 0; csh < vValueChild.size(); ++csh){
					vValueChild[csh] = 0.0;
					for(size_t psh = 0; psh < vValueParent.size(); ++psh)
						vValueChild[csh] += vvShape[csh][psh] * vValueParent[psh];
				}
			}
		}


		template <typename TElem>
		void prolongate(TElem* parent,
		                TransferValueAccessor& vValueChild,
		                TransferValueAccessor& vValueParent)
		{
			static constexpr int pdim = TElem::dim;

		//	get reference object ids
			const ReferenceObjectID parentRoid = parent->reference_object_id();

		//	access the values
			vValueParent.access_closure(parent);

		//	get vertices if required by some prolongation
			std::vector<MathVector<dim> > vCornerParent;
			CollectCornerCoordinates(vCornerParent, *parent, *IElemProlongation<TDomain>::m_spDomain);

		//	get local finite element trial spaces
			const LocalShapeFunctionSet<pdim>& lsfs =
					LocalFiniteElementProvider::get<pdim>(parentRoid, m_lfeid);

		//	prolongate to all children with dimension pdim and lower
			switch(pdim){
				case 3: prolongate<TElem,Volume>(parent, vValueChild, vValueParent, lsfs, vCornerParent, parentRoid);
				case 2: prolongate<TElem,Face>(parent, vValueChild, vValueParent, lsfs, vCornerParent, parentRoid);
				case 1: prolongate<TElem,Edge>(parent, vValueChild, vValueParent, lsfs, vCornerParent, parentRoid);
						prolongate<TElem,Vertex>(parent, vValueChild, vValueParent, lsfs, vCornerParent, parentRoid);
				break;
				default: UG_THROW("Dimension "<<pdim<<" not supported");
			}
		}

		void prolongate(Vertex* parent) override
		{
			const MultiGrid* mg = IElemProlongation<TDomain>::m_spGrid.get();
			TransferValueAccessor& vValueChild = *IElemProlongation<TDomain>::m_vValueChild;
            TransferValueAccessor& vValueParent = *IElemProlongation<TDomain>::m_vValueParent;

		//	access the values
			vValueParent.access_closure(parent);

		//	 number of children of the base type
			const int numChild = mg->num_children<Vertex>(parent);
			if(numChild == 0) return;
			if(numChild != 1) UG_THROW("Num child Vertex must be 1");

		//	get child
			Vertex* child = mg->get_child<Vertex>(parent, 0);

		//	access the values
			vValueChild.access_inner(child);

		//	copy value
			vValueChild[0] = vValueParent[0];
		}


		bool perform_restriction_on(GridBaseObjectId gbo) override
		{
			if(m_lfeid.dim() < gbo) return false;
			return true;
		}

		template <typename TElem>
		void do_restrict(TElem* parent,
		                TransferValueAccessor& vValueChild,
		                TransferValueAccessor& vValueParent)
		{
		//	access the values
			vValueParent.access_inner(parent);

		//	global positions of fine dofs
			std::vector<MathVector<dim> > vDoFPos;
			InnerDoFPosition(vDoFPos, parent, *IElemRestriction<TDomain>::m_spDomain, m_lfeid);

			const MultiGrid* mg = IElemRestriction<TDomain>::m_spGrid.get();

		//	 number of children of the base type
			const int numChild = mg->num_children<TElem>(parent);

		//	loop all children
			for(size_t ip = 0; ip < vDoFPos.size(); ++ip)
			{
				int c = 0;
				for(; c < numChild; ++c)
				{
				//	get child
					TElem* child = mg->get_child<TElem>(parent, c);

				//	global positions of fine dofs
					std::vector<MathVector<dim> > vChildDoFPos;
					DoFPosition(vChildDoFPos, child, *IElemRestriction<TDomain>::m_spDomain, m_lfeid);

				//	check if pos matches
					int cip = -1;
					for(size_t i = 0; i < vChildDoFPos.size(); ++i)
						if(VecDistanceSq(vChildDoFPos[i], vDoFPos[ip]) < 1e-8)
							cip = i;

				//	if not found in this child, continue
					if(cip == -1) continue;

				//	access the values
					vValueChild.access_closure(child);

				//	sum up interpolated values
					vValueParent[ip] = vValueChild[cip];

					break;
				}
			//	check that DoF found
				if(c == numChild)
					UG_THROW("Lagrange-Element-Coarsening: Cannot find DoF position "
							"in child elems for "<<ip<<"'s Pos: "<<vDoFPos[ip])
			}

		}


		void do_restrict(Vertex* parent) override
		{
			const MultiGrid* mg = IElemRestriction<TDomain>::m_spGrid.get();
			const int numChild = mg->num_children<Vertex>(parent);
			if(numChild != 1) UG_THROW("Num child Vertex must be 1");

			TransferValueAccessor& vValueChild = *IElemRestriction<TDomain>::m_vValueChild;
            TransferValueAccessor& vValueParent = *IElemRestriction<TDomain>::m_vValueParent;

			Vertex* child = mg->get_child<Vertex>(parent, 0);

			vValueChild.access_inner(child);
			vValueParent.access_inner(parent);

			vValueParent[0] = vValueChild[0];
		}


	protected:
		LFEID m_lfeid;
};


template <typename TDomain>
class CrouzeixRaviartElemTransfer
	: public ElemProlongationBase<TDomain, CrouzeixRaviartElemTransfer<TDomain> >
	, public ElemRestrictionBase<TDomain, CrouzeixRaviartElemTransfer<TDomain> >
{
	public:
	///	world dimension
		static constexpr int dim = TDomain::dim;

	public:
		CrouzeixRaviartElemTransfer(const LFEID& lfeid) : m_lfeid(lfeid) {}

		~CrouzeixRaviartElemTransfer() override = default;

		bool perform_prolongation_on(GridBaseObjectId gbo) override
		{
		//	prolongation from elems that have a side as child
			if(m_lfeid.dim()-1 == gbo) return true;
			if(m_lfeid.dim()   == gbo) return true;
			return false;
		}

		// the following line silences -Woverloaded-virtual
		using ElemProlongationBase<TDomain, CrouzeixRaviartElemTransfer<TDomain> >::prolongate;
		void prolongate(Vertex* parent,
						TransferValueAccessor& vValueChild,
						TransferValueAccessor& vValueParent)
		{
			UG_THROW("No implementation for Vertex and Crouzeix-Raviart.")
		}

				template <typename TSide>
		void prolongate(const std::vector<typename TSide::sideof*>& vParentElem,
						const std::vector<TSide*>& vChildSide,
						TransferValueAccessor& vValueChild,
						TransferValueAccessor& vValueParent)
		{
		// 	get associated macro elements, this elem is a side of
			using TElem = typename TSide::sideof;
			if(vParentElem.size() > 2)
				UG_THROW("More than 2 Elements share a side.")

		//	reset values to zero
			for(size_t c = 0; c < vChildSide.size(); ++c){
				vValueChild.access_inner(vChildSide[c]);
				for(size_t csh = 0; csh < vValueChild.size(); ++csh)
					vValueChild[csh] = 0.0;
			}

		//	loop parent Elems
			for(size_t p = 0; p < vParentElem.size(); ++p)
			{
			//	get parent elem
				TElem* parentElem = vParentElem[p];

			//	get reference object ids
				static constexpr int pdim = TElem::dim;
				const ReferenceObjectID parentRoid = parentElem->reference_object_id();

			//	access the values
				vValueParent.access_closure(parentElem);

			//	get vertices if required by some prolongation
				std::vector<MathVector<dim> > vCornerParent;
				CollectCornerCoordinates(vCornerParent, *parentElem, *IElemProlongation<TDomain>::m_spDomain);

			//	get local finite element trial spaces
				const LocalShapeFunctionSet<pdim>& lsfs =
						LocalFiniteElementProvider::get<pdim>(parentRoid, m_lfeid);

			//	loop child sides
				for(size_t c = 0; c < vChildSide.size(); ++c)
				{
				//	get child side
					TSide* childSide = vChildSide[c];

				//	access the value
					vValueChild.access_inner(childSide);

				//	global positions of fine dofs
					std::vector<MathVector<dim> > vDoFPos;
					InnerDoFPosition(vDoFPos, childSide, *IElemProlongation<TDomain>::m_spDomain, m_lfeid);

				//	get Reference Mapping
					DimReferenceMapping<pdim, dim>& map =
						ReferenceMappingProvider::get<pdim, dim>(parentRoid, vCornerParent);

				//	get local position of DoF
					std::vector<MathVector<pdim> > vLocPos(vDoFPos.size(), 0.0);
					map.global_to_local(vLocPos, vDoFPos);

				//	evaluate coarse shape fct at fine local point
					std::vector<std::vector<number> > vvShape;
					lsfs.shapes(vvShape, vLocPos);

				//	sum up interpolated values
					for(size_t csh = 0; csh < vValueChild.size(); ++csh){
						for(size_t psh = 0; psh < vValueParent.size(); ++psh)
							vValueChild[csh] += (1./vParentElem.size()) *
												vvShape[csh][psh] * vValueParent[psh];
					}
				}
			}
		}
		
		template <typename TParent>
		void prolongate(TParent* parent,
						TransferValueAccessor& vValueChild,
						TransferValueAccessor& vValueParent)
		{
			const MultiGrid* mg = IElemProlongation<TDomain>::m_spGrid.get();
		//	a) prolongation from a side
			if(TParent::dim == m_lfeid.dim()-1){
				using TElem = typename TParent::sideof;
				using TSide = TParent;

			//	get child sides
				const int numChild = mg->num_children<TParent>(parent);
				if(numChild == 0) return;

			//	get children
				std::vector<TSide*> vChildSide(numChild);
				for(int c = 0; c < numChild; ++c)
					vChildSide[c] = mg->get_child<TParent>(parent, c);

			// 	get associated macro elements, this elem is a side of
				typename Grid::traits<TElem>::secure_container vElem;
				const_cast<MultiGrid*>(mg)->associated_elements(vElem, parent);
				std::vector<TElem*> vParentElem(vElem.size());
				for(size_t p = 0; p < vElem.size(); ++p)
					vParentElem[p] = vElem[p];

			//	call prolongation
				this->template prolongate<TSide>(vParentElem, vChildSide, vValueChild, vValueParent);
			}

		//	b) prolongation from an element
			if(TParent::dim == m_lfeid.dim()){
				using TElem = TParent;
				using TSide = typename TParent::side;

			//	get child sides
				const int numChild = mg->num_children<TSide>(parent);
				if(numChild == 0) return;

			//	get childs
				std::vector<TSide*> vChildSide(numChild);
				for(int c = 0; c < numChild; ++c)
					vChildSide[c] = mg->get_child<TSide>(parent, c);

			// 	get associated macro elements, this elem is a side of
				std::vector<TElem*> vParentElem(1);
				vParentElem[0] = parent;

			//	call prolongation
				this->template prolongate<TSide>(vParentElem, vChildSide, vValueChild, vValueParent);
			}
		};

		bool perform_restriction_on(GridBaseObjectId gbo) override
		{
		//	restriction only on sides
			if(m_lfeid.dim()-1 == gbo) return true;
			return false;
		}

		template <typename TParent>
		void do_restrict(TParent* parent,
		                 TransferValueAccessor& vValueChild,
		                 TransferValueAccessor& vValueParent)
		{
			const MultiGrid* mg = IElemRestriction<TDomain>::m_spGrid.get();

			if(TParent::dim != m_lfeid.dim()-1)
				UG_THROW("Only restrict on sides.")

		//	get child sides
			const int numChild = mg->num_children<TParent>(parent);
			if(numChild == 0) return;

		//	access the values
			vValueParent.access_inner(parent);
			UG_ASSERT(vValueParent.size() == 1, "Num DoFs per Side must be 1 for CR");

		//	reset values
			vValueParent[0] = 0.0;

		//	get children
			for(int c = 0; c < numChild; ++c)
			{
			//	get child side
				TParent* childSide = mg->get_child<TParent>(parent, c);

			//	access the value
				vValueChild.access_inner(childSide);
				UG_ASSERT(vValueChild.size() == 1, "Num DoFs per Side must be 1 for CR");

			//	sum up values
				vValueParent[0] += vValueChild[0];
			}

			vValueParent[0] *= (1./numChild);
		}

	protected:
		LFEID m_lfeid;
};


} // end namespace ug

#endif