/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter, Andreas Vogel
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

#ifndef __H__UG__LIB_DISC__FUNCTION_SPACE__GRID_FUNCTION_GLOBAL_USER_DATA__
#define __H__UG__LIB_DISC__FUNCTION_SPACE__GRID_FUNCTION_GLOBAL_USER_DATA__

#include "common/common.h"

#include "lib_grid/tools/subset_group.h"

#include "lib_disc/common/function_group.h"
#include "lib_disc/common/groups_util.h"
#include "lib_disc/quadrature/quadrature.h"
#include "lib_disc/domain_util.h"
#include "lib_disc/local_finite_element/local_finite_element_provider.h"
#include "lib_disc/spatial_disc/user_data/std_glob_pos_data.h"
#include "lib_disc/reference_element/reference_mapping_provider.h"
#include "lib_grid/algorithms/space_partitioning/lg_ntree.h"

#include <math.h>       /* fabs */

namespace ug{



template <typename TGridFunction, int elemDim = TGridFunction::dim>
class GlobalGridFunctionNumberData
	: public StdGlobPosData<GlobalGridFunctionNumberData<TGridFunction, elemDim>, number, TGridFunction::dim>
{
	public:
	///	world dimension of grid function
		static const int dim = TGridFunction::dim;
		typedef typename TGridFunction::template dim_traits<elemDim>::grid_base_object element_t;

		private:
	/// grid function
		SmartPtr<TGridFunction> m_spGridFct;

	///	component of function
		size_t m_fct;

	///	local finite element id
		LFEID m_lfeID;

		typedef lg_ntree<dim, dim, element_t>	tree_t;
		tree_t	m_tree;

	public:
	/// constructor
		GlobalGridFunctionNumberData(SmartPtr<TGridFunction> spGridFct, const char* cmp)
		: m_spGridFct(spGridFct),
		  m_tree(*spGridFct->domain()->grid(), spGridFct->domain()->position_attachment())
		{
			//this->set_functions(cmp);

			//	get function id of name
			m_fct = spGridFct->fct_id_by_name(cmp);

			//	check that function exists
			if(m_fct >= spGridFct->num_fct())
				UG_THROW("GridFunctionNumberData: Function space does not contain"
						" a function with name " << cmp << ".");

			//	local finite element id
			m_lfeID = spGridFct->local_finite_element_id(m_fct);



//			const MGSubsetHandler& ssh = *m_spGridFct->domain()->subset_handler();

			SubsetGroup ssGrp(m_spGridFct->domain()->subset_handler());
			ssGrp.add_all();

			std::vector<size_t> subsetsOfGridFunction;
			std::vector<element_t*> elemsWithGridFunctions;

			typename TGridFunction::template dim_traits<elemDim>::const_iterator iterEnd, iter;


			for(size_t si = 0; si < ssGrp.size(); si++){
				if( spGridFct->is_def_in_subset(m_fct, si) )
					subsetsOfGridFunction.push_back(si);
			}


			for(size_t i = 0; i<subsetsOfGridFunction.size(); i++){
				size_t si = subsetsOfGridFunction[i];
				iter = spGridFct->template begin<element_t>(si);
				iterEnd = spGridFct->template end<element_t>(si);

				for(;iter!=iterEnd; ++iter){
					element_t *elem = *iter;
					elemsWithGridFunctions.push_back(elem);
				}
			}

//			m_tree.create_tree(spGridFct->template begin<element_t>(si),
//														spGridFct->template end<element_t>(si));

			m_tree.create_tree(elemsWithGridFunctions.begin(), elemsWithGridFunctions.end());

		};

		virtual ~GlobalGridFunctionNumberData() {}

		virtual bool continuous() const
		{
			return LocalFiniteElementProvider::continuous(m_lfeID);
		}

		///	to full-fill UserData-Interface
		inline void evaluate(number& value, const MathVector<dim>& x, number time, int si) const
		{
			if(!evaluate_global(value, x))
				UG_THROW("For function "<<m_fct<<" couldn't find an element containing the specified point: " << x);
		}

		inline number evaluate(const MathVector<dim>& x) const
		{
			number value;
			if(!evaluate(value, x))
				UG_THROW("For function "<<m_fct<<" couldn't find an element containing the specified point: " << x);
			return value;
		}

		///	evaluates the data at a given point, returns false if point not found
		inline bool evaluate(number& value, const MathVector<dim>& x) const
		{
			element_t* elem = NULL;
			//try{

				if(!FindContainingElement(elem, m_tree, x)){
					return false;
				}

			//	get corners of element
				std::vector<MathVector<dim> > vCornerCoords;
				CollectCornerCoordinates(vCornerCoords, *elem, *m_spGridFct->domain());

			//	reference object id
				const ReferenceObjectID roid = elem->reference_object_id();

			//	get local position of DoF
				DimReferenceMapping<elemDim, dim>& map
					= ReferenceMappingProvider::get<elemDim, dim>(roid, vCornerCoords);
				MathVector<elemDim> locPos;
				VecSet(locPos, 0.5);
				map.global_to_local(locPos, x);

			//	evaluate at shapes at ip
				const LocalShapeFunctionSet<elemDim>& rTrialSpace =
						LocalFiniteElementProvider::get<elemDim>(roid, m_lfeID);
				std::vector<number> vShape;
				rTrialSpace.shapes(vShape, locPos);

			//	get multiindices of element
				std::vector<DoFIndex> ind;
				m_spGridFct->dof_indices(elem, m_fct, ind);

			// 	compute solution at integration point
				value = 0.0;
				for(size_t sh = 0; sh < vShape.size(); ++sh)
				{
					const number valSH = DoFRef(*m_spGridFct, ind[sh]);
					value += valSH * vShape[sh];
				}

			//	point is found
				return true;
			//}
			//UG_CATCH_THROW("GlobalGridFunctionNumberData: Evaluation failed."
			//			   << "Point: " << x << ", Element: "
			//			   << ElementDebugInfo(*m_spGridFct->domain()->grid(), elem));
		}

		/// evaluate value on all procs
		inline bool evaluate_global(number& value, const MathVector<dim>& x) const
		{
			// evaluate at this proc
			bool bFound = this->evaluate(value, x);

#ifdef UG_PARALLEL
			// share value between all procs
			pcl::ProcessCommunicator com;
			int numFound = (bFound ? 1 : 0);
			numFound = com.allreduce(numFound, PCL_RO_SUM);

			// get overall value
			if(!bFound) value = 0.0;
			number globValue = com.allreduce(value, PCL_RO_SUM) / numFound;

			if(numFound == 0)
				UG_THROW("Point "<<x<<" not found on all "<<pcl::NumProcs()<<" procs.");

			// check correctness for continuous spaces
			// note: if the point is found more than one it is located on the
			// boundary of some element. thus, if the space is continuous, those
			// values should match on all procs.
			if(bFound)
				if(LocalFiniteElementProvider::continuous(m_lfeID)){
					if( fabs(value) > 1e-10 && fabs((globValue - value) / value) > 1e-8)
						UG_THROW("For point " << x << " Global mean "<<globValue<<" != local value "<<value);
				}

			// set as global value
			value = globValue;
			bFound = true;
#endif

			if(!bFound)
				UG_THROW("Couldn't find an element containing the specified point: " << x);
			return bFound;
		}

		// evaluates at given position
		number evaluate_global(std::vector<number> vPos)
		{
			if((int)vPos.size() != dim)
				UG_THROW("Expected "<<dim<<" components, but given "<<vPos.size());

			MathVector<dim> x;
			for(int i = 0; i < dim; i++) x[i] = vPos[i];

			number value;
			evaluate_global(value, x);

			return value;
		}
};



template <typename TGridFunction>
class GlobalGridFunctionGradientData
	: public StdGlobPosData<GlobalGridFunctionGradientData<TGridFunction>, MathVector<TGridFunction::dim>, TGridFunction::dim>
{
	public:
	///	world dimension of grid function
		static const int dim = TGridFunction::dim;
		typedef typename TGridFunction::element_type  element_t;

		private:
	/// grid function
		SmartPtr<TGridFunction> m_spGridFct;

	///	component of function
		size_t m_fct;

	///	local finite element id
		LFEID m_lfeID;

		typedef lg_ntree<dim, dim, element_t>	tree_t;
		tree_t	m_tree;

	public:
	/// constructor
		GlobalGridFunctionGradientData(SmartPtr<TGridFunction> spGridFct, const char* cmp)
		: m_spGridFct(spGridFct),
		  m_tree(*spGridFct->domain()->grid(), spGridFct->domain()->position_attachment())
		{
			//this->set_functions(cmp);

			//	get function id of name
			m_fct = spGridFct->fct_id_by_name(cmp);

			//	check that function exists
			if(m_fct >= spGridFct->num_fct())
				UG_THROW("GridFunctionNumberData: Function space does not contain"
						" a function with name " << cmp << ".");

			//	local finite element id
			m_lfeID = spGridFct->local_finite_element_id(m_fct);



//			const MGSubsetHandler& ssh = *m_spGridFct->domain()->subset_handler();

			SubsetGroup ssGrp(m_spGridFct->domain()->subset_handler());
			ssGrp.add_all();

			std::vector<size_t> subsetsOfGridFunction;
			std::vector<element_t*> elemsWithGridFunctions;

			typename TGridFunction::const_element_iterator iterEnd, iter;


			for(size_t si = 0; si < ssGrp.size(); si++){
				if( spGridFct->is_def_in_subset(m_fct, si) )
					subsetsOfGridFunction.push_back(si);
			}


			for(size_t i = 0; i<subsetsOfGridFunction.size(); i++){
				size_t si = subsetsOfGridFunction[i];
				iter = spGridFct->template begin<element_t>(si);
				iterEnd = spGridFct->template end<element_t>(si);

				for(;iter!=iterEnd; ++iter){
					element_t *elem = *iter;
					elemsWithGridFunctions.push_back(elem);
				}
			}

//			m_tree.create_tree(spGridFct->template begin<element_t>(si),
//														spGridFct->template end<element_t>(si));

			m_tree.create_tree(elemsWithGridFunctions.begin(), elemsWithGridFunctions.end());

		};

		virtual ~GlobalGridFunctionGradientData() {}

		virtual bool continuous() const
		{
			return false;
		}

		///	to full-fill UserData-Interface
		inline void evaluate(MathVector<dim>& value, const MathVector<dim>& x, number time, int si) const
		{
			if(!evaluate(value, x))
				UG_THROW("For function "<<m_fct<<" couldn't find an element containing the specified point: " << x);
		}

		///	evaluates the data at a given point, returns false if point not found
		inline bool evaluate(MathVector<dim>& value, const MathVector<dim>& x) const
		{
			static const int refDim = dim;

			element_t* elem = NULL;
			try{

				if(!FindContainingElement(elem, m_tree, x)){
					return false;
				}

			//	get corners of element
				std::vector<MathVector<dim> > vCornerCoords;
				CollectCornerCoordinates(vCornerCoords, *elem, *m_spGridFct->domain());

			//	reference object id
				const ReferenceObjectID roid = elem->reference_object_id();

			//	get local position of DoF
				DimReferenceMapping<dim, dim>& map
					= ReferenceMappingProvider::get<dim, dim>(roid, vCornerCoords);
				MathVector<dim> locPos;
				VecSet(locPos, 0.5);
				map.global_to_local(locPos, x);

				
				MathMatrix<refDim, dim> JT;
				try{
					DimReferenceMapping<refDim, dim>& mapping
					= ReferenceMappingProvider::get<refDim, dim>(roid, vCornerCoords);

					//	compute transformation matrices
					mapping.jacobian_transposed(JT, x);

				}UG_CATCH_THROW("GlobalGridFunctionGradientData: failed.");

			//	evaluate at shapes at ip
				const LocalShapeFunctionSet<dim>& rTrialSpace =
						LocalFiniteElementProvider::get<dim>(roid, m_lfeID);
				std::vector<MathVector<refDim> > vLocGrad;
				rTrialSpace.grads(vLocGrad, locPos);


			//	Reference Mapping
				MathMatrix<dim, refDim> JTInv;
				RightInverse (JTInv, JT);

			//	get multiindices of element
				std::vector<DoFIndex> ind;
				m_spGridFct->dof_indices(elem, m_fct, ind);

			//	storage for shape function at ip
				MathVector<refDim> locGrad;

			//	compute grad at ip
				VecSet(locGrad, 0.0);
				for(size_t sh = 0; sh < vLocGrad.size(); ++sh)
				{
					const number valSH = DoFRef( *m_spGridFct, ind[sh]);
					VecScaleAppend(locGrad, valSH, vLocGrad[sh]);
				}

			// 	transform to global space
				MatVecMult(value, JTInv, locGrad);

			//	point is found
				return true;
			}
			UG_CATCH_THROW("GlobalGridFunctionGradientData: Evaluation failed."
						   << "Point: " << x << ", Element: "
						   << ElementDebugInfo(*m_spGridFct->domain()->grid(), elem));
		}

		/// evaluate value on all procs
		inline void evaluate_global(MathVector<dim>& value, const MathVector<dim>& x) const
		{
			// evaluate at this proc
			bool bFound = this->evaluate(value, x);

			// \todo: (optinal) check for evaluation on other procs

			if(!bFound)
				UG_THROW("Couldn't find an element containing the specified point: " << x);
		}

		// evaluates at given position
		std::vector<number> evaluate_global(std::vector<number> vPos)
		{
			if((int)vPos.size() != dim)
				UG_THROW("Expected "<<dim<<" components, but given "<<vPos.size());

			MathVector<dim> x;
			for(int i = 0; i < dim; i++) x[i] = vPos[i];

			MathVector<dim> value;
			evaluate_global(value, x);

			for(int i = 0; i < dim; i++) vPos[i] = value[i];
			return vPos;
		}
};

} // end namespace ug

#endif /* __H__UG__LIB_DISC__FUNCTION_SPACE__GRID_FUNCTION_GLOBAL_USER_DATA__ */
