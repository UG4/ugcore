/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
 * Author: Ivo Muha
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

// extern headers
#include <iostream>
#include <sstream>
#include <string>
#include <limits>
#include <set>

// include bridge
#include "bridge/bridge.h"
#include "bridge/util.h"
#include "bridge/util_domain_algebra_dependent.h"

// lib_disc includes
#include "lib_disc/dof_manager/dof_distribution.h"
#include "lib_disc/function_spaces/grid_function.h"

// user data
#include "lib_disc/spatial_disc/user_data/user_data.h"

#include "lib_disc/reference_element/reference_element_traits.h"
#include "lib_disc/reference_element/reference_mapping_provider.h"
#include "lib_grid/algorithms/space_partitioning/lg_ntree.h"
#include "common/util/provider.h"
#include "lib_disc/domain_util.h"
#include "lib_disc/time_disc/time_integrator_observers/time_integrator_observer_interface.h"

// pcl includes
#ifdef UG_PARALLEL
	#include "pcl/pcl_util.h"
	#include "pcl/pcl_process_communicator.h"
#endif

#ifdef UG_PARALLEL
	#include "lib_grid/parallelization/distributed_grid.h"
#endif

using namespace std;

namespace ug{
namespace bridge{
namespace Evaluate{


template <typename TDomain, typename TAlgebra>
class NumberValuedUserDataEvaluator
{
	static const int dim = TDomain::dim;
	typedef GridFunction<TDomain, TAlgebra> TGridFunction;
	typedef typename TGridFunction::template dim_traits<dim>::grid_base_object TElem;	

	typedef lg_ntree<dim, dim, TElem>	tree_t;

	public:

		NumberValuedUserDataEvaluator(SmartPtr<UserData<number, dim> > userData) : m_userData(userData)
		{

		}

		number evaluateLua(const std::vector<number>& pos, SmartPtr<TGridFunction> u, number time)
		{
			number result = 0;

			this->evaluate(pos, result, u, time);

			return result;
		}


		bool evaluate(const std::vector<number>& pos,
									number& result,
									SmartPtr<TGridFunction> u,
									number time)
		{
			if(!m_initialized) 
				initialize(u);

			bool found = evaluateOnThisProcess(pos, result, u, time);
		

		#ifndef UG_PARALLEL
			return found;
		#endif

		#ifdef UG_PARALLEL
			// share value between all procs
			pcl::ProcessCommunicator com;
			int numFound = (found ? 1 : 0);
			numFound = com.allreduce(numFound, PCL_RO_SUM);

			// get overall value
			// if found on more than one processor, the data will be the same
			number globalResult = com.allreduce(result, PCL_RO_SUM);

			if(numFound == 0)
				return false;

			// set as result			
			result = globalResult / numFound;
				

			return true;
		#endif
		}


	private:

		bool evaluateOnThisProcess(const std::vector<number>& pos,
									number& result,
									SmartPtr<TGridFunction> u,
									number time)
		{
			TElem* elem = NULL;

			MathVector<dim> globalPosition;

			for(int i = 0; i < dim; i++)
			{
				globalPosition[i] = pos[i];
			}

			if(!FindContainingElement(elem, *m_tree, globalPosition))
			{
				return false;
			}

			//	get corners of element
			std::vector<MathVector<dim> > vCornerCoords;
			CollectCornerCoordinates(vCornerCoords, *elem, *u->domain());

			//	get subset
			int si = u->domain()->subset_handler()->get_subset_index(elem);

			//	reference object id
			const ReferenceObjectID roid = elem->reference_object_id();

			//	get local position of DoF
			DimReferenceMapping<dim, dim>& map
				= ReferenceMappingProvider::get<dim, dim>(roid, vCornerCoords);
			MathVector<dim> locPos;
			VecSet(locPos, 0.5);
			map.global_to_local(locPos, globalPosition);

			// storage for the result
			number value;

			//	get local solution if needed
			if(m_userData->requires_grid_fct())
			{
				//	create storage
				LocalIndices ind;
				LocalVector localU;

				// 	get global indices
				u->indices(elem, ind);

				// 	adapt local algebra
				localU.resize(ind);

				// 	read local values of u
				GetLocalVector(localU, *u);
							
				try
				{
					(*m_userData)(&value, &globalPosition, time, si, elem,
							&vCornerCoords[0], &locPos, 1, &localU, NULL);
				}
				UG_CATCH_THROW("NumberValuedUserDataEvaluator: Cannot evaluate data.");
			}
			else
			{
				try
				{
					(*m_userData)(&value, &globalPosition, time, si, elem,
							&vCornerCoords[0], &locPos, 1, NULL, NULL);
				}
				UG_CATCH_THROW("NumberValuedUserDataEvaluator: Cannot evaluate data.");
			}

			result = value;			

			return true;
		}

		void initialize(SmartPtr<TGridFunction> u)
		{
			m_initialized = true;

			m_tree = make_sp(new tree_t(*u->domain()->grid(), u->domain()->position_attachment()));
			m_tree->create_tree(u->template begin<TElem>(), u->template end<TElem>());
		}

		bool m_initialized = false;
		SmartPtr<tree_t> m_tree;
		SmartPtr<UserData<number, dim> > m_userData;

};

template <typename TDomain, typename TAlgebra>
class VectorValuedUserDataEvaluator
{
	static const int dim = TDomain::dim;
	typedef GridFunction<TDomain, TAlgebra> TGridFunction;
	typedef typename TGridFunction::template dim_traits<dim>::grid_base_object TElem;	

	typedef lg_ntree<dim, dim, TElem>	tree_t;

	public:

		VectorValuedUserDataEvaluator(SmartPtr<UserData<MathVector<dim>, dim> > userData) : m_userData(userData)
		{

		}

		std::vector<number> evaluateLua(const std::vector<number>& pos, SmartPtr<TGridFunction> u, number time, bool interpolateOverNeighbouringTriangles)
		{
			std::vector<number> result;

			this->evaluate(pos, result, u, time, interpolateOverNeighbouringTriangles);

			return result;
		}


		bool evaluate(const std::vector<number>& pos,
									std::vector<number>& result,
									SmartPtr<TGridFunction> u,
									number time,
									bool interpolateOverNeighbouringTriangles)
		{
			if(!m_initialized) 
				initialize(u);

			bool found = false;

			if(interpolateOverNeighbouringTriangles)
				found = evaluateOnThisProcessNeighbouring(pos, result, u, time);
			else
				found = evaluateOnThisProcess(pos, result, u, time);
		

		#ifndef UG_PARALLEL
			return found;
		#endif

		#ifdef UG_PARALLEL
			// share value between all procs
			pcl::ProcessCommunicator com;
			int numFound = (found ? 1 : 0);
			numFound = com.allreduce(numFound, PCL_RO_SUM);

			// get overall value
			// if found on more than one processor, the data will be the same
			std::vector<number> globalResult;
			com.allreduce(result, globalResult, PCL_RO_SUM);

			if(numFound == 0)
				return false;

			// set as result
			for(int i = 0; i < dim; i++)
			{
				result[i] = globalResult[i] / numFound;
			}	

			return true;
		#endif
		}


	private:

		bool evaluateOnThisProcessNeighbouring(const std::vector<number>& pos,
									std::vector<number>& result,
									SmartPtr<TGridFunction> u,
									number time)
		{
			typedef typename TDomain::grid_type TGrid;
			typedef typename Grid::traits<TElem>::secure_container TElemContainer;

			result.resize(TDomain::dim);

			TElem* elem = NULL;

			MathVector<dim> globalPosition;

			for(int i = 0; i < dim; i++)
			{
				globalPosition[i] = pos[i];
			}

			if(!FindContainingElement(elem, *m_tree, globalPosition))
			{
				return false;
			}

			// find all adjacent faces/volumes/edges
			std::set<TElem* > elements;
			elements.insert(elem);

			TGrid* grid = u->domain()->grid().get();

			for(size_t i = 0; i < elem->num_vertices(); i++)
			{
				// for each vertex of the found element: get associated faces/volumes/edges
				Vertex* vrt = elem->vertex(i);

				TElemContainer associatedElements;
				grid->associated_elements(associatedElements, vrt);

				for(size_t j = 0; j < associatedElements.size(); j++)
				{
					elements.insert(associatedElements[j]);
				}
			}

			double sumWeights = 0;

			for(TElem * element : elements)
			{
				//	get corners of element
				std::vector<MathVector<dim> > vCornerCoords;
				CollectCornerCoordinates(vCornerCoords, *element, *u->domain());

				// find center of element
				MathVector<dim> centerGlobal;
				for(size_t i = 0; i < vCornerCoords.size(); i++)
				{
					for(int d = 0; d < dim; d++)
						centerGlobal[d] += vCornerCoords[i][d] / vCornerCoords.size();					
				}

				// calculate distance to evaluated position
				double distance = 0;
				for(int d = 0; d < dim; d++)
					distance += (pos[d]-centerGlobal[d])*(pos[d]-centerGlobal[d]);				
				distance = sqrt(distance);

				// calculate weight
				double weight = exp(-distance);
				sumWeights += weight;
				
				//	reference object id
				const ReferenceObjectID roid = element->reference_object_id();
				DimReferenceElement<dim> refElem = ReferenceElementProvider::get<dim>(roid);

				// calculate the elements center in local coordinates
				MathVector<dim> centerLocal;
				for(size_t i = 0; i < refElem.corners()->size(); i++)
				{
					for(int d = 0; d < dim; d++)
						centerLocal[d] += refElem.corners()[i][d] / refElem.corners()->size();
				}

				//	get subset
				int si = u->domain()->subset_handler()->get_subset_index(element);

				// storage for the result
				MathVector<dim> value;

				//	get local solution if needed
				if(m_userData->requires_grid_fct())
				{
					//	create storage
					LocalIndices ind;
					LocalVector localU;

					// 	get global indices
					u->indices(element, ind);

					// 	adapt local algebra
					localU.resize(ind);

					// 	read local values of u
					GetLocalVector(localU, *u);
								
					try
					{
						(*m_userData)(&value, &centerGlobal, time, si, element,
								&vCornerCoords[0], &centerLocal, 1, &localU, NULL);
					}
					UG_CATCH_THROW("VectorValuedUserDataEvaluator: Cannot evaluate data.");
				}
				else
				{
					try
					{
						(*m_userData)(&value, &centerGlobal, time, si, element,
								&vCornerCoords[0], &centerLocal, 1, NULL, NULL);
					}
					UG_CATCH_THROW("VectorValuedUserDataEvaluator: Cannot evaluate data.");
				}

				// add the vector with the calculated weight
				for(int i = 0; i < dim; i++)
				{
					result[i] = value[i] * weight;
				}
			}

			// scale vector by sumWeights
			for(int i = 0; i < dim; i++)
			{
				result[i] = result[i] / sumWeights;
			}			

			return true;
		}

		bool evaluateOnThisProcess(const std::vector<number>& pos,
									std::vector<number>& result,
									SmartPtr<TGridFunction> u,
									number time)
		{
			result.resize(TDomain::dim);

			TElem* elem = NULL;

			MathVector<dim> globalPosition;

			for(int i = 0; i < dim; i++)
			{
				globalPosition[i] = pos[i];
			}

			if(!FindContainingElement(elem, *m_tree, globalPosition))
			{
				return false;
			}

			//	get corners of element
			std::vector<MathVector<dim> > vCornerCoords;
			CollectCornerCoordinates(vCornerCoords, *elem, *u->domain());

			//	get subset
			int si = u->domain()->subset_handler()->get_subset_index(elem);

			//	reference object id
			const ReferenceObjectID roid = elem->reference_object_id();

			//	get local position of DoF
			DimReferenceMapping<dim, dim>& map
				= ReferenceMappingProvider::get<dim, dim>(roid, vCornerCoords);
			MathVector<dim> locPos;
			VecSet(locPos, 0.5);
			map.global_to_local(locPos, globalPosition);

			// storage for the result
			MathVector<dim> value;

			//	get local solution if needed
			if(m_userData->requires_grid_fct())
			{
				//	create storage
				LocalIndices ind;
				LocalVector localU;

				// 	get global indices
				u->indices(elem, ind);

				// 	adapt local algebra
				localU.resize(ind);

				// 	read local values of u
				GetLocalVector(localU, *u);
							
				try
				{
					(*m_userData)(&value, &globalPosition, time, si, elem,
							&vCornerCoords[0], &locPos, 1, &localU, NULL);
				}
				UG_CATCH_THROW("VectorValuedUserDataEvaluator: Cannot evaluate data.");
			}
			else
			{
				try
				{
					(*m_userData)(&value, &globalPosition, time, si, elem,
							&vCornerCoords[0], &locPos, 1, NULL, NULL);
				}
				UG_CATCH_THROW("VectorValuedUserDataEvaluator: Cannot evaluate data.");
			}

			for(int i = 0; i < dim; i++)
			{
				result[i] = value[i];
			}

			return true;
		}

		void initialize(SmartPtr<TGridFunction> u)
		{
			m_initialized = true;

			m_tree = make_sp(new tree_t(*u->domain()->grid(), u->domain()->position_attachment()));
			m_tree->create_tree(u->template begin<TElem>(), u->template end<TElem>());
		}

		bool m_initialized = false;
		SmartPtr<tree_t> m_tree;
		SmartPtr<UserData<MathVector<dim>, dim> > m_userData;

};

template <typename TDomain, typename TAlgebra>
class PointEvaluatorBase : public ITimeIntegratorObserver<TDomain, TAlgebra>
{
	typedef std::vector<number> TPoint;	
	static const int dim = TDomain::dim;
	typedef GridFunction<TDomain, TAlgebra> TGridFunction;

	public:
		void add_evaluation_point(TPoint point)
		{
			UG_COND_THROW(m_initialized, "PointEvaluator: Can't add point, Evaluator in use.");

			m_evaluationPoints.push_back(point);
		}		
		
		void set_filename(std::string filename)
		{
			UG_COND_THROW(m_initialized, "PointEvaluator: Can't change filename, Evaluator in use.");

			m_filename = filename;
		}	

		void set_separator(std::string separator)
		{
			UG_COND_THROW(m_initialized, "PointEvaluator: Can't change separator, Evaluator in use.");

			m_separator = separator;
		}

		virtual bool step_process(SmartPtr<TGridFunction> uNew, int step, number time, number dt) override
		{
			if(!m_initialized)
				initialize();

			std::fstream output_file;

			#ifdef UG_PARALLEL
            if(pcl::ProcRank() != 0)     
			{
				// passing an unopenend file object. nothing will be written!
				write_evaluations(output_file, uNew, step, time, dt);
				return true;
			}
            #endif

			output_file.open(m_filename.c_str(), std::fstream::app | std::fstream::out);
			
			UG_COND_THROW(output_file.fail(), "PointEvaluator: Can't open output file.");
			
            output_file.precision(15);
			write_evaluations(output_file, uNew, step, time, dt);
			
			output_file.close();

			return true;
		}

	protected:

		void initialize()
		{
			m_initialized = true;

			#ifdef UG_PARALLEL
            if(pcl::ProcRank() != 0)     
			{
				return;
			}
            #endif

			std::fstream output_file;
			output_file.open(m_filename.c_str(), std::fstream::trunc | std::fstream::out);
			
			UG_COND_THROW(output_file.fail(), "PointEvaluator: Can't open output file.");

			write_header(output_file);

			output_file.close();
		}

		virtual void write_header(std::ostream& output) {}
		virtual void write_evaluations(std::ostream& output, SmartPtr<TGridFunction> uNew, int step, number time, number dt) {}

		void write_point(std::ostream& output, const TPoint& point)
		{
			output << "{";
			for(int i = 0; i < dim; i++)
			{
				if(i > 0)
					output << ",";
				output << point[i];
			}
			output << "}";
		}

		bool m_initialized = false;

		std::vector<TPoint> m_evaluationPoints;
		std::string m_filename;
		std::string m_separator = "\t";
};

template <typename TDomain, typename TAlgebra>
class VectorValuedUserDataPointEvaluator : public PointEvaluatorBase<TDomain, TAlgebra>
{
	static const int dim = TDomain::dim;
	typedef GridFunction<TDomain, TAlgebra> TGridFunction;

	public:
		VectorValuedUserDataPointEvaluator(SmartPtr<UserData<MathVector<dim>, dim> > userData) 
		: m_evaluator(userData)
		{}	

		void set_interpolate_over_neighbouring_elements(bool interpolate)
		{
			m_interpolateOverNeighbouringTriangles = interpolate;
		}	

	protected:

		virtual void write_evaluations(std::ostream& output, SmartPtr<TGridFunction> uNew, int step, number time, number dt) override
		{
			UG_LOG(" * Write Vector-valued Position Data to '" << this->m_filename << "' ... \n");
			output << time << this->m_separator;
			std::vector<number> result;
			for (auto point : this->m_evaluationPoints)
			{
				result.clear();
				if(m_evaluator.evaluate(point, result, uNew, time, m_interpolateOverNeighbouringTriangles))
				{
					for(int d = 0; d < dim; d++)
					{
						output << result[d] << this->m_separator;
					}
				}
				else
				{
					for(int d = 0; d < dim; d++)
					{
						output << "NaN" << this->m_separator;
					}
				}

			}
			output << "\n";
		}

		virtual void write_header(std::ostream& output) override
		{
			output << "# Position Evaluating file - vector valued\n";
			char axis[3] = { 'x', 'y', 'z'};

			output << "time" << this->m_separator;

			for (auto point : this->m_evaluationPoints)
			{
				for(int d = 0; d < dim; d++)
				{
					this->write_point(output, point);
					output << "-" << axis[d] << this->m_separator;
				}
			}

			output << "\n";
		}

	private:
		VectorValuedUserDataEvaluator<TDomain, TAlgebra> m_evaluator;
		bool m_interpolateOverNeighbouringTriangles = false;

};

template <typename TDomain, typename TAlgebra>
class NumberValuedUserDataPointEvaluator : public PointEvaluatorBase<TDomain, TAlgebra>
{
	static const int dim = TDomain::dim;
	typedef GridFunction<TDomain, TAlgebra> TGridFunction;

	public:
		NumberValuedUserDataPointEvaluator(SmartPtr<UserData<number, dim> > userData) : m_evaluator(userData)
		{}		

	protected:

		virtual void write_evaluations(std::ostream& output, SmartPtr<TGridFunction> uNew, int step, number time, number dt) override
		{
			UG_LOG(" * Write Number-valued Position Data to '" << this->m_filename << "' ... \n");
			output << time << this->m_separator;
			number result;
			for (auto point : this->m_evaluationPoints)
			{
				if(m_evaluator.evaluate(point, result, uNew, time))
				{					
					output << result << this->m_separator;					
				}
				else
				{					
					output << "NaN" << this->m_separator;					
				}

			}
			output << "\n";
		}

		virtual void write_header(std::ostream& output) override
		{
			output << "# Position Evaluating file - number valued\n";
			output << "time" << this->m_separator;

			for (auto point : this->m_evaluationPoints)
			{
				this->write_point(output, point);					
			}

			output << "\n";
		}

	private:
		NumberValuedUserDataEvaluator<TDomain, TAlgebra> m_evaluator;
};

template <typename TDomain, typename TAlgebra>
SmartPtr<PointEvaluatorBase<TDomain,TAlgebra> > GetUserDataPointEvaluator(SmartPtr<UserData<MathVector<TDomain::dim>, TDomain::dim> > userData)
{
	return make_sp(new VectorValuedUserDataPointEvaluator<TDomain, TAlgebra>(userData));
}

template <typename TDomain, typename TAlgebra>
SmartPtr<PointEvaluatorBase<TDomain,TAlgebra> > GetUserDataPointEvaluator(SmartPtr<UserData<number, TDomain::dim> > userData)
{
	return make_sp(new NumberValuedUserDataPointEvaluator<TDomain, TAlgebra>(userData));
}

template <typename TDomain>
bool CloseVertexExists(const MathVector<TDomain::dim>& globPos,
					   TDomain* dom,
					   const char* subsets,
					   SmartPtr<typename TDomain::subset_handler_type> sh,
					   number maxDist)
{
//	domain type
	typedef TDomain domain_type;
	typedef typename domain_type::grid_type grid_type;
	typedef typename domain_type::subset_handler_type subset_handler_type;
// get position accessor
	grid_type* grid = dom->grid().get();

	const typename domain_type::position_accessor_type& aaPos
										= dom->position_accessor();

	typename subset_handler_type::template traits<Vertex>::const_iterator iterEnd, iter;
	number minDistanceSq = numeric_limits<number>::max();

	#ifdef UG_PARALLEL
		DistributedGridManager* dgm = grid->distributed_grid_manager();
	#endif

	SubsetGroup ssGrp(sh);
	if(subsets != NULL)
		ssGrp.add(TokenizeString(subsets));
	else
		ssGrp.add_all();

	for(size_t i = 0; i < ssGrp.size(); ++i)
	{
	//	get subset index
		const int si = ssGrp[i];
	// 	iterate over all elements
		for(size_t lvl = 0; lvl < sh->num_levels(); ++lvl){
			iterEnd = sh->template end<Vertex>(si, lvl);
			iter = sh->template begin<Vertex>(si, lvl);
			for(; iter != iterEnd; ++iter)
			{
			//	get element
			//todo: replace most of the following checks by a spGridFct->contains(...)

				Vertex* vrt = *iter;
				if(grid->has_children(vrt)) continue;

				#ifdef UG_PARALLEL
					if(dgm->is_ghost(vrt))	continue;
					if(dgm->contains_status(vrt, INT_H_SLAVE)) continue;
				#endif
			//	global position
				number buffer = VecDistanceSq(globPos, aaPos[vrt]);
				if(buffer < minDistanceSq)
				{
					minDistanceSq = buffer;
				}
			}
		}
	}
	return 	minDistanceSq < sq(maxDist);
}

/**
 * \defgroup interpolate_bridge Interpolation Bridge
 * \ingroup disc_bridge
 * \{
 */

template <typename TGridFunction>
number EvaluateAtVertex(const MathVector<TGridFunction::dim>& globPos,
						SmartPtr<TGridFunction> spGridFct,
						size_t fct,
						const SubsetGroup& ssGrp,
						typename TGridFunction::domain_type::subset_handler_type* sh,
						bool minimizeOverAllProcs = false)
{

//	domain type
	typedef typename TGridFunction::domain_type domain_type;
	typedef typename domain_type::grid_type grid_type;
	typedef typename domain_type::subset_handler_type subset_handler_type;
// get position accessor
	domain_type* dom = spGridFct->domain().get();
	grid_type* grid = dom->grid().get();

	subset_handler_type* domSH = dom->subset_handler().get();

	const typename domain_type::position_accessor_type& aaPos
										= dom->position_accessor();

	std::vector<DoFIndex> ind;
	Vertex* chosen = NULL;

	typename subset_handler_type::template traits<Vertex>::const_iterator iterEnd, iter;
	number minDistanceSq = std::numeric_limits<double>::max();

	#ifdef UG_PARALLEL
		DistributedGridManager* dgm = grid->distributed_grid_manager();
	#endif

	for(size_t i = 0; i < ssGrp.size(); ++i)
	{
	//	get subset index
		const int si = ssGrp[i];

	// 	iterate over all elements
		for(size_t lvl = 0; lvl < sh->num_levels(); ++lvl){
			iterEnd = sh->template end<Vertex>(si, lvl);
			iter = sh->template begin<Vertex>(si, lvl);
			for(; iter != iterEnd; ++iter)
			{
			//	get element
			//todo: replace most of the following checks by a spGridFct->contains(...)
				Vertex* vrt = *iter;
				if(grid->has_children(vrt)) continue;

				#ifdef UG_PARALLEL
					if(dgm->is_ghost(vrt))	continue;
					if(dgm->contains_status(vrt, INT_H_SLAVE)) continue;
				#endif

				int domSI = domSH->get_subset_index(vrt);

			//	skip if function is not defined in subset
				if(!spGridFct->is_def_in_subset(fct, domSI)) continue;

			//	global position
				number buffer = VecDistanceSq(globPos, aaPos[vrt]);
				if(buffer < minDistanceSq)
				{
					minDistanceSq = buffer;
					chosen = *iter;
				}
			}
		}
	}

	// get corresponding value (if vertex found, otherwise take 0)
	number value = 0.0;
	if (chosen)
	{
		spGridFct->inner_dof_indices(chosen, fct, ind);
		value = DoFRef(*spGridFct, ind[0]);
	}

	// in parallel environment, find global minimal distance and corresponding value
#ifdef UG_PARALLEL
	if (minimizeOverAllProcs)
	{
		pcl::MinimalKeyValuePairAcrossAllProcs<number, number, std::less<number> >(minDistanceSq, value);

		// check that a vertex has been found
		UG_COND_THROW(minDistanceSq == std::numeric_limits<double>::max(),
		"No vertex of given subsets could be located on any process.");

		return value;
	}
#endif

	// check that a vertex has been found
	UG_COND_THROW(!chosen, "No vertex of given subsets could be located.");

	return value;
}


template <typename TGridFunction>
number EvaluateAtClosestVertex(const MathVector<TGridFunction::dim>& pos,
							  SmartPtr<TGridFunction> spGridFct, const char* cmp,
							  const char* subsets,
							  SmartPtr<typename TGridFunction::domain_type::subset_handler_type> sh)
{
//	get function id of name
	const size_t fct = spGridFct->fct_id_by_name(cmp);

//	check that function found
	if(fct > spGridFct->num_fct())
		UG_THROW("Evaluate: Name of component '"<<cmp<<"' not found.");


//	create subset group
	SubsetGroup ssGrp(sh);
	if(subsets != NULL)
	{
		ssGrp.add(TokenizeString(subsets));
	}
	else
	{
	//	add all subsets and remove lower dim subsets afterwards
		ssGrp.add_all();
	}

	return EvaluateAtVertex<TGridFunction>(pos, spGridFct, fct, ssGrp, sh.get());
}


template <typename TGridFunction>
number EvaluateAtClosestVertexAllProcs
(
	const MathVector<TGridFunction::dim>& pos,
	SmartPtr<TGridFunction> spGridFct,
	const char* cmp,
	const char* subsets,
	SmartPtr<typename TGridFunction::domain_type::subset_handler_type> sh
)
{
	// get function id of name
	const size_t fct = spGridFct->fct_id_by_name(cmp);

	// check that function found
	if (fct > spGridFct->num_fct())
		UG_THROW("Evaluate: Name of component '"<<cmp<<"' not found.");


	// create subset group
	SubsetGroup ssGrp(sh);
	if (subsets != NULL)
		ssGrp.add(TokenizeString(subsets));
	else
	{
	//	add all subsets and remove lower dim subsets afterwards
		ssGrp.add_all();
	}

	return EvaluateAtVertex<TGridFunction>(pos, spGridFct, fct, ssGrp, sh.get(), true);
}


/**
 * Class exporting the functionality. All functionality that is to
 * be used in scripts or visualization must be registered here.
 */
struct Functionality
{

/**
 * Function called for the registration of Domain and Algebra dependent parts.
 * All Functions and Classes depending on both Domain and Algebra
 * are to be placed here when registering. The method is called for all
 * available Domain and Algebra types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <typename TDomain, typename TAlgebra>
static void DomainAlgebra(Registry& reg, string grp)
{
	string suffix = GetDomainAlgebraSuffix<TDomain,TAlgebra>();
	string tag = GetDomainAlgebraTag<TDomain,TAlgebra>();

//	typedef
	typedef ug::GridFunction<TDomain, TAlgebra> TFct;

	{
	//	reg.add_function("Integral", static_cast<number (*)(SmartPtr<UserData<number,dim> >, SmartPtr<TFct>, const char*, number, int)>(&Integral<TFct>), grp, "Integral", "Data#GridFunction#Subsets#Time#QuadOrder");

		//reg.add_function("EvaluateAtClosestVertex", static_cast<number (*)(const std::vector<number>&, SmartPtr<TFct>, const char*, const char*)>(&EvaluateAtClosestVertex<TFct>),grp, "Evaluate_at_closest_vertex", "Position#GridFunction#Component#Subsets");
		reg.add_function("EvaluateAtClosestVertex",
						 &EvaluateAtClosestVertex<TFct>,
						 grp, "Evaluate_at_closest_vertex", "Position#GridFunction#Component#Subsets#SubsetHandler");
		reg.add_function("EvaluateAtClosestVertexAllProcs",
						 &EvaluateAtClosestVertexAllProcs<TFct>,
						 grp, "Evaluate_at_closest_vertex", "Position#GridFunction#Component#Subsets#SubsetHandler");

	}


	{
		typedef PointEvaluatorBase<TDomain, TAlgebra> T;
		typedef ITimeIntegratorObserver<TDomain, TAlgebra> TBase;
		
		string name = string("PointEvaluatorBase").append(suffix);

		reg.add_class_<T, TBase>(name, grp)
					   	.add_method("add_evaluation_point", &T::add_evaluation_point, "point", "")
					   	.add_method("set_filename", &T::set_filename, "filename", "")
					   	.add_method("set_separator", &T::set_separator, "separator", "")
						.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "PointEvaluatorBase", tag);
	}

	{
		typedef VectorValuedUserDataPointEvaluator<TDomain, TAlgebra> T;
		typedef PointEvaluatorBase<TDomain, TAlgebra> TBase;

		string name = string("VectorValuedUserDataPointEvaluator").append(suffix);

		reg.add_class_<T, TBase>(name, grp)
					   	.template add_constructor<void (*)(SmartPtr<UserData<MathVector<TDomain::dim>, TDomain::dim> >) >("")					   
						.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "VectorValuedUserDataPointEvaluator", tag);
	}

	{
		typedef VectorValuedUserDataEvaluator<TDomain, TAlgebra> T;
		string name = string("VectorValuedUserDataEvaluator").append(suffix);

		reg.add_class_<T>(name, grp)
					   	.template add_constructor<void (*)(SmartPtr<UserData<MathVector<TDomain::dim>, TDomain::dim> >) >("")
					   	.add_method("evaluate", &T::evaluateLua, "point#result#solution#time", "")
						.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "VectorValuedUserDataEvaluator", tag);
	}

	{
		typedef NumberValuedUserDataPointEvaluator<TDomain, TAlgebra> T;
		typedef PointEvaluatorBase<TDomain, TAlgebra> TBase;

		string name = string("NumberValuedUserDataPointEvaluator").append(suffix);

		reg.add_class_<T, TBase>(name, grp)
					   	.template add_constructor<void (*)(SmartPtr<UserData<number, TDomain::dim> >) >("")					   
						.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "NumberValuedUserDataPointEvaluator", tag);
	}

	{
		typedef NumberValuedUserDataEvaluator<TDomain, TAlgebra> T;
		string name = string("NumberValuedUserDataEvaluator").append(suffix);

		reg.add_class_<T>(name, grp)
					   	.template add_constructor<void (*)(SmartPtr<UserData<number, TDomain::dim> >) >("")
					   	.add_method("evaluate", &T::evaluateLua, "point#result#solution#time", "")
						.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "NumberValuedUserDataEvaluator", tag);
	}
	{

		reg.add_function("GetUserDataPointEvaluator", static_cast<SmartPtr<PointEvaluatorBase<TDomain,TAlgebra> > (*)(SmartPtr<UserData<MathVector<TDomain::dim>, TDomain::dim> >)>(&GetUserDataPointEvaluator<TDomain,TAlgebra>), 
			grp, 
			"GetUserDataPointEvaluator", 
			"UserDataObject");
		
		reg.add_function("GetUserDataPointEvaluator", static_cast<SmartPtr<PointEvaluatorBase<TDomain,TAlgebra> > (*)(SmartPtr<UserData<number, TDomain::dim> >)>(&GetUserDataPointEvaluator<TDomain,TAlgebra>), 
			grp, 
			"GetUserDataPointEvaluator", 
			"UserDataObject");
		
	}
}

/**
 * Function called for the registration of Domain dependent parts.
 * All Functions and Classes depending on the Domain
 * are to be placed here when registering. The method is called for all
 * available Domain types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <typename TDomain>
static void Domain(Registry& reg, string grp)
{
	string suffix = GetDomainSuffix<TDomain>();
	string tag = GetDomainTag<TDomain>();
	reg.add_function("CloseVertexExists", &CloseVertexExists<TDomain>, grp);
}

/**
 * Function called for the registration of Dimension dependent parts.
 * All Functions and Classes depending on the Dimension
 * are to be placed here when registering. The method is called for all
 * available Dimension types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <int dim>
static void Dimension(Registry& reg, string grp)
{
	string suffix = GetDimensionSuffix<dim>();
	string tag = GetDimensionTag<dim>();

}

/**
 * Function called for the registration of Algebra dependent parts.
 * All Functions and Classes depending on Algebra
 * are to be placed here when registering. The method is called for all
 * available Algebra types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <typename TAlgebra>
static void Algebra(Registry& reg, string grp)
{
	string suffix = GetAlgebraSuffix<TAlgebra>();
	string tag = GetAlgebraTag<TAlgebra>();

}

/**
 * Function called for the registration of Domain and Algebra independent parts.
 * All Functions and Classes not depending on Domain and Algebra
 * are to be placed here when registering.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
static void Common(Registry& reg, string grp)
{
}

}; // end Functionality

// end group
/// \}

}// namespace Evaluate

///
void RegisterBridge_Evaluate(Registry& reg, string grp)
{
	grp.append("/Evaluate");
	typedef Evaluate::Functionality Functionality;

	try{
//		RegisterCommon<Functionality>(reg,grp);
//		RegisterDimensionDependent<Functionality>(reg,grp);
		RegisterDomainDependent<Functionality>(reg,grp);
//		RegisterAlgebraDependent<Functionality>(reg,grp);
		RegisterDomainAlgebraDependent<Functionality>(reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

}//	end of namespace bridge
}//	end of namespace ug
