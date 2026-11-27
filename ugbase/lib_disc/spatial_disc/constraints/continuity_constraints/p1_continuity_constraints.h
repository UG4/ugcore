/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__CONSTRAINTS__CONTINUITY_CONSTRAINTS__P1_CONTINUITY_CONSTRAINTS__
#define __H__UG__LIB_DISC__SPATIAL_DISC__CONSTRAINTS__CONTINUITY_CONSTRAINTS__P1_CONTINUITY_CONSTRAINTS__

#include "lib_disc/assemble_interface.h"
#include "lib_disc/spatial_disc/constraints/constraint_interface.h"
#include "lib_grid/algorithms/geom_obj_util/vertex_util.h"

namespace ug {

/// returns the vertices of the object constraining a hanging vertex
void CollectConstraining(std::vector<Vertex*>& vConstrainingVrt,
						 const Grid& grid,
                         ConstrainedVertex* hgVrt,
                         bool bClearContainer = true);


template <typename TDomain, typename TAlgebra>
class SymP1Constraints
	: public IDomainConstraint<TDomain, TAlgebra>
{
	public:
	// 	Algebra type
		typedef TAlgebra algebra_type;

	// 	Type of algebra matrix
		typedef typename algebra_type::matrix_type matrix_type;

	// 	Type of algebra vector
		typedef typename algebra_type::vector_type vector_type;

	public:
		SymP1Constraints() : IDomainConstraint<TDomain, TAlgebra>() {}
		virtual ~SymP1Constraints() {}

		virtual int type() const {return CT_HANGING;}

		void adjust_jacobian(matrix_type& J, const vector_type& u,
		                     ConstSmartPtr<DoFDistribution> dd, int type, number time = 0.0,
                             ConstSmartPtr<VectorTimeSeries<vector_type> > vSol = NULL,
   							 const number s_a0 = 1.0);

		void adjust_defect(vector_type& d, const vector_type& u,
		                   ConstSmartPtr<DoFDistribution> dd, int type, number time = 0.0,
                           ConstSmartPtr<VectorTimeSeries<vector_type> > vSol = NULL,
						   const std::vector<number>* vScaleMass = NULL,
                           const std::vector<number>* vScaleStiff = NULL);

		void adjust_rhs(vector_type& rhs, const vector_type& u,
		                ConstSmartPtr<DoFDistribution> dd, int type, number time = 0.0);

		void adjust_linear(matrix_type& mat, vector_type& rhs, const vector_type& u,
		                   ConstSmartPtr<DoFDistribution> dd, int type, number time = 0.0);

		void adjust_solution(vector_type& u, ConstSmartPtr<DoFDistribution> dd,
							 int type, number time = 0.0);

		void adjust_prolongation(matrix_type& P,
								 ConstSmartPtr<DoFDistribution> ddFine,
								 ConstSmartPtr<DoFDistribution> ddCoarse,
								 int type,
								 number time = 0.0);

		void adjust_restriction(matrix_type& R,
								ConstSmartPtr<DoFDistribution> ddCoarse,
								ConstSmartPtr<DoFDistribution> ddFine,
								int type,
								number time = 0.0);

		virtual void adjust_correction
		(	vector_type& u,
			ConstSmartPtr<DoFDistribution> dd,
			int type,
			number time = 0.0
		);
};



template <typename TDomain, typename TAlgebra>
class OneSideP1Constraints
	: public IDomainConstraint<TDomain, TAlgebra>
{
	public:
	// 	Algebra type
		typedef TAlgebra algebra_type;

	// 	Type of algebra matrix
		typedef typename algebra_type::matrix_type matrix_type;

	// 	Type of algebra vector
		typedef typename algebra_type::vector_type vector_type;

	protected:
		typedef IDomainConstraint<TDomain, TAlgebra> base_type;

	public:
		OneSideP1Constraints() : IDomainConstraint<TDomain, TAlgebra>() {}
		virtual ~OneSideP1Constraints() {}

		virtual int type() const {return CT_HANGING;}

		void adjust_jacobian(matrix_type& J, const vector_type& u,
		                     ConstSmartPtr<DoFDistribution> dd, int type, number time = 0.0,
                             ConstSmartPtr<VectorTimeSeries<vector_type> > vSol = NULL,
							 const number s_a0 = 1.0);

		void adjust_defect(vector_type& d, const vector_type& u,
		                   ConstSmartPtr<DoFDistribution> dd, int type, number time = 0.0,
                           ConstSmartPtr<VectorTimeSeries<vector_type> > vSol = NULL,
						   const std::vector<number>* vScaleMass = NULL,
                           const std::vector<number>* vScaleStiff = NULL);

		void adjust_rhs(vector_type& rhs, const vector_type& u,
		                ConstSmartPtr<DoFDistribution> dd, int type, number time = 0.0);

		void adjust_linear(matrix_type& mat, vector_type& rhs, const vector_type& u,
		                   ConstSmartPtr<DoFDistribution> dd, int type, number time = 0.0);

		void adjust_solution(vector_type& u, ConstSmartPtr<DoFDistribution> dd,
							 int type, number time = 0.0);

		void adjust_prolongation(matrix_type& P,
								 ConstSmartPtr<DoFDistribution> ddFine,
								 ConstSmartPtr<DoFDistribution> ddCoarse,
								 int type,
								 number time = 0.0);

		void adjust_restriction(matrix_type& R,
								ConstSmartPtr<DoFDistribution> ddCoarse,
								ConstSmartPtr<DoFDistribution> ddFine,
								int type,
								number time = 0.0);

		virtual void adjust_correction
		(	vector_type& u,
			ConstSmartPtr<DoFDistribution> dd,
			int type,
			number time = 0.0
		);
};

}; // namespace ug

#include "p1_continuity_constraints_impl.h"

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__CONSTRAINTS__CONTINUITY_CONSTRAINTS__P1_CONTINUITY_CONSTRAINTS__ */
