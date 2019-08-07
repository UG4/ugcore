/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
 * Author: Susanne Höllbacher
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

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__LOCAL_TO_GLOBAL_MAPPER__
#define __H__UG__LIB_DISC__SPATIAL_DISC__LOCAL_TO_GLOBAL_MAPPER__

// extern headers
#include <vector>

// intern headers
#include "lib_disc/dof_manager/dof_distribution.h"
#include "lib_disc/time_disc/solution_time_series.h"


namespace ug{

/// interface for definition of special LocalToGlobal mappings
/**
 * \tparam	TAlgebra			type of Algebra
 */
template <typename TAlgebra>
class ILocalToGlobalMapper
{
	public:
	///	Algebra type
		typedef TAlgebra algebra_type;

	///	Type of algebra matrix
		typedef typename algebra_type::matrix_type matrix_type;

	///	Type of algebra vector
		typedef typename algebra_type::vector_type vector_type;

	public:
	///	default Constructor
		ILocalToGlobalMapper() {};

	///	send local entries to global matrix
		virtual void add_local_vec_to_global(vector_type& vec, const LocalVector& lvec,
				ConstSmartPtr<DoFDistribution> dd) = 0;

	///	send local entries to global rhs
		virtual void add_local_mat_to_global(matrix_type& mat, const LocalMatrix& lmat,
				ConstSmartPtr<DoFDistribution> dd) = 0;

    /// sets certain entries to Dirichlet DoFs (in elem_disc_assembla_util)
        virtual void adjust_mat_global(matrix_type& mat, const LocalMatrix& lmat,
                                   ConstSmartPtr<DoFDistribution> dd) = 0;
    
    ///	modifies local solution vector for adapted defect computation ( in elem_disc_assemble_util)
        virtual void modify_LocalData(LocalMatrix& locJ, LocalVector& locU,
                                  ConstSmartPtr<DoFDistribution> dd) = 0;
        virtual void modify_LocalData(LocalVectorTimeSeries& uT, LocalMatrix& locJ, LocalVector& locU,
                                  ConstSmartPtr<DoFDistribution> dd) = 0;
    
        virtual void modify_LocalData(LocalVector& locD, LocalVector& tmpLocD, LocalVector& locU,
                                  ConstSmartPtr<DoFDistribution> dd) = 0;
        virtual void modify_LocalData(LocalVectorTimeSeries& uT, LocalVector& locD, LocalVector& tmpLocD, LocalVector& locU,
                                  ConstSmartPtr<DoFDistribution> dd, size_t t) = 0;
    
    ///	modifies global solution vector for adapted defect computation (in domain_disc_impl)
        virtual void modify_GlobalSol(vector_type& vecMod, const vector_type& vec, ConstSmartPtr<DoFDistribution> dd) = 0;
    
        virtual void modify_GlobalSol(SmartPtr<VectorTimeSeries<vector_type> > vSolMod,
                                  ConstSmartPtr<VectorTimeSeries<vector_type> > vSol, ConstSmartPtr<DoFDistribution> dd) = 0;

	///	virtual destructor
		virtual ~ILocalToGlobalMapper() {};

};


} // end namespace ug

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__LOCAL_TO_GLOBAL_MAPPER__*/
