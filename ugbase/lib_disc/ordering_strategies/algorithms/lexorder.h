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

#ifndef __H__UG__LIB_DISC__ORDERING_STRATEGIES_ALGORITHMS__LEXORDER__
#define __H__UG__LIB_DISC__ORDERING_STRATEGIES_ALGORITHMS__LEXORDER__

#include <vector>
#include <utility> // for pair
#include <cstring> //strlen

#include "lib_disc/domain.h"
#include "lib_disc/function_spaces/grid_function.h"

#include "lib_algebra/ordering_strategies/algorithms/IOrderingAlgorithm.h"
#include "lib_algebra/ordering_strategies/algorithms/util.cpp"
#include "lib_disc/function_spaces/dof_position_util.h"

#include "common/error.h"

namespace ug{

template<int dim>
void ComputeLexicographicOrder(std::vector<size_t>& vNewIndex,
                               std::vector<std::pair<MathVector<dim>, size_t> >& vPos,
							   size_t orderDim = 0, bool increasing = true);

/// orders the dof distribution using lexicographic order
template <typename TDomain>
void OrderLexForDofDist(SmartPtr<DoFDistribution> dd, ConstSmartPtr<TDomain> domain,
	   size_t orderDim = 0, bool increasing = true);

/// orders the all DofDistributions of the ApproximationSpace using lexicographic order
template <typename TDomain>
void OrderLex(ApproximationSpace<TDomain>& approxSpace, const char *order);


template <typename TAlgebra, typename TDomain, typename O_t>
class LexOrdering : public IOrderingAlgorithm<TAlgebra, O_t>
{
public:
	typedef typename TAlgebra::matrix_type M_t;
	typedef typename TAlgebra::vector_type V_t;
	typedef IOrderingAlgorithm<TAlgebra, O_t> baseclass;

	/// Grid function type for the solution
	typedef GridFunction<TDomain, TAlgebra> GridFunc_t;

	/// Position attachment type
	typedef typename std::pair<MathVector<TDomain::dim>, size_t> Position_t;

	LexOrdering(){}

	/// clone constructor
	LexOrdering( const LexOrdering<TAlgebra, TDomain, O_t> &parent )
			: baseclass(){}

	SmartPtr<IOrderingAlgorithm<TAlgebra, O_t> > clone()
	{
		return make_sp(new LexOrdering<TAlgebra, TDomain, O_t>(*this));
	}

	void parse(){
		if(sizeof(m_dir) == 0){
			UG_THROW(name() << "::parse: empty string\n");
		}

		std::cout << name() << "::parse()" << std::endl;
		std::cout << "length: " << strlen(m_dir) << std::endl;
	}

	void compute(){
		size_t len = strlen(m_dir);

		if(len == 0){
			UG_THROW(name() << "::compute: Empty direction!");
		}

		size_t pos = 0;
		bool increasing = true;
		char sign;
		while(pos < len){
			if(increasing){
				sign = '+';
			}
			else{
				sign = '-';
			}

			switch(m_dir[pos]){
				case '+':
					//UG_LOG(name() << "::compute: LexOrdering in + direction.\n")
					++pos;
					increasing = true;
					break;
				case '-':
					//UG_LOG(name() << "::compute: LexOrdering in - direction.\n")
					++pos;
					increasing = false;
					break;
				case 'x':
					UG_LOG(name() << "::compute: LexOrdering in " << sign << "x direction.\n")
					ComputeLexicographicOrder<TDomain::dim>(o, m_vPositions, 0, increasing);
					++pos;
					increasing = true;
					break;
				case 'y':
					UG_LOG(name() << "::compute: LexOrdering in " << sign << "y direction.\n")
					ComputeLexicographicOrder<TDomain::dim>(o, m_vPositions, 1, increasing);
					++pos;
					increasing = true;
					break;
				case 'z':
					UG_LOG(name() << "::compute: LexOrdering in " << sign << "z direction.\n")
					ComputeLexicographicOrder<TDomain::dim>(o, m_vPositions, 2, increasing);
					++pos;
					increasing = true;
					break;
				default:
					UG_THROW(name() << "::compute: Invalid token in direction string, valid tokens: +, -, x, y, z");
			}
		}

		mat = NULL;
	}

	void check(){
		if(!is_permutation(o)){
			UG_THROW(name() << "::check: Not a permutation!");
		}
	}

	O_t& ordering(){
		return o;
	}

	void init(M_t* A, const V_t& V){
		if(strcmp(m_dir, "") == 0){
			UG_THROW(name() << "::init: no direction chosen!");
		}

		try{
			const GridFunc_t* pGridF;
			if((pGridF = dynamic_cast<const GridFunc_t*>(&V)) == 0){
				UG_THROW(name() << "::init: No DoFDistribution specified.");
			}

			SmartPtr<DoFDistribution> dd = ((GridFunc_t*) pGridF)->dof_distribution();

			size_t indN = dd->num_indices();

			if(indN != A->num_rows ()){
				UG_THROW(name() << "::init: #indices != #rows");
			}

			o.resize(indN);
			ExtractPositions(pGridF->domain(), dd, m_vPositions);
		}
		catch(...){
			throw;
		}

		#ifdef UG_ENABLE_DEBUG_LOGS
		UG_LOG("Using " << name() << " in " << m_dir << " direction\n");
		#endif

		mat = A;
	}

	void init(M_t*){
		UG_THROW(name() << "::init: Cannot initialize smoother without a geometry. Specify the 2nd argument for init!");
	}

	void init(M_t*, const V_t&, const O_t&){
		UG_THROW(name() << "::init: induced subgraph version not implemented yet!");
	}

	void init(M_t*, const O_t&){
		UG_THROW(name() << "::init: induced subgraph version not implemented yet!");
	}

	virtual const char* name() const {return "LexOrdering";}

	void set_direction(const char *dir){
		m_dir = dir;
	}

private:
	O_t o;
	M_t* mat;

	const char *m_dir;
	std::vector<Position_t> m_vPositions;
};

} // end namespace ug

#endif
