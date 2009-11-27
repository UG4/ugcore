/*
 * spacialdiscretization.h
 *
 *  Created on: 22.10.2009
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__SPACIALDISCRETIZATION__
#define __H__LIBDISCRETIZATION__SPACIALDISCRETIZATION__

#include "subsetdiscretization.h"
#include "numericalsolution.h"
#include "lib_algebra/lib_algebra.h"
#include <vector>
#include <string>
#include "newton.h"

namespace ug{

class SpacialDiscretization : public NonLinearAssemble {

	public:
		SpacialDiscretization();
		SpacialDiscretization(std::string name);

		void set_name(std::string name);
		std::string name();

		void add_SubsetDiscretization(SubsetDiscretization& psd);
		void delete_SubsetDiscretization(int nr);
		void clear_SubsetDiscretizations();
		int numberOfSubsetDiscretizations();

		void print_info();

		bool assemble_defect(Vector& vec, NumericalSolution& u, number time, number s_m, number s_a)
		{
			bool b = true;

			for(uint i=0; i < m_SubsetDiscretization.size(); ++i)
			{
				b = b && (m_SubsetDiscretization[i]->assemble_defect(vec, u, time, s_m, s_a));
			}
			return b;
		}

		bool assemble_jacobian(Matrix& mat, NumericalSolution& u, number time, number s_m, number s_a)
		{
			bool b = true;

			for(uint i=0; i < m_SubsetDiscretization.size(); ++i)
			{
				b = b && (m_SubsetDiscretization[i]->assemble_jacobian(mat, u, time, s_m, s_a));
			}
			return b;
		}

		// only now, remove later.
		bool assemble_linear(Matrix& mat, Vector& rhs, NumericalSolution& u)
		{
			DirichletValues dirVal;
			bool b = true;

			for(uint i=0; i < m_SubsetDiscretization.size(); ++i)
			{
				b = b && (m_SubsetDiscretization[i]->assemble_linear(mat, rhs, u));
			}

			if(get_dirichlet_values(u, dirVal) == false) return false;
			dirVal.set_values(rhs);
			dirVal.set_rows(mat);

			return b;
		}

		bool assemble_defect(Vector& vec, NumericalSolution& u)
		{
			bool b = true;
			for(uint i=0; i < m_SubsetDiscretization.size(); ++i)
			{
				b = b && (m_SubsetDiscretization[i]->assemble_defect(vec, u));
			}

			DirichletValues dirVal;
			if(get_dirichlet_values(u, dirVal) == false) return false;
			dirVal.set_zero_values(vec);
			return b;
		}

		bool assemble_jacobian(Matrix& mat, NumericalSolution& u)
		{
			bool b = true;
			for(uint i=0; i < m_SubsetDiscretization.size(); ++i)
			{
				b = b && (m_SubsetDiscretization[i]->assemble_jacobian(mat, u));
			}

			DirichletValues dirVal;
			if(get_dirichlet_values(u, dirVal) == false) return false;
			dirVal.set_rows(mat);
			return b;
		}

		bool assemble_solution(NumericalSolution& u)
		{
			DirichletValues dirVal;

			if(get_dirichlet_values(u, dirVal) == false) return false;
			dirVal.set_values(*u.GridVector());
			return true;
		}

		bool get_dirichlet_values(NumericalSolution& u, DirichletValues& dirVal)
		{
			bool b = true;

			for(uint i=0; i < m_SubsetDiscretization.size(); ++i)
			{
				b = b && (m_SubsetDiscretization[i]->get_dirichlet_values(u, dirVal));
			}
			return b;
		}

	protected:
		typedef std::vector<SubsetDiscretization*> SubsetDiscretizationContainer;

	protected:
		std::string m_name;
		SubsetDiscretizationContainer m_SubsetDiscretization;

};






}



#endif /* __H__LIBDISCRETIZATION__SPACIALDISCRETIZATION__ */
