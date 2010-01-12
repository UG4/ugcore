/*
 * problem.h
 *
 *  Created on: 25.06.2009
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__PROBLEM__
#define __H__LIBDISCRETIZATION__PROBLEM__

#include "systemdiscretization.h"
#include "numericalsolution.h"
#include "lib_algebra/lib_algebra.h"
#include <vector>
#include <string>

namespace ug {

template <int d>
class SubsetDiscretization {

	public:
		SubsetDiscretization();
		SubsetDiscretization(std::string name);

		void set_name(std::string name);
		std::string name();

		void add_Subset(SubsetHandler& sh, uint subsetIndex);
		void clear_Subset();

		void add_System(SystemDiscretization<d>& sys);
		void delete_System(int nr);
		void clear_Systems();
		int numberOfSystems();

		bool assemble_defect(Vector& vec, NumericalSolution<d>& u, number time, number s_m, number s_a)
		{
			if(m_sh == NULL) return false;
			bool b = true;
			b = b && assemble_defect<Triangle>(vec, u, time, s_m, s_a);
			b = b && assemble_defect<Quadrilateral>(vec, u, time, s_m, s_a);
			return b;
		}

		template <typename TElem>
		bool assemble_defect(Vector& vec, NumericalSolution<d>& u, number time, number s_m, number s_a)
		{
			typename geometry_traits<TElem>::iterator iterBegin, iterEnd, iter;

			iterBegin = m_sh->begin<TElem>(m_subsetIndex);
			iterEnd = m_sh->end<TElem>(m_subsetIndex);

			bool b = true;

			/* loop over all Geometric Objects of type "TElem" */
			for(iter = iterBegin; iter != iterEnd; iter++)
			{
				TElem *elem = *iter;

				// Create local Matrix

				for(unsigned int i=0; i<m_SystemDiscretization.size();i++)
				{
					b = b && (m_SystemDiscretization[i]->prepare_element(elem, u, *m_sh, m_subsetIndex));
					b = b && (m_SystemDiscretization[i]->assemble_defect(elem, vec, u, time, s_m, s_a));
				}
			}
			return b;
		}

		bool assemble_jacobian(Matrix& mat, NumericalSolution<d>& u, number time, number s_m, number s_a)
		{
			bool b = true;
			return b;
		}
		template <typename TElem>
		bool assemble_jacobian(Matrix& mat, NumericalSolution<d>& u, number time, number s_m, number s_a)
		{
			typename geometry_traits<TElem>::iterator iterBegin, iterEnd, iter;

			iterBegin = m_sh->begin<TElem>(m_subsetIndex);
			iterEnd = m_sh->end<TElem>(m_subsetIndex);

			bool b = true;

			/* loop over all Geometric Objects of type "TElem" */
			for(iter = iterBegin; iter != iterEnd; iter++)
			{
				TElem *elem = *iter;

				// Create local Matrix

				for(unsigned int i=0; i<m_SystemDiscretization.size();i++)
				{
					b = b && (m_SystemDiscretization[i]->prepare_element(elem, u, *m_sh));
					b = b && (m_SystemDiscretization[i]->assemble_jacobian(elem, mat, u, time, s_m, s_a));
				}
			}
			return b;
		}

		bool assemble_linear(Matrix& mat, Vector& rhs, NumericalSolution<d>& u)
		{
			if(m_sh == NULL) return false;
			bool b = true;
			b = b && assemble_linear<Triangle>(mat, rhs, u);
			b = b && assemble_linear<Quadrilateral>(mat, rhs, u);
			return b;
		}

		template <typename TElem>
		bool assemble_linear(Matrix& mat, Vector& rhs, NumericalSolution<d>& u)
		{
			typename geometry_traits<TElem>::iterator iterBegin, iterEnd, iter;

			iterBegin = m_sh->begin<TElem>(m_subsetIndex);
			iterEnd = m_sh->end<TElem>(m_subsetIndex);

			bool b = true;

			/* loop over all Geometric Objects of type "TElem" */
			for(iter = iterBegin; iter != iterEnd; iter++)
			{
				TElem *elem = *iter;

				// Create local Matrix

				for(unsigned int i=0; i<m_SystemDiscretization.size();i++)
				{
					b = b && (m_SystemDiscretization[i]->prepare_element(elem, u, *m_sh, m_subsetIndex));
					b = b && (m_SystemDiscretization[i]->assemble_linear(elem, mat, rhs, u));
				}
			}
			return b;
		}

		bool assemble_defect(Vector& vec, NumericalSolution<d>& u)
		{
			if(m_sh == NULL) return false;
			bool b = true;
			b = b && assemble_defect<Triangle>(vec, u);
			b = b && assemble_defect<Quadrilateral>(vec, u);
			return b;
		}

		template <typename TElem>
		bool assemble_defect(Vector& vec, NumericalSolution<d>& u)
		{
			typename geometry_traits<TElem>::iterator iterBegin, iterEnd, iter;

			iterBegin = m_sh->begin<TElem>(m_subsetIndex);
			iterEnd = m_sh->end<TElem>(m_subsetIndex);

			bool b = true;

			/* loop over all Geometric Objects of type "TElem" */
			for(iter = iterBegin; iter != iterEnd; iter++)
			{
				TElem *elem = *iter;

				// Create local Matrix

				for(unsigned int i=0; i<m_SystemDiscretization.size();i++)
				{
					b = b && (m_SystemDiscretization[i]->prepare_element(elem, u, *m_sh, m_subsetIndex));
					b = b && (m_SystemDiscretization[i]->assemble_defect(elem, vec, u));
				}
			}
			return b;
		}

		bool assemble_jacobian(Matrix& mat, NumericalSolution<d>& u)
		{
			if(m_sh == NULL) return false;
			bool b = true;
			b = b && assemble_jacobian<Triangle>(mat, u);
			b = b && assemble_jacobian<Quadrilateral>(mat, u);
			return b;
		}

		template <typename TElem>
		bool assemble_jacobian(Matrix& mat, NumericalSolution<d>& u)
		{
			typename geometry_traits<TElem>::iterator iterBegin, iterEnd, iter;

			iterBegin = m_sh->begin<TElem>(m_subsetIndex);
			iterEnd = m_sh->end<TElem>(m_subsetIndex);

			bool b = true;

			/* loop over all Geometric Objects of type "TElem" */
			for(iter = iterBegin; iter != iterEnd; iter++)
			{
				TElem *elem = *iter;

				// Create local Matrix

				for(unsigned int i=0; i<m_SystemDiscretization.size();i++)
				{
					b = b && (m_SystemDiscretization[i]->prepare_element(elem, u, *m_sh, m_subsetIndex));
					b = b && (m_SystemDiscretization[i]->assemble_jacobian(elem, mat, u));
				}
			}
			return b;
		}

		bool get_dirichlet_values(NumericalSolution<d>& u, DirichletValues<d>& dirVal)
		{
			bool b = true;

			for(uint i=0; i < m_SystemDiscretization.size(); ++i)
			{
				b = b && (m_SystemDiscretization[i]->get_dirichlet_values(*m_sh, m_subsetIndex, u, dirVal));
			}

			return b;
		}


		void print_info();

	protected:
		typedef std::vector<SystemDiscretization<d>*> SystemDiscretizationContainer;

	protected:
		std::string m_name;
		SystemDiscretizationContainer m_SystemDiscretization;
		SubsetHandler* m_sh;
		uint m_subsetIndex;
};



}

#endif /* __H__LIBDISCRETIZATION__PROBLEM__ */
