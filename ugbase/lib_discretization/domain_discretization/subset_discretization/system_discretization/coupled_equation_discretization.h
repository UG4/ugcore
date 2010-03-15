/*
 * systemdiscretization.h
 *
 *  Created on: 16.11.2009
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__SYSTEMDISCRETIZATION__
#define __H__LIBDISCRETIZATION__SYSTEMDISCRETIZATION__

#include "equation.h"
//#include "../lib_grid.h"
#include <typeinfo>
#include <string>

namespace ug {


// general implementation for system of equations
template <int d>
class EquationSystem : public SystemDiscretization<d> {

	public:
		bool prepare_element(GeometricObject* elem, NumericalSolution<d>& u, SubsetHandler& sh, int SubsetIndex)
		{
			if(typeid(*elem) == typeid(Triangle)) return prepare_element<Triangle>(dynamic_cast<Triangle*>(elem), u, sh, SubsetIndex);
			if(typeid(*elem) == typeid(Quadrilateral)) return prepare_element<Quadrilateral>(dynamic_cast<Quadrilateral*>(elem), u, sh, SubsetIndex);
			return false;
		}

		bool assemble_defect(GeometricObject* elem, Vector& vec, NumericalSolution<d>& u, number time, number s_m, number s_a)
		{
			if(typeid(*elem) == typeid(Triangle)) return assemble_defect<Triangle>(dynamic_cast<Triangle*>(elem), vec, u, time, s_m, s_a);
			if(typeid(*elem) == typeid(Quadrilateral)) return assemble_defect<Quadrilateral>(dynamic_cast<Quadrilateral*>(elem), vec, u, time, s_m, s_a);
			return false;
		}

		bool assemble_jacobian(GeometricObject* elem, Matrix& mat, NumericalSolution<d>& u, number time, number s_m, number s_a)
		{
			if(typeid(*elem) == typeid(Triangle)) return assemble_jacobian<Triangle>(dynamic_cast<Triangle*>(elem), mat, u, time, s_m, s_a);
			if(typeid(*elem) == typeid(Quadrilateral)) return assemble_jacobian<Quadrilateral>(dynamic_cast<Quadrilateral*>(elem), mat, u, time, s_m, s_a);
			return false;
		}

		bool assemble_defect(GeometricObject* elem, Vector& vec, NumericalSolution<d>& u)
		{
			if(typeid(*elem) == typeid(Triangle)) return assemble_defect<Triangle>(dynamic_cast<Triangle*>(elem), vec, u);
			if(typeid(*elem) == typeid(Quadrilateral)) return assemble_defect<Quadrilateral>(dynamic_cast<Quadrilateral*>(elem), vec, u);
			return false;
		}

		bool assemble_jacobian(GeometricObject* elem, Matrix& mat, NumericalSolution<d>& u)
		{
			if(typeid(*elem) == typeid(Triangle)) return assemble_jacobian<Triangle>(dynamic_cast<Triangle*>(elem), mat, u);
			if(typeid(*elem) == typeid(Quadrilateral)) return assemble_jacobian<Quadrilateral>(dynamic_cast<Quadrilateral*>(elem), mat, u);
			return false;
		}

		bool assemble_linear(GeometricObject* elem, Matrix& mat, Vector& rhs, NumericalSolution<d>& u)
		{
			if(typeid(*elem) == typeid(Triangle)) return assemble_linear<Triangle>(dynamic_cast<Triangle*>(elem), mat, rhs, u);
			if(typeid(*elem) == typeid(Quadrilateral)) return assemble_linear<Quadrilateral>(dynamic_cast<Quadrilateral*>(elem), mat, rhs, u);
			return false;
		}

		bool get_dirichlet_values(SubsetHandler& sh, uint subsetIndex, NumericalSolution<d>& u, DirichletValues<d>& dirVal)
		{
			bool b = true;
			for(uint i=0; i < m_Equation.size(); ++i)
			{
				b = b && (m_Equation[i]->get_dirichlet_values(sh, subsetIndex, u, dirVal));
			}
			return b;
		}

	public:
		EquationSystem();
		EquationSystem(std::string name);

		void set_name(std::string name);
		std::string name();

		void add_Equation(Equation<d>& eq);
		void delete_Equation(int nr);
		void clear_Equations();
		int numberOfEquations();

		void print_info();

	protected:
		template <typename TElem>
		inline bool prepare_element(TElem* elem, NumericalSolution<d>& u, SubsetHandler& sh, int SubsetIndex)
		{
			bool b = true;
			for(unsigned int i=0; i<m_Equation.size();i++)
			{
				b = b && (m_Equation[i]->prepare_element<TElem>(elem, u, sh, SubsetIndex));
			}
			return b;
		}

		template <typename TElem>
		inline bool assemble_defect(TElem* elem, Vector& vec, NumericalSolution<d>& u, number time, number s_m, number s_a)
		{
			bool b = true;
			for(unsigned int i=0; i<m_Equation.size();i++)
			{
				b = b && (m_Equation[i]->assemble_defect<TElem>(elem, vec, u, time, s_m, s_a));
			}
			return b;
		}

		template <typename TElem>
		inline bool assemble_jacobian(TElem* elem, Matrix& mat, NumericalSolution<d>& u, number time, number s_m, number s_a)
		{
			bool b = true;
			for(unsigned int i=0; i<m_Equation.size();i++)
			{
				b = b && (m_Equation[i]->assemble_jacobian<TElem>(elem, mat, u, time, s_m, s_a));
			}
			return b;
		}

		template <typename TElem>
		inline bool assemble_defect(TElem* elem, Vector& vec, NumericalSolution<d>& u)
		{
			bool b = true;
			for(unsigned int i=0; i<m_Equation.size();i++)
			{
				b = b && (m_Equation[i]->assemble_defect<TElem>(elem, vec, u));
			}
			return b;
		}

		template <typename TElem>
		inline bool assemble_jacobian(TElem* elem, Matrix& mat, NumericalSolution<d>& u)
		{
			bool b = true;
			for(unsigned int i=0; i<m_Equation.size();i++)
			{
				b = b && (m_Equation[i]->assemble_jacobian<TElem>(elem, mat, u));
			}
			return b;
		}

		template <typename TElem>
		inline bool assemble_linear(TElem* elem, Matrix& mat, Vector& rhs, NumericalSolution<d>& u)
		{
			bool b = true;
			for(unsigned int i=0; i<m_Equation.size();i++)
			{
				b = b && (m_Equation[i]->assemble_linear<TElem>(elem, mat, rhs, u));
			}
			return b;
		}

	protected:
		typedef std::vector<Equation<d>*> EquationContainer;

	protected:
		std::string m_name;
		EquationContainer m_Equation;

};








}



#endif /* __H__LIBDISCRETIZATION__SYSTEMDISCRETIZATION__ */
