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

// interface
class SystemDiscretization {
public:
	virtual bool prepare_element(GeometricObject* elem, NumericalSolution& u, SubsetHandler& sh, int SubsetIndex) = 0;

	virtual bool assemble_defect(GeometricObject* elem, Vector& vec, NumericalSolution& u, number time, number s_m, number s_a) = 0;
	virtual bool assemble_jacobian(GeometricObject* elem, Matrix& mat, NumericalSolution& u, number time, number s_m, number s_a) = 0;
	virtual bool assemble_defect(GeometricObject* elem, Vector& vec, NumericalSolution& u) = 0;
	virtual bool assemble_jacobian(GeometricObject* elem, Matrix& mat, NumericalSolution& u) = 0;
	virtual bool assemble_linear(GeometricObject* elem, Matrix& mat, Vector& rhs, NumericalSolution& u) = 0;

	virtual bool get_dirichlet_values(SubsetHandler& sh, uint subsetIndex, NumericalSolution& u, DirichletValues& dirVal) = 0;

	virtual ~SystemDiscretization(){}

};

// general implementation for system of equations
class EquationSystem : public SystemDiscretization {

	public:
		bool prepare_element(GeometricObject* elem, NumericalSolution& u, SubsetHandler& sh, int SubsetIndex)
		{
			if(typeid(*elem) == typeid(Triangle)) return prepare_element<Triangle>(dynamic_cast<Triangle*>(elem), u, sh, SubsetIndex);
			if(typeid(*elem) == typeid(Quadrilateral)) return prepare_element<Quadrilateral>(dynamic_cast<Quadrilateral*>(elem), u, sh, SubsetIndex);
			return false;
		}

		bool assemble_defect(GeometricObject* elem, Vector& vec, NumericalSolution& u, number time, number s_m, number s_a)
		{
			if(typeid(*elem) == typeid(Triangle)) return assemble_defect<Triangle>(dynamic_cast<Triangle*>(elem), vec, u, time, s_m, s_a);
			if(typeid(*elem) == typeid(Quadrilateral)) return assemble_defect<Quadrilateral>(dynamic_cast<Quadrilateral*>(elem), vec, u, time, s_m, s_a);
			return false;
		}

		bool assemble_jacobian(GeometricObject* elem, Matrix& mat, NumericalSolution& u, number time, number s_m, number s_a)
		{
			if(typeid(*elem) == typeid(Triangle)) return assemble_jacobian<Triangle>(dynamic_cast<Triangle*>(elem), mat, u, time, s_m, s_a);
			if(typeid(*elem) == typeid(Quadrilateral)) return assemble_jacobian<Quadrilateral>(dynamic_cast<Quadrilateral*>(elem), mat, u, time, s_m, s_a);
			return false;
		}

		bool assemble_defect(GeometricObject* elem, Vector& vec, NumericalSolution& u)
		{
			if(typeid(*elem) == typeid(Triangle)) return assemble_defect<Triangle>(dynamic_cast<Triangle*>(elem), vec, u);
			if(typeid(*elem) == typeid(Quadrilateral)) return assemble_defect<Quadrilateral>(dynamic_cast<Quadrilateral*>(elem), vec, u);
			return false;
		}

		bool assemble_jacobian(GeometricObject* elem, Matrix& mat, NumericalSolution& u)
		{
			if(typeid(*elem) == typeid(Triangle)) return assemble_jacobian<Triangle>(dynamic_cast<Triangle*>(elem), mat, u);
			if(typeid(*elem) == typeid(Quadrilateral)) return assemble_jacobian<Quadrilateral>(dynamic_cast<Quadrilateral*>(elem), mat, u);
			return false;
		}

		bool assemble_linear(GeometricObject* elem, Matrix& mat, Vector& rhs, NumericalSolution& u)
		{
			if(typeid(*elem) == typeid(Triangle)) return assemble_linear<Triangle>(dynamic_cast<Triangle*>(elem), mat, rhs, u);
			if(typeid(*elem) == typeid(Quadrilateral)) return assemble_linear<Quadrilateral>(dynamic_cast<Quadrilateral*>(elem), mat, rhs, u);
			return false;
		}

		bool get_dirichlet_values(SubsetHandler& sh, uint subsetIndex, NumericalSolution& u, DirichletValues& dirVal)
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

		void add_Equation(Equation& eq);
		void delete_Equation(int nr);
		void clear_Equations();
		int numberOfEquations();

		void print_info();

	protected:
		template <typename TElem>
		inline bool prepare_element(TElem* elem, NumericalSolution& u, SubsetHandler& sh, int SubsetIndex)
		{
			bool b = true;
			for(unsigned int i=0; i<m_Equation.size();i++)
			{
				b = b && (m_Equation[i]->prepare_element<TElem, APosition>(elem, u, sh, SubsetIndex));
			}
			return b;
		}

		template <typename TElem>
		inline bool assemble_defect(TElem* elem, Vector& vec, NumericalSolution& u, number time, number s_m, number s_a)
		{
			bool b = true;
			for(unsigned int i=0; i<m_Equation.size();i++)
			{
				b = b && (m_Equation[i]->assemble_defect<TElem, APosition>(elem, vec, u, time, s_m, s_a));
			}
			return b;
		}

		template <typename TElem>
		inline bool assemble_jacobian(TElem* elem, Matrix& mat, NumericalSolution& u, number time, number s_m, number s_a)
		{
			bool b = true;
			for(unsigned int i=0; i<m_Equation.size();i++)
			{
				b = b && (m_Equation[i]->assemble_jacobian<TElem, APosition>(elem, mat, u, time, s_m, s_a));
			}
			return b;
		}

		template <typename TElem>
		inline bool assemble_defect(TElem* elem, Vector& vec, NumericalSolution& u)
		{
			bool b = true;
			for(unsigned int i=0; i<m_Equation.size();i++)
			{
				b = b && (m_Equation[i]->assemble_defect<TElem, APosition>(elem, vec, u));
			}
			return b;
		}

		template <typename TElem>
		inline bool assemble_jacobian(TElem* elem, Matrix& mat, NumericalSolution& u)
		{
			bool b = true;
			for(unsigned int i=0; i<m_Equation.size();i++)
			{
				b = b && (m_Equation[i]->assemble_jacobian<TElem, APosition>(elem, mat, u));
			}
			return b;
		}

		template <typename TElem>
		inline bool assemble_linear(TElem* elem, Matrix& mat, Vector& rhs, NumericalSolution& u)
		{
			bool b = true;
			for(unsigned int i=0; i<m_Equation.size();i++)
			{
				b = b && (m_Equation[i]->assemble_linear<TElem, APosition>(elem, mat, rhs, u));
			}
			return b;
		}

	protected:
		typedef std::vector<Equation*> EquationContainer;

	protected:
		std::string m_name;
		EquationContainer m_Equation;

};








}



#endif /* __H__LIBDISCRETIZATION__SYSTEMDISCRETIZATION__ */
