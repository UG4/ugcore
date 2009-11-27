/*
 * globaldiscretization.cpp
 *
 *  Created on: 04.05.2009
 *      Author: andreasvogel
 */



#include <vector>
#include "equation.h"
#include "differentialoperator.h"
#include "lib_grid/lib_grid.h"
#include "lib_algebra/lib_algebra.h"

namespace ug {

Equation::Equation()
{
	m_name = "No Name";
}

Equation::Equation(std::string name, int nr_func)
{
	m_name = name;
	_nr_func = nr_func;
}

void Equation::set_name(std::string name)
{
	m_name = name;
}

std::string Equation::name()
{
	return m_name;
}

bool Equation::add_differentialoperator(TimeOperator& op)
{
	m_DifferentialOperatorVector.push_back(&op);
	m_TimeOperatorVector.push_back(&op);

	return true;
};


bool Equation::add_differentialoperator(DivergenzDifferentialOperator& op)
{
	m_DifferentialOperatorVector.push_back(&op);
	m_DivergenzDifferentialOperatorVector.push_back(&op);

	return true;
};

bool Equation::add_differentialoperator(ScalarDifferentialOperator& op)
{
	m_DifferentialOperatorVector.push_back(&op);
	m_ScalarDifferentialOperatorVector.push_back(&op);

	return true;
};

bool Equation::delete_differentialoperator(const int nr)
{
	DifferentialOperator* toDelete = m_DifferentialOperatorVector[nr];

	m_DifferentialOperatorVector.erase(m_DifferentialOperatorVector.begin() + nr);

	for(int i; i < m_TimeOperatorVector.size(); ++i)
	{
		if(toDelete == m_TimeOperatorVector[i])
			m_TimeOperatorVector.erase(m_TimeOperatorVector.begin() + i);
	}
	for(int i; i < m_ScalarDifferentialOperatorVector.size(); ++i)
	{
		if(toDelete == m_ScalarDifferentialOperatorVector[i])
			m_ScalarDifferentialOperatorVector.erase(m_ScalarDifferentialOperatorVector.begin() + i);
	}
	for(int i; i < m_DivergenzDifferentialOperatorVector.size(); ++i)
	{
		if(toDelete == m_DivergenzDifferentialOperatorVector[i])
			m_DivergenzDifferentialOperatorVector.erase(m_DivergenzDifferentialOperatorVector.begin() + i);
	}

	return true;
}


bool Equation::clear_differentialoperator()
{
	m_DifferentialOperatorVector.clear();
	m_TimeOperatorVector.clear();
	m_ScalarDifferentialOperatorVector.clear();
	m_DivergenzDifferentialOperatorVector.clear();

	return true;
}

bool Equation::add_rhs(RHS& rhs)
{
	m_RHSVector.push_back(&rhs);

	return true;
};

bool Equation::delete_rhs(const int nr)
{
	m_RHSVector.erase(m_RHSVector.begin() + nr);

	return true;
}


bool Equation::clear_rhs()
{
	m_RHSVector.clear();

	return true;
}

void Equation::print_info()
{
	std::cout << "Equation \""<< this->name() << "\" contains the following Terms:" << std::endl;

	std::cout << "Differentialoperator:" << std::endl;
	for(unsigned int i = 0; i < m_DifferentialOperatorVector.size(); i++)
	{
		std::cout << i << ":" << m_DifferentialOperatorVector[i]->name() << std::endl;
	}
	std::cout << std::endl;
	std::cout << "RHS:" << std::endl;
	for(unsigned int i = 0; i < m_RHSVector.size(); i++)
	{
		std::cout << i << ":" << m_RHSVector[i]->name() << std::endl;
	}
	std::cout << std::endl;
}

bool Equation::set_discretzationscheme(DiscretizationSchemeID type)
{
	m_DiscretizationSchemeID = type;

	return true;
}


bool Equation::get_dirichlet_values(SubsetHandler& sh, uint subsetIndex, NumericalSolution& u, DirichletValues& dirVal)
{
	switch(m_DiscretizationSchemeID)
	{
	case DISC_SCHEME_FE1:
		if(dirVal.add_dirichlet_nodes(u, _nr_func, m_DirichletBND, sh, subsetIndex) == false) return false;
		break;
	default:
		printf("ERROR in Equation::assemble: Can not find Discretization Scheme\n");
		return false;
	}

	return true;
}

};
