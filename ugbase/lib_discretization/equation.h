/*
 * globaldiscretization.h
 *
 *  Created on: 04.05.2009
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__EQUATION__
#define __H__LIBDISCRETIZATION__EQUATION__

#include "differentialoperator.h"
#include "rhs.h"
#include "lib_grid/subset_handler.h"
#include "lib_algebra/lib_algebra.h"
#include "dofhandler.h"
#include "discretizationscheme.h"
#include <vector>
#include <string>

namespace ug {

template <int d>
class Equation
{
	public:
		Equation();
		Equation(std::string name, int _nr_func);

		void set_name(std::string name);
		std::string name();

		bool add_differentialoperator(TimeOperator<d>& op);
		bool add_differentialoperator(ScalarDifferentialOperator<d>& op);
		bool add_differentialoperator(DivergenzDifferentialOperator<d>& op);
		bool delete_differentialoperator(const int nr);
		bool clear_differentialoperator();

		bool add_rhs(RHS<d>& rhs);
		bool delete_rhs(const int nr);
		bool clear_rhs();

		bool add_dirichletBND(DirichletBNDCond<d>& cond)
		{
			m_DirichletBND = &cond;
			return true;
		}
		DirichletBNDCond<d>* get_dirichletBND()
		{
			return m_DirichletBND;
		}
		bool clear_dirichletBND()
		{
			m_DirichletBND = NULL;
			return true;
		}

		void print_info();

		bool set_discretzationscheme(DiscretizationSchemeID type);

		template <typename TElem>
		bool prepare_element(TElem* elem, NumericalSolution<d>& u, SubsetHandler& sh, int SubsetIndex)
		{
			static DiscretizationScheme<TElem, d>* DiscScheme = &DiscretizationSchemes<TElem, d>::DiscretizationScheme(m_DiscretizationSchemeID);

			DiscScheme->prepareElement(elem, u, _nr_func, sh, SubsetIndex);

			return true;
		}

		/**** instationary assemblings ****/
		template <typename TElem>
		bool assemble_defect(TElem* elem,Vector& vec, NumericalSolution<d>& u, number time, number s_m, number s_a);

		template <typename TElem>
		bool assemble_jacobian(TElem* elem, Matrix& mat, NumericalSolution<d>& u, number time, number s_m, number s_a);

		/**** stationary assemblings ****/
		template <typename TElem>
		bool assemble_defect(TElem* elem, Vector& vec, NumericalSolution<d>& u);

		template <typename TElem>
		bool assemble_jacobian(TElem* elem, Matrix& mat, NumericalSolution<d>& u);

		template <typename TElem>
		bool assemble_linear(TElem* elem, Matrix& mat, Vector& vec, NumericalSolution<d>& u);

		bool get_dirichlet_values(SubsetHandler& sh, uint subsetIndex, NumericalSolution<d>& u, DirichletValues<d>& dirVal);

	protected:
		typedef std::vector<DifferentialOperator*> DifferentialOperatorContainer;
		typedef std::vector<TimeOperator<d>*> TimeOperatorContainer;
		typedef std::vector<DivergenzDifferentialOperator<d>*> DivergenzDifferentialOperatorContainer;
		typedef std::vector<ScalarDifferentialOperator<d>*> ScalarDifferentialOperatorContainer;
		typedef std::vector<RHS<d>*> RHSContainer;

	protected:
		std::string m_name;
		DifferentialOperatorContainer m_DifferentialOperatorVector;
		TimeOperatorContainer m_TimeOperatorVector;
		ScalarDifferentialOperatorContainer m_ScalarDifferentialOperatorVector;
		DivergenzDifferentialOperatorContainer m_DivergenzDifferentialOperatorVector;
		RHSContainer m_RHSVector;
		DirichletBNDCond<d>* m_DirichletBND;
		DiscretizationSchemeID m_DiscretizationSchemeID;

		int _nr_func;
};


} /* end of namespace ug */

#include "equation_impl.h"

#endif /* __H__LIBDISCRETIZATION__EQUATION__ */

