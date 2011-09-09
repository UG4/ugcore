/*
 * navier_stokes_bnd.h
 *
 *  Created on: 01.09.2011
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__NAVIER_STOKES__NAVIER_STOKES_BND__
#define __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__NAVIER_STOKES__NAVIER_STOKES_BND__

#include "lib_discretization/spatial_discretization/disc_item.h"

#include "lib_discretization/spatial_discretization/elem_disc/neumann_boundary/neumann_boundary.h"
#include "lib_discretization/spatial_discretization/constraints/dirichlet_boundary/lagrange_dirichlet_boundary.h"

namespace ug{

template <	typename TDomain,
			typename TDoFDistribution,
			typename TAlgebra>
class NavierStokesInflow
	: public IDiscretizationItem<TDomain,TDoFDistribution,TAlgebra>
{
	private:
		const static int dim = TDomain::dim;
		typedef boost::function<void (MathVector<dim>& value, const MathVector<dim>& x, number time)> VectorFunctor;

	public:
	///	returns the number of element discs
		virtual size_t num_elem_disc() const {return 1;}

	///	returns the element disc
		virtual IDomainElemDisc<TDomain>* get_elem_disc(size_t i) {return &m_NeumannDisc;}

	///	returns the number of constraints
		virtual size_t num_constraint() const {return 1;}

	///	returns an element disc
		virtual IConstraint<TDoFDistribution, TAlgebra>* get_constraint(size_t i) {return &m_DirichletConstraint;}

	public:
	///	sets the symbolic names for velocity and pressure
		bool set_functions(const char* functions)
		{
			std::string strFunctions(functions);
			std::vector<std::string> tokens;

			TokenizeString(strFunctions, tokens, ',');

			if(tokens.size() != TDomain::dim + 1)
			{
				UG_LOG("ERROR in 'NavierStokesInflow::set_functions': This Boundary "
						"Condition works on exactly dim+1 (velocity+pressure) "
						"components, but "<<tokens.size()<<"components given.\n");
				return false;
			}

			m_velNames.clear();
			for(int i=0;i<TDomain::dim; ++i)
			{
				if(i>0) m_velNames.append(",");
				m_velNames.append(tokens[i]);
			}

			m_pressureName = tokens[TDomain::dim];

			return true;
		}

	///	sets the velocity to a given value
		bool add(VectorFunctor& user, const char* subsetsBND)
		{
			if(m_velNames.empty() || m_pressureName.empty())
			{
				UG_LOG("ERROR in 'NavierStokesInflow::add': Symbolic names for"
						" velocity and pressure not set. Please set them first.\n");
				return false;
			}
			m_NeumannDisc.add(user, m_pressureName.c_str(), subsetsBND);
			m_DirichletConstraint.add(user, m_velNames.c_str(), subsetsBND);

			return true;
		}

	///	sets the subsets reguarded as inner
		void set_subsets(const char* subsetInner){m_NeumannDisc.set_subsets(subsetInner);}

	///	sets the approximation space
		void set_approximation_space(IApproximationSpace<TDomain>& approxSpace)
		{
			m_DirichletConstraint.set_approximation_space(approxSpace);
		}

	protected:
	///	neumann disc for pressure equation
		FV1NeumannBoundaryElemDisc<TDomain> m_NeumannDisc;

	///	dirichlet disc for velocity components
		LagrangeDirichletBoundary<TDomain,TDoFDistribution,TAlgebra> m_DirichletConstraint;

	///	name of velocity components
		std::string m_velNames;

	///	name of pressure components
		std::string m_pressureName;
};

template <	typename TDomain,
			typename TDoFDistribution,
			typename TAlgebra>
class NavierStokesWall
	: public IDiscretizationItem<TDomain,TDoFDistribution,TAlgebra>
{
	public:
	///	returns the number of element discs
		virtual size_t num_elem_disc() const {return 0;}

	///	returns the element disc
		virtual IDomainElemDisc<TDomain>* get_elem_disc(size_t i) {return NULL;}

	///	returns the number of constraints
		virtual size_t num_constraint() const {return 1;}

	///	returns an element disc
		virtual IConstraint<TDoFDistribution, TAlgebra>* get_constraint(size_t i) {return &m_DirichletConstraint;}

	///	virtual destructor
		~NavierStokesWall() {}

	public:
	///	sets the symbolic names for velocity and pressure
		bool set_functions(const char* functions)
		{
			std::string strFunctions(functions);

			m_velNames.clear();
			TokenizeString(strFunctions, m_velNames, ',');

			if(m_velNames.size() != TDomain::dim + 1)
			{
				UG_LOG("ERROR in 'NavierStokesWall::set_functions': This Boundary "
						"Condition works on exactly dim+1 (velocity+pressure) "
						"components, but "<<m_velNames.size()<<"components given.\n");
				return false;
			}

			return true;
		}

	///	sets the velocity to a given value
		bool add(const char* subsetsBND)
		{
			if(m_velNames.empty())
			{
				UG_LOG("ERROR in 'NavierStokesWall::add': Symbolic names for"
						" velocity and pressure not set. Please set them first.\n");
				return false;
			}

			for(int i = 0; i < TDomain::dim; ++i)
			{
				m_DirichletConstraint.add(0.0, m_velNames[i].c_str(), subsetsBND);
			}

			return true;
		}

	///	sets the approximation space
		void set_approximation_space(IApproximationSpace<TDomain>& approxSpace)
		{
			m_DirichletConstraint.set_approximation_space(approxSpace);
		}

	protected:
	///	dirichlet disc for velocity components
		LagrangeDirichletBoundary<TDomain,TDoFDistribution,TAlgebra> m_DirichletConstraint;

	///	name of velocity components
		std::vector<std::string> m_velNames;
};

} // end namespace ug

#endif /* __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__NAVIER_STOKES__NAVIER_STOKES_BND__ */
