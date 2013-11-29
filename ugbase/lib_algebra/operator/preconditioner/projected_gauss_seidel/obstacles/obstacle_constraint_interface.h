/*
 * obstacle_constraint_interface.h
 *
 *  Created on: 25.11.2013
 *      Author: raphaelprohl
 */

#ifndef __H__UG__LIB_ALGEBRA__OPERATOR__PRECONDITIONER__PROJECTED_GAUSS_SEIDEL__OBSTACLE_CONSTRAINT_INTERFACE__
#define __H__UG__LIB_ALGEBRA__OPERATOR__PRECONDITIONER__PROJECTED_GAUSS_SEIDEL__OBSTACLE_CONSTRAINT_INTERFACE__

#include "lib_disc/common/multi_index.h"
#include "lib_disc/function_spaces/grid_function.h"

using namespace std;

namespace ug{

/// Interface for Obstacle Constraints
/**
 *  The Interface for Obstacle Constraints provides the framework to define obstacle constraints
 *  of the form
 *
 * 			c(u) <= upObs 		(cf. 'set_upper_obstacle')
 *
 * 	and
 *
 * 			c(u) >= lowObs 		(cf. 'set_lower_obstacle')
 *
 * 	where u is the solution vector. Here, 'upObs' and 'lowObs' are user-defined functions,
 * 	which need to be of the same size as the function of unknowns u. The obstacle function
 * 	c is defined by creating a derived class of this interface and by specializing the way
 * 	the correction should be computed (cf. correction_for_lower_obs, correction_for_upper_obs,
 * 	correction_for_lower_and_upper_obs).
 *
 * 	One simple example is the ScalarObstacle-class. Such an obstacle functions can be used in
 * 	combination with projected preconditioners. They should be passed to the preconditioner
 * 	by 'IProjPreconditioner::set_obstacle_constraint'.
 */
template <typename TAlgebra>
class IObstacleConstraint
{
	public:
	///	Algebra type
		typedef TAlgebra algebra_type;

	///	Matrix type
		typedef typename algebra_type::matrix_type matrix_type;

	///	Vector type
		typedef typename algebra_type::vector_type vector_type;

	///	Value type
		typedef typename vector_type::value_type value_type;

	public:
	/// constructor for an obstacle defined on some subset(s)
		template <typename TDomain>
		IObstacleConstraint(const GridFunction<TDomain, TAlgebra>& u, const char* subsets):
			m_bLowerObs(false), m_bUpperObs(false)
		{
			m_vActiveIndicesLow.resize(0);
			m_vActiveIndicesUp.resize(0);

			m_spDD = u.dof_distribution();
			m_ssName = subsets;
			init(u);
		};

	/// constructor
		IObstacleConstraint(): m_bLowerObs(false), m_bUpperObs(false)
		{
			m_vActiveIndicesLow.resize(0);
			m_vActiveIndicesUp.resize(0);
			m_vIndicesOfObsSubsets.resize(0);
		};

	///	set constraint/obstacle g_low (for c(u) >= g_low)
		void set_lower_obstacle(SmartPtr<vector_type> lowObs){
			m_spVecOfLowObsValues = lowObs; m_bLowerObs = true;}

	///	set constraint/obstacle g_up (for c(u) <= g_up)
		void set_upper_obstacle(SmartPtr<vector_type> upObs){
			m_spVecOfUpObsValues = upObs; m_bUpperObs = true;}


	///	is lower obstacle set
		bool lower_obs_set(){return m_bLowerObs;}
	///	is upper obstacle set
		bool upper_obs_set(){return m_bUpperObs;}

	///	resets the vectors storing the active indices
		void reset_active_indices(){m_vActiveIndicesLow.resize(0);
			m_vActiveIndicesUp.resize(0);}

	///	access to the vector of active indices wrt the lower obstacle constraint
		SmartPtr<std::vector<MultiIndex<2> > > lower_active_indices(){
			return m_spLowerActiveInd;}

	///	access to the vector of active indices wrt the upper obstacle constraint
		SmartPtr<std::vector<MultiIndex<2> > > upper_active_indices(){
			return m_spUpperActiveInd;}

	///	returns the number of indices (lying in the obstacle surface(s)),
	///	which are stored in 'm_vIndicesOfObsSubsets'
		size_t nrOfIndicesOfObsSubsets() {return m_vIndicesOfObsSubsets.size();};

	///	checks if a given index is in an obstacle subset
		bool index_is_in_obs_subset(size_t index){return true;}



	///	init function checks the obstacle constraints with respect to consistency
		void init(const vector_type& u);

	///	computes the correction for the case that only a lower obstacle is set, i.e. u >= g_low
		virtual void correction_for_lower_obs(vector_type& c, vector_type& lastSol,
				const size_t index) = 0;

	///	computes the correction for the case that only an upper obstacle is set, i.e. u <= g_up
		virtual void correction_for_upper_obs(vector_type& c, vector_type& lastSol,
				const size_t index) = 0;

	///	computes the correction for the case that a lower and an upper obstacle is set
		virtual void correction_for_lower_and_upper_obs(vector_type& c, vector_type& lastSol,
				const size_t index) = 0;

	///	Destructor
		virtual ~IObstacleConstraint(){};

	private:
	///	stores all indices of obstacle subset(s)
		template <typename TElem>
		void obstacle_indices_on_subset(const size_t si);

	protected:
	///	pointer to constraint/obstacle values
		SmartPtr<vector_type> m_spVecOfLowObsValues, m_spVecOfUpObsValues;

	///	pointer to vector of active indices
		SmartPtr<std::vector<MultiIndex<2> > > m_spLowerActiveInd;
		SmartPtr<std::vector<MultiIndex<2> > > m_spUpperActiveInd;

	private:
	/// flag indicating if an obstacle is set
		bool m_bLowerObs, m_bUpperObs;

	///	vector of the algebra indices, which are in the obstacle subsets
		vector<size_t> m_vIndicesOfObsSubsets;

	///	store the indices, which satisfy the constraints (lower resp. upper constraint)
	/// with equality in m_vActiveIndices.
		vector<MultiIndex<2> > m_vActiveIndicesLow, m_vActiveIndicesUp;

	///	name of subset(s), on which the obstacle is defined
		string m_ssName;

	/// subsetGroup
		SubsetGroup m_ssGrp;

	///	pointer to the DofDistribution on the whole domain
		ConstSmartPtr<DoFDistribution> m_spDD;

};

} // end namespace ug

// include implementation
#include "obstacle_constraint_interface_impl.h"

#endif /* __H__UG__LIB_ALGEBRA__OPERATOR__PRECONDITIONER__PROJECTED_GAUSS_SEIDEL__OBSTACLE_CONSTRAINT_INTERFACE__ */

