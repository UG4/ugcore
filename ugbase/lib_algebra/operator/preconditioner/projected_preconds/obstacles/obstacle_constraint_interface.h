/*
 * obstacle_constraint_interface.h
 *
 *  Created on: 25.11.2013
 *      Author: raphaelprohl
 */

#ifndef __H__UG__LIB_ALGEBRA__OPERATOR__PRECONDITIONER__PROJECTED_PRECONDS__OBSTACLE_CONSTRAINT_INTERFACE__
#define __H__UG__LIB_ALGEBRA__OPERATOR__PRECONDITIONER__PROJECTED_PRECONDS__OBSTACLE_CONSTRAINT_INTERFACE__

#include "lib_disc/common/multi_index.h"

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
 * 	One simple example is the ScalarObstacle-class.
 *
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
	/// constructor
		IObstacleConstraint(): m_bLowerObs(false), m_bUpperObs(false){
			m_vActiveIndicesLow.resize(0); m_vActiveIndicesUp.resize(0);};

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
		void reset_active_indices(){m_vActiveIndicesLow.resize(0); m_vActiveIndicesUp.resize(0);}


	///	computes the correction for the case that only a lower obstacle is set, i.e. u >= g_low
		virtual void correction_for_lower_obs(vector_type& c, vector_type& lastSol, const size_t index, const value_type& tmpSol) = 0;

	///	computes the correction for the case that only an upper obstacle is set, i.e. u <= g_up
		virtual void correction_for_upper_obs(vector_type& c, vector_type& lastSol, const size_t index, const value_type& tmpSol) = 0;

	///	computes the correction for the case that a lower and an upper obstacle is set
		virtual void correction_for_lower_and_upper_obs(vector_type& c, vector_type& lastSol, const size_t index, const value_type& tmpSol) = 0;

	///	Destructor
		virtual ~IObstacleConstraint(){};

	public:
	///	store the indices, which satisfy the constraints (lower resp upper constraint)
	/// with equality in m_vActiveIndices.
		std::vector<MultiIndex<2> > m_vActiveIndicesLow, m_vActiveIndicesUp;

	protected:
	///	pointer to constraint/obstacle values
		SmartPtr<vector_type> m_spVecOfLowObsValues, m_spVecOfUpObsValues;

	private:
	/// flag indicating if an obstacle is set
		bool m_bLowerObs, m_bUpperObs;
};

} // end namespace ug

#endif /* __H__UG__LIB_ALGEBRA__OPERATOR__PRECONDITIONER__PROJECTED_PRECONDS__OBSTACLE_CONSTRAINT_INTERFACE__ */

