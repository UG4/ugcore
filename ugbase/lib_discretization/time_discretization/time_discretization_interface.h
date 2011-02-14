/*
 * time_discretization_interface.h
 *
 *  Created on: 23.10.2009
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISCRETIZATION__TIME_DISCRETIZATION__TIME_DISCRETIZATION_INTERFACE__
#define __H__UG__LIB_DISCRETIZATION__TIME_DISCRETIZATION__TIME_DISCRETIZATION_INTERFACE__

// extern libraries
#include <deque>

// other ug libraries
#include "common/common.h"

// module intern libraries
#include "lib_discretization/assemble.h"
#include "lib_discretization/spatial_discretization/domain_discretization_interface.h"

namespace ug{

/// \ingroup lib_disc_time_assemble
/// @{

/// old solutions and time steps
/**
 * This class holds solutions and corresponding points in time. It is
 * intended to group previous computed solutions for a time stepping scheme, such
 * that previous steps can be passed to a time stepping scheme at once. Internally,
 * this object is basically a deque of old solutions, and adding a newly
 * computed solution lets the object pop the oldest stored solution.
 */
template <typename TVector>
class PreviousSolutions
{
	public:
	///	vector type of solutions
		typedef TVector vector_type;

	public:

	///	returns number of previous time steps handled
		size_t size() const {return m_vPreviousSolution.size();}

	///	adds new time point, oldest solution is discarded and returned
		vector_type* push_discard_oldest(vector_type& vec, number time)
		{
			vector_type* discardVec = m_vPreviousSolution.back().solution();
			m_vPreviousSolution.pop_back();
			m_vPreviousSolution.push_front(PreviousSolution(vec, time));
			return discardVec;
		}

	///	adds new time point, not discarding the oldest
		void push(vector_type& vec, number time)
		{
			m_vPreviousSolution.push_front(PreviousSolution(vec, time));
		}

	///	returns previous solution
		const vector_type& solution(size_t i) const {return *(m_vPreviousSolution.at(i).solution());}

	///	returns previous solution
		vector_type& solution(size_t i) {return *(m_vPreviousSolution.at(i).solution());}

	///	returns point in time for previous solution
		number time(size_t i) const {return m_vPreviousSolution.at(i).time();}

	///	returns oldest solution
		vector_type& oldest_solution() {return *(m_vPreviousSolution.back().solution());}

	protected:
		class PreviousSolution
		{
			public:
				PreviousSolution() : vec(NULL), t(0.0) {}

				PreviousSolution(vector_type& vec_, number t_)
					: vec(&vec_), t(t_) {}

			///	access solution
				vector_type* solution() {return vec;}

			///	const access solution
				const vector_type* solution() const {return vec;}

			///	access time
				number& time() {return t;}

			///	const access time
				const number& time() const {return t;}

			protected:
			//	solution vector at time point
				vector_type* vec;

			//	point in time
				number t;
		};

	//	deque of previous solutions
		std::deque<PreviousSolution> m_vPreviousSolution;
};



/// Time Discretization Interface
/**
 * Defines the time discretization interface.
 *
 * This class uses a ISpatialDiscratization in order to implement the
 * IAssemble interface.
 *
 * After the method prepare step has been called, Jacobian/Defect can be computed.
 * \tparam 	TDoFDistribution	DoF Distribution Type
 * \tparam	TAlgebra			Algebra Type
 */
template <	typename TDoFDistribution,
			typename TAlgebra>
class ITimeDiscretization
	: public IAssemble<TDoFDistribution, TAlgebra>
{
	public:
	//	DoF Distribution type
		typedef IDoFDistribution<TDoFDistribution> dof_distribution_type;

	//	Algebra type
		typedef TAlgebra algebra_type;

	//	Vector type
		typedef typename algebra_type::vector_type vector_type;

	// 	Domain Discretization type
		typedef IDomainDiscretization<TDoFDistribution, algebra_type>
			domain_discretization_type;

	public:
	/// Default constructor
	/**
	 * Creates an empty Time Discretization. In order to use this class a
	 * spatial Discretization has to be set.
	 */
		ITimeDiscretization()
			: m_pDomDisc(NULL)	{}

	/// create and set domain discretization
	/**
	 * \param[in] 	dd	Domain Discretization
	 */
		ITimeDiscretization(domain_discretization_type& dd)
			: m_pDomDisc(&dd)
		{}

	///	set the domain discretization
		void set_domain_discretization(domain_discretization_type& dd)
		{
			m_pDomDisc = &dd;
		}

	///	get the domain discretization
		domain_discretization_type* get_domain_discretization()
		{
			return m_pDomDisc;
		}

	/// prepares the assembling of Defect/Jacobian for a time step
	/**
	 *	This function supplies the TimeDiscretization with previous time
	 *	steps and step size before the assembling routines can be called.
	 *
	 * \param[in] u_old 	the solution at the previous time steps
	 * \param[in] time_old	the time at the previous time steps
	 * \param[in] dt		size of time step
	 */
		virtual bool prepare_step(const PreviousSolutions<vector_type>& prevSol,
		                          number dt) = 0;

	/// returns number of previous time steps needed
		virtual size_t num_prev_steps() = 0;

		// todo: Remove
		size_t num_fct() const{return m_pDomDisc->num_fct();}

	protected:
		domain_discretization_type* m_pDomDisc; ///< Domain Discretization
};

/// @}

} // end namespace ug

#endif /* __H__UG__LIB_DISCRETIZATION__TIME_DISCRETIZATION__TIME_DISCRETIZATION_INTERFACE__ */
