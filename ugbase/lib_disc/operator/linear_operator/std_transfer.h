/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#ifndef __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__STD_TRANSFER__
#define __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__STD_TRANSFER__

// extern headers
#include <iostream>

// other ug4 modules
#include "common/common.h"
#include "transfer_interface.h"
#include "lib_algebra/operator/debug_writer.h"

#ifdef UG_PARALLEL
#include "lib_disc/parallelization/parallelization_util.h"
#endif

namespace ug{

///	Standard Prolongation Operator
/**	By default a special optimization is performed for p1-lagrange-elements.
 * This optimization is only valid if all elements have been refined with
 * standard refinement rules. If closure elements are generated, this optimization
 * has to be deactivated (use StdTransfer::enable_p1_lagrange_optimization(false)).
 */
template <typename TDomain, typename TAlgebra>
class StdTransfer :
	virtual public ITransferOperator<TDomain, TAlgebra>
{
	public:
	///	Type of base class
		typedef ITransferOperator<TDomain, TAlgebra> base_type;

	///	Type of Algebra
		typedef TAlgebra algebra_type;

	///	Type of Vector
		typedef typename TAlgebra::vector_type vector_type;

	///	Type of Matrix
		typedef typename TAlgebra::matrix_type matrix_type;

	///	Type of Domain
		typedef TDomain domain_type;

	///	Type of GridFunction
		typedef GridFunction<TDomain, TAlgebra> GF;

	public:
	/// Default constructor
		StdTransfer() : ITransferOperator<TDomain, TAlgebra>(),
						m_p1LagrangeOptimizationEnabled(true),
						m_dampRes(1.0), m_dampProl(1.0),
						bCached(true), m_bUseTransposed(true),
						m_spDebugWriter(NULL)
		{};

	/// virtual destructor
		virtual ~StdTransfer(){};

	///	set restriction damping (only applied on vector operation, not (!!) in assembled matrices)
		void set_restriction_damping(number damp) {m_dampRes = damp;}

	///	set prolongation damping (only applied on vector operation, not (!!) in assembled matrices)
		void set_prolongation_damping(number damp) {m_dampProl = damp;}

	///	set debug writer
		void set_debug(SmartPtr<IDebugWriter<TAlgebra> > spDebugWriter) {
			m_spDebugWriter = spDebugWriter;
		}

	///	enables/disables an assembling optimization for p1-lagrange elements
	/**	The optimization is enabled by default. It can however only be used,
	 * if all elements are refined with their standard refinement rule. If one
	 * uses anisotropic refinement or refinement with closure, the optimization
	 * should be disabled.
	 * \todo	The normal assembling strategy should be optimized in such a way
	 * 			that the p1-lagrange optimization is no longer required. This
	 * 			however involves something like a ref-type-hash for each element,
	 * 			which returns a unique number based on the types and order of children.*/
		void enable_p1_lagrange_optimization(bool enable)	{m_p1LagrangeOptimizationEnabled = enable;}
		bool p1_lagrange_optimization_enabled() const		{return m_p1LagrangeOptimizationEnabled;}

	///	sets if restriction and prolongation are transposed
		void set_use_transposed(bool bTransposed) {m_bUseTransposed = bTransposed;}

	public:
	///	Set levels
		virtual void set_levels(GridLevel coarseLevel, GridLevel fineLevel) {}

	///	initialize the operator
		virtual void init() {}

	///	returns new instance with same setting
		virtual SmartPtr<ITransferOperator<TDomain, TAlgebra> > clone();

	/// apply Operator, interpolate function
		virtual void prolongate(vector_type& uFine, const vector_type& uCoarse){
			GF* pFine = dynamic_cast<GF*>(&uFine);
			const GF* pCoarse = dynamic_cast<const GF*>(&uCoarse);
			if(!pFine || !pCoarse)
				UG_THROW("StdTransfer: fine and coarse vectors expected to be "
						"a grid function.");
			prolongate(*pFine, *pCoarse);
		}

	/// apply transposed Operator, restrict function
		virtual void do_restrict(vector_type& uCoarse, const vector_type& uFine){
			const GF* pFine = dynamic_cast<const GF*>(&uFine);
			GF* pCoarse = dynamic_cast<GF*>(&uCoarse);
			if(!pFine || !pCoarse)
				UG_THROW("StdTransfer: fine and coarse vectors expected to be "
						"a grid function.");
			do_restrict(*pCoarse, *pFine);
		}

	public:
	///	returns prolongation as a matrix
		virtual SmartPtr<matrix_type>
		prolongation(const GridLevel& fineGL, const GridLevel& coarseGL,
		             ConstSmartPtr<ApproximationSpace<TDomain> > spApproxSpace);

	///	returns restriction as a matrix
		virtual SmartPtr<matrix_type>
		restriction(const GridLevel& coarseGL, const GridLevel& fineGL,
		            ConstSmartPtr<ApproximationSpace<TDomain> > spApproxSpace);

	///	apply operator to a grid function
		void prolongate(GF& uFine, const GF& uCoarse);

	///	apply operator to a grid function
		void do_restrict(GF& uCoarse, const GF& uFine);

	protected:
	///	debug writing of matrix
		void write_debug(const matrix_type& mat, std::string name,
		                 const GridLevel& glTo, const GridLevel& glFrom);

		template <typename TChild>
		void assemble_restriction(matrix_type& mat,
		                          const DoFDistribution& coarseDD,
		                          const DoFDistribution& fineDD,
		                          ConstSmartPtr<TDomain> spDomain);
		void assemble_restriction(matrix_type& mat,
		                          const DoFDistribution& coarseDD,
		                          const DoFDistribution& fineDD,
		                          ConstSmartPtr<TDomain> spDomain);

		template <typename TChild>
		void assemble_prolongation(matrix_type& mat,
                                   const DoFDistribution& fineDD,
                                   const DoFDistribution& coarseDD,
                                   ConstSmartPtr<TDomain> spDomain);
		void assemble_prolongation(matrix_type& mat,
                                   const DoFDistribution& fineDD,
                                   const DoFDistribution& coarseDD,
                                   ConstSmartPtr<TDomain> spDomain);

		void assemble_prolongation_p1(matrix_type& mat,
		                             const DoFDistribution& fineDD,
		                             const DoFDistribution& coarseDD);

	protected:
	///	struct to distinguish already assembled operators
		struct TransferKey{
			TransferKey(const GridLevel& toGL_, const GridLevel& fromGL_,
			            const RevisionCounter& revCnt_)
			: toGL(toGL_), fromGL(fromGL_), revCnt(revCnt_) {}
			GridLevel toGL, fromGL;
			RevisionCounter revCnt;

			bool operator<(const TransferKey& other) const {
				if(revCnt != other.revCnt) return revCnt < other.revCnt;
				if(toGL != other.toGL) return toGL < other.toGL;
				return fromGL < other.fromGL;
			}
		};

		typedef std::map<TransferKey, SmartPtr<matrix_type> > TransferMap;
		TransferMap m_mRestriction;
		TransferMap m_mProlongation;

		void remove_outdated(TransferMap& map, const RevisionCounter& revCnt) {
			typedef typename TransferMap::iterator iterator;
			for(iterator iter = map.begin(); iter != map.end();)
			{
				const RevisionCounter& cnt = iter->first.revCnt;
				if((cnt.obj() == revCnt.obj()) && (cnt != revCnt)){
					map.erase(iter++);
				} else {
					++iter;
				}
			}
		}

	protected:
	///	list of post processes
		using base_type::m_vConstraint;

	///	flag for p1-lagrange-optimization
		bool m_p1LagrangeOptimizationEnabled;

	///	damping parameter
		number m_dampRes;
		number m_dampProl;

	///	flag if cached (matrix) transfer used
		bool bCached;

	///	flag if transposed is used
		bool m_bUseTransposed;

	///	debug writer
		SmartPtr<IDebugWriter<TAlgebra> > m_spDebugWriter;
};

} // end namespace ug

#include "std_transfer_impl.h"

#endif /* __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__STD_TRANSFER__ */
