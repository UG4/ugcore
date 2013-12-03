/*
 * std_transfer.h
 *
 *  Created on: 04.12.2009
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__STD_TRANSFER__
#define __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__STD_TRANSFER__

// extern headers
#include <iostream>

// other ug4 modules
#include "common/common.h"
#include "transfer_interface.h"
#include "lib_disc/spatial_disc/constraints/constraint_interface.h"
#include "lib_algebra/operator/debug_writer.h"

#ifdef UG_PARALLEL
#include "lib_disc/parallelization/parallelization_util.h"
#endif

namespace ug{

///	Standard Prologation Operator
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
		StdTransfer() : m_p1LagrangeOptimizationEnabled(true),
						m_dampRes(1.0), bCached(true), m_spDebugWriter(NULL)
		{clear_constraints();};

	/// virtual destructor
		virtual ~StdTransfer(){};

	///	set interpolation damping
		void set_restriction_damping(number damp) {m_dampRes = damp;}

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

	///	clears dirichlet post processes
		void clear_constraints() {m_vConstraint.clear();}

	///	adds a dirichlet post process (not added if already registered)
		void add_constraint(SmartPtr<IConstraint<TAlgebra> > pp);

	///	removes a post process
		void remove_constraint(SmartPtr<IConstraint<TAlgebra> > pp);

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

		void assemble_restriction_elemwise(matrix_type& mat,
		                                    const DoFDistribution& coarseDD, const DoFDistribution& fineDD,
		                                    ConstSmartPtr<TDomain> spDomain);

		void assemble_restriction_p1(matrix_type& mat,
		                              const DoFDistribution& coarseDD, const DoFDistribution& fineDD);

		template <typename TElem>
		void set_identity_on_pure_surface(matrix_type& mat,
		                                  const DoFDistribution& coarseDD, const DoFDistribution& fineDD);

		void set_identity_on_pure_surface(matrix_type& mat,
		                                  const DoFDistribution& coarseDD, const DoFDistribution& fineDD);

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
		std::vector<SmartPtr<IConstraint<TAlgebra> > > m_vConstraint;

	///	flag for p1-lagrange-optimization
		bool m_p1LagrangeOptimizationEnabled;

	///	damping parameter
		number m_dampRes;

	///	flag if cached (matrix) transfer used
		bool bCached;

	///	debug writer
		SmartPtr<IDebugWriter<TAlgebra> > m_spDebugWriter;
};

} // end namespace ug

#include "std_transfer_impl.h"

#endif /* __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__STD_TRANSFER__ */
