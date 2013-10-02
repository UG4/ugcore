/*
 * nl_gauss_seidel.h
 *
 *  Created on: 07.01.2013
 *  (main parts are based on the structure of
 *  	newton.h and some ideas of Sebastian Reiter & Andreas Vogel)
 *
 *      Author: raphaelprohl
 */

#ifndef NL_GAUSS_SEIDEL_H_
#define NL_GAUSS_SEIDEL_H_

#include "lib_algebra/operator/interface/operator_inverse.h"

// modul intern headers
#include "lib_disc/assemble_interface.h"
#include "lib_disc/operator/non_linear_operator/assembled_non_linear_operator.h"
#include "lib_disc/spatial_disc/local_to_global/local_to_global_mapper.h"

namespace ug {

template <typename TAlgebra>
class LocalToGlobalMapper_NL_GS : public ILocalToGlobalMapper<TAlgebra>
{
	public:
	///	Algebra type
		typedef TAlgebra algebra_type;

	///	Type of algebra matrix
		typedef typename algebra_type::matrix_type matrix_type;

	///	Type of algebra vector
		typedef typename algebra_type::vector_type vector_type;

	public:
	///	default constructor
		LocalToGlobalMapper_NL_GS() {}//m_assemblingDoFIndex = 0;}

	///	adds a local vector to the global rhs
		void AddLocalVec(vector_type& vec, const LocalVector& lvec, ConstSmartPtr<DoFDistribution> dd);

	///	adds a local matrix to the global matrix
		void AddLocalMatToGlobal(matrix_type& mat, const LocalMatrix& lmat, ConstSmartPtr<DoFDistribution> dd);

	/// sets assembling index
		void set_assemblingDoFindex(const DoFIndex assIndex){ m_assemblingDoFIndex = assIndex;}

	///	destructor
		~LocalToGlobalMapper_NL_GS() {};

	private:
		DoFIndex m_assemblingDoFIndex;
};

template <typename TDomain, typename TAlgebra>
class NLGaussSeidelSolver
	: public IOperatorInverse<typename TAlgebra::vector_type>,
	  public DebugWritingObject<TAlgebra>
{
	private:
	///	own type
		typedef NLGaussSeidelSolver<TDomain, TAlgebra> this_type;

	public:
	///	Algebra type
		typedef TAlgebra algebra_type;

	///	Matrix type
		typedef typename TAlgebra::matrix_type matrix_type;

	///	Vector type
		typedef typename TAlgebra::vector_type vector_type;

	///	Value type
		typedef typename vector_type::value_type value_type;

	///	Domain type
		typedef TDomain domain_type;

	///	Type of approximation space
		typedef ApproximationSpace<domain_type>	approx_space_type;

	///	Type of geometric base object
		typedef typename domain_traits<TDomain::dim>::geometric_base_object geometric_base_object;

	protected:
		typedef DebugWritingObject<TAlgebra> base_writer_type;

	public:
	///	constructor
		NLGaussSeidelSolver();

	///	constructor
		NLGaussSeidelSolver(SmartPtr<approx_space_type> spApproxSpace,
					SmartPtr<IConvergenceCheck<vector_type> > spConvCheck);

		void set_approximation_space(SmartPtr<approx_space_type> spApproxSpace)
		{m_spApproxSpace = spApproxSpace;}
		void set_convergence_check(SmartPtr<IConvergenceCheck<vector_type> > spConvCheck);
		void set_damp(const number damp) {m_damp = damp;}
	///	sets constraint/obstacle
		void set_constraint(const vector_type& cons) {m_ConsVec = cons; m_bProjectedGS = true;}

	/// This operator inverts the Operator N: Y -> X
		virtual bool init(SmartPtr<IOperator<vector_type> > N);

	/// prepare Operator
		virtual bool prepare(vector_type& u);

	/// apply Operator, i.e. N^{-1}(0) = u
		virtual bool apply(vector_type& u);

	private:
	///	help functions for debug output
	///	\{
		void write_debug(const vector_type& vec, const char* filename);
		void write_debug(const matrix_type& mat, const char* filename);
	/// \}

	private:
		///	Approximation Space
		SmartPtr<approx_space_type> m_spApproxSpace;

		///	DoF distribution pointer
		ConstSmartPtr<DoFDistribution> m_spLevDD;
		ConstSmartPtr<DoFDistribution> m_spSurfDD;

		/// DoF Distribution used
		GridLevel m_gridLevel;

		SmartPtr<IConvergenceCheck<vector_type> > m_spConvCheck;

		SmartPtr<AssembledOperator<algebra_type> > m_N;
		SmartPtr<AssembledLinearOperator<algebra_type> > m_J_block;
		SmartPtr<IAssemble<TAlgebra> > m_spAss;

		///	damping factor
		number m_damp;

		///	vector describing a constraint
		vector_type m_ConsVec;
		bool m_bProjectedGS;

		///	call counter
		int m_dgbCall;

		/* TODO: hier alle typen von geometric_base_object zulassen?
		typedef Attachment<vector<geometric_base_object*> > AElemList; 	//attachment type: attachment of ElemDatas
		AElemList m_aElemList;						//the instance of the attachment type
		typedef Grid::VertexAttachmentAccessor<AElemList>	ElemListAccessor;
		ElemListAccessor m_aaElemList;
		///	use this method to make sure that all required attachments are attached
		void attach_attachments()
		{
			typename TDomain::grid_type& grid = *this->domain().grid();
			grid.attach_to_vertices(m_aElemList);
			m_aaElemList.access(grid, m_aElemList);
		}*/

		// TODO: hier alle typen von geometric_base_object zulassen?
		typedef std::vector<geometric_base_object*> elemList;
		std::vector<elemList> m_vElemList;

		///	selector of elements with contributions to a specific DoF
		Selector m_sel;

		///	LocalToGlobal-Mapper which accounts for index-wise assembling
		LocalToGlobalMapper_NL_GS<TAlgebra> m_map;
};

}

#include "nl_gauss_seidel_impl.h"

#endif /* NL_GAUSS_SEIDEL_H_ */
