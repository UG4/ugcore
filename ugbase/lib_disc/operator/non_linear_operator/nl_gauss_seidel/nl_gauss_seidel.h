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

//#include "lib_grid/lg_base.h"
#include "lib_algebra/operator/interface/operator_inverse.h"

// modul intern headers
#include "lib_disc/assemble_interface.h"
#include "lib_disc/operator/non_linear_operator/assembled_non_linear_operator.h"

namespace ug {

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

	///	Vector type
		typedef typename TAlgebra::vector_type vector_type;

	///	Matrix type
		typedef typename TAlgebra::matrix_type matrix_type;

	///	Value type
		typedef typename TAlgebra::vector_type::value_type val_type;

	///	Domain type
		typedef TDomain domain_type;

	///	Type of approximation space
		typedef ApproximationSpace<domain_type>	approx_space_type;

	///	Type of geometric base object
		typedef typename domain_traits<TDomain::dim>::geometric_base_object geometric_base_object;

	protected:
		typedef DebugWritingObject<TAlgebra> base_writer_type;
		using base_writer_type::write_debug;

	public:
	///	default constructor
		NLGaussSeidelSolver();

	///	constructor
		NLGaussSeidelSolver(SmartPtr<approx_space_type> spApproxSpace,
					SmartPtr<IConvergenceCheck<vector_type> > spConvCheck);

		void set_approximation_space(SmartPtr<approx_space_type> spApproxSpace)
		{m_spApproxSpace = spApproxSpace;}
		void set_convergence_check(SmartPtr<IConvergenceCheck<vector_type> > spConvCheck);
		void set_damp(number damp) {m_damp = damp;}

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

		vector_type m_d;
		vector_type m_d_block;
		vector_type m_c_block;

		SmartPtr<AssembledOperator<algebra_type> > m_N;
		SmartPtr<AssembledLinearOperator<algebra_type> > m_J_block;
		IAssemble<algebra_type>* m_pAss;

		number m_damp;

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
		typedef vector<geometric_base_object*> elemList;
		vector<elemList> m_vElemList;

		///	selector of elements with contributions to a specific DoF
		Selector m_sel;
};

}

#include "nl_gauss_seidel_impl.h"

#endif /* NL_GAUSS_SEIDEL_H_ */
