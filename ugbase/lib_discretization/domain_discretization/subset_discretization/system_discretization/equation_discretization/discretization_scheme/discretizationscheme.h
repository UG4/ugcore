/*
 * discretizationscheme.h
 *
 *  Created on: 26.06.2009
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__DISCRETIZATIONSCHEME__
#define __H__LIBDISCRETIZATION__DISCRETIZATIONSCHEME__

#include "differentialoperator.h"
#include "dirichletbndcond.h"
#include "rhs.h"
#include "referenceelement.h"
#include "dofhandler.h"
#include "trialspace.h"

#include "../lib_grid/lib_grid.h"
#include "../lib_algebra/lib_algebra.h"
#include "../common/common.h"

#include <iostream>
#include <vector>

namespace ug {

enum DiscretizationSchemeID
{
	DISC_SCHEME_INVALID = -1,
	DISC_SCHEME_FE1,
	DISC_SCHEME_FVE1lump,
	NUM_DISC_SCHEME
};

template <int d>
class DirichletValues{

	public:

		bool add_dirichlet_nodes(NumericalSolution<d>& u, int nr_func, DirichletBNDCond<d>* dirichbnd, SubsetHandler& sh, uint subsetIndex);

		bool set_values(Vector& vec);

		bool set_zero_values(Vector& vec);

		bool set_rows(Matrix& mat);

	private:
		std::vector<int> m_vector_indices;
		std::vector<int> m_matrixrow_indices;
		std::vector<number> m_vector_values;
};


template<typename TElem, int d>
class DiscretizationScheme{

	public:
		virtual void prepareElement(TElem* elem, NumericalSolution<d>& u, int nr_func, SubsetHandler& sh, int SubsetIndex) = 0;

		virtual void reset_local_jacobian() = 0;
		virtual void reset_local_defect() = 0;

		virtual void add_op_to_local_jacobian(TimeOperator<d>* op, TElem* elem, NumericalSolution<d>& u) = 0;
		virtual void add_op_to_local_jacobian(DivergenzDifferentialOperator<d>* op, TElem* elem, NumericalSolution<d>& u) = 0;
		virtual void add_op_to_local_jacobian(ScalarDifferentialOperator<d>* op, TElem* elem, NumericalSolution<d>& u) = 0;

		virtual void add_rhs_to_local_defect(RHS<d>* rhs, TElem* elem) = 0;
		virtual void add_op_to_local_defect(TimeOperator<d>* op, TElem* elem, NumericalSolution<d>& u) = 0;
		virtual void add_op_to_local_defect(DivergenzDifferentialOperator<d>* op, TElem* elem, NumericalSolution<d>& u) = 0;
		virtual void add_op_to_local_defect(ScalarDifferentialOperator<d>* op, TElem* elem, NumericalSolution<d>& u) = 0;

		virtual bool send_local_jacobian_to_global_jacobian(TElem* elem, Matrix& mat) = 0;
		virtual bool send_local_jacobian_to_global_jacobian(TElem* elem, Matrix& mat, number scale) = 0;
		virtual bool send_local_defect_to_global_defect(TElem* elem, Vector& vec) = 0;
		virtual bool send_local_defect_to_global_defect(TElem* elem, Vector& vec, number scale) = 0;

		virtual ~DiscretizationScheme()
		{}

	protected:
};


template <typename TElem, int d>
class FVE1lumpDiscretization : public DiscretizationScheme<TElem, d>{

	enum{DiscretizationSchemeID = DISC_SCHEME_FVE1lump};
	static const int RefDim = reference_element_traits<TElem>::RefDim;

	//! TODO: change this to some trial_space_traits<TElem>::nsh;
	static const int m_nsh = reference_element_traits<TElem>::NumberCorners;

	public:
		FVE1lumpDiscretization();

		void prepareElement(TElem* elem, NumericalSolution<d>& u, int nr_func, SubsetHandler& sh, int SubsetIndex);

		void add_op_to_local_jacobian(TimeOperator<d>* op, TElem* elem, NumericalSolution<d>& u);
		void add_op_to_local_jacobian(ScalarDifferentialOperator<d>* op, TElem* elem, NumericalSolution<d>& u);
		void add_op_to_local_jacobian(DivergenzDifferentialOperator<d>* op, TElem* elem, NumericalSolution<d>& u);

		void add_rhs_to_local_defect(RHS<d>* rhs, TElem* elem);
		void add_op_to_local_defect(TimeOperator<d>* op, TElem* elem, NumericalSolution<d>& u);
		void add_op_to_local_defect(ScalarDifferentialOperator<d>* op, TElem* elem, NumericalSolution<d>& u);
		void add_op_to_local_defect(DivergenzDifferentialOperator<d>* op, TElem* elem, NumericalSolution<d>& u);

		void reset_local_jacobian();
		void reset_local_defect();

		bool send_local_jacobian_to_global_jacobian(TElem* tr, Matrix& mat);
		bool send_local_jacobian_to_global_jacobian(TElem* tr, Matrix& mat, number scale);
		bool send_local_defect_to_global_defect(TElem* tr, Vector& vec);
		bool send_local_defect_to_global_defect(TElem* tr, Vector& vec, number scale);

		virtual ~FVE1lumpDiscretization();

	private:
		MathVector<d> m_corners[m_nsh];
		int m_rows[m_nsh];

		ReferenceElement<TElem> m_RefElem;
		const TrialSpace<TElem>* m_TrialSpace;

		const static int m_nip  = 1;

		number m_Shape[m_nip][m_nsh];
		MathVector<RefDim> m_LocShapeGrad[m_nip][m_nsh];
		MathVector<d> m_GlobShapeGrad[m_nip][m_nsh];

		MathVector<RefDim>* m_LocalIP;
		MathVector<RefDim>* m_LocalNormal;
		number* m_SCVF_Size;

		MathMatrix<d,RefDim> m_IPTrafo[m_nip];
		number m_det[m_nip];
		MathVector<d> m_GlobalIP[m_nip];

		number m_loc_u[m_nsh];

		MathMatrix<m_nsh, m_nsh> m_LocalJacobian;
		MathVector<m_nsh> m_LocalDefect;

	private:

		void InitializeIntegrationPoints();

};

template <typename TElem, int d>
class FE1Discretization : public DiscretizationScheme<TElem, d>{

	enum{DiscretizationSchemeID = DISC_SCHEME_FE1};
	static const int RefDim = reference_element_traits<TElem>::RefDim;

	//! TODO: change this to some trial_space_traits<TElem>::nsh;
	static const int m_nsh = reference_element_traits<TElem>::NumberCorners;

	public:
		FE1Discretization();

		void prepareElement(TElem* elem, NumericalSolution<d>& u, int nr_func, SubsetHandler& sh, int SubsetIndex);

		void add_op_to_local_jacobian(TimeOperator<d>* op, TElem* elem, NumericalSolution<d>& u);
		void add_op_to_local_jacobian(ScalarDifferentialOperator<d>* op, TElem* elem, NumericalSolution<d>& u);
		void add_op_to_local_jacobian(DivergenzDifferentialOperator<d>* op, TElem* elem, NumericalSolution<d>& u);

		void add_rhs_to_local_defect(RHS<d>* rhs, TElem* elem);
		void add_op_to_local_defect(TimeOperator<d>* op, TElem* elem, NumericalSolution<d>& u);
		void add_op_to_local_defect(ScalarDifferentialOperator<d>* op, TElem* elem, NumericalSolution<d>& u);
		void add_op_to_local_defect(DivergenzDifferentialOperator<d>* op, TElem* elem, NumericalSolution<d>& u);

		void reset_local_jacobian();
		void reset_local_defect();

		bool send_local_jacobian_to_global_jacobian(TElem* tr, Matrix& mat);
		bool send_local_jacobian_to_global_jacobian(TElem* tr, Matrix& mat, number scale);
		bool send_local_defect_to_global_defect(TElem* tr, Vector& vec);
		bool send_local_defect_to_global_defect(TElem* tr, Vector& vec, number scale);

		virtual ~FE1Discretization();

	private:
		MathVector<d> m_corners[m_nsh];
		int m_rows[m_nsh];

		ReferenceElement<TElem> m_RefElem;
		const TrialSpace<TElem>* m_TrialSpace;

		const static int m_nip  = 1;

		number m_Shape[m_nip][m_nsh];
		MathVector<RefDim> m_LocShapeGrad[m_nip][m_nsh];
		MathVector<d> m_GlobShapeGrad[m_nip][m_nsh];

		MathVector<RefDim>* m_LocalIP;

		MathMatrix<d,RefDim> m_IPTrafo[m_nip];
		number m_det[m_nip];
		MathVector<d> m_GlobalIP[m_nip];

		number m_loc_u[m_nsh];

		MathMatrix<m_nsh, m_nsh> m_LocalJacobian;
		MathVector<m_nsh> m_LocalDefect;

	private:

		void InitializeIntegrationPoints();

};

// Singleton, holding all Discretization Schemes available
template <typename TElem, int d>
class DiscretizationSchemes {

	private:

		static DiscretizationSchemes& inst()
		{
			static DiscretizationSchemes myInst;
			return myInst;
		};

		// private constructor
		DiscretizationSchemes()
		{};

		inline static DiscretizationScheme<TElem, d>& get_DiscretizationScheme(DiscretizationSchemeID type)
		{
			static FVE1lumpDiscretization<TElem, d> FVE1lumpDiscretization;
			static FE1Discretization<TElem, d> FE1Discretization;

			if(type == DISC_SCHEME_FE1)
				return FE1Discretization;
			if(type == DISC_SCHEME_FVE1lump)
				return FVE1lumpDiscretization;

			assert(0 && "ERROR in get_DiscretizationScheme: Not available");
		}

	public:
		static DiscretizationScheme<TElem, d>& DiscretizationScheme(DiscretizationSchemeID type)
		{
			return inst().get_DiscretizationScheme(type);
		}
};

} // End of namespace ug

#include "discretizationscheme_impl.h"


#endif /* __H__LIBDISCRETIZATION__DISCRETIZATIONSCHEME__ */
