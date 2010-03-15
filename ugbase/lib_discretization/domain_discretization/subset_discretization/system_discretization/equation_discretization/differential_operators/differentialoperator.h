/*
 * differentialoperator.h
 *
 *  Created on: 04.05.2009
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__DIFFERENTIALOPERATOR__
#define __H__LIBDISCRETIZATION__DIFFERENTIALOPERATOR__

#include "lib_grid/lib_grid.h"
#include "numericalsolution.h"
#include "../common/common.h"
#include <string>

namespace ug {

class DifferentialOperator {

	public:
		DifferentialOperator(std::string name);

		void set_name(std::string name);
		std::string name();

		virtual ~DifferentialOperator()
		{}

	protected:
		std::string m_name;
};

template <int d>
class DivergenzDifferentialOperator : public DifferentialOperator {

	public:
		DivergenzDifferentialOperator(std::string str):DifferentialOperator(str)
		{};

		virtual void compute_jacobian_at_ip(MathVector<d> GlobIP, number Shape[], MathVector<d> ShapeGrad[], MathVector<d> JacValue[], int nsh) = 0;
		virtual void compute_defect_at_ip(MathVector<d> GlobIP, number Shape, MathVector<d> GlobShapeGrad, MathVector<d>& DefectValue) = 0;

		virtual ~DivergenzDifferentialOperator()
		{};

};

template <int d>
class ScalarDifferentialOperator : public DifferentialOperator {

	public:
		ScalarDifferentialOperator(std::string str):DifferentialOperator(str)
		{};

		virtual void compute_jacobian_at_ip(MathVector<d> GlobIP, number Shape[], MathVector<d> ShapeGrad[], number JacValue[], int nsh) = 0;
		virtual void compute_defect_at_ip(MathVector<d> GlobIP, number Shape, MathVector<d> GlobShapeGrad, number& DefectValue) = 0;

		virtual ~ScalarDifferentialOperator()
		{};

};

template <int d>
class TimeOperator : public DifferentialOperator {

	public:
		TimeOperator(std::string str):DifferentialOperator(str)
		{};

		virtual void compute_jacobian_at_ip(MathVector<d> GlobIP, number Shape[], MathVector<d> ShapeGrad[], number JacValue[], int nsh) = 0;
		virtual void compute_defect_at_ip(MathVector<d> GlobIP, number Shape, MathVector<d> GlobShapeGrad, number& DefectValue) = 0;

		virtual ~TimeOperator()
		{};

};

template <int d>
class TimeIdentity : public TimeOperator<d> {

	public:
		TimeIdentity(std::string str, NumericalSolution<d>& NumSol):TimeOperator<d>(str)
		{
			m_NumericalSolution = &NumSol;
		};

		void compute_jacobian_at_ip(MathVector<d> GlobIP, number Shape[], MathVector<d> ShapeGrad[], number JacValue[], int nsh)
		{
			for(int j=0; j < nsh; j++)
			{
				JacValue[j] = Shape[j];
			}
		}

		void compute_defect_at_ip(MathVector<d> GlobIP, number Shape, MathVector<d> GlobShapeGrad, number& DefectValue)
		{
			DefectValue = Shape;
		}


		~TimeIdentity()
		{}

	protected:
		NumericalSolution<d>* m_NumericalSolution;

};

template <int d>
class ReactionOp : public ScalarDifferentialOperator<d> {

	protected:
		typedef bool (*ReactionFunction)(MathVector<d>,number&);

	public:
		ReactionOp(std::string str, ReactionFunction ReactionFct, NumericalSolution<d>& NumSol) : ScalarDifferentialOperator<d>(str)
		{
			m_NumericalSolution = &NumSol;
			m_ReactionFct = ReactionFct;
		};

		void compute_jacobian_at_ip(MathVector<d> GlobIP, number Shape[], MathVector<d> GlobShapeGrad[], number JacValue[], int nsh)
		{
			number f;

			m_ReactionFct(GlobIP, f);

			for(int j=0; j < nsh; j++)
			{
				JacValue[j] = f * Shape[j];
			}
		}

		void compute_defect_at_ip(MathVector<d> GlobIP, number Shape, MathVector<d> GlobShapeGrad, number& DefectValue)
		{
			number f;

			m_ReactionFct(GlobIP, f);

			DefectValue = f* Shape;
		}

		~ReactionOp()
		{}

	protected:
		NumericalSolution<d>* m_NumericalSolution;
		ReactionFunction m_ReactionFct;

};

template <int d>
class DivDGradOp : public DivergenzDifferentialOperator<d> {

	protected:
		typedef bool (*DiffTensorFunction)(MathVector<d>,MathMatrix<d,d>&);

	public:
		DivDGradOp(std::string str, DiffTensorFunction DiffTensor, NumericalSolution<d>& NumSol):DivergenzDifferentialOperator<d>(str)
		{
			m_NumericalSolution = &NumSol;
			m_DiffTensor = DiffTensor;
		};

		void compute_jacobian_at_ip(MathVector<d> GlobIP, number Shape[], MathVector<d> GlobShapeGrad[], MathVector<d> JacValue[], int nsh)
		{
			MathMatrix<d,d> D;

			m_DiffTensor(GlobIP, D);

			for(int j=0; j < nsh; j++)
			{
				MatVecMult(JacValue[j], D, GlobShapeGrad[j]);
				VecMultiply(JacValue[j], JacValue[j], -1.0);
			}
		}

		void compute_defect_at_ip(MathVector<d> GlobIP, number Shape, MathVector<d> GlobShapeGrad, MathVector<d>& DefectValue)
		{
			MathMatrix<d,d> D;

			m_DiffTensor(GlobIP, D);

			MatVecMult(DefectValue, D, GlobShapeGrad);
			VecMultiply(DefectValue, DefectValue, -1.0);
		}

		~DivDGradOp()
		{}

	protected:
		NumericalSolution<d>* m_NumericalSolution;
		DiffTensorFunction m_DiffTensor;

};

} /* end namespace libDiscretization */



#endif /* __H__LIBDISCRETIZATION__DIFFERENTIALOPERATOR__ */
