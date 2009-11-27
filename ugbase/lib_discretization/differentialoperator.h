/*
 * differentialoperator.h
 *
 *  Created on: 04.05.2009
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__DIFFERENTIALOPERATOR__
#define __H__LIBDISCRETIZATION__DIFFERENTIALOPERATOR__

#include "lib_grid/lib_grid.h"
#include "lib_grid/geometric_objects.h"
#include "numericalsolution.h"
#include "../common/math/ugmath.h"
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

class DivergenzDifferentialOperator : public DifferentialOperator {

	public:
		DivergenzDifferentialOperator(std::string str):DifferentialOperator(str)
		{};

		virtual void compute_jacobian_at_ip(MathVector<3> GlobIP, number Shape[], MathVector<2> ShapeGrad[], MathVector<2> JacValue[], int nsh) = 0;
		virtual void compute_defect_at_ip(MathVector<3> GlobIP, number Shape, MathVector<2> GlobShapeGrad, MathVector<2>& DefectValue) = 0;

		virtual ~DivergenzDifferentialOperator()
		{};

};

class ScalarDifferentialOperator : public DifferentialOperator {

	public:
		ScalarDifferentialOperator(std::string str):DifferentialOperator(str)
		{};

		virtual void compute_jacobian_at_ip(MathVector<3> GlobIP, number Shape[], MathVector<2> ShapeGrad[], number JacValue[], int nsh) = 0;
		virtual void compute_defect_at_ip(MathVector<3> GlobIP, number Shape, MathVector<2> GlobShapeGrad, number& DefectValue) = 0;

		virtual ~ScalarDifferentialOperator()
		{};

};

class TimeOperator : public DifferentialOperator {

	public:
		TimeOperator(std::string str):DifferentialOperator(str)
		{};

		virtual void compute_jacobian_at_ip(MathVector<3> GlobIP, number Shape[], MathVector<2> ShapeGrad[], number JacValue[], int nsh) = 0;
		virtual void compute_defect_at_ip(MathVector<3> GlobIP, number Shape, MathVector<2> GlobShapeGrad, number& DefectValue) = 0;

		virtual ~TimeOperator()
		{};

};

class TimeIdentity : public TimeOperator {

	public:
		TimeIdentity(std::string str, NumericalSolution& NumSol):TimeOperator(str)
		{
			m_NumericalSolution = &NumSol;
		};

		void compute_jacobian_at_ip(MathVector<3> GlobIP, number Shape[], MathVector<2> ShapeGrad[], number JacValue[], int nsh)
		{
			for(int j=0; j < nsh; j++)
			{
				JacValue[j] = Shape[j];
			}
		}

		void compute_defect_at_ip(MathVector<3> GlobIP, number Shape, MathVector<2> GlobShapeGrad, number& DefectValue)
		{
			DefectValue = Shape;
		}


		~TimeIdentity()
		{}

	protected:
		NumericalSolution* m_NumericalSolution;

};


class ReactionOp : public ScalarDifferentialOperator {

	protected:
		typedef bool (*ReactionFunction)(MathVector<3>,number&);

	public:
		ReactionOp(std::string str, ReactionFunction ReactionFct, NumericalSolution& NumSol):ScalarDifferentialOperator(str)
		{
			m_NumericalSolution = &NumSol;
			m_ReactionFct = ReactionFct;
		};

		void compute_jacobian_at_ip(MathVector<3> GlobIP, number Shape[], MathVector<2> GlobShapeGrad[], number JacValue[], int nsh)
		{
			number f;

			m_ReactionFct(GlobIP, f);

			for(int j=0; j < nsh; j++)
			{
				JacValue[j] = f * Shape[j];
			}
		}

		void compute_defect_at_ip(MathVector<3> GlobIP, number Shape, MathVector<2> GlobShapeGrad, number& DefectValue)
		{
			number f;

			m_ReactionFct(GlobIP, f);

			DefectValue = f* Shape;
		}

		~ReactionOp()
		{}

	protected:
		NumericalSolution* m_NumericalSolution;
		ReactionFunction m_ReactionFct;

};


class DivDGradOp : public DivergenzDifferentialOperator {

	protected:
		typedef bool (*DiffTensorFunction)(MathVector<3>,MathMatrix<2,2>&);

	public:
		DivDGradOp(std::string str, DiffTensorFunction DiffTensor, NumericalSolution& NumSol):DivergenzDifferentialOperator(str)
		{
			m_NumericalSolution = &NumSol;
			m_DiffTensor = DiffTensor;
		};

		void compute_jacobian_at_ip(MathVector<3> GlobIP, number Shape[], MathVector<2> GlobShapeGrad[], MathVector<2> JacValue[], int nsh)
		{
			MathMatrix<2,2> D;

			m_DiffTensor(GlobIP, D);

			for(int j=0; j < nsh; j++)
			{
				MatVecMult(JacValue[j], D, GlobShapeGrad[j]);
				VecMultiply(JacValue[j], JacValue[j], -1.0);
			}
		}

		void compute_defect_at_ip(MathVector<3> GlobIP, number Shape, MathVector<2> GlobShapeGrad, MathVector<2>& DefectValue)
		{
			MathMatrix<2,2> D;

			m_DiffTensor(GlobIP, D);

			MatVecMult(DefectValue, D, GlobShapeGrad);
			VecMultiply(DefectValue, DefectValue, -1.0);
		}

		~DivDGradOp()
		{}

	protected:
		NumericalSolution* m_NumericalSolution;
		DiffTensorFunction m_DiffTensor;

};

} /* end namespace libDiscretization */



#endif /* __H__LIBDISCRETIZATION__DIFFERENTIALOPERATOR__ */
