/**
 * \file amg_rs_prolongation.h
 *
 * \author Martin Rupp
 *
 * \date 16.06.2010
 *
 * Goethe-Center for Scientific Computing 2009-2010.
 */


#ifndef __H__LIB_DISCRETIZATION__AMG_SOLVER__AMG_SOLVER__
#define __H__LIB_DISCRETIZATION__AMG_SOLVER__AMG_SOLVER__

#include "amg.h"

namespace ug{


template <typename TAlgebra>
class AMGPreconditioner : public IPreconditioner<TAlgebra>
{
	public:
	//	Algebra type
		typedef TAlgebra algebra_type;

	//	Vector type
		typedef typename TAlgebra::vector_type vector_type;

	//	Matrix type
		typedef typename TAlgebra::matrix_type matrix_type;

	public:
	//	Constructor
		AMGPreconditioner() : m_amg() {};

	// 	Clone
		virtual ILinearIterator<vector_type,vector_type>* clone()
		{
			AMGPreconditioner<algebra_type>* clone = new AMGPreconditioner<algebra_type>();
			return dynamic_cast<ILinearIterator<vector_type,vector_type>* >(clone);
		}

	protected:
	//	Name of preconditioner
		virtual const char* name() const {return "AMGPreconditioner";}

	//	Preprocess routine
		virtual bool preprocess(matrix_type& mat)
		{
		//	todo : error check
			m_amg.init(mat);

			return true;
		}

	//	Postprocess routine
		virtual bool postprocess() {return true;}

	//	Stepping routine
		virtual bool step(matrix_type& mat, vector_type& c, const vector_type& d)
		{
			return 		m_amg.get_correction_and_update_defect(d, c);
		}

		void set_debug_positions(const MathVector<2> *positions2d, size_t size)
		{
			m_amg.set_debug_positions(positions2d, size);
		}

		void set_debug_positions(const MathVector<3> *positions3d, size_t size)
		{
			m_amg.set_debug_positions(positions3d, size);
		}


	protected:
		amg<matrix_type, vector_type> m_amg;
		int m_nu1, m_nu2;

};



} // namespace ug


#endif /* __H__LIB_DISCRETIZATION__AMG_SOLVER__AMG_SOLVER__ */
