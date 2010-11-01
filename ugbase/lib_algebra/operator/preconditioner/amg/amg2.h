/**
 * \file amg.h
 *
 * \author Martin Rupp
 *
 * \date 06.08.2010
 *
 * class declaration for amg
 *
 * Goethe-Center for Scientific Computing 2009-2010.
 */


#ifndef __H__LIB_DISCRETIZATION__AMG_SOLVER__AMG_H__
#define __H__LIB_DISCRETIZATION__AMG_SOLVER__AMG_H__


namespace ug{

template <typename TAlgebra>
class amg:
	public IPreconditioner<	TAlgebra >
{
public:
//	Algebra type
	typedef TAlgebra algebra_type;

//	Vector type
	typedef typename TAlgebra::vector_type vector_type;

//	Matrix type
	typedef typename TAlgebra::matrix_type matrix_type;

	typedef typename matrix_type::value_type value_type;
	
//  functions
	amg() {}
	virtual ILinearIterator<vector_type,vector_type>* clone()
	{
		amg<algebra_type>* clone = new amg<algebra_type>();
		return dynamic_cast<ILinearIterator<vector_type,vector_type>* >(clone);
	}
	//	Name of preconditioner
	virtual ~amg() {}
	void cleanup() {}

protected:
	virtual const char* name() const {return "AMGPreconditioner";}

//	Preprocess routine
	virtual bool preprocess(matrix_type& mat) { return true; }

//	Postprocess routine
	virtual bool postprocess() {return true;}

//	Stepping routine
	virtual bool step(matrix_type& mat, vector_type& c, const vector_type& d)
	{
		return true;
	}

//  data
	void set_nu1(int new_nu1) { }
	void set_nu2(int new_nu2) { }
	void set_gamma(int new_gamma) {  }
	void set_theta(double new_theta) {  }
	void set_sigma(double new_sigma) {  }
	void set_max_levels(int new_max_levels)
	{
	}
	void set_aggressive_coarsening_A_2() { }
	void set_aggressive_coarsening_A_1() { }
	void set_presmoother(ILinearIterator<vector_type, vector_type> *presmoother) {	 }
	void set_postsmoother(ILinearIterator<vector_type, vector_type> *postsmoother) { }
	void set_base_solver(ILinearOperatorInverse<vector_type, vector_type> *basesolver) {  }

	void tostring() const { }

};
	
}
#endif // __H__LIB_DISCRETIZATION__AMG_SOLVER__AMG_H__
