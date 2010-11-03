/*
 * user_data.h
 *
 *  Created on: 13.10.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__USER_DATA__
#define __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__USER_DATA__

#include <boost/function.hpp>

namespace ug {


template <int dim>
class IUserNumberProvider
{
	public:
	//	Functor Type
		typedef boost::function<void (number& n, const MathVector<dim>& x, number& time)> functor_type;

	/// provides the functor
		virtual functor_type get_functor() const = 0;

	///	virtual destructor
		virtual ~IUserNumberProvider(){}
};



template <int dim>
class IUserVectorProvider
{
	public:
	//	Functor Type
		typedef boost::function<void (MathVector<dim>& v, const MathVector<dim>& x, number& time)> functor_type;

	/// provides the functor
		virtual functor_type get_functor() const = 0;

	///	virtual destructor
		virtual ~IUserVectorProvider(){}
};


template <int dim>
class IUserMatrixProvider
{
	public:
	//	Functor Type
		typedef boost::function<void (MathMatrix<dim,dim>& D, const MathVector<dim>& x, number& time)> functor_type;

	/// provides the functor
		virtual functor_type get_functor() const = 0;

	///	virtual destructor
		virtual ~IUserMatrixProvider(){}
};



template <int dim>
class IBoundaryNumberProvider
{
	public:
	//	Functor Type
		typedef boost::function<bool (number& n, const MathVector<dim>& x, number& time)> functor_type;

	/// provides the functor
		virtual functor_type get_functor() const = 0;

	///	virtual destructor
		virtual ~IBoundaryNumberProvider(){}
};


///////////////////////////////////////
///////////////////////////////////////
// Const Data
///////////////////////////////////////
///////////////////////////////////////

template <int dim>
class ConstUserNumber : public IUserNumberProvider<dim>
{
	public:
	//	Functor Type
		typedef typename IUserNumberProvider<dim>::functor_type functor_type;

	//	return functor
		virtual functor_type get_functor() const {return boost::ref(*this);}

	public:
		ConstUserNumber() {m_Number = 0.0;}

		void set(number val)
		{
			m_Number = val;
		}

		void print() const
		{
			UG_LOG("ConstUserNumber:" << m_Number << "\n");
		}

		void operator() (number& c, const MathVector<dim>& x, number time = 0.0) const
		{
			c = m_Number;
		}

	protected:
		number m_Number;
};

template <int dim>
class ConstUserVector : public IUserVectorProvider<dim>
{
	public:
	//	Functor Type
		typedef typename IUserVectorProvider<dim>::functor_type functor_type;

	//	return functor
		virtual functor_type get_functor() const {return boost::ref(*this);}

	public:
		ConstUserVector() {set_all_entries(0.0);}

		void set_all_entries(number val) { m_Vector = val;}

		void set_entry(size_t i, number val)
		{
			m_Vector[i] = val;
		}

		void print() const
		{
			UG_LOG("ConstUserVector:" << m_Vector << "\n");
		}

		void operator() (MathVector<dim>& v, const MathVector<dim>& x, number time = 0.0) const
		{
			v = m_Vector;
		}

	protected:
		MathVector<dim> m_Vector;
};

template <int dim>
class ConstUserMatrix : public IUserMatrixProvider<dim>
{
	public:
	//	Functor Type
		typedef typename IUserMatrixProvider<dim>::functor_type functor_type;

	//	return functor
		virtual functor_type get_functor() const {return boost::ref(*this);}

	public:
		ConstUserMatrix() {set_diag_tensor(1.0);}

		void set_diag_tensor(number val)
		{
			for(size_t i = 0; i < dim; ++i){
				for(size_t j = 0; j < dim; ++j){
					m_Tensor[i][j] = 0;
				}
				m_Tensor[i][i] = val;
			}
		}

		void set_all_entries(number val)
		{
			for(size_t i = 0; i < dim; ++i){
				for(size_t j = 0; j < dim; ++j){
					m_Tensor[i][j] = val;
				}
			}
		}

		void set_entry(size_t i, size_t j, number val)
		{
			m_Tensor[i][j] = val;
		}

		void print() const
		{
			UG_LOG("ConstUserMatrix:\n" << m_Tensor << "\n");
		}

		void operator() (MathMatrix<dim, dim>& D, const MathVector<dim>& x, number time = 0.0) const
		{
			D = m_Tensor;
		}

	protected:
		MathMatrix<dim, dim> m_Tensor;
};

template <int dim>
class ConstBoundaryNumber : public IBoundaryNumberProvider<dim>
{
	public:
	//	Functor Type
		typedef typename IBoundaryNumberProvider<dim>::functor_type functor_type;

	//	return functor
		virtual functor_type get_functor() const {return boost::ref(*this);}

	public:
		ConstBoundaryNumber() {m_Number = 0.0;}

		void set(number val)
		{
			m_Number = val;
		}

		void print() const
		{
			UG_LOG("ConstBoundaryNumber:" << m_Number << "\n");
		}

		bool operator() (number& c, const MathVector<dim>& x, number time = 0.0) const
		{
			c = m_Number;
			return true;
		}

	protected:
		number m_Number;
};

}

#endif /* __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__USER_DATA__ */
