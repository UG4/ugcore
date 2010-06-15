#ifndef __BOOST_BLOCK_HH__
#define __BOOST_BLOCK_HH__

#include <boost/timer.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/vector_sparse.hpp>
#include <boost/numeric/ublas/vector_of_vector.hpp>
#include <boost/numeric/ublas/vector_expression.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/lu.hpp>
 
 #include <boost/numeric/ublas/expression_types.hpp>
 #include <boost/type_traits.hpp>
namespace ublas = boost::numeric::ublas;


//! A fixed size vector
template <typename T, int N>
class block_vector : public ublas::vector<T, ublas::bounded_array<T,N> >
{
  typedef ublas::vector<T, ublas::bounded_array<T,N> > Base_vector;
public:
  // Default construction
  block_vector() : Base_vector(N) 
  {}
  // Construction and assignment from a uBLAS vector expression or copy assignment
  template <class R> block_vector (const ublas::vector_expression<R>& r) : Base_vector(r)
  {}
  template <class R> void operator=(const ublas::vector_expression<R>& r)
  {
    Base_vector::operator=(r);
  }
  template <class R> void operator=(const Base_vector& r)
  {
    Base_vector::operator=(r);
  }
};




//! A fixed size matrix
template <typename T, int M, int N=M>
class block_matrix : public ublas::compressed_matrix<T, ublas::row_major, 0, ublas::bounded_array<std::size_t,M*N>, ublas::bounded_array<T,M*N> >
{ 
  typedef ublas::compressed_matrix<T, ublas::row_major, 0, ublas::bounded_array<std::size_t,M*N>, ublas::bounded_array<T,M*N> > base_matrix;

public:
  // Default construction
  block_matrix() : base_matrix(M,N) 
  {}
  // Construction and assignment from a uBLAS vector expression or copy assignment
  template <class R> block_matrix (const ublas::matrix_expression<R>& r) : base_matrix(r)
  {}
  template <class R> void operator=(const ublas::matrix_expression<R>& r)
  {
    base_matrix::operator=(r);
  }
  template <class R> void operator=(const base_matrix& r)
  {
    base_matrix::operator=(r);
  }
};



//
// determine level of recursion
// todo: implement as a metafunction!
//
/*

// empty dummy policy
template <typename VE, bool k>
struct recursive_traits_dummy{
  enum{depth = 0};
};

//! general type to define the level of recursion
template <typename T> 
struct recursive_traits : public
recursive_traits_dummy<T, boost::is_arithmetic<typename T::value_type>::value > {};

//dummy: recursive definition
template <typename VE>
struct recursive_traits_dummy<VE, false> 
{
  enum{depth = (recursive_traits<typename VE::value_type>::depth) + 1};
};

//dummy: specialization for leave
template <typename VE>
struct recursive_traits_dummy<VE, true> 
{
  enum{depth = 1};
};

//! general type to define the level of recursion
template <typename T> 
struct recursive_depth {
  enum {value = 
	recursive_traits_dummy<T, (boost::is_arithmetic<typename T::value_type>::value) >::depth };
};
*/

/*
//
// inner products (recursive)
//


// meta function for block recursive inner product 
// general implementation 
template <typename X, int L>
struct frec_inner_prod 
{
  template <typename E1, typename E2> 
  X apply(const E1 &x, const E2 &y)
  {
    // typedef typename E1::value_type result_type; 
    typedef typename E1::value_type value_type;
    typedef typename E1::size_type vector_size_type;
    
    vector_size_type size = x.size();
    X t=X(0);
    frec_inner_prod<X, L-1> frec;
    for (vector_size_type i=0; i<size; ++i)
      t += frec.apply (x(i), y(i) );
    return X(t);
    
  }
};

// meta function for block recursive inner product 
// specialized implementation 
template <typename X>
struct frec_inner_prod<X,1> 
{
  template <typename E1, typename E2> 
  X apply(const E1 &x, const E2 &y)
  {return inner_prod (x, y);}
};


// actual function call
template <typename X, typename E1, typename E2>
X r_inner_prod(const E1 &x, const E2 &y)
{
  frec_inner_prod<double, recursive_depth<E1>::value > frec;
  return frec.apply(x,y);
};

*/









  
#endif
