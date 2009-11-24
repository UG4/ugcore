#ifndef __RECURSIVE_VECTOR_HH__
#define __RECURSIVE_VECTOR_HH__

#include "RecursiveTraits.hh"

// general recursive definition
template<typename X, int L=recursive_depth<X>::value >
struct r_innerproduct
{
  typedef typename r_innerproduct<typename X::value_type, L-1>::result_type result_type;
  
  static result_type apply(const X &x, const X&y)
  {   
    typedef typename X::size_type vector_size_type;
    
    //   r_innerproduct<typename X::value_type, L-1> rip; 
    vector_size_type size = x.size();
    result_type s(0);
    
    for (vector_size_type i=0; i<size; ++i)
      {
	s +=  r_innerproduct<typename X::value_type, L-1>::apply(x(i), y(i));
      }
    
    return s;
  }
};


// specialization for leaf object
template<typename X>
struct r_innerproduct<X, 1>
{
  typedef double result_type;
  
  static result_type apply(const X &x, const X &y)
  { return inner_prod (x, y); }
};




template <typename X>
double GenInnerProduct(const X &x, const X &y)
{
  //  r_innerproduct<X, recursive_traits<X>::depth > frec;
  //return frec.compute(x,y);
  return r_innerproduct<X, recursive_traits<X>::depth >::apply(x,y);
};



#endif
