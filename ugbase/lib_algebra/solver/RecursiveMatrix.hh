#ifndef __RECURSIVE_MATRIX_HH__
#define __RECURSIVE_MATRIX_HH__


#include <boost/type_traits.hpp>
#include "RecursiveTraits.hh"
#include "BoostBlock.hh"
//
// matrix vector operations
//

// meta function for block recursive matrix vector product 
// general implementation 
template <int L>
struct frec_mv_assign 
{
  template <typename M, typename X, typename Y> 
  static void apply(const M &A, const X &x, Y& y)
  {
    frec_mv_assign<L-1> frec_mv;

    typedef typename Y::value_type value_type;
    
    typedef typename M::const_iterator1 row_iter_type;
    typedef typename M::const_iterator2 entry_iter_type;
    
    for (row_iter_type rowi = A.begin1(); rowi != A.end1(); ++rowi) 
      {
	entry_iter_type iter_ij = rowi.begin();	
	const entry_iter_type iter_end = rowi.end();

	// todo: avoid using a temporary
	value_type yi; 
	// todo: avoid stupid call
	y(iter_ij.index1()) *= 0.0;

	while (iter_ij != iter_end)
	  {
	    //std::cout << iter_ij.index1() << " "<< iter_ij.index2()<< std::endl;
	    frec_mv_assign<L-1>::apply((*iter_ij),  x(iter_ij.index2()), yi);
	    y(iter_ij.index1()) += yi;
	    ++iter_ij;
	  }
      }
    
    return;
  }
};



// meta function for block recursive matrix vector product
// specialized implementation 
template <>
struct frec_mv_assign<1> 
{
  template <typename M, typename X, typename Y> 
  static void apply(const M &A, const X &x, Y &y)
  {axpy_prod(A, x, y, true); return;}
};


// meta function for block recursive matrix vector product 
// general implementation 
template <int L>
struct frec_mvsub 
{
  template <typename M, typename X, typename Y> 
  static void apply(const M &A, const X &x, Y& y)
  {
    typedef typename Y::value_type value_type;
    
    typedef typename M::const_iterator1 row_iter_type;
    typedef typename M::const_iterator2 entry_iter_type;
    
    for (row_iter_type rowi = A.begin1(); rowi != A.end1(); ++rowi) 
      {
	entry_iter_type iter_ij = rowi.begin();	
	const entry_iter_type iter_end = rowi.end();

	// todo: avoid using a temporary
	value_type &yi = y(iter_ij.index1());

	while (iter_ij != iter_end)
	  {
	    //  std::cout << iter_ij.index1() << " "<< iter_ij.index2()<< std::endl;
	    frec_mvsub<L-1>::apply((*iter_ij),  x(iter_ij.index2()), yi);
	    ++iter_ij;
	  }
      }
    
    return;
  }
};

// specialized implementation 
template <>
struct frec_mvsub<1> 
{
  template <typename M, typename X, typename Y> 
  static void apply(const M &A, const X &x, Y &y)
  {axpy_prod (A, -x, y, false); return;}
};

template <>
struct frec_mvsub<0> 
{
  template <typename M, typename X, typename Y> 
  static void apply(const M &A, const X &x, Y &y)
  {axpy_prod (A, -x, y, false); return;}

  static void apply(const double A, const double x, double &y)
  {y=A*x; return;}
};
/*
template <>
struct frec_mvsub<0> 
{
    static void apply(const double A, const double x, double &y)
  {y -=A*x; return;}
};

*/

///////////////////////////////////////////////////////
//
// recursive diagonal solver
//
///////////////////////////////////////////////////////


// general implementation 
template <int L>
struct frec_dsolve
{
  template <typename M, typename X, typename Y> 
  static void dsolve(const M &A, X &x, const Y& y)
  {
    
    BOOST_STATIC_ASSERT((boost::is_same<typename X::value_type, 
			 typename Y::value_type>::value));

    typedef typename M::const_iterator1 row_iterator;    

    row_iterator end = A.end1();
    for (row_iterator rowi = A.begin1(); rowi != end; ++rowi) 
      {	
	const size_t i =rowi.index1();
	frec_dsolve<L-1>::dsolve(A(i,i), x(i), y(i));
      }
    return;
  }
};



// here we solve up to level zero!
template <>
struct frec_dsolve<0> 
{
  template <typename M, typename X, typename Y> 
  static void dsolve(const M &A, X &x, const Y &y)
  { mvsolve (A, x, y); return;}
};


///////////////////////////////////////////////////////
//
// recursive lower part solver
//
///////////////////////////////////////////////////////

// general implementation 
template <int L>
struct frec_lsolve
{
  template <typename M, typename X, typename Y> 
  static void lsolve(const M &A, X &x, const Y& y)
  {
    
    BOOST_STATIC_ASSERT((boost::is_same<typename X::value_type, 
			 typename Y::value_type>::value));
    BOOST_STATIC_ASSERT((L>=1));

    typedef typename M::const_iterator1 row_iterator;    
    typedef typename M::const_iterator2 entry_iterator;

    row_iterator end = A.end1();
    for (row_iterator rowi = A.begin1(); rowi != end; ++rowi) 
      {	
	const size_t ind =rowi.index1();
	typedef typename Y::value_type vvalue_type;
	vvalue_type  temp(y(ind));
	
	entry_iterator iter_ij = rowi.begin();
	int jind=iter_ij.index2();
	for (jind=iter_ij.index2(); jind<ind; )
	  {
	    mvsub(*iter_ij,  x(jind), temp);
	    iter_ij++;
	      //frec_mvsub<recursive_depth<vvalue_type>::value>::apply(*iter_ij, x(jind), temp);
	  }
	
	frec_lsolve<L-1>::lsolve(*iter_ij, x(ind), temp);
      }
    return;
  }

};

// here we solve up to level zero!
template <>
struct frec_lsolve<0> 
{
  template <typename M, typename X, typename Y> 
  static void lsolve(const M &A, X &x, const Y &y)
  { mvsolve (A, x, y); return;}
};




template <typename M, typename X, typename Y>
void mvsolve(const M &A, X &x, const Y &y)
{
  // must specify a specialized routine for your data types!
  // BOOST_STATIC_ASSERT((false));
  assert(0);
  return;
};

void mvsolve(const double A, double &x, const double y)
{
  x = y/A; 
};



template<typename X, int N>
void mvsolve(const block_matrix<X,N,N> &A, block_vector<X,N> &x, const block_vector<X,N> &y)
{
  
};

template<typename T>
void mvsolve(const block_matrix<T,2,2> &A, block_vector<T,2> &x, const block_vector<T,2> &y)
{
  T detA (A(0,0)*A(1,1)- A(0,1)*A(1,0) ); 
  x(0) = ( A(1,1) * y(0) - A(0,1) * y(1)) / detA;
  x(1) = ( -A(1,0) * y(0) + A(0,0) * y(1)) / detA;
};

template<typename T>
void mvsolve(const block_matrix<T,3,3> &A, block_vector<T,3> &x, const block_vector<T,3> &y)
{
  //BOOST_STATIC_ASSERT((false));
  T detA = 0;
  assert(0);
};


// user accessible functions
template <typename M, typename X, typename Y>
void mv_prod(const M &A, const X &x, Y& y) 
{
  BOOST_STATIC_ASSERT(((recursive_depth<X>::value) == (recursive_depth<Y>::value) ));
  frec_mv_assign<recursive_depth<X>::value>::apply(A, x, y);
  return;
} 


template <typename M, typename X, typename Y>
void mvsub(const M &A, const X &x, Y& y) 
{
  BOOST_STATIC_ASSERT((recursive_depth<X>::value == recursive_depth<Y>::value ));
  frec_mvsub<recursive_depth<X>::value>::apply(A,x,y);
  return;
}  

void mvsub(const double A, const double x, double& y) 
{
  y-=A*x;
 }  




template <typename M, typename X, typename Y>
void mv_dsolve(const M &A, X &x, const Y& y) 
{
  BOOST_STATIC_ASSERT(((recursive_depth<X>::value) == (recursive_depth<Y>::value) ));
  frec_dsolve<recursive_depth<X>::value>::dsolve(A,x,y);
  return;
}  

template <int L, typename M, typename X, typename Y>
void mv_dsolve_rec(const M &A, X &x, const Y& y) 
{  
  BOOST_STATIC_ASSERT(((recursive_depth<X>::value) >= L ));
  frec_dsolve<L>::dsolve(A,x,y);
  return;
}  

template <typename M, typename X, typename Y>
void mvlsolve(const M &A, X &x, const Y& y) 
{
  BOOST_STATIC_ASSERT(((recursive_depth<X>::value) == (recursive_depth<Y>::value) ));
  frec_lsolve<recursive_depth<X>::value>::lsolve(A,x,y);
  return;
}  

template <int L, typename M, typename X, typename Y>
void mvlsolveRec(const M &A, X &x, const Y& y) 
{  
  BOOST_STATIC_ASSERT(((recursive_depth<X>::value) >= L ));
  frec_lsolve<L>::lsolve(A,x,y);
  return;
}  

#endif
