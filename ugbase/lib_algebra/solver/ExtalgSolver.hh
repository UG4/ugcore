#ifndef __EXTALG_SOLVER_HH__
#define __EXTALG_SOLVER_HH__

// LinOp:              A: X->Y,
// Preconditioner:     M: Y->X,
// LinOpInverse:       M: Y->X

#include "RecursiveVector.hh"
#include "RecursiveMatrix.hh"

#include <boost/numeric/ublas/banded.hpp>


//
// Linear operators
// A: X -> Y
//


//! matrix free base class
template <class X, class Y>
class LinOp
{
  typedef X domain_type;
  typedef Y range_type;

  void apply(const X &x, Y &y);
  void addscaled(const X &x, Y &y);
};


template <class M, class X, class Y>
class MatrixLinOp
{
public:
  typedef M matrix_type;
  typedef X domain_type;
  typedef Y range_type;

  MatrixLinOp(const M &A) : m_mat (A)
  {}

  //! apply a matrix based lin op
  void apply(const X &x, Y& y)
  { mv_prod(m_mat, x, y); }

  void apply_sub(const X &x, Y& y)
  { mvsub(m_mat, x, y); }

private:
  const M &m_mat;
};


//
// preconditioner
// M: Y -> X
//


//! matrix dependent preconditioner
template <class M, class X, class Y>
class PJacobi
{
public:
  typedef M matrix_type;
  typedef X domain_type;

  typedef Y range_type;

  PJacobi(const M &A): m_mat(A) {}
  virtual ~PJacobi() {}

  void pre(X &c, Y &d);
  void apply(X &c, const Y &d);
  void post(X &c, Y &d);

private:
  const M &m_mat;
};


template <class M, class X, class Y>
void PJacobi<M,X,Y>::pre(X &c, Y &d)
{

}

template <class M, class X, class Y>
void PJacobi<M,X,Y>::apply(X &c, const Y &d)
{
  mv_dsolve(m_mat, c, d);
}

template <class M, class X, class Y>
void PJacobi<M,X,Y>::post(X &c, Y &d)
{

}

/*

//! matrix dependent preconditioner
template <class M, class X, class Y>
class PTransform
{
public:
  typedef M matrix_type;
  typedef X domain_type;

  typedef Y range_type;

  PTransform(const M &A): m_mat(A) {}
  virtual ~PTransform() {}

  void pre(X &c, Y &d);
  void apply(X &c, const Y &d);
  void post(X &c, Y &d);

private:
  const M &m_mat;
  const P1 &m_p11;
  const P2 &m_p22;
};


template <class M, class X, class Y>
void PTransform<M,X,Y>::pre(X &c, Y &d)
{

}

template <class M, class X, class Y>
void PTransform<M,X,Y>::apply(X &c, const Y &d)
{

}

template <class M, class X, class Y>
void PTransform<M,X,Y>::post(X &c, Y &d)
{

}


*/

//
// inverse linear operator
// \todo: must also provide a Preconditioner proxy!
//

template <class L, class P, class C>
class LinOpInverse
{
public:

  typedef L op_type;
  typedef P invop_type;


  LinOpInverse (L &op, P &inv, C &check)
    : m_lin_op(op), m_inv_op(inv), m_conv_check(check)
  {

  }

  template <class X, class Y>
  void apply(X &x, Y &y);


  L& LinOp() {return m_lin_op;}
  P& InvOp() {return m_inv_op;}
  //V& VectorSpace() {return m_vspace;}

  C& ConvCheck() {return m_conv_check;}

private:
  L &m_lin_op;
  P &m_inv_op;
  C &m_conv_check;
};


template <class L, class P, class C>
template <class X, class Y>
void LinOpInverse<L,P,C>::apply(X &x, Y &d)
{
  // compute defect (move somewhere else ?)
  LinOp().apply_sub(x,d);

  // create correction vector
  X c(x);

  // initialize iteration
  InvOp().pre(c, d);


  // status output (pre)
  ConvCheck().init(d);
  ConvCheck().status(std::cout);

  while (!ConvCheck().has_converged())
    {
      // approximate solution (one step)
      InvOp().apply(c, d);

      // update vector and rhs
      x+=c;
      LinOp().apply_sub(c,d);


      // next step and
      // status output (intermediate)
      ConvCheck().update(d);
      ConvCheck().status(std::cout);
    }



  // clean up
  InvOp().post(c, d);

  // status output (post)
  //std::cout << sqrt(GenInnerProduct(d,d)) << std::endl;
  return;
}

//
// default convergence measure
//

class StdConvCheck
{
public:
  StdConvCheck(int nsteps, double red, bool verbose)
    : m_max_steps(nsteps), m_reduction(red), m_verbose(verbose)
  {}

  template <typename Y>
  void init(const Y &y) {
    m_initial_defect = sqrt(GenInnerProduct(y,y));
    m_current_defect = m_initial_defect;
    m_step=0;
  };

  template <typename Y>
  void update(const Y &y){

    m_current_defect = sqrt(GenInnerProduct(y,y));
    m_step++;
  }

  bool has_converged() {
    return ((step() >=m_max_steps) || ( reduction() < m_reduction));
  }

  double reduction() {return m_current_defect / m_initial_defect;}
  int step() {return m_step;}

  template <typename OS>
  void status (OS &os){
   if(m_verbose) os << step() << ":" << m_current_defect << " " << reduction() << std::endl;
  }


private:
  double m_initial_defect;
  double m_current_defect;
  int m_step;

  int m_max_steps;
  double m_reduction;

  bool m_verbose;
};

#endif
