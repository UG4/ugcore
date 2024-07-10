#ifndef __SUPPORT_H__
#define __SUPPORT_H__

#include <iostream>
#include <stdexcept>
#include <cmath>
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <assert.h>
#include <string>
#include <sstream>
#include <utility>
#include <vector>
#include <type_traits>


 //////////////////////////////////////////////////////////////////////////////

// class to help count and store a bool and a number of templete type
// comparable to std::pair<bool,int> but more dedicated to the specific aim
 
 template< typename T >
 class VertexFractureProperties
 {
 public:

 	VertexFractureProperties( bool isBndFracVertex,  T numberCrossingFracsInVertex )
 	: m_isBndFracVertex(isBndFracVertex), m_numberCountedFracsInVertex(numberCrossingFracsInVertex)
 	{
 	};


 	VertexFractureProperties()
 	: VertexFractureProperties( false, 0 )
 	{
 	};

 	void setIsBndFracVertex( bool iBDV = true )
 	{
 		m_isBndFracVertex = iBDV;
 	}

 	void setNumberCrossingFracsInVertex( T const & nCFIV )
 	{
 		m_numberCountedFracsInVertex = nCFIV;
 	}

 	bool getIsBndFracVertex()
 	{
 		return m_isBndFracVertex;
 	}

// 	T getCountedNumberFracsInVertex()
// 	{
// 		return m_numberCountedFracsInVertex;
// 	}


 	T getNumberFracEdgesInVertex()
 	{
 		return m_numberCountedFracsInVertex;
 	}

// 	T getNumberCrossingFracsInVertex()
//  {
// 		if( m_isBndFracVertex )
// 			return m_numberCountedFracsInVertex;
//
// 		// for inner vertices, each edge passed when
// 		// fractures are counted along their edges
// 		// that the vertizes get hit twice for each fracture run
// 		// only for boundary vertices, this happens only once per fracture
// 		T multipeInnerHits = 2;
//
// 		T rest = m_numberCountedFracsInVertex % multipeInnerHits;
//
// 		if( rest != 0 )
// 		{
//// 			UG_THROW("Expand layers: rest division frac counting not zero " << m_numberCountedFracsInVertex << std::endl);
//
// 			throw std::runtime_error("error");
//
// 			return 0;
// 		}
//
// 		return m_numberCountedFracsInVertex / multipeInnerHits;
//  }

 	VertexFractureProperties & operator++( int a )
 	{
 		m_numberCountedFracsInVertex++;
 		return *this;
 	}


 private:
 	bool m_isBndFracVertex;
 	T m_numberCountedFracsInVertex;
 };
 
 //////////////////////////////////////////////////////////////////////////////
 
// a class to store a matrix with two indices
 
template< typename I, typename D,
typename std::enable_if<std::is_integral<I>::value, int>::type = 0,
typename std::enable_if<std::is_arithmetic<D>::value, int>::type = 0
//, typename std::conditional_t<std::is_same_v<D, bool>, char, D>
	>
class MatrixTwoIndices
{
public:

   // standard constructor
   MatrixTwoIndices() : x_degree(-1), y_degree(-1) {};

   // constructor, wants to know the degrees
   MatrixTwoIndices( I _x_degree_, I _y_degree_, D defVal = 0 )
          : x_degree(_x_degree_), y_degree(_y_degree_)
 	{ values = std::vector<D>(  (_x_degree_)*(_y_degree_),  defVal  ); }

   // asking for a special element, cout << object(i,j) ...
   D const operator()( I i, I j ) const
   {
 	  assert( x_degree > 0 &&  y_degree > 0 );
      return values[ j*(x_degree) + i ];
   }

   // giving a special element a special value , object(i,j) = xx -> values[...] = xx
   D & operator()( I i, I j )
   {
 	  assert( x_degree > 0 &&  y_degree > 0 );
      return values[ j*(x_degree) + i ];
   }

private:

   std::vector<D> values;
   I x_degree, y_degree; // the degree of the polynom in x and y

 };

//	std::vector<double> minDist2Center(fracInfos.size(), std::numeric_limits<double>);
// matrix, wo auch index drin ist der subdom

//	MatrixTwoIndices<IndexType,double> mat_fracInd_minFacePep( fracInfos.size(), 100 );
//	MatrixTwoIndices<IndexType,double> mat_fracInd_minFacePep( fracInfos.size(), std::numeric_limits<double>::max() );
	//MatrixTwoIndices mat_fracInd_minFacePep( fracInfos.size(), std::numeric_limits<double> );

//	MatrixTwoIndices<IndexType,double> mat_fracInd_minFacePep( 30, 20, std::numeric_limits<double>::max() );
//
//	class bla{
//
//	};
//
//	bla blubb;
//
//	MatrixTwoIndices<IndexType, bla> mat_by( 30, 20, blubb );


/////////////////////////////////////////////////////////////////////////////


template <class T>
class T_min
{
 
public:
  // constructor, initializes minval (cause the first told is in this
  // moment the maximum) 
  T_min( T val ) : minval( val ) {}; 

  // tells the minimal value
  T const operator()() const  { return minval; };
 
  // wants to know values, saves the minimum
  void operator()(T val) 
  { if( minval > val )  minval = val; };


protected:

  T minval;
  
private:
  
  T_min() {};
};

//////////////////////////////////////////////////////////////////////


template <
typename ECKENTYP,
typename GESICHTSTYP, 
typename SENKRECHTENTYP 
>
class VertexFractureTriple
{

public:
	
	VertexFractureTriple( ECKENTYP const & edge, GESICHTSTYP const & face, SENKRECHTENTYP const & normal   )
	: m_edge(edge), m_face(face), m_normal(normal)
	{	
	};

	ECKENTYP const getEdge() const { return m_edge; } 
	
	GESICHTSTYP const getFace() const { return m_face; } 
	
	SENKRECHTENTYP const getNormal() const { return m_normal; } 

private:
	
	ECKENTYP m_edge;
	GESICHTSTYP m_face;
	SENKRECHTENTYP  m_normal;
	
	VertexFractureTriple()
	{};
	
};


//////////////////////////////////////////////////////////////////

template < typename VRT, typename IndTyp > //, typename EDG >
class CrossingVertexInfo
{
public:

	CrossingVertexInfo( VRT const & crossVrt, IndTyp numbCrossFracs )
	: m_crossVrt(crossVrt), m_numbCrossFracs( numbCrossFracs ), m_vecShiftedVrts(std::vector<VRT>())
	//,m_vecOrigEdges(std::vector<EDG>())
	{

	}

	VRT getCrossVertex() const { return m_crossVrt; }

	IndTyp getNumbCrossFracs() const { return m_numbCrossFracs; }

	void addShiftVrtx( VRT const & vrt )
	{
		m_vecShiftedVrts.push_back(vrt);
	}

//	void addOriginalFracEdge( EDG const & edg ) { m_vecOrigEdges.push_back(edg); }

	std::vector<VRT> getVecShiftedVrts() const { return m_vecShiftedVrts; }

//	std::vector<EDG> getVecOrigFracEdges() const { return m_vecOrigEdges; }

private:

	VRT m_crossVrt;
	IndTyp m_numbCrossFracs;
	std::vector<VRT> m_vecShiftedVrts;
//	std::vector<EDG> m_vecOrigEdges;
};


//////////////////////////////////////////////////////////////////



#endif

