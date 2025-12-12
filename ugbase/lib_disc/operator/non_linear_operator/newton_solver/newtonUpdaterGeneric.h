/*
 * newtonUpdaterGeneric.h
 *
 *  Created on: 28.07.2021
 *      Author: Markus Knodel
 */

#ifndef UGCORE_UGBASE_LIB_DISC_OPERATOR_NON_LINEAR_OPERATOR_NEWTON_SOLVER_NEWTONUPDATERGENERIC_H_
#define UGCORE_UGBASE_LIB_DISC_OPERATOR_NON_LINEAR_OPERATOR_NEWTON_SOLVER_NEWTONUPDATERGENERIC_H_

//#include <ostream>
//#include <string>
//#include <vector>
//#include <cmath>
//#include <sstream>
//#include <type_traits>

//#include "common/common.h"

namespace ug {

template <typename TVector>
class NewtonUpdaterGeneric
{

public:

	virtual ~NewtonUpdaterGeneric() = default;


	using vector_type = TVector;

	virtual bool updateSolution( vector_type & sol, vector_type const & corr, bool signNegative = true )
	{
		if( signNegative )
			sol -= corr;
		else
			sol += corr;

		return true;
//		return (*this)(sol, 1, sol, - 1, corr );
	}


	virtual bool updateSolution( vector_type & solNew, number scaleOldSol,  vector_type const & solOld,
					             number scaleCorr, vector_type const & corr )
	{
		VecScaleAdd(solNew, scaleOldSol, solOld, scaleCorr, corr);

		return true;

	}

	virtual bool resetSolution( vector_type & resettedSol, vector_type const & oldSol )
	{
		resettedSol = oldSol;

		// reset not only u, but also internal variables if available!!!!!

		return true;
	}

	virtual bool tellAndFixUpdateEvents( vector_type const & sol )
	{
		return true;
	}

};

}


#endif

// muss noch im Newton und line search irgendwie automatisch so initialisiert werden, aber auf User-Wunsch ueberschrieben
// mit Hilfe von gettern und settern, damit entweder der line search oder der Newton selber sich das vom anderen holen koennen
// und dann muss es soweit registriert werden, wie es noetig ist, damit man von aussen zu greifen kann
// und die abgeleiteten Klassen brauchen ggf einen Konstruktor und so weiter
///	Vector type
// using vector_type = typename TAlgebra::vector_type;

//	using vector_type = TVector; // typename TAlgebra::vector_type;
//	using vec_typ_val_typ = typename vector_type::value_type;

// using vector_type = TVector; // typename TAlgebra::vector_type;
// using vec_typ_val_typ = typename vector_type::value_type;
// using value_type = typename vector_type::value_type;

//	template <typename VecTyp> //, typename ScalTyp >
//	template <typename ScalTyp> //, typename ScalTyp >
//	bool updateNewton( vector_type & solNew, value_type scaleOldSol,  vector_type const & solOld,
//			 	 	   value_type scaleCorr, vector_type const & corr )
//	bool updateNewton( VecTyp & solNew, ScalTyp scaleOldSol,  VecTyp const & solOld,
//			 	 	   ScalTyp scaleCorr, VecTyp const & corr )
//	bool updateNewton( VecTyp & solNew, typename VecTyp::value_type scaleOldSol,  VecTyp const & solOld,
//					   typename VecTyp::value_type scaleCorr, VecTyp const & corr )
//	bool updateNewton( vector_type & solNew, ScalTyp scaleOldSol,  vector_type const & solOld,
//					   ScalTyp scaleCorr, vector_type const & corr )
//		static_assert( std::is_same<ScalTyp, double>::value, "is ScalType double" );
//
//		static_assert( !std::is_same<ScalTyp, number>::value, "! is ScalType number" );
//
//		static_assert( std::is_same<ScalTyp, value_type>::value, "is it value_type" );
//
//		static_assert( std::is_same<double, value_type>::value, "is value_type double" );
//
//		static_assert( std::is_same<number, value_type>::value, "is value_type number" );

		// 	Standard: try on line u := u - lambda*p
		//	VecScaleAdd(u, 1.0, u, (-1)*lambda, p);

// using size_type = typename vector_type::size_type;
////	    using size_type = typename vector_type::size_type;
//
//		for(size_type i = 0; i < solNew.size(); ++i)
//		{
//			solNew[i] = scaleOldSol * solOld[i] + scaleCorr * corr[i];
//		}
