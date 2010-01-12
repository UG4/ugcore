/*
 * dirichletbndcond.h
 *
 *  Created on: 03.07.2009
 *      Author: andreasvogel
 */

#ifndef DIRICHLETBNDCOND_H_
#define DIRICHLETBNDCOND_H_

namespace ug{

template <int d>
class DirichletBNDCond{

public:
	typedef bool (*BNDCond)(MathVector<d>, number& val);

public:
	DirichletBNDCond(BNDCond func)
	{
		m_BNDFunc = func;
	}

	bool BNDValueFunction(MathVector<d> coord, number& val)
	{
		return m_BNDFunc(coord, val);
	}

private:
	BNDCond m_BNDFunc;
};

}


#endif /* DIRICHLETBNDCOND_H_ */
