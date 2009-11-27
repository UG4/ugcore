/*
 * dirichletbndcond.h
 *
 *  Created on: 03.07.2009
 *      Author: andreasvogel
 */

#ifndef DIRICHLETBNDCOND_H_
#define DIRICHLETBNDCOND_H_

namespace ug{

class DirichletBNDCond{

public:
	typedef bool (*BNDCond)(MathVector<3>, number& val);

public:
	DirichletBNDCond(BNDCond func)
	{
		m_BNDFunc = func;
	}

	bool BNDValueFunction(MathVector<3> coord, number& val)
	{
		return m_BNDFunc(coord, val);
	}

private:
	BNDCond m_BNDFunc;
};

}


#endif /* DIRICHLETBNDCOND_H_ */
