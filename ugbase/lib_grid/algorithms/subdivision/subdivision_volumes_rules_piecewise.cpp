// created by mstepnie
// martin.stepniewski@gcsc.uni-frankfurt.de
// Juli 14, 2014

#include <cmath>
#include "subdivision_rules_piecewise_loop.h"

using namespace std;

namespace ug
{

SubdivRules_PLoop::
SubdivRules_PLoop()
{
//	precalculate betas
//	the number is quite arbitrary here.
	size_t numPrecals = 16;
	m_betas.resize(numPrecals);
	for(size_t i = 0; i < numPrecals; ++i)
		m_betas[i] = calculate_beta(i);
}

SubdivRules_PLoop::
SubdivRules_PLoop(const SubdivRules_PLoop& src)
{
//	since this method won't ever be executed it can stay empty
}

SubdivRules_PLoop& SubdivRules_PLoop::
operator=(const SubdivRules_PLoop& src)
{
//	since this method won't ever be executed it can stay empty
	return *this;
}


number SubdivRules_PLoop::
get_beta(size_t valency)
{
	if(valency < m_betas.size())
		return m_betas[valency];
		
	return calculate_beta(valency);
}

number SubdivRules_PLoop::
calculate_beta(size_t valency)
{
	if(valency == 6)
		return 0.0625;
		
	if(valency > 0){
		const number tmp = 0.375 + 0.25 * cos((2.0*PI)/(number)valency);
		return (0.625 - tmp*tmp)/(number)valency;
	}

	return 0;
}


}//	end of namespace
