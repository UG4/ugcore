#include "domain.h"
#include "common/util/table.h"

using namespace std;

namespace ug{

std::string DomainInfo::
to_string() const
{
	if((m_numElems.size() != m_numLocalElems.size())
	   || (m_numElems.size() != m_numLocalGhosts.size()))
	{
		UG_THROW("elem-arrays have to have the same number of entries!");
	}

	StringStreamTable t;

	t(0, 0) << "lvl";
	t(0, 1) << "#total-elems";
	t(0, 2) << "#local-elems";
	t(0, 3) << "(% of total)";
	t(0, 4) << "#local-ghosts";

	for(size_t i = 0; i < m_numElems.size(); ++i){
		int r = i+1;
		t(r, 0) << i;
		t(r, 1) << m_numElems[i];
		t(r, 2) << m_numLocalElems[i];
		if(m_numElems[i] > 0)
			t(r, 3) << (float)m_numLocalElems[i] / (float)m_numElems[i];
		else
			t(r, 3) << "-";
		t(r, 4) << m_numLocalGhosts[i];
	}

	return t.to_string();
}

}//	end of namespace
