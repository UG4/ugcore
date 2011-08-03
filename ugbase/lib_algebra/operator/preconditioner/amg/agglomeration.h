#include <vector>
#include <map>

namespace ug
{

void EasyAgglomeration(const std::vector<size_t> sizes,
		const std::vector<std::map<int, size_t> > connections,
		std::vector<std::vector<int> > &mergeWith,
		size_t minimalSize, size_t preferredSize);

}
