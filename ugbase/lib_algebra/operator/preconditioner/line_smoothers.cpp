#include "line_smoothers.h"
namespace ug{
template<int dim>
bool ComparePosDimYDir(const std::pair<MathVector<dim>, size_t> &p1,
					   const std::pair<MathVector<dim>, size_t> &p2)
{return false;}

template<>
bool ComparePosDimYDir<1>(const std::pair<MathVector<1>, size_t> &p1,
						  const std::pair<MathVector<1>, size_t> &p2)
{
	return (p1.first[0]<p2.first[0]);
};

// Order for 2D
template<>
bool ComparePosDimYDir<2>(const std::pair<MathVector<2>, size_t> &p1,
						  const std::pair<MathVector<2>, size_t> &p2)
{
	if (p1.first[0]==p2.first[0]) {
		return (p1.first[1] < p2.first[1]);
	}
	else {
		return (p1.first[0] < p2.first[0]);
	}
};

// Order for 3D
template<>
bool ComparePosDimYDir<3>(const std::pair<MathVector<3>, size_t> &p1,
						  const std::pair<MathVector<3>, size_t> &p2)
{
	if (p1.first[2]==p2.first[2]){
		if (p1.first[0]==p2.first[0]) {
			return (p1.first[1] < p2.first[1]);
		}
		else {
			return (p1.first[0] < p2.first[0]);
		}
	}
	else{
		return (p1.first[2] < p2.first[2]);
	}
};




template<int dim>
bool ComparePosDimZDir(const std::pair<MathVector<dim>, size_t> &p1,
					   const std::pair<MathVector<dim>, size_t> &p2)
{return false;}

// Order for 1D
template<>
bool ComparePosDimZDir<1>(const std::pair<MathVector<1>, size_t> &p1,
						  const std::pair<MathVector<1>, size_t> &p2)
{
	return (p1.first[0]<p2.first[0]);
};

// Order for 2D
template<>
bool ComparePosDimZDir<2>(const std::pair<MathVector<2>, size_t> &p1,
						  const std::pair<MathVector<2>, size_t> &p2)
{
	if (p1.first[0]==p2.first[0]) {
		return (p1.first[1] < p2.first[1]);
	}
	else {
		return (p1.first[0] < p2.first[0]);
	}
};

// Order for 3D
template<>
bool ComparePosDimZDir<3>(const std::pair<MathVector<3>, size_t> &p1,
						  const std::pair<MathVector<3>, size_t> &p2)
{
	if (p1.first[1]==p2.first[1]){
		if (p1.first[0]==p2.first[0]) {
			return (p1.first[2] < p2.first[2]);
		}
		else {
			return (p1.first[0] < p2.first[0]);
		}
	}
	else{
		return (p1.first[1] < p2.first[1]);
	}
};

} // namespace ug
