#ifndef __H__LIB_ALGEBRA__AMG__FAMG_STRUCTS_H__
#define __H__LIB_ALGEBRA__AMG__FAMG_STRUCTS_H__

namespace ug{

template<typename value_type>
struct s_interpolation
{
	size_t from;
	value_type value;
};

struct neighborstruct2
{
	void print()
	{
		UG_LOG(parents.size() << " parents, F = " << F << ": ");
		for(size_t i=0; i<parents.size(); i++)
			UG_LOG((i>0 ? "," : "") << "p" << i << ": [" << GetOriginalIndex(parents[i].from) << "] -> " << parents[i].value);
		UG_LOG(std::endl);
	}
	FixedArray1<s_interpolation<double>, 2> parents;
	double F;
};

struct neighborstruct_var
{
	void print()
	{
		UG_LOG(parents.size() << " parents, F = " << F << ": ");
		for(size_t i=0; i<parents.size(); i++)
			UG_LOG((i>0 ? "," : "") << "p" << i << ": [" << GetOriginalIndex(parents[i].from) << "] -> " << parents[i].value);
		UG_LOG(std::endl);
	}
	stdvector<s_interpolation<double> > parents;
	double F;
};


}

#endif
