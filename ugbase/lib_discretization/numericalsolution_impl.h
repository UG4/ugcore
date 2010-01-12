/*
 * numericalsolution_impl.h
 *
 *  Created on: 03.12.2009
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__NUMERICALSOLUTION_IMPL__
#define __H__LIBDISCRETIZATION__NUMERICALSOLUTION_IMPL__

namespace ug{

template<typename TElem>
bool DoFPattern::get_indices(TElem* elem, unsigned int nr_solution, int* ind)
{
	assert(nr_solution < _SingleSolutionInfoVec.size());
	int subsetIndex = _sh->get_subset_index(elem);
	int i;
	for(i = 0; i < _SingleSolutionInfoVec[nr_solution]->num_SubsetIndex; ++i)
	{
		if(subsetIndex == _SingleSolutionInfoVec[nr_solution]->SubsetIndex[i]) break;
	}
	if(i ==  _SingleSolutionInfoVec[nr_solution]->num_SubsetIndex) assert(0 && "ERROR in get_index. Solution not defined for SubsetIndex.");

	int ncomp = _SingleSolutionInfoVec[nr_solution]->group_Comp[subsetIndex];

	typename geometry_traits<TElem>::Descriptor TDesc;
	for(unsigned int i=0; i< TDesc.num_vertices(); i++)
	{
		VertexBase* vert = elem->vertex(i);
		ind[i] = _aaDoFVRT[vert].nr + ncomp;
	}
	return true;
}



////////////////////////////
// numerical solution
////////////////////////////

template <int d>
template <typename TElem>
const TrialSpace<TElem>& NumericalSolution<d>::get_TrialSpace(int nr_func)
{
	return TrialSpaces<TElem>::TrialSpace(_pattern->get_TrialSpaceType(nr_func));
}

template <int d>
template<typename TElem>
bool NumericalSolution<d>::get_indices(TElem* elem, int nr_solution, int* ind)
{
	_pattern->get_indices(elem, nr_solution, ind);
	return true;
}

template <int d>
template <typename TElem>
bool NumericalSolution<d>::get_local_DoFValues(TElem* elem, int nr_func, number* DoFValues)
{
	typename geometry_traits<TElem>::Descriptor TDesc;

	int nvalues = TDesc.num_vertices();

	int* indices = new int[nvalues];

	for(int i=0; i< nvalues; i++)
	{
		VertexBase* vert = elem->vertex(i);
		indices[i] = (int) (_pattern->get_index(vert, nr_func));
	}

	int level = 0;
	if(_pattern->num_levels() >= 2)
	{
		MultiGrid* mg = dynamic_cast<MultiGrid*>(_pattern->get_assigned_subset()->get_assigned_grid());
		assert(mg != NULL && "ERROR in get_local_DoFValues: num_levels > 0 but not a MultiGrid. Aborting.");
		level = mg->get_level(elem);
	}

	double *doubleDoFValues = new double[nvalues];
	if(_GridVector[level]->get_values(nvalues, indices, doubleDoFValues) == false) return false;

	for(int i=0; i< nvalues; i++)
	{
		DoFValues[i] = (number) doubleDoFValues[i];
	}


	delete[] indices;
	delete[] doubleDoFValues;
	return true;
}

template <int d>
NumericalSolution<d>::NumericalSolution(std::string shortname, std::string longname, std::string description, DoFPattern& pattern, Domain<d>& domain)
{
 	_shortname = shortname;
 	_longname = longname;
 	_description = description;
 	_domain = &domain;

	_pattern = &pattern;

	_GridVector.resize(_pattern->num_levels());
	for(uint i = 0; i < _pattern->num_levels(); ++i)
	{
		_GridVector[i] = new Vector();
		_GridVector[i]->create_vector(_pattern->num_doubles(i));
	}

	return;
}

template <int d>
bool NumericalSolution<d>::assign(const NumericalSolution<d>& v)
{
	assert(v.get_pattern()->num_levels() == _pattern->num_levels());

	for(uint lev = 0; lev < _pattern->num_levels(); ++lev)
	{
		*(this->_GridVector[lev]) = *(v._GridVector[lev]);
	}
	return true;
}


template <int d>
bool NumericalSolution<d>::set_values(bool (*fct)(MathVector<d>, number&), int nr_func, ISubsetHandler& sh, uint subsetIndex)
{
	return set_values(fct, nr_func, sh, subsetIndex, 0);
}

template <int d>
bool NumericalSolution<d>::set_values(bool (*fct)(MathVector<d>, number&), int nr_func, ISubsetHandler& sh, uint subsetIndex, int level)
{
	VertexBaseIterator iterBegin, iterEnd, iter;

	Grid* grid = sh.get_assigned_grid();
	Grid::VertexAttachmentAccessor<Attachment<MathVector<d> > > aaPos(*grid, *(_domain->get_position_attachment()));

	MultiGrid* mg = dynamic_cast<MultiGrid*>(grid);

	if(mg == NULL)
	{
		if(level != 0)
		{
			std::cout << "Numerical Solution is only defined on a Grid, no MultiGrid assignment possible." << std::endl;
			return false;
		}
		iterBegin = grid->begin<VertexBase>();
		iterEnd = grid->end<VertexBase>();
	}
	else
	{
		iterBegin = mg->begin<VertexBase>(level);
		iterEnd = mg->end<VertexBase>(level);
	}

	int index;
	number val;
	MathVector<d> corner;
	std::vector<int> indices;
	std::vector<number> values;

	/* loop over all Vertices */
	int nvertices = 0;
	for(iter = iterBegin; iter != iterEnd; iter++)
	{
		VertexBase *vert = *iter;
		if((uint)sh.get_subset_index(vert) != subsetIndex) continue;

		corner = aaPos[vert];
		index = _pattern->get_index(vert, nr_func);
		fct(corner, val);
		values.push_back(val);
		indices.push_back(index);

		nvertices++;
	}

	double* valueArray = new double[indices.size()];
	int* indexArray = new int[indices.size()];

	for(uint i=0; i< indices.size(); i++)
	{
		valueArray[i] = values[i];
		indexArray[i] = indices[i];
	}

	if(_GridVector[level]->set_values((int) indices.size(), indexArray, valueArray) != true)
		return false;
	delete valueArray;
	delete indexArray;

	return true;
}

template <int d>
Vector* NumericalSolution<d>::GridVector(int level)
{
	return _GridVector[level];
}

template <int d>
Vector* NumericalSolution<d>::GridVector()
{
	return _GridVector[0];
}

template <int d>
DoFPattern* NumericalSolution<d>::get_pattern()
{
	return _pattern;
}

template <int d>
TrialSpaceType NumericalSolution<d>::get_TrialSpaceType(int nr_func)
{
	return _pattern->get_TrialSpaceType(nr_func);
}

template <int d>
bool NumericalSolution<d>::print()
{
	std::cout << "shortname: \"" << _shortname << "\", ";
	std::cout << "longname: \"" << _longname <<"\", ";
	std::cout << "description: \"" << _description << "\"" << std::endl;
	return true;
}

template <int d>
NumericalSolution<d>::~NumericalSolution()
{
	for(uint i = 0; i < _pattern->num_levels(); ++i)
	{
		delete _GridVector[i];
	}

}

template <int d>
std::string NumericalSolution<d>::get_name(int comp)
{
	return _pattern->get_name(comp);
}



}


#endif /* __H__LIBDISCRETIZATION__NUMERICALSOLUTION_IMPL__ */
