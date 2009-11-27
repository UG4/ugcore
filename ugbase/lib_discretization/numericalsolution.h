/*
 * numericalsolution.h
 *
 *  Created on: 12.05.2009
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__NUMERICALSOLUTION__
#define __H__LIBDISCRETIZATION__NUMERICALSOLUTION__

#include <string>
#include <algorithm>
#include "../lib_grid/lib_grid.h"
#include "../lib_algebra/lib_algebra.h"
#include "dofhandler.h"
#include "trialspace.h"
#include "referenceelement.h"

namespace ug {

class NumericalSolutionPattern : public GridObserver{
	protected:
		struct SingleSolutionInfo {
			std::string name; // name of solution
			TrialSpaceType Ansatzspace; // associated trial space
			int num_SubsetIndex; // number of subsets where solution lives
			int* SubsetIndex; // subset indices where solution lives
			int* group_Comp; // nr of component in grouping per subset
		};

		struct GroupSolutionInfo {
			std::string name; // name of grouping
			int num_SingleSolutions; // number of solutions grouped here
			int* SingleSolutions; // index of single solutions contained in this grouping
			int num_SubsetIndex; // number of subsets where grouping lives
			int* SubsetIndex; // subset index where grouping lives
			int* num_Comp;    // number of Components per Subset
		};

	public:
		NumericalSolutionPattern(std::string name, SubsetHandler& sh)
		{
			_name = name;
			_sh = &sh;
			_num_sh = _sh->num_subsets();
			_lock = false;
		};

		bool add_solution(std::string name, TrialSpaceType Ansatzspace, int* SubsetIndex, int n)
		{
			if(_lock == true)
			{
				std::cout << "Pattern already locked. Can not add solution" << std::endl;
				return false;
			}

			// Create SingleSolutionInfo
			SingleSolutionInfo* info = new SingleSolutionInfo();
			info->name = name;
			info->Ansatzspace = Ansatzspace;
			info->num_SubsetIndex = n;
			info->SubsetIndex = new int[n];
			info->group_Comp = new int[n];
			for(int i = 0; i < n; ++i)
			{
				info->SubsetIndex[i] = SubsetIndex[i];
				info->group_Comp[i] = 0;
			}
			_SingleSolutionInfoVec.push_back(info);

			// Create Grouping with only on component: the SingleSolution
			GroupSolutionInfo* groupInfo = new GroupSolutionInfo();
			groupInfo->name = name;
			groupInfo->num_SingleSolutions = 1;
			groupInfo->SingleSolutions = new int[1];
			for(int i = 0; i < 1; ++i)
			{
				groupInfo->SingleSolutions[i] = _SingleSolutionInfoVec.size() - 1;
			}
			groupInfo->num_SubsetIndex = n;
			groupInfo->SubsetIndex = new int[n];
			groupInfo->num_Comp = new int[n];
			for(int i = 0; i < n; ++i)
			{
				groupInfo->SubsetIndex[i] = SubsetIndex[i];
				groupInfo->num_Comp[i] = 1;
			}
			_GroupSolutionInfoVec.push_back(groupInfo);

			return true;
		}

		bool group_solution(std::string name, int* GroupSolutions, int n)
		{
			if(_lock == true)
			{
				std::cout << "Pattern already locked. Can not add solution" << std::endl;
				return false;
			}

			// Sort GroupSolutions to be grouped
			std::sort(GroupSolutions, GroupSolutions + n);

			// collect subsets and num_Comp of all solutions that are to be grouped
			std::vector<int> SubsetIndex;
			std::vector<int> num_Comp;
			std::vector<int> SingleSolutions;
			for(int i = 0; i < n; ++i)
			{
				if(GroupSolutions[i] >= _GroupSolutionInfoVec.size() || GroupSolutions[i] < 0)
				{
					std::cout << "ERROR in group_solution: Trying to read Group " << GroupSolutions[i] << ", but grouping does not exist. Aborting" << std::endl;
					return false;
				}

				int nr = GroupSolutions[i];

				for(int j = 0; j < _GroupSolutionInfoVec[nr]->num_SubsetIndex; ++j)
				{
					int k;
					for(k = 0; k < SubsetIndex.size(); ++k)
					{
						if(SubsetIndex[k] == _GroupSolutionInfoVec[nr]->SubsetIndex[j])
						{
							num_Comp[k] += _GroupSolutionInfoVec[nr]->num_Comp[j];
							break;
						}
					}
					if(k == SubsetIndex.size())
					{
						SubsetIndex.push_back(_GroupSolutionInfoVec[nr]->SubsetIndex[j]);
						num_Comp.push_back(_GroupSolutionInfoVec[nr]->num_Comp[j]);
					}
				}
				for(int j = 0; j < _GroupSolutionInfoVec[nr]->num_SingleSolutions; ++j)
				{
					for(int k = 0; k < SingleSolutions.size(); ++k)
					{
						if(SingleSolutions[k] == _GroupSolutionInfoVec[nr]->SingleSolutions[j])
						{
							// This should never happen
							std::cout << "ERROR in group_solutions: Multiple occurance of a single solution in grouping. Internal error. Aborting" << std::endl;
							return false;
						}
					}
					SingleSolutions.push_back(_GroupSolutionInfoVec[nr]->SingleSolutions[j]);
				}
			}

			std::sort(SingleSolutions.begin(), SingleSolutions.end());

			// create new grouping
			GroupSolutionInfo* info = new GroupSolutionInfo();
			info->name = name;
			info->num_SubsetIndex = SubsetIndex.size();
			info->SubsetIndex = new int[SubsetIndex.size()];
			info->num_Comp = new int[SubsetIndex.size()];
			for(int i = 0; i < SubsetIndex.size(); ++i)
			{
				info->SubsetIndex[i] = SubsetIndex[i];
				info->num_Comp[i] = num_Comp[i];
			}

			info->num_SingleSolutions = SingleSolutions.size();
			info->SingleSolutions = new int[SingleSolutions.size()];
			for(int i = 0; i < SingleSolutions.size(); ++i)
			{
				int nr = SingleSolutions[i];
				info->SingleSolutions[i] = nr;
			}

			for(int i = 0; i < SubsetIndex.size(); ++i)
			{
				num_Comp[i] = 0;
				for(int j = 0; j < SingleSolutions.size(); ++j)
				{
					int nr = SingleSolutions[j];
					int k;
					for(int k = 0; k < _SingleSolutionInfoVec[nr]->num_SubsetIndex; ++k)
					{
						//_SingleSolutionInfoVec[nr]->group_Comp[k] = -100;
						if(_SingleSolutionInfoVec[nr]->SubsetIndex[k] == SubsetIndex[i])
						{
							_SingleSolutionInfoVec[nr]->group_Comp[k] = (num_Comp[i])++;
						}
					}
				}
			}
			for(int i = 0; i < SubsetIndex.size(); ++i)
			{
				if(info->num_Comp[i] != num_Comp[i])
				{
					std::cout << "ERROR in group_solutions: Something wrong with num_Comp. Aborting" << std::endl;
					return false;
				}
			}

			_GroupSolutionInfoVec.push_back(info);

			for(int i = 0; i < n; ++i)
			{
				int nr = GroupSolutions[i];
				delete _GroupSolutionInfoVec[nr]->SubsetIndex;
				delete _GroupSolutionInfoVec[nr]->num_Comp;
				delete _GroupSolutionInfoVec[nr]->SingleSolutions;
				delete _GroupSolutionInfoVec[nr];
			}

			for(int i = n-1; i >= 0; --i)
			{
				int nr = GroupSolutions[i];
				_GroupSolutionInfoVec.erase(_GroupSolutionInfoVec.begin() + nr);
			}

			return true;
		}

		bool clear_grouping()
		{
			if(_lock == true)
			{
				std::cout << "Pattern already locked. Can not clear grouping" << std::endl;
				return false;
			}

			// Clear old grouping
			for(int i = 0; i< _GroupSolutionInfoVec.size();++i)
			{
				delete _GroupSolutionInfoVec[i]->SubsetIndex;
				delete _GroupSolutionInfoVec[i]->num_Comp;
				delete _GroupSolutionInfoVec[i]->SingleSolutions;
				delete _GroupSolutionInfoVec[i];
			}
			_GroupSolutionInfoVec.clear();


			// restore single solution grouping
			for(int i = 0; i < _SingleSolutionInfoVec.size(); ++i)
			{
				GroupSolutionInfo* groupInfo = new GroupSolutionInfo();
				groupInfo->name = _SingleSolutionInfoVec[i]->name;
				groupInfo->num_SingleSolutions = 1;
				groupInfo->SingleSolutions = new int[1];
				for(int j = 0; j < 1; ++j)
				{
					groupInfo->SingleSolutions[j] = j;
				}
				groupInfo->num_SubsetIndex = _SingleSolutionInfoVec[i]->num_SubsetIndex;
				groupInfo->SubsetIndex = new int[groupInfo->num_SubsetIndex];
				groupInfo->num_Comp = new int[groupInfo->num_SubsetIndex];
				for(int j = 0; j < groupInfo->num_SubsetIndex; ++j)
				{
					groupInfo->SubsetIndex[j] = _SingleSolutionInfoVec[i]->SubsetIndex[j];
					groupInfo->num_Comp[j] = 1;
				}
				_GroupSolutionInfoVec.push_back(groupInfo);
			}

			return true;
		}

		bool print_info()
		{
			std::cout << "Pattern contains the following single Solutions:" << std::endl;
			for(int i = 0; i < _SingleSolutionInfoVec.size(); ++i)
			{
				std::cout << i << ": '"<< _SingleSolutionInfoVec[i]->name << "' on Subset(s) ";
				for(int j = 0; j < _SingleSolutionInfoVec[i]->num_SubsetIndex; ++j)
				{
					std::cout << _SingleSolutionInfoVec[i]->SubsetIndex[j] << " ";
				}
				std::cout << std::endl;
			}
			std::cout << "Grouped in the following components:" << std::endl;
			for(int i = 0; i < _GroupSolutionInfoVec.size(); ++i)
			{
				std::cout << i << ": '"<< _GroupSolutionInfoVec[i]->name << "' grouping single Solution(s) ";
				for(int j = 0; j < _GroupSolutionInfoVec[i]->num_SingleSolutions; ++j)
				{
					std::cout << _GroupSolutionInfoVec[i]->SingleSolutions[j] << " ";
				}
				//std::cout << std::endl;
				std::cout << "with DoF (Subset id, Nr. of Comps): ";
				for(int j = 0; j < _GroupSolutionInfoVec[i]->num_SubsetIndex; ++j)
				{
					std::cout << "(" << _GroupSolutionInfoVec[i]->SubsetIndex[j] << ",";
					std::cout << _GroupSolutionInfoVec[i]->num_Comp[j] << ") ";
				}
				std::cout << std::endl;
			}

			return true;
		}

		bool finalize()
		{
			// Check if only one Grouping present
			if(_GroupSolutionInfoVec.size() > 1)
			{
				std::cout << "ERROR in finalize: Currently only one Group Solution (with node values) allowed. Aborting." << std::endl;
				return false;
			}

			_lock = true;
			Grid* grid = _sh->get_assigned_grid();
			grid->register_observer(this, OT_GRID_OBSERVER | OT_VERTEX_OBSERVER);

			grid->attach_to_vertices(_aDoF);
			_aaDoFVRT.access(*grid, _aDoF);

			// assign DoFs
			int n = 0;

			VertexIterator iterBegin, iterEnd, iter;

			for(int i=0; i<_GroupSolutionInfoVec.size(); ++i)
			{
				for(int j = 0; j < _sh->num_subsets(); ++j)
				{
					int ncomp = _GroupSolutionInfoVec[i]->num_Comp[j];

					std::cout << "Assigning "<< ncomp << " Dof(s) per Vertex in Subset " << j << std::endl;

					iterBegin = _sh->begin<Vertex>(j);
					iterEnd = _sh->end<Vertex>(j);

					int count = 0;
					for(iter = iterBegin; iter != iterEnd; iter++)
					{
						Vertex* vrt = *iter;
						_aaDoFVRT[vrt].nr = n;
						_aaDoFVRT[vrt].ncomp = ncomp;
						n += ncomp;
						count += ncomp;
					}
					std::cout << count << " DoFs assigned to subset" << j << std::endl;
				}
			}

			std::cout << n << " DoFs assigned to grid" << std::endl;
			_num_doubles = n;
			return true;
		}

		int num_doubles()
		{
			return _num_doubles;
		}

		TrialSpaceType get_trial_space_type(int nr_solution)
		{
			return _SingleSolutionInfoVec[nr_solution]->Ansatzspace;
		}

		template<typename TElem>
		bool get_indices(TElem* elem, int nr_solution, int* ind)
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

			for(int i=0; i< geometry_traits<TElem>::Descriptor.num_vertices(); i++)
			{
				VertexBase* vert = elem->vertex(i);
				ind[i] = _aaDoFVRT[vert].nr + ncomp;
				assert(ind[i] < _num_doubles);
				std::cout << ind[i];
			}
		}

		int get_index(VertexBase* vrt, int nr_solution)
		{
			assert(nr_solution < _SingleSolutionInfoVec.size());
			int subsetIndex = _sh->get_subset_index(vrt);
			int i;
			for(i = 0; i < _SingleSolutionInfoVec[nr_solution]->num_SubsetIndex; ++i)
			{
				if(subsetIndex == _SingleSolutionInfoVec[nr_solution]->SubsetIndex[i]) break;
			}
			if(i ==  _SingleSolutionInfoVec[nr_solution]->num_SubsetIndex) assert(0 && "ERROR in get_index. Solution not defined for SubsetIndex.");

			int ncomp = _SingleSolutionInfoVec[nr_solution]->group_Comp[i];

			if(_aaDoFVRT[vrt].nr + ncomp >= _num_doubles || _aaDoFVRT[vrt].nr + ncomp < 0)
			{
				std::cout << "Nr.: " << _aaDoFVRT[vrt].nr << ", Ncomp: " << ncomp << " for Solution " << nr_solution << " on Subset: "<< subsetIndex<< std::endl;
				assert(0);
			}

			//std::cout << _aaDoFVRT[vrt].nr + ncomp << std::endl;

			return _aaDoFVRT[vrt].nr + ncomp;
		}

		~NumericalSolutionPattern()
		{
			for(int i = 0; i < _SingleSolutionInfoVec.size(); ++i)
			{
				delete _SingleSolutionInfoVec[i]->SubsetIndex;
				delete _SingleSolutionInfoVec[i]->group_Comp;
				delete _SingleSolutionInfoVec[i];
			}

			for(int i = 0; i < _GroupSolutionInfoVec.size(); ++i)
			{
				delete _GroupSolutionInfoVec[i]->SubsetIndex;
				delete _GroupSolutionInfoVec[i]->num_Comp;
				delete _GroupSolutionInfoVec[i]->SingleSolutions;
				delete _GroupSolutionInfoVec[i];
			}

			Grid* grid = _sh->get_assigned_grid();
			grid->detach_from_vertices(_aDoF);
			_aaDoFVRT.invalidate();
			grid->unregister_observer(this);

		}

		std::string get_name(int comp)
		{
			assert(comp < _SingleSolutionInfoVec.size());
			return _SingleSolutionInfoVec[comp]->name;
		}
		std::string get_name()
		{
			return _name;
		}
		int num_comp()
		{
			return _SingleSolutionInfoVec.size();
		}
		bool comp_def_in_subset(int nrComp, int subsetIndex)
		{
			for(int i = 0; i < _SingleSolutionInfoVec[nrComp]->num_SubsetIndex; ++i)
			{
				if(subsetIndex == _SingleSolutionInfoVec[nrComp]->SubsetIndex[i])
					return true;
			}
			return false;
		}

		SubsetHandler* get_assigned_subset()
		{
			return _sh;
		}

	protected:
		struct DoF {
			int nr;
			int ncomp;
		};

		typedef ug::Attachment<DoF> ADoF;

	protected:
		ADoF _aDoF;
		Grid::VertexAttachmentAccessor<ADoF> _aaDoFVRT;

		bool _lock;
		int* _neededDoFperSubset;

		std::vector<SingleSolutionInfo*> _SingleSolutionInfoVec;
		std::vector<GroupSolutionInfo*>  _GroupSolutionInfoVec;

		uint _num_sh;
		SubsetHandler* _sh;

		int _num_doubles;

		std::string _name;
};

class NumericalSolution {
	protected:
		typedef bool (*ValueFunction)(MathVector<3>, number&);

	public:
		NumericalSolution(std::string shortname, std::string longname, std::string description, NumericalSolutionPattern& pattern, Grid& grid);

		bool set_values(ValueFunction fct, int nr_func, SubsetHandler& sh, uint subsetIndex);

		Vector* GridVector();

		NumericalSolutionPattern* get_pattern();

		template <class TElem>
		ug::TrialSpace<TElem>& TrialSpace(int nr_func, int SubsetIndex)
		{
			return TrialSpaces<TElem>::TrialSpace(_pattern->get_trial_space_type(nr_func));
		}

		template <class TElem>
		bool get_local_DoFValues(TElem* elem, int nr_func, number* DoFValues)
		{
			typename geometry_traits<TElem>::Descriptor TDesc;

			int nvalues = TDesc.num_vertices();

			int* indices = new int[nvalues];

			for(int i=0; i< nvalues; i++)
			{
				VertexBase* vert = elem->vertex(i);
				indices[i] = (int) (_pattern->get_index(vert, nr_func));
			}

			double *doubleDoFValues = new double[nvalues];
			if(_GridVector->get_values(nvalues, indices, doubleDoFValues) == false) return false;

			for(int i=0; i< nvalues; i++)
			{
				DoFValues[i] = (number) doubleDoFValues[i];
			}


			return true;
		}

		bool print();

		~NumericalSolution();

		std::string get_name(int comp);

	protected:
		std::string m_shortname;
		std::string m_longname;
		std::string m_description;

		TrialSpaceType m_TrialSpaceType;

		Vector* _GridVector;
		NumericalSolutionPattern* _pattern;
		Grid* m_grid;

};




} /* end namespace ug */


#endif /* __H__LIBDISCRETIZATION__NUMERICALSOLUTION__ */
