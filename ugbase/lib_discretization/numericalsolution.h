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
		NumericalSolutionPattern(std::string name, SubsetHandler& sh);

		bool add_solution(std::string name, TrialSpaceType Ansatzspace, int* SubsetIndex, int n);

		bool group_solution(std::string name, int* GroupSolutions, int n);

		bool clear_grouping();

		bool print_info();

		bool finalize();

		int num_doubles(int i);
		int num_doubles();

		int num_levels();

		TrialSpaceType get_trial_space_type(int nr_solution);

		template<typename TElem>
		bool get_indices(TElem* elem, int nr_solution, int* ind);

		int get_index(VertexBase* vrt, int nr_solution);

		~NumericalSolutionPattern();

		std::string get_name(int comp);

		std::string get_name();

		int num_comp();

		bool comp_def_in_subset(int nrComp, int subsetIndex);

		SubsetHandler* get_assigned_subset();

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

		int _num_levels;
		int* _num_doubles;

		std::string _name;
};

class NumericalSolution {
	protected:
		typedef bool (*ValueFunction)(MathVector<3>, number&);

	public:
		NumericalSolution(std::string shortname, std::string longname, std::string description, NumericalSolutionPattern& pattern, Grid& grid);

		bool set_values(ValueFunction fct, int nr_func, SubsetHandler& sh, uint subsetIndex);
		bool set_values(ValueFunction fct, int nr_func, SubsetHandler& sh, uint subsetIndex, int level);

		Vector* GridVector(int level);
		Vector* GridVector();

		NumericalSolutionPattern* get_pattern();

		template <class TElem>
		ug::TrialSpace<TElem>& TrialSpace(int nr_func, int SubsetIndex);

		template <class TElem>
		bool get_local_DoFValues(TElem* elem, int nr_func, number* DoFValues);

		bool print();

		~NumericalSolution();

		std::string get_name(int comp);

	protected:
		std::string m_shortname;
		std::string m_longname;
		std::string m_description;

		TrialSpaceType m_TrialSpaceType;

		std::vector<Vector*> _GridVector;
		NumericalSolutionPattern* _pattern;
		Grid* m_grid;

};


} /* end namespace ug */

#include "numericalsolution_impl.h"

#endif /* __H__LIBDISCRETIZATION__NUMERICALSOLUTION__ */
