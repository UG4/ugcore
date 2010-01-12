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
#include "domain.h"

namespace ug {

class DoFPattern : public GridObserver{
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
		DoFPattern(std::string name, ISubsetHandler& sh);

		bool add_solution(std::string name, TrialSpaceType Ansatzspace, int* SubsetIndex, int n);

		bool group_solution(std::string name, int* GroupSolutions, int n);

		bool clear_grouping();

		bool print_info();

		bool finalize();

		int num_doubles(int i);
		int num_doubles();

		uint num_levels();

		TrialSpaceType get_TrialSpaceType(int nr_solution);

		template<typename TElem>
		bool get_indices(TElem* elem, unsigned int nr_solution, int* ind);

		int get_index(VertexBase* vrt, int nr_solution);

		~DoFPattern();

		std::string get_name(int comp);

		std::string get_name();

		int num_comp();

		bool comp_def_in_subset(int nrComp, int subsetIndex);

		ISubsetHandler* get_assigned_subset();

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
		ISubsetHandler* _sh;

		uint _num_levels;
		int* _num_doubles;

		std::string _name;
};

template <int d>
class NumericalSolution {
	public:
		NumericalSolution(std::string shortname, std::string longname, std::string description, DoFPattern& pattern, Domain<d>& domain);

		bool set_values(bool (*fct)(MathVector<d>, number&), int nr_func, ISubsetHandler& sh, uint subsetIndex);
		bool set_values(bool (*fct)(MathVector<d>, number&), int nr_func, ISubsetHandler& sh, uint subsetIndex, int level);

		Vector* GridVector(int level);
		Vector* GridVector();

		DoFPattern* get_pattern();

		template <class TElem>
		const TrialSpace<TElem>& get_TrialSpace(int nr_func);

		TrialSpaceType get_TrialSpaceType(int nr_func);

		Domain<d>* get_domain()
		{
			return _domain;
		}

		template<typename TElem>
		bool get_indices(TElem* elem, int nr_solution, int* ind);

		template <class TElem>
		bool get_local_DoFValues(TElem* elem, int nr_func, number* DoFValues);

		bool print();

		bool assign(const NumericalSolution<d>& v);

		~NumericalSolution();

		std::string get_name(int comp);

	protected:
		std::string _shortname;
		std::string _longname;
		std::string _description;

		std::vector<Vector*> _GridVector;
		DoFPattern* _pattern;
		Domain<d>* _domain;
};


} /* end namespace ug */

#include "numericalsolution_impl.h"

#endif /* __H__LIBDISCRETIZATION__NUMERICALSOLUTION__ */
