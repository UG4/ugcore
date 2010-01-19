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
////////////////////////////////////////////////
// DoFPattern
/// manages the distribution of degrees of freedom
/**
 * Manages the distribution of DoFs on a Domain, that has been separated into subsets.
 * Given a SubsetHandler, single solutions can be added to selected subsets. The chosen
 * TrialSpaceType determines the number of DoFs on the Geometric Objects in the Subset.
 * The user can group several solutions to a multivalued solution.
 * Invoking the finalize command let the DoFPattern distribute the DoFs.
 */
class DoFPattern : public GridObserver{
	protected:
	typedef unsigned short comp_type;
	typedef uint index_type;

	protected:
		struct SingleSolutionInfo {
			std::string name; 					// name of solution
			TrialSpaceType TrialSpace; 			// associated trial space
			std::vector<int> SubsetIndex; 		// subset indices where solution lives
			std::vector<int> group_comp; 		// nr of component in grouping per subset
		};

		struct GroupSolutionInfo {
			std::string name; 					// name of grouping
			std::vector<int> SingleSolutions; 	// index of single solutions contained in this grouping
			std::vector<int> SubsetIndex; 		// subset index where grouping lives
			std::vector<uint> num_comp;    		// number of Components per Subset
		};

	public:
		/// constructor
		/**
		 * \param[in] name	Name of this DoFPattern
		 * \param[in] sh	SubsetHandler
		 */
		DoFPattern(std::string name, ISubsetHandler& sh);

		/// add a single solution of TrialSpaceType to selected subsets
		/**
		 * \param[in] name			Name of this Single Solution
		 * \param[in] TrialSpace	Trial Space for this function
		 * \param[in] SubsetIndecex	Std::Vector of subset indeces, where this solution lives
		 */
		bool add_solution(std::string name, TrialSpaceType TrialSpace, std::vector<int>& SubsetIndices);

		/// group single solutions
		/**
		 * By this function a user can group single solutions to a new one. The single solutions will be
		 * removed from the pattern and a new group solution containing those will be added. It is
		 * also possible to group 'Group Solutions'.
		 * The Grouping has an effect on the dof numbering. While single solutions get an own index with only one component for
		 * each dof, a Grouped Solution has on index with several components.
		 *
		 */
		bool group_solution(std::string name, std::vector<uint>& GroupSolutions);

		/// clears grouping
		bool clear_grouping();

		/// gives informations about the current status
		bool print_info();

		/// performs a finalizing step. The Pattern can not be altered after finishing
		bool finalize();

		/// returns the number of dofs on level 'level'
		uint num_dofs(uint level = 0);

		/// returns the number of dofs on level 'level', subset s
		//uint num_dofs(int s, uint level = 0);

		/// returns the number of levels
		uint num_levels();

		/// returns the trial space of the discrete function nr_fct
		TrialSpaceType get_TrialSpaceType(uint nr_fct);

		/// returns the indices of the dofs on the Element elem for the discrete function 'nr_fct'
		template<typename TElem>
		bool get_indices(TElem* elem, uint nr_fct, uint* ind);

		/// returns the index of the dofs on the Element elem for the discrete function 'nr_fct'
		index_type get_index(VertexBase* vrt, uint nr_fct);

		/// returns the name of the discrete function nr_fct
		std::string get_name(uint nr_fct);

		/// returns the name of this dof pattern
		std::string get_name();

		/// returns the number of discrete functions in this dof pattern
		int num_comp();

		/// returns true if the discrete function nr_fct is defined on subset s
		bool comp_def_in_subset(uint nr_fct, int s);

		/// returns the assigned SubsetHandler
		ISubsetHandler* get_assigned_subset();

		/// destructor
		~DoFPattern();

	protected:

		/// general dof structure
		/**
		 * index - processor global consecutive index for dofs
		 * ncomp - number of subcomponents of dof
		 */
		struct DoF {
			uint index; // index
			unsigned short ncomp;  // number of components
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
		std::vector<uint> _num_dofs;

		std::string _name;
};

template <int d>
class NumericalSolution {
	public:
		/// constructor
		NumericalSolution(std::string shortname, std::string longname, std::string description, DoFPattern& pattern, Domain<d>& domain);

		/// sets the value of the function
		bool set_values(bool (*fct)(MathVector<d>, number&), uint nr_fct, ISubsetHandler& sh, uint subsetIndex, uint level = 0);

		/// returns the associated GridVector
		Vector* GridVector(uint level = 0);

		/// returns the associated pattern
		DoFPattern* get_pattern();

		/// returns the trial space type of numerical solution 'nr_fct'
		TrialSpaceType get_TrialSpaceType(uint nr_func);

		/// returns the trial space of numerical solution 'nr_fct'
		template <class TElem>
		const TrialSpace<TElem>& get_TrialSpace(uint nr_fct);

		/// returns the associated domain
		Domain<d>* get_domain();

		/// returns the dof indices of a given element for the solution 'nr_fct'
		template<typename TElem>
		bool get_indices(TElem* elem, uint nr_fct, int* ind);

		/// returns the dof values of a given element for the solution 'nr_fct'
		template <class TElem>
		bool get_local_DoFValues(TElem* elem, uint nr_func, number* DoFValues);

		/// print information
		bool print();

		/// assignment
		bool assign(NumericalSolution<d>& v);

		/// destructor
		~NumericalSolution();

		/// name of numerical solution 'nr_fct'
		std::string get_name(uint nr_fct);

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
