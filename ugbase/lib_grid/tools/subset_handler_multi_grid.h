//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y08 m11 d24   (reworked y09 m12 d15)

#ifndef __H__LIBGRID__SUBSET_HANDLER_MULTI_GRID__
#define __H__LIBGRID__SUBSET_HANDLER_MULTI_GRID__

#include <vector>
#include "lib_grid/multi_grid.h"
#include "common/util/section_container.h"
#include "subset_handler_interface.h"

namespace ug
{
////////////////////////////////////////////////////////////////////////
//	MultiGridSubsetHandlerBase
class MultiGridSubsetHandler : public ISubsetHandler
{
	public:
		using ISubsetHandler::assign_subset;
		
	public:
		MultiGridSubsetHandler(uint supportedElements = SHE_ALL);
		MultiGridSubsetHandler(MultiGrid& mg, uint supportedElements = SHE_ALL);
	/**	WARNING: Don't call the copy-constructor from derived classes,
	  *	Since it calls virtual methods.*/
		MultiGridSubsetHandler(const MultiGridSubsetHandler& sh);
		~MultiGridSubsetHandler();
		
		inline void assign_grid(MultiGrid& mg)	{m_pMG = &mg; ISubsetHandler::assign_grid(mg);}
		
	///	Makes sure that the subset with the given index exists.
	/**	If required the subsets between num_subsets() and index will be created.*/
		inline void subset_required(int index);
		
	///	returns the number of subsets
		inline uint num_subsets() const	{return m_numSubsets;}
		
	///	returns the number of levels
		inline uint num_levels() const	{return (uint)m_levels.size();}
		
	////////////////////////////////////////////////
	//	implementation of public virtual methdos of ISubsetHandler.
	///	assigns a vertex to a subset.
	/**	If the subset doesn't already exist, it will be created.*/
		void assign_subset(VertexBase* elem, int subsetIndex);

	///	assigns an edge to a subset.
	/**	If the subset doesn't already exist, it will be created.*/
		void assign_subset(EdgeBase* elem, int subsetIndex);

	///	assigns a face to a subset.
	/**	If the subset doesn't already exist, it will be created.*/
		void assign_subset(Face* elem, int subsetIndex);

	///	assigns a volume to a subset.
	/**	If the subset doesn't already exist, it will be created.*/
		void assign_subset(Volume* elem, int subsetIndex);

	////////////////////////////////////////////////
	//	element-access
	///	returns the begin-iterator for the elements of type TElem in the given subset.
	/**	e.g. begin<Triangle>(0)*/
		template <class TElem>
		typename geometry_traits<TElem>::iterator
		begin(int subsetIndex, int level);

	///	returns the end-iterator for the elements of type TElem in the given subset.
	/**	e.g. end<Triangle>(0)*/
		template <class TElem>
		typename geometry_traits<TElem>::iterator
		end(int subsetIndex, int level);

	///	returns the number of elements in the given subset
		template <class TElem>
		uint num(int subsetIndex);

	///	returns the number of elements in the given subset on the given level
		template <class TElem>
		uint num(int subsetIndex, int level);

	///	removes all elements of type TElem from the specified subset.
		template <class TElem>
		void clear_subset_elements(int subsetIndex);

	///	removes all elements of type TElem from the specified subset on the given level.
		template <class TElem>
		void clear_subset_elements(int subsetIndex, int level);

	///	returns a GeometricObjectCollection
	/**	the returned GeometricObjectCollection hold the elements of the
	 *	specified subset on the given level.*/
		GeometricObjectCollection
		get_goc(int subsetIndex, int level);
		
	///	returns a MultiLevelGeometricObjectCollection
	/**	the returned MultiLevelGeometricObjectCollection hold the
	 *	elements of the specified subset.*/
		MultiLevelGeometricObjectCollection
		get_mlgoc_by_subset(int subsetIndex);

	///	returns a MultiLevelGeometricObjectCollection
	/**	the returned MultiLevelGeometricObjectCollection hold the
	 *	elements of the specified level.*/
		MultiLevelGeometricObjectCollection
		get_mlgoc_by_level(int level);

	//	derived from GridObserver
		virtual void unregistered_from_grid(Grid* grid);
		
	protected:
	////////////////////////////////////////////////
	//	implementation of protected virtual methdos of ISubsetHandler.
	///	erases the subsets. Doesn't alter any indices.
		void erase_subset_lists();
		
	///	clears the element lists in the given subset. Does not alter any indices.
		void clear_subset_lists(int index);

	///	changes the subset-indices of all elements int the subset.
	/**	WARNING: subsets are not automatically changed accordingly.
	 *	After termination Subset-Indices and Subset-Infos/iterators are asynchronous.
	 *	Make sure to change subset-infos and iterators accordingly.*/		
		void change_subset_indices(int indOld, int indNew);		

		
	///	add a subset
		void add_required_subset_lists(int maxIndex);
		
	///	erases the subset but does not touch the subset-indices.
		void erase_subset_lists(int index);

	///	swaps the subsets but does not touch the subset-indices.
		void swap_subset_lists(int ind1, int ind2);

	///	moves the subset but does not touch the subset-indices.
		void move_subset_lists(int indexFrom, int indexTo);

	///	this method is called by ISubsetHandler when attachment_support has been enabled.
		void register_subset_elements_at_pipe();

	////////////////////////////////////////////////
	//	protected helper methods
	///	a helper method for the public assign_subset methods.
		template<class TElemPtr>
		void assign_subset(TElemPtr elem, int subsetIndex, int elemType);

	///	helper for change_subset_indices
		template<class TElem>
		void change_elem_subset_indices(int indOld, int indNew);
		
		inline void level_required(int level)		{while(m_levels.size() <= level) add_level();}

		void add_level();
		void add_subset_to_all_levels();///< increases m_numSubsets.

	protected:
		typedef ISubsetHandler::SectionContainer SectionContainer;
		
		struct Subset
		{
			SectionContainer 	m_elements[NUM_GEOMETRIC_BASE_OBJECTS];	/// holds pointers to elements.
			//attachment_pipe
		};

		typedef std::vector<Subset*>	SubsetVec;
		typedef std::vector<SubsetVec>	LevelVec;

	///	returns the subset with index si on the given level
		inline Subset* subset(int si, int level)	{return m_levels[level][si];}

	protected:
		MultiGrid*		m_pMG;
		LevelVec		m_levels;
		int				m_numSubsets;
};

inline void
MultiGridSubsetHandler::
subset_required(int index)
{
	if(index >= m_numSubsets)
	{
		add_required_subset_lists(index);
		ISubsetHandler::subset_info_required(index);
	}
}

typedef MultiGridSubsetHandler MGSubsetHandler;

}//	end of namespace

//	include implementation
#include "subset_handler_multi_grid_impl.hpp"

#endif
