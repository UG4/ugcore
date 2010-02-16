//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y08 m11 d24   (reworked y09 m12 d15)

#ifndef __H__LIBGRID__SUBSET_HANDLER_GRID__
#define __H__LIBGRID__SUBSET_HANDLER_GRID__

#include <vector>
#include "lib_grid/grid/grid.h"
#include "common/util/section_container.h"
#include "subset_handler_interface.h"

namespace ug
{
////////////////////////////////////////////////////////////////////////
//	GridSubsetHandler
class GridSubsetHandler : public ISubsetHandler
{
	public:
		using ISubsetHandler::assign_subset;
		
	public:
		GridSubsetHandler(uint supportedElements = SHE_ALL);
		GridSubsetHandler(Grid& grid, uint supportedElements = SHE_ALL);
		GridSubsetHandler(const GridSubsetHandler& sh);
		~GridSubsetHandler();
		
		GridSubsetHandler& operator = (const GridSubsetHandler& sh);

		inline void assign_grid(Grid& grid)		{ISubsetHandler::assign_grid(grid);}
		
	///	Makes sure that the subset with the given index exists.
	/**	If required the subsets between num_subsets() and index will be created.
	 *	ISubsetHandler::subset_info_required is called automatically.*/
		inline void subset_required(int index);
		
	///	returns the number of subsets
		inline uint num_subsets() const	{return m_subsets.size();}
		
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
		begin(int subsetIndex) const;

	///	returns the end-iterator for the elements of type TElem in the given subset.
	/**	e.g. end<Triangle>(0)*/
		template <class TElem>
		typename geometry_traits<TElem>::iterator
		end(int subsetIndex) const;

	///	returns the number of elements in the given subset
		template <class TElem>
		uint num_elements(int subsetIndex) const;

	///	returns the number of elements in the given subset
		template <class TElem>
		uint num(int subsetIndex) const;

	///	removes all elements of type TElem from the specified subset.
		template <class TElem>
		void clear_subset_elements(int subsetIndex);

	//	geometric-object-collection
		GeometricObjectCollection
		get_geometric_object_collection(int subsetIndex);
		
	//	multi-level-geometric-object-collection
		MultiLevelGeometricObjectCollection
		get_multi_level_geometric_object_collection();

	//	derived from GridObserver
		virtual void unregistered_from_grid(Grid* grid);
		
	///	only for debug purposes
		template <class TElem>
		bool perform_self_tests();
		
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

	protected:
		typedef ISubsetHandler::SectionContainer SectionContainer;
		
		struct Subset
		{
			SectionContainer 	m_elements[NUM_GEOMETRIC_BASE_OBJECTS];	/// holds pointers to elements.
			//attachment_pipe
		};

		typedef std::vector<Subset*>	SubsetVec;

	protected:
		SubsetVec		m_subsets;
};


inline void GridSubsetHandler::
subset_required(int index)
{
	if(index >= m_subsets.size())
	{
		add_required_subset_lists(index);
		ISubsetHandler::subset_info_required(index);
	}
}


typedef GridSubsetHandler SubsetHandler;

}//	end of namespace

//	include implementation
#include "subset_handler_grid_impl.hpp"

#endif
