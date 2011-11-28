// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 24.11.2011 (m,d,y)

#ifndef __H__UG__surface_view_TMP__
#define __H__UG__surface_view_TMP__

#include <boost/iterator/iterator_facade.hpp>
#include "lib_grid/multi_grid.h"
#include "../subset_handler_multi_grid.h"

namespace ug{
namespace tmp{

class SurfaceView;

template <class TElem>
class SurfaceViewElementIterator
	: public boost::iterator_facade<SurfaceViewElementIterator<TElem>, TElem,
									boost::forward_traversal_tag>
{
	private:
		friend class boost::iterator_core_access;

		inline bool equal(SurfaceViewElementIterator<TElem> const& other) const;
		void increment();
		inline TElem& dereference() const;

	private:
		typename geometry_traits<TElem>::iterator m_elemIter;
		SurfaceView* m_surfView;
		int m_si;
		short int m_lvl;
		short int m_topLvl;// should be min(topLvl, maxFilledLvl)
};

class SurfaceView
{
	public:
		template <class TElem>
		struct traits{
			typedef SurfaceViewElementIterator<TElem>		iterator;
			typedef SurfaceViewElementIterator<TElem const>	const_iterator;
		};

	public:
		SurfaceView(MGSubsetHandler& sh) : m_sh(sh) {}

		MGSubsetHandler& subset_handler()	{return m_sh;}


	///ATTENTION!!!
		bool is_shadow() const	{return false;}

	protected:
		MGSubsetHandler& m_sh;
};


#ifdef OUT_OF_ORDER
class SurfaceView : public ISubsetHandler
{
	public:
		using ISubsetHandler::assign_subset;

	public:
		SurfaceView();
		SurfaceView(MultiGrid& mg);

		void assign_grid(MultiGrid& mg);
		inline MultiGrid* grid()	{return m_pMG;}

	///	perform cleanup
		virtual void grid_to_be_destroyed(Grid* grid);

	protected:
	////////////////////////////////////////////////
	//	protected virtual methdos of ISubsetHandler.
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
		void cleanup();

	protected:
		typedef ISubsetHandler::SectionContainer 	SectionContainer;
		typedef ISubsetHandler::AttachedElemList	AttachedElemList;

		struct Subset
		{
		/// holds pointers to elements in the surface grid.
			SectionContainer 	m_surf_elements[NUM_GEOMETRIC_BASE_OBJECTS];
		/// holds pointers to shadow elements.
			SectionContainer 	m_shadows[NUM_GEOMETRIC_BASE_OBJECTS];
		};

		typedef std::vector<Subset*>	SubsetVec;
		typedef std::vector<SubsetVec>	LevelVec;

	protected:
		MultiGrid*		m_pMG;
		AttachedElemList::AEntry	m_aSharedEntry;

};
#endif
}//	end of namespace
}//	end of namespace


////////////////////////////////////////
//	include implementation
#include "surface_view_impl.hpp"

#endif
