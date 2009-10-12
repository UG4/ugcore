//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y08 m11 d24

#ifndef __H__LIBGRID__SUBSET_HANDLER__
#define __H__LIBGRID__SUBSET_HANDLER__

#include <list>
#include <string>
#include <vector>
#include "grid/grid.h"
#include "common_attachments.h"
#include "common/util/section_container.h"

namespace ug
{

////////////////////////////////////////////////////////////////////////
//	SubsetInfo
///	a struct that holds information associated with subsets.
struct SubsetInfo
{
	SubsetInfo();
	std::string	name;
	int			materialIndex;
};

////////////////////////////////////////////////////////////////////////
//	SubsetHandler
///	groups the grids elements to subsets.
/**
 * Using the subset handler you can organize your grid in subsets.
 * Each element can be assigned to exactly one subset.
 * New elements are assigned to no subset by default (index -1).
 */
class SubsetHandler : public GridObserver
{
	public:
		SubsetHandler();
		SubsetHandler(Grid& grid);
		SubsetHandler(const SubsetHandler& sh);

		virtual ~SubsetHandler();

		void assign_grid(Grid& grid);
		Grid* get_assigned_grid();

	/**	new elements will be automatically assigned to this subset.
	 * 	set this to a negative value to avoid automatic assignment (-1 by default).
	 *	only used if subset_inheritance is disabled or if no parent is specified.*/
		void set_default_subset_index(int subsetIndex);
		inline int get_default_subset_index()	{return m_defaultSubsetIndex;}

	/**	if enabled, newly created elements derive their subset-index from their parents.
	 * Enabled by default.
	 * If enabled, the default subset index will be ignored if a parent is specified
	 * on element creation.*/
		void enable_subset_inheritance(bool bEnable);
		inline bool subset_inheritance_enabled()	{return m_bSubsetInheritanceEnabled;}

		void clear();
		void clear_subset(int subsetIndex);
		void clear_subsets();

		template <class TElem>
		void clear_subsets_elements(int subsetIndex);

	///	inserts a subset at the given index. Moves all other subsets 1 index higher.
		void insert_subset(int subsetIndex);///< changes subset-indices of other subsets.
	///	erases the subset at the given index. Assigns -1 to all entries. Moves all other subsets 1 index up.
		void erase_subset(int subsetIndex);///< changes subset-indices of other subsets.
	///	Swaps the given subsets,
		void swap_subsets(int subsetIndex1, int subsetIndex2);
	///	Moves the subset from index From to index To. Moves all subsets between indexFrom+1 and indexTo in the opposite direction.
		void move_subset(int indexFrom, int indexTo);///< changes subset indices of other subsets.

		void set_subset_info(int subsetIndex, const SubsetInfo& subsetInfo);
		SubsetInfo& subset_info(int subsetIndex);
		const SubsetInfo& subset_info(int subsetIndex) const;

		void assign_subset(VertexBase* elem, int subsetIndex);
		void assign_subset(EdgeBase* elem, int subsetIndex);
		void assign_subset(Face* elem, int subsetIndex);
		void assign_subset(Volume* elem, int subsetIndex);

		template <class TIterator>
		void assign_subset(TIterator iterBegin, TIterator iterEnd, int subsetIndex);

		int get_subset_index(GeometricObject* elem);
		int get_subset_index(VertexBase* elem);
		int get_subset_index(EdgeBase* elem);
		int get_subset_index(Face* elem);
		int get_subset_index(Volume* elem);

		inline uint num_subsets() const	{return m_subsets.size();}

		template <class TElem>
		uint num_elements(int subsetIndex);

		template <class TElem>
		uint num(int subsetIndex);

		template <class TElem>
		typename geometry_traits<TElem>::iterator
		begin(int subsetIndex);

		template <class TElem>
		typename geometry_traits<TElem>::iterator
		end(int subsetIndex);

	//	geometric-object-collection
		GeometricObjectCollection get_geometric_object_collection(int subsetIndex);

	//	grid callbacks
		virtual void registered_at_grid(Grid* grid);
		virtual void unregistered_from_grid(Grid* grid);
		virtual void elements_to_be_cleared(Grid* grid);

	//	vertex callbacks
		virtual void vertex_created(Grid* grid, VertexBase* vrt, GeometricObject* pParent = NULL);
		virtual void vertex_to_be_erased(Grid* grid, VertexBase* vrt);

	//	edge callbacks
		virtual void edge_created(Grid* grid, EdgeBase* edge, GeometricObject* pParent = NULL);
		virtual void edge_to_be_erased(Grid* grid, EdgeBase* edge);

	//	face callbacks
		virtual void face_created(Grid* grid, Face* face, GeometricObject* pParent = NULL);
		virtual void face_to_be_erased(Grid* grid, Face* face);

	//	volume callbacks
		virtual void volume_created(Grid* grid, Volume* vol, GeometricObject* pParent = NULL);
		virtual void volume_to_be_erased(Grid* grid, Volume* vol);

	protected:
		template <class TElem>
		void reset_subset_indices(typename geometry_traits<TElem>::iterator iterBegin,
									typename geometry_traits<TElem>::iterator iterEnd);

		template <class TElem>
		void set_subset_indices(typename geometry_traits<TElem>::iterator iterBegin,
								typename geometry_traits<TElem>::iterator iterEnd,
								int subsetIndex);///<	does not alter subset iterators. Use only if you know what you're doing!

		void resize_subset_vec(int newSize);

	protected:
		typedef SectionContainer<GeometricObject*, std::list<GeometricObject*> >	SectionContainer;
		struct Subset
		{
			SectionContainer 	m_elements[NUM_GEOMETRIC_BASE_OBJECTS];	/// holds pointers to elements.
			//attachment_pipe
		};

	protected:
		typedef AInt	ASubsetIndex;
		typedef Attachment<std::list<GeometricObject*>::iterator>	AIterator;
		typedef std::vector<Subset*>	SubsetVec;
		typedef std::vector<SubsetInfo>	SubsetInfoVec;

	protected:
		Grid*			m_pGrid;
		ASubsetIndex	m_aSubsetIndex;
		AIterator		m_aIterator;
		SubsetVec		m_subsets;
		SubsetInfoVec	m_subsetInfos;
		int				m_defaultSubsetIndex;
		bool			m_bSubsetInheritanceEnabled;

		Grid::VertexAttachmentAccessor<ASubsetIndex>	m_aaSubsetIndexVRT;
		Grid::EdgeAttachmentAccessor<ASubsetIndex>		m_aaSubsetIndexEDGE;
		Grid::FaceAttachmentAccessor<ASubsetIndex>		m_aaSubsetIndexFACE;
		Grid::VolumeAttachmentAccessor<ASubsetIndex>	m_aaSubsetIndexVOL;

		Grid::VertexAttachmentAccessor<AIterator>		m_aaIteratorVRT;
		Grid::EdgeAttachmentAccessor<AIterator>			m_aaIteratorEDGE;
		Grid::FaceAttachmentAccessor<AIterator>			m_aaIteratorFACE;
		Grid::VolumeAttachmentAccessor<AIterator>		m_aaIteratorVOL;
};

}//	end of namespace

//	include implementation
#include "subset_handler_impl.hpp"

#endif
