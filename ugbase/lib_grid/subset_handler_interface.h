//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y08 m11 d24

#ifndef __H__LIBGRID__SUBSET_HANDLER_INTERFACE__
#define __H__LIBGRID__SUBSET_HANDLER_INTERFACE__

#include <list>
#include <string>
#include <vector>
#include "grid/grid.h"
#include "common_attachments.h"

namespace ug
{
////////////////////////////////////////////////////////////////////////
//	predeclarations
class ISubsetHandler;

////////////////////////////////////////////////////////////////////////
//	SubsetHandlerElements
///	Use these constants to specify which elements shall be supported by a SubsetHandler.
/**
 * You may combine the constants using or-operations.
 */
enum SubsetHandlerElements
{
	SHE_NONE = 0,
	SHE_VERTEX = 1,
	SHE_EDGE = 1<<1,
	SHE_FACE = 1<<2,
	SHE_VOLUME = 1 << 3,
	SHE_ALL = SHE_VERTEX | SHE_EDGE | SHE_FACE | SHE_VOLUME
};

////////////////////////////////////////////////////////////////////////
//	SubsetState
///	The SubsetState is not yet really used inside of libGrid.
/**
 * The main reason why a SubsetState is introduced, is that
 * applications that use libGrid need a mechanism to store
 * information in a subset.
 * It would be a good idea to think about an attachment-like system
 * for subsets.
 */
enum SubsetState
{
	SS_NONE = 0,
	SS_USER_STATE = 1 << 16
};

////////////////////////////////////////////////////////////////////////
//	SubsetInfo
///	a struct that holds information associated with subsets.
/**
 * In the moment a SubsetInfo is a collection of various types.
 * None of them are really required for libGrid (indeed only name and
 * materialIndex are used in the moment).
 * The other variables are introduced mainly for applications that use
 * libGrid. This is not the best way to do this!
 * It would be a good idea to think about an attachment-like system
 * for subsets.
 */
struct SubsetInfo
{
	SubsetInfo();
	std::string	name;
	int			materialIndex;
	vector4		color;
	uint		subsetState;///< an or-combination of SubsetState flags.
};


////////////////////////////////////////////////////////////////////////
//	specialization of attachment_traits for VertexBase
template<>
class attachment_traits<VertexBase*, ISubsetHandler>
{
	public:
		typedef VertexBase*&			ElemRef;
		typedef VertexBase*				ElemPtr;
		typedef const VertexBase*		ConstElemPtr;
		typedef ISubsetHandler*			ElemHandlerPtr;
		typedef const ISubsetHandler*	ConstElemHandlerPtr;

		static inline void invalidate_entry(ElemHandlerPtr pHandler, ElemRef elem)				{elem = NULL;}
		static inline bool entry_is_invalid(ElemHandlerPtr pHandler, ElemRef elem)				{return elem != NULL;}
		static inline uint get_data_index(ElemHandlerPtr pHandler, ConstElemPtr elem);
		static inline void set_data_index(ElemHandlerPtr pHandler, ElemPtr elem, uint index);
};

//	specialization of attachment_traits for EdgeBase
template<>
class attachment_traits<EdgeBase*, ISubsetHandler>
{
	public:
		typedef EdgeBase*&				ElemRef;
		typedef EdgeBase*				ElemPtr;
		typedef const EdgeBase*			ConstElemPtr;
		typedef ISubsetHandler*			ElemHandlerPtr;
		typedef const ISubsetHandler*	ConstElemHandlerPtr;

		static inline void invalidate_entry(ElemHandlerPtr pHandler, ElemRef elem)				{elem = NULL;}
		static inline bool entry_is_invalid(ElemHandlerPtr pHandler, ElemRef elem)				{return elem != NULL;}
		static inline uint get_data_index(ElemHandlerPtr pHandler, ConstElemPtr elem);
		static inline void set_data_index(ElemHandlerPtr pHandler, ElemPtr elem, uint index);
};

//	specialization of attachment_traits for Face
template<>
class attachment_traits<Face*, ISubsetHandler>
{
	public:
		typedef Face*&					ElemRef;
		typedef Face*					ElemPtr;
		typedef const Face*				ConstElemPtr;
		typedef ISubsetHandler*			ElemHandlerPtr;
		typedef const ISubsetHandler*	ConstElemHandlerPtr;

		static inline void invalidate_entry(ElemHandlerPtr pHandler, ElemRef elem)				{elem = NULL;}
		static inline bool entry_is_invalid(ElemHandlerPtr pHandler, ElemRef elem)				{return elem != NULL;}
		static inline uint get_data_index(ElemHandlerPtr pHandler, ConstElemPtr elem);
		static inline void set_data_index(ElemHandlerPtr pHandler, ElemPtr elem, uint index);
};

//	specialization of attachment_traits for Volume
template<>
class attachment_traits<Volume*, ISubsetHandler>
{
	public:
		typedef Volume*&				ElemRef;
		typedef Volume*					ElemPtr;
		typedef const Volume*			ConstElemPtr;
		typedef ISubsetHandler*			ElemHandlerPtr;
		typedef const ISubsetHandler*	ConstElemHandlerPtr;

		static inline void invalidate_entry(ElemHandlerPtr pHandler, ElemRef elem)				{elem = NULL;}
		static inline bool entry_is_invalid(ElemHandlerPtr pHandler, ElemRef elem)				{return elem != NULL;}
		static inline uint get_data_index(ElemHandlerPtr pHandler, ConstElemPtr elem);
		static inline void set_data_index(ElemHandlerPtr pHandler, ElemPtr elem, uint index);
};


////////////////////////////////////////////////////////////////////////
//	ISubsetHandler
/**
 * A derived class has to implement the following public methods:
 * <code>
 * virtual void assign_subset(VertexBase* elem, int subsetIndex)
 * virtual void assign_subset(EdgeBase* elem, int subsetIndex)
 * virtual void assign_subset(Face* elem, int subsetIndex)
 * virtual void assign_subset(Volume* elem, int subsetIndex)
 * </code>
 *
 * In those methods 
 * Derived classes have to store the objects that are selected for
 * a subset in a ISubsetHandler::SectionContainer. Note that multiple
 * SectionContainers per subset may be used.
 *
 * A SubsetHandler also supports subset-attachments.
 * That means that you can attach data to the elements of a subset.
 * Use attach_to_vertices, attach_to_edges, ... to attach data to
 * the elements of the subset.
 * use ISubsetHandler::VertexAttachmentAccessor, ... to access the 
 * attached data.
 * Please note that you may only use a subset-attachment-accessor
 * with elements that are contained in the subset for which the
 * accessor has been created.
 *
 * Subset-attachments currently do not support pass-on behaviours.
 */
class ISubsetHandler : public GridObserver
{
	friend class attachment_traits<VertexBase*, ISubsetHandler>;
	friend class attachment_traits<EdgeBase*, ISubsetHandler>;
	friend class attachment_traits<Face*, ISubsetHandler>;
	friend class attachment_traits<Volume*, ISubsetHandler>;

	public:
		typedef AttachmentPipe<VertexBase*, ISubsetHandler>	VertexAttachmentPipe;
		typedef AttachmentPipe<EdgeBase*, ISubsetHandler>	EdgeAttachmentPipe;
		typedef AttachmentPipe<Face*, ISubsetHandler>		FaceAttachmentPipe;
		typedef AttachmentPipe<Volume*, ISubsetHandler>		VolumeAttachmentPipe;

	public:
	///	pass an or-combination of SubsetHandlerElements to supportedElements.
	/**	supportedElements define the elements on which the SubsetHandler works.
	 *	Default is SHE_ALL (all element-types).*/
		ISubsetHandler(uint supportedElements = SHE_ALL);

	///	assigns subsets based on the subsets in the given subset-handler
	/**	The constructed subset-handler will work on the same grid as the
	 *	subset handler that was passed to the constructor.
	 *	Elements will be assigned to the same subsets as in the given handler.
	 *	All properties are copied too.
	 *
	 *	Please note, that attachments are not copied in the current version.*/
		ISubsetHandler(const ISubsetHandler& sh);

	/**	The destructor automatically unregisters the subset-handler from the grid.
	 *	on deregistration erase_subset_lists of the derived class will be called.*/
		virtual ~ISubsetHandler();

	///	assigns subsets based on the subsets in the given subset-handler
	/**	Elements of this handler will be assigned to subsets based on their
	 *	order in the underlying grids.
	 *	The underlying grid of this handler will not be changed. This is particularly
	 *	useful if you just copied a grid and if you now want to copy the subsets
	 *	in the associated subset-handlers.
	 *	
	 *	Please note, that attachments are not copied in the current version.*/
		virtual ISubsetHandler& operator = (const ISubsetHandler& sh);
		
	///	returns a pointer to the grid on which the subset-handler works.
	/**	returns NULL if no grid is assigned.*/
		Grid* get_assigned_grid() const;
		
	///	returns true if the given element-types are supported.
	/**	pass an or-combination of constants enumerated in SubsetHandlerElements.*/
		bool elements_are_supported(uint shElements);

	///	set the type of elements that shall be handled by the SubsetHandler.
	/**	Pass an or-combination of constants enumerated in SubsetHandlerElements.
	 *	\sa SubsetHandler::enable_element_support*/
		void set_supported_elements(uint shElements);

	///	enable support for element-types. Does not invalidate previous settings.
	/**	pass an or-combination of constants enumerated in SubsetHandlerElements.*/
		void enable_element_support(uint shElements);

	///	disable support for element-types.
	/**	pass an or-combination of constants enumerated in SubsetHandlerElements.*/
		void disable_element_support(uint shElements);
		
	/**	new elements will be automatically assigned to this subset.
	 * 	set this to a negative value to avoid automatic assignment (-1 by default).
	 *	only used if subset_inheritance is disabled or if no parent is specified.*/
		void set_default_subset_index(int subsetIndex);
		inline int get_default_subset_index()	{return m_defaultSubsetIndex;}

	/**	if enabled, newly created elements derive their subset-index from their parents.
	 *	Enabled by default.
	 *	If enabled, the default subset index will be ignored if a parent is specified
	 *	on element creation.*/
		void enable_subset_inheritance(bool bEnable);
		inline bool subset_inheritance_enabled()	{return m_bSubsetInheritanceEnabled;}
		
	///	if the subset with the given index does not yet exist, it will be created.
	/**	All subsets in between num_subsets and index will be created, too.*/
		inline void subset_info_required(int index);

	///	returns the number of subset-infos
		inline uint num_subset_infos() const		{return (uint)m_subsetInfos.size();}
		
	/** if the subset at subsetIndex does not yet exist, it will be created.*/
		void set_subset_info(int subsetIndex, const SubsetInfo& subsetInfo);

	/** if the subset at subsetIndex does not yet exist, it will be created.*/
		SubsetInfo& subset_info(int subsetIndex);

	/** if the subset at subsetIndex does not yet exist, it will be created.*/
		const SubsetInfo& subset_info(int subsetIndex) const;
		
		void clear();
		void clear_subset(int subsetIndex);
		void clear_subsets();


	///	inserts a subset at the given index. Moves all other subsets 1 index higher.
		void insert_subset(int subsetIndex);///< changes subset-indices of other subsets.
	///	erases the subset at the given index. Assigns -1 to all entries. Moves all other subsets 1 index up.
		void erase_subset(int subsetIndex);///< changes subset-indices of other subsets.
	///	Swaps the given subsets,
		void swap_subsets(int subsetIndex1, int subsetIndex2);
	///	Moves the subset from index From to index To. Moves all subsets between indexFrom+1 and indexTo in the opposite direction.
		void move_subset(int indexFrom, int indexTo);///< changes subset indices of other subsets.

		template <class TIterator>
		void assign_subset(TIterator iterBegin, TIterator iterEnd, int subsetIndex);

		int get_subset_index(GeometricObject* elem) const;
		inline int get_subset_index(VertexBase* elem) const	{return m_aaSubsetIndexVRT[elem];}
		inline int get_subset_index(EdgeBase* elem) const	{return m_aaSubsetIndexEDGE[elem];}
		inline int get_subset_index(Face* elem) const		{return m_aaSubsetIndexFACE[elem];}
		inline int get_subset_index(Volume* elem) const		{return m_aaSubsetIndexVOL[elem];}

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
		
	/**	The implementation in a derived class should store the element in a list
	 *	and call subset_assigned with the iterators position and the subset-index.
	 *	The iterator can later be retrieved with get_list_iterator(...).
	 *	The index can be retrieved with get_subset_index(...).*/
		virtual void assign_subset(VertexBase* elem, int subsetIndex) = 0;

	/**	The implementation in a derived class should store the element in a list
	 *	and call subset_assigned with the iterators position and the subset-index.
	 *	The iterator can later be retrieved with get_list_iterator(...).
	 *	The index can be retrieved with get_subset_index(...).*/
		virtual void assign_subset(EdgeBase* elem, int subsetIndex) = 0;

	/**	The implementation in a derived class should store the element in a list
	 *	and call subset_assigned with the iterators position and the subset-index.
	 *	The iterator can later be retrieved with get_list_iterator(...).
	 *	The index can be retrieved with get_subset_index(...).*/
		virtual void assign_subset(Face* elem, int subsetIndex) = 0;

	/**	The implementation in a derived class should store the element in a list
	 *	and call subset_assigned with the iterators position and the subset-index.
	 *	The iterator can later be retrieved with get_list_iterator(...).
	 *	The index can be retrieved with get_subset_index(...).*/
		virtual void assign_subset(Volume* elem, int subsetIndex) = 0;
	
	
	////////////////////////////////
	//	attachments
	
	///	enable subset-attachment support
	/**	if subset-attachments are enabled you may attach data to the elements
	 *	of a subset. This is useful if you want to store different data in the
	 *	elements of different subsets.*/
		void enable_subset_attachments(bool bEnable);
		
	///	returns true if subset-attachments are enabled.
		inline bool subset_attachments_are_enabled()	{return m_bSubsetAttachmentsEnabled;};

		inline uint get_attachment_data_index(const VertexBase* v) const	{return m_aaDataIndVRT[v];}
		inline uint get_attachment_data_index(const EdgeBase* e) const		{return m_aaDataIndEDGE[e];}
		inline uint get_attachment_data_index(const Face* f) const			{return m_aaDataIndFACE[f];}
		inline uint get_attachment_data_index(const Volume* v) const		{return m_aaDataIndVOL[v];}

	///	attach with unspecified default value.
	/**	Pass either VertexBase, EdgeBase, Face or Volume as TGeomObjClass.*/
		template <class TGeomObjClass>
		inline void attach_to(IAttachment& attachment, int subsetIndex);

	///	attach with specified default value
	/**	Pass either VertexBase, EdgeBase, Face or Volume as TGeomObjClass.*/
		template <class TGeomObjClass, class TAttachment>
		void attach_to_dv(TAttachment& attachment, int subsetIndex,
						const typename TAttachment::ValueType& defaultValue);

	//	detach
	/**	Pass either VertexBase, EdgeBase, Face or Volume as TGeomObjClass.*/
		template <class TGeomObjClass>
		void detach_from(IAttachment& attachment, int subsetIndex);

	////////////////////////////////
	//	attachments helper
	//	attach with attachments default pass-on behaviour
		inline void attach_to_vertices(IAttachment& attachment, int subsetIndex)	{attach_to<VertexBase>(attachment, subsetIndex);}
		inline void attach_to_edges(IAttachment& attachment, int subsetIndex)		{attach_to<EdgeBase>(attachment, subsetIndex);}
		inline void attach_to_faces(IAttachment& attachment, int subsetIndex)		{attach_to<Face>(attachment, subsetIndex);}
		inline void attach_to_volumes(IAttachment& attachment, int subsetIndex)		{attach_to<Volume>(attachment, subsetIndex);}

	//	attach with default value and attachments default pass-on behaviour
		template <class TAttachment>
		inline void attach_to_vertices_dv(TAttachment& attachment, int subsetIndex, const typename TAttachment::ValueType& defaultValue)	{attach_to_dv<VertexBase>(attachment, subsetIndex, defaultValue);}
		template <class TAttachment>
		inline void attach_to_edges_dv(TAttachment& attachment, int subsetIndex, const typename TAttachment::ValueType& defaultValue)		{attach_to_dv<EdgeBase>(attachment, subsetIndex, defaultValue);}
		template <class TAttachment>
		inline void attach_to_faces_dv(TAttachment& attachment, int subsetIndex, const typename TAttachment::ValueType& defaultValue)		{attach_to_dv<Face>(attachment, subsetIndex, defaultValue);}
		template <class TAttachment>
		inline void attach_to_volumes_dv(TAttachment& attachment, int subsetIndex, const typename TAttachment::ValueType& defaultValue)		{attach_to_dv<Volume>(attachment, subsetIndex, defaultValue);}

	//	detach
		inline void detach_from_vertices(IAttachment& attachment, int subsetIndex)	{detach_from<VertexBase>(attachment, subsetIndex);}
		inline void detach_from_edges(IAttachment& attachment, int subsetIndex)		{detach_from<EdgeBase>(attachment, subsetIndex);}
		inline void detach_from_faces(IAttachment& attachment, int subsetIndex)		{detach_from<Face>(attachment, subsetIndex);}
		inline void detach_from_volumes(IAttachment& attachment, int subsetIndex)	{detach_from<Volume>(attachment, subsetIndex);}

	///	returns the attachment data container for elements of type TGeomObj for the given subset.
	/**	Use the data-container with care! You should never clear or resize it.
	 *
	 *	Valid types for TGeomObj are VertexBase, EdgeBase, Face and Volume.
	 *	call it like this (let sh be an instance of ISubsetHandler):
	 *	sh.get_attachment_data_container<VertexBase>(aSomeAttachment, someSubsetIndex);*/
		template <class TGeomObj, class TAttachment>
		inline typename TAttachment::ContainerType*
		get_attachment_data_container(TAttachment& attachment, int subsetIndex);

	protected:
		typedef SectionContainer<GeometricObject*, std::list<GeometricObject*> >	SectionContainer;
		typedef SectionContainer::iterator iterator;
		
	protected:
	///	selects elements based on the selection in the srcHandler
		void assign_subset_handler(const ISubsetHandler& sh);
		
	///	set the grid on which the subset-handler shall work.
	/**	The subset-handler can only work on one grid at a time.
	 *	It is cruicial that assign_grid methods of derived classes call
	 *	this method.*/
		void assign_grid(Grid& grid);
		
	///	sets the subset-indices of all elements of m_pGrid to -1.
	/**	Use with care! Only indices are affected. The elements are not
	 *	removed from any lists.
	 *	pass an or-combination of constants enumerated in SubsetHandlerElements.*/
		void reset_subset_indices(uint shElements = SHE_ALL);
		
	///	creates all required infos (and pipes) up to the given index.
		void create_required_subset_infos(int index);
		
		inline void subset_assigned(VertexBase* v, iterator iter, int subsetIndex);
		inline void subset_assigned(EdgeBase* e, iterator iter, int subsetIndex);
		inline void subset_assigned(Face* f, iterator iter, int subsetIndex);
		inline void subset_assigned(Volume* v, iterator iter, int subsetIndex);
		
		inline iterator get_list_iterator(VertexBase* v)	{return m_aaIteratorVRT[v];}
		inline iterator get_list_iterator(EdgeBase* e)		{return m_aaIteratorEDGE[e];}
		inline iterator get_list_iterator(Face* f)			{return m_aaIteratorFACE[f];}
		inline iterator get_list_iterator(Volume* v)		{return m_aaIteratorVOL[v];}
		
	/**	alters the subset index only. Suited as a helper for methods like 
	 *	change_subset_indices or reset_subset_indices.
	 *	WARNING: This method only alters the index but does not actually
	 *	move the element to another subset. Use assign_subset instead for this task.*/	
		inline void alter_subset_index(VertexBase* v, int subsetIndex)	{m_aaSubsetIndexVRT[v] = subsetIndex;}
	/**	alters the subset index only. Suited as a helper for methods like 
	 *	change_subset_indices or reset_subset_indices.
	 *	WARNING: This method only alters the index but does not actually
	 *	move the element to another subset. Use assign_subset instead for this task.*/
		inline void alter_subset_index(EdgeBase* e, int subsetIndex)	{m_aaSubsetIndexEDGE[e] = subsetIndex;}
	/**	alters the subset index only. Suited as a helper for methods like 
	 *	change_subset_indices or reset_subset_indices.
	 *	WARNING: This method only alters the index but does not actually
	 *	move the element to another subset. Use assign_subset instead for this task.*/
		inline void alter_subset_index(Face* f, int subsetIndex)		{m_aaSubsetIndexFACE[f] = subsetIndex;}
	/**	alters the subset index only. Suited as a helper for methods like 
	 *	change_subset_indices or reset_subset_indices.
	 *	WARNING: This method only alters the index but does not actually
	 *	move the element to another subset. Use assign_subset instead for this task.*/
		inline void alter_subset_index(Volume* v, int subsetIndex)		{m_aaSubsetIndexVOL[v] = subsetIndex;}
		
		virtual void erase_subset_lists() = 0;
		
		virtual void clear_subset_lists(int index) = 0;
		
		virtual void change_subset_indices(int indOld, int indNew) = 0;
		
		
	///	add a subset if requiered - so that the subset with maxIndex exists.
		virtual void add_required_subset_lists(int maxIndex) = 0;
		
	///	erase the subset-lists but do not touch the subset-indices.
		virtual void erase_subset_lists(int index) = 0;

	///	swap the subset-lists but do not touch the subset-indices.
		virtual void swap_subset_lists(int ind1, int ind2) = 0;

	///	move the subset-lists but do not touch the subset-indices.
		virtual void move_subset_lists(int indexFrom, int indexTo) = 0;

	
	////////////////////////////////
	//	attachments
		inline void set_attachment_data_index(VertexBase* v, uint index)	{m_aaDataIndVRT[v] = index;}
		inline void set_attachment_data_index(EdgeBase* e, uint index)		{m_aaDataIndEDGE[e] = index;}
		inline void set_attachment_data_index(Face* f, uint index)			{m_aaDataIndFACE[f] = index;}
		inline void set_attachment_data_index(Volume* v, uint index)		{m_aaDataIndVOL[v] = index;}

		void resize_attachment_pipes(size_t newSize);
		void clear_attachment_pipes(int subsetIndex);
		void clear_attachment_pipes();
	
		template <class TGeomObj>
		inline AttachmentPipe<TGeomObj*, ISubsetHandler>&
		get_attachment_pipe(int subsetIndex);

	////////////////////////////////
	//	virtual methods for attachments
	///	this method is called by ISubsetHandler when attachment_support has been enabled.
	/**	derived classes have to implement this method.
	 *	during this method they have to call register_at_pipe for all elements
	 *	that are contained in one of the subsets.
	 *	WARNING: This method is cruical for the attachment system.
	 *	You should never call it yourself.*/
		virtual void register_subset_elements_at_pipe() = 0;
		
	///	this method should be called during \sa register_subset_elements_at_pipe.
	/**	WARNING: This method is cruical for the attachment system.
	 *	You should only call it during \sa register_subset_elements_at_pipe
	 *	Only call this method for elements that are contained in a subset.*/
		inline void register_at_pipe(VertexBase* elem)	{m_vertexAttachmentPipes[get_subset_index(elem)]->register_element(elem);}

	///	this method should be called x \sa register_subset_elements_at_pipe.
	/**	WARNING: This method is cruical for the attachment system.
	 *	You should only call it during \sa register_subset_elements_at_pipe
	 *	Only call this method for elements that are contained in a subset.*/
		inline void register_at_pipe(EdgeBase* elem)	{m_edgeAttachmentPipes[get_subset_index(elem)]->register_element(elem);}

	///	this method should be called during \sa register_subset_elements_at_pipe.
	/**	WARNING: This method is cruical for the attachment system.
	 *	You should only call it during \sa register_subset_elements_at_pipe
	 *	Only call this method for elements that are contained in a subset.*/
		inline void register_at_pipe(Face* elem)		{m_faceAttachmentPipes[get_subset_index(elem)]->register_element(elem);}

	///	this method should be called during \sa register_subset_elements_at_pipe.
	/**	WARNING: This method is cruical for the attachment system.
	 *	You should only call it during \sa register_subset_elements_at_pipe
	 *	Only call this method for elements that are contained in a subset.*/
		inline void register_at_pipe(Volume* elem)		{m_volumeAttachmentPipes[get_subset_index(elem)]->register_element(elem);}

	public:
	///	attachment accessor grants access to data associated with elements of a subset.
	/**	Valid types for TGeomObj are VertexBase, EdgeBase, Face and Volume*/
		template <class TGeomObj, class TAttachment>
		class AttachmentAccessor : public ug::AttachmentAccessor<TGeomObj*, TAttachment, ISubsetHandler>
		{
			protected:
				typedef ug::AttachmentAccessor<TGeomObj*, TAttachment, ISubsetHandler>	BaseClass;
				
			public:
				AttachmentAccessor()	{}
				AttachmentAccessor(const AttachmentAccessor& aa) : BaseClass(aa)	{}
				AttachmentAccessor(ISubsetHandler& sh, TAttachment& a, int subsetIndex) : BaseClass(sh.get_attachment_pipe<TGeomObj>(subsetIndex), a)	{}

				inline void access(ISubsetHandler& sh, TAttachment& a, int subsetIndex)
					{BaseClass::access(sh.get_attachment_pipe<TGeomObj>(subsetIndex), a);}
		};


	///	the multi-subset-attachment-accessor allows to access an attachment that has been attached to multiple subsets.
	/**	Note that this accessor is slower than the normal attachment-accessors.
	 *	If subsets are added to the subset-handler or attachments are added/removed,
	 *	you should not rely on the accessor to still work properly.
	 *	Call access() in that case to re-validate him.*/
/*
		template <class TAttachmet, class TGeomBaseObj>
		class MultiSubsetAttachmentAccessor
		{
		//TODO: implement this!
		};
*/
	protected:
		typedef AInt					ASubsetIndex;
		typedef Attachment<iterator>	AIterator;
		typedef Attachment<uint>		ADataIndex;
		typedef std::vector<SubsetInfo>	SubsetInfoVec;
		typedef std::vector<VertexAttachmentPipe*>	VertexAttachmentPipeVec;
		typedef std::vector<EdgeAttachmentPipe*>		EdgeAttachmentPipeVec;
		typedef std::vector<FaceAttachmentPipe*>		FaceAttachmentPipeVec;
		typedef std::vector<VolumeAttachmentPipe*>	VolumeAttachmentPipeVec;

	protected:
		Grid*				m_pGrid;
		VertexAttachmentPipeVec	m_vertexAttachmentPipes;
		EdgeAttachmentPipeVec	m_edgeAttachmentPipes;
		FaceAttachmentPipeVec	m_faceAttachmentPipes;
		VolumeAttachmentPipeVec	m_volumeAttachmentPipes;

		SubsetInfoVec	m_subsetInfos;
		uint			m_supportedElements;

		ASubsetIndex	m_aSubsetIndex;
		AIterator		m_aIterator;
		ADataIndex		m_aDataIndex;
		
		int				m_defaultSubsetIndex;
		bool			m_bSubsetInheritanceEnabled;
		bool			m_bSubsetAttachmentsEnabled;
				
		Grid::VertexAttachmentAccessor<ASubsetIndex>	m_aaSubsetIndexVRT;
		Grid::EdgeAttachmentAccessor<ASubsetIndex>		m_aaSubsetIndexEDGE;
		Grid::FaceAttachmentAccessor<ASubsetIndex>		m_aaSubsetIndexFACE;
		Grid::VolumeAttachmentAccessor<ASubsetIndex>	m_aaSubsetIndexVOL;

		Grid::VertexAttachmentAccessor<AIterator>		m_aaIteratorVRT;
		Grid::EdgeAttachmentAccessor<AIterator>			m_aaIteratorEDGE;
		Grid::FaceAttachmentAccessor<AIterator>			m_aaIteratorFACE;
		Grid::VolumeAttachmentAccessor<AIterator>		m_aaIteratorVOL;

		Grid::VertexAttachmentAccessor<ADataIndex>		m_aaDataIndVRT;
		Grid::EdgeAttachmentAccessor<ADataIndex>		m_aaDataIndEDGE;
		Grid::FaceAttachmentAccessor<ADataIndex>		m_aaDataIndFACE;
		Grid::VolumeAttachmentAccessor<ADataIndex>		m_aaDataIndVOL;
};

}//	end of namespace

////////////////////////////////////////////////
//	include implementation
#include "subset_handler_interface_impl.hpp"

#endif
