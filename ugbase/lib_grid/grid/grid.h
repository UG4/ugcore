//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y08 m10 d09

#ifndef __H__LIB_GRID__GRID__
#define __H__LIB_GRID__GRID__

#include <vector>
#include <list>
#include "grid_constants.h"
#include "geometric_base_objects.h"
#include "grid_observer.h"
#include "common/util/section_container.h"
#include "geometric_object_collection.h"

namespace ug
{

/// \addtogroup lib_grid
///	@{

////////////////////////////////////////////////////////////////////////////////////////////////
//	Grid
///	Manages the elements of a grid and their interconnection.
/**
 * The Grid class is the heart of libGrid. It can be used to create elements of
 * custom types, for example Vertices, Triangles or Tetrahedrons.
 * All elements have to be derived from either VertexBase, EdgeBase, Face or Volume
 * and have to specialize the template-class geometry_traits
 * (both defined in geometric_base_objects.h).
 * The grid can automatically create information about connected geometric objects.
 * You could for example query a face for associated volumes.
 * Associated elements are however only stored, if the according options are set.
 * This reduces storage space if only particular associations are required.
 * Since these options can be changed dynamically at runtime on could in practice always
 * start with a grid holding minimal connection informations. Options are then enabled
 * on the fly, as required by called methods. Note however that some algorithms do
 * run faster if particular options are enabled.
 * On the other hand a grid serves as the connection between geometric objects and their
 * associated data. This connection is realized by so called attachments.
 * You can attach data to vertices, edges, faces, and volumes. Note, that data is attached
 * to all elements of one type at once. No sub-types are  considered. This allows for
 * efficient storage and data-access.
 * thanks to this way of storing data, methods and classes can temporarily attach data
 * to the grids elements and thus store and access data of given elements in an
 * efficient way. The data-access is handled by so called attachment-accessors.
 * Grid defines one for each base-type as well as a generic one.
 *
 * The Grid class itself only knows about the concepts of Vertex, Edge, Face
 * and Volume (classes VertexBase, EdgeBase, Face, Volume).
 * In order to implement algorithms one needs more specialized element-types,
 * like Triangles, Quadrilaterals or Tetrahedrons.
 * Those special types are supported by the use of templates. All that a custom type
 * has to satisfy is that it is derived from one of the four basic elements and
 * that it specializes the class geometry_traits.
 * By calling Grids methods with those specialized types as template parameters,
 * you can create, erase and iterate over specialized elements.
 * You could for example iterate over all faces using the calls
 * Grid::begin<Face>() and Grid::end<Face>()  (or by calling Grid::faces_begin() and Grid::faces_end()).
 * if you want to iterate only over elements of type Triangle (given this type has
 * been defined), you could call
 * Grid::begin<Triangle>() and Grid::end<Triangle>().
 * This is applicable to all methods that take an TGeomObj as their template parameter.
 * If a template function takes the parameter TGeomObjClass, then only one of the four
 * base objects should be passed (VertexBase, EdgeBase, Face and Volume).
 * This use if templates allows for a arbitrary number of custom elements, without
 * requiring any changes to the Grid-class itself.
 */
class Grid
{
	public:
	///	the attachment-pipe used by Grid
		typedef ug::AttachmentPipe<GeometricObject*, Grid>	AttachmentPipe;

	///	the generic attachment-accessor for access to grids attachment pipes.
		template <class TElem, class TAttachment>
		class AttachmentAccessor : public ug::AttachmentAccessor<GeometricObject*, TAttachment, Grid>
		{
			public:
				AttachmentAccessor();
				AttachmentAccessor(const AttachmentAccessor& aa);
				AttachmentAccessor(Grid& grid, TAttachment& a);

				inline void access(Grid& grid, TAttachment& a)
					{ug::AttachmentAccessor<GeometricObject*, TAttachment, Grid>::access(grid.get_attachment_pipe<TElem>(), a);}
		};

	//	half-specialized AttachmentAccessors:
		template <class TAttachment>
		class VertexAttachmentAccessor : public AttachmentAccessor<VertexBase, TAttachment>
		{
			public:
				VertexAttachmentAccessor();
				VertexAttachmentAccessor(const VertexAttachmentAccessor& aa);
				VertexAttachmentAccessor(Grid& grid, TAttachment& a);
		};

		template <class TAttachment>
		class EdgeAttachmentAccessor : public AttachmentAccessor<EdgeBase, TAttachment>
		{
			public:
				EdgeAttachmentAccessor();
				EdgeAttachmentAccessor(const EdgeAttachmentAccessor& aa);
				EdgeAttachmentAccessor(Grid& grid, TAttachment& a);
		};

		template <class TAttachment>
		class FaceAttachmentAccessor : public AttachmentAccessor<Face, TAttachment>
		{
			public:
				FaceAttachmentAccessor();
				FaceAttachmentAccessor(const FaceAttachmentAccessor& aa);
				FaceAttachmentAccessor(Grid& grid, TAttachment& a);
		};

		template <class TAttachment>
		class VolumeAttachmentAccessor : public AttachmentAccessor<Volume, TAttachment>
		{
			public:
				VolumeAttachmentAccessor();
				VolumeAttachmentAccessor(const VolumeAttachmentAccessor& aa);
				VolumeAttachmentAccessor(Grid& grid, TAttachment& a);
		};

		typedef std::list<VertexBase*>	VertexContainer;
		typedef std::list<EdgeBase*>	EdgeContainer;
		typedef std::list<Face*>		FaceContainer;
		typedef std::list<Volume*>		VolumeContainer;
	///	used to iterate over associated edges of vertices, faces and volumes
		typedef EdgeContainer::iterator 	AssociatedEdgeIterator;
	///	used to iterate over associated faces of vertices, edges and volumes
		typedef FaceContainer::iterator 	AssociatedFaceIterator;
	///	used to iterate over associated volumes of vertices, edges and faces
		typedef VolumeContainer::iterator 	AssociatedVolumeIterator;

	public:
	////////////////////////////////////////////////
	//	constructors and destructor
	///	initialises the grid and sets the option GRIDOPT_DEFAULT.
	/**	Pass a custom option to the alternative constructor in order
	 *  to initialise the grid with your options.
	 *	\sa Grid(uint options)*/
		Grid();

	///	initialises the grid with the given option.
	/**	pass an or-combination of constants enumerated in
	 *  VertexOptions, EdgeOptions, FaceOptions, VolumeOptions and GridOptions.*/
		Grid(uint options);

	///	copies all elements and some attachments from the passed grid to this grid.
	/**	While all elements and the options are copied completly from the source-grid,
	 *	the attachments are only copied if their pass-on behaviour is set to true.*/
		Grid(const Grid& grid);

		virtual ~Grid();

	///	copies all elements and some attachments from the passed grid to this grid.
	/**	While all elements and the options are copied completly from the source-grid,
	 *	the attachments are only copied if their pass-on behaviour is set to true.
	 *	Attachments that were already attached to this grid are removed prior to
	 *	copying if their pass-on behaviour was set to true. They will be kept otherwise.
	 *	This is relevant to ensure that observers like GridSubsetHandler
	 *	will work after the assignment.*/	
		Grid& operator = (const Grid& grid);
		
	////////////////////////////////////////////////
	//	grid options
	/**	you can pass any option enumerated in GridOptions or specify a custom option
	 *	using an or-combination of the constants enumerated in
	 *	VertexOptions, EdgeOptions, FaceOptions and VolumeOptions.
	 *	See GridOptions for possible combinations.*/
		void set_options(uint options);
		uint get_options() const;
		void enable_options(uint options);	///< see set_options for a description of valid parameters.
		void disable_options(uint options);	///< see set_options for a description of valid parameters.
		bool option_is_enabled(uint option) const;///< see set_options for a description of valid parameters.

		void clear();
		void clear_geometry();
		void clear_attachments();

	////////////////////////////////////////////////
	//	element creation
	///	create a custom element.
	/**
	 * TGeomObj has to be a geometric object type as described in geometric_base_objects.h.
	 * You may optionally specify a GeometricObject pParent.
	 * pParent may be used by observers to initialize the created object in a specific way.
	 */
		template<class TGeomObj>
		typename geometry_traits<TGeomObj>::iterator
		create(GeometricObject* pParent = NULL);

	///	create a custom element from a descriptor.
	/**
	 * TGeomObj has to be a geometric object type as described in geometric_base_objects.h.
	 * You may optionally specify a GeometricObject pParent.
	 * pParent may be used by observers to initialize the created object in a specific way.
	 */
		template <class TGeomObj>
		typename geometry_traits<TGeomObj>::iterator
		create(const typename geometry_traits<TGeomObj>::Descriptor& descriptor,
				GeometricObject* pParent = NULL);

	///	create a custom element and replaces an old one.
	/**
	 * Be sure, that TGeomObj has the same reference-element as pReplaceMe.
	 * Similar to create. With pReplaceMe you may specify an element of the same base-type
	 * as TGeomObj, which shall be replaced. The replaced element will be deleted from grid.
	 * pReplaceMe will be treated as pParent for GridObservers.
	 *
	 * Calls pass_on_values.
	 *
	 * Notes for GridObservers:
	 * create_and_replace will call in the given order (replace elem with the appropriate name).
	 * - elem_created(newElem)
	 * - elem_to_be_replaced(oldElem, newElem)
	 * - elem_to_be_deleted(oldElem)
	 */
		template <class TGeomObj>
		typename geometry_traits<TGeomObj>::iterator
		create_and_replace(typename geometry_traits<TGeomObj>::geometric_base_object* pReplaceMe);

	///	this method creates a new vertex, which has the same type as pCloneMe.
		VertexBaseIterator create_by_cloning(VertexBase* pCloneMe, GeometricObject* pParent = NULL);

	///	this method creates a new edge, which has the same type as pCloneMe.
		EdgeBaseIterator create_by_cloning(EdgeBase* pCloneMe, const EdgeVertices& ev, GeometricObject* pParent = NULL);

	///	this method creates a new face, which has the same type as pCloneMe.
		FaceIterator create_by_cloning(Face* pCloneMe, const FaceVertices& fv, GeometricObject* pParent = NULL);

	///	this method creates a new volume, which has the same type as pCloneMe.
		VolumeIterator create_by_cloning(Volume* pCloneMe, const VolumeVertices& vv, GeometricObject* pParent = NULL);

	////////////////////////////////////////////////
	//	element deletion
		void erase(GeometricObject* geomObj);
		void erase(VertexBase* vrt);
		void erase(EdgeBase* edge);
		void erase(Face* face);
		void erase(Volume* vol);

		template <class GeomObjIter>
		void erase(const GeomObjIter& iterBegin, const GeomObjIter& iterEnd);

	////////////////////////////////////////////////
	//	replace
	///	Replace vrtOld with vrtNew.
	/**
	 * WARNING: USE THIS METHOD WITH CARE!
	 * vrtOld and vrtNew have both to be registered vertices of this grid.
	 * vrtOld will be erased during this method.
	 * Make sure that no geometric object in the grid references both
	 * vrtOld and vrtNew.
	 * This method iterates through all geometric objects that are
	 * connected with vrtOld and replaces vrtOld by vrtNew in each.
	 * Connectivity information will be updated on the fly.
	 * Elements that reference both vrtOld and vrtNew will reference
	 * vrtNew two times after the completion of replace_vertex.
	 * This leads to degenerate elements and most likely to bad
	 * runtime behavior.
	 *
	 * requires options in GRIDOPT_VERTEXCENTRIC_INTERCONNECTION.
	 */
		bool replace_vertex(VertexBase* vrtOld, VertexBase* vrtNew);

	///	checks if replace_vertex would be a valid operation
	/**
	 * Checks if a call of replace_vertex with vertices vrtOld and vrtNew
	 * would lead to degenerate elements. If so false is returned.
	 * If not true is returned.
	 *
	 * requires options in GRIDOPT_VERTEXCENTRIC_INTERCONNECTION.
	 */
		bool replace_vertex_is_valid(VertexBase* vrtOld, VertexBase* vrtNew);

	/*
	///	VrtPairIterator has to be an iterator with value-type std::pair<VertexBase*, VertexBase*>
		template <class VrtPairIter>
		void replace_vertices(VrtPairIter& iterBegin, VrtPairIter& iterEnd);
	*/

	////////////////////////////////////////////////
	//	geometric-object-collection
	///	returns the the GeometricObjectCollection of the grid:
		virtual GeometricObjectCollection get_geometric_object_collection();

	////////////////////////////////////////////////
	///	flips the orientation of a face.
		void flip_orientation(Face* f);

	////////////////////////////////////////////////
	///	flips the orientation of a volume.
		void flip_orientation(Volume* vol);
		
	////////////////////////////////////////////////
	//	Iterators
		template <class TGeomObj>
		typename geometry_traits<TGeomObj>::iterator
		begin();

		template <class TGeomObj>
		typename geometry_traits<TGeomObj>::iterator
		end();

		inline VertexBaseIterator	vertices_begin()	{return begin<VertexBase>();}
		inline VertexBaseIterator	vertices_end()		{return end<VertexBase>();}
		inline EdgeBaseIterator		edges_begin()		{return begin<EdgeBase>();}
		inline EdgeBaseIterator		edges_end()			{return end<EdgeBase>();}
		inline FaceIterator			faces_begin()		{return begin<Face>();}
		inline FaceIterator			faces_end()			{return end<Face>();}
		inline VolumeIterator		volumes_begin()		{return begin<Volume>();}
		inline VolumeIterator		volumes_end()		{return end<Volume>();}
	
		template <class TGeomObj>
		typename geometry_traits<TGeomObj>::const_iterator
		begin() const;

		template <class TGeomObj>
		typename geometry_traits<TGeomObj>::const_iterator
		end() const;

	//	element manipulation
	/*
		void set_vertices(EdgeBase* e, EdgeDesctiptor& ed);
		void set_vertices(Face* f, FaceDescriptor& fd);
		void set_vertices(Volume* v, VolumeDescriptor& vd);
	*/

	////////////////////////////////////////////////
	//	element numbers
		template <class TGeomObj>
		size_t num() const;
		inline size_t num_vertices() const	{return num<VertexBase>();}
		inline size_t num_edges() const		{return num<EdgeBase>();}
		inline size_t num_faces()	const		{return num<Face>();}
		inline size_t num_volumes()const		{return num<Volume>();}

		size_t vertex_fragmentation();	///< returns the number of unused vertex-data-entries.
		size_t edge_fragmentation();		///< returns the number of unused edge-data-entries.
		size_t face_fragmentation();		///< returns the number of unused face-data-entries.
		size_t volume_fragmentation();	///< returns the number of unused volume-data-entries.

	///	returns the size of the associated attachment containers.
	/**	valid types for TGeomObj are VertexBase, EdgeBase, Face, Volume.*/
		template <class TGeomObj>
		size_t attachment_container_size() const;

	////////////////////////////////////////////////
	//	connectivity-information
	///	returns the edge between v1 and v2, if it exists. Returns NULL if not.
		EdgeBase* get_edge(VertexBase* v1, VertexBase* v2);
	///	returns the edge that is described by ev.
	/**	Note that you may pass an EdgeDescriptor to this method.*/
		EdgeBase* get_edge(EdgeVertices& ev);
	///	If it exists, this method returns the i-th edge of the given Face. If not NULL is returned.
	/**	To make sure that associated edges always exist, enable the grid-option
	 *	FACEOPT_AUTOGENERATE_EDGES.
	 *	For maximal performance, the option FACEOPT_STORE_ASSOCIATED_EDGES
	 *	should be enabled.*/
		EdgeBase* get_edge(Face* f, int ind);
	///	If it exists, this method returns the i-th edge of the given Volume. If not NULL is returned.
	/**	To make sure that associated edges always exist, enable the grid-option
	 *	VOLOPT_AUTOGENERATE_EDGES.
	 *	For maximal performance, the option VOLOPT_STORE_ASSOCIATED_EDGES
	 *	should be enabled.*/
		EdgeBase* get_edge(Volume* v, int ind);
	///	returns the face that is described by fv.
	/**	Note that you may pass a FaceDescriptor to this method.*/
		Face* get_face(FaceVertices& fv);
	///	If it exists, this method returns the i-th face of the given Volume. If not NULL is returned.
	/**	To make sure that associated faces always exist, enable the grid-option
	 *	VOLOPT_AUTOGENERATE_FACES.
	 *	For maximal performance, the option VOLOPT_STORE_ASSOCIATED_FACES
	 *	should be enabled.*/
		Face* get_face(Volume* v, int ind);
	///	returns the volume that is described by ev.
	/**	Note that you may pass an VolumeDescriptor to this method.*/
		Volume* get_volume(VolumeVertices& vv);
		
	////////////////////////////////////////////////
	//	access to the sides of an geometric object
	///	This method returns the i-th side of an EdgeBase, Face or Volume.
	/**	If obj has dimension d, then all associated elements of dimension d-1
	 *	are regarded as sides. (Face -> EdgeBase). Only derivates of Volume,
	 *	Face or EdgeBase may be queried for their sides. You may not pass a general
	 * 	GeometricObject nor objects of a type not mentioned above.
	 *
	 *	It is not in all cases guaranteed that an object has sides.
	 *	If i.e. the FACEOPT_AUTOGENERATE_EDGES is not enabled in the grids options,
	 *	then it is not guaranteed that all side-edges of each face exist (in most
	 *	cases they wont exist). The method returns NULL in this case.
	 *	If however FACEOPT_AUTOGENERATE_EDGES is enabled then this method always
	 *	will return an edge if queried for the side of a face.
	 *	For volumes the appropriate option is called VOLOPT_AUTOGENERATE_FACES.
	 *
	 *	The method will be faster if elements store associated lower dimensional
	 *	elements. Options VOLOPT_STORE_ASSOCIATED_FACES and
	 *	FACEOPT_STORE_ASSOCIATED_EDGES have to be enabled for this.
	 */
		template <class TGeomObj>
		typename TGeomObj::lower_dim_base_object*
		get_side(TGeomObj* obj, size_t side);

	////////////////////////////////////////////////
	//	attachments
	////////////////////////////////////////////////////////////////////////
	///	attach with custom pass-on-behaviour and unspecified default value.
		template <class TGeomObjClass>
		void attach_to(IAttachment& attachment, bool passOnValues);

		inline void attach_to_vertices(IAttachment& attachment, bool passOnValues)	{attach_to<VertexBase>(attachment, passOnValues);}
		inline void attach_to_edges(IAttachment& attachment, bool passOnValues)		{attach_to<EdgeBase>(attachment, passOnValues);}
		inline void attach_to_faces(IAttachment& attachment, bool passOnValues)		{attach_to<Face>(attachment, passOnValues);}
		inline void attach_to_volumes(IAttachment& attachment, bool passOnValues)		{attach_to<Volume>(attachment, passOnValues);}

	////////////////////////////////////////////////////////////////////////
	///	attach with default pass-on behaviour and unspecified default value.
		template <class TGeomObjClass>
		inline void attach_to(IAttachment& attachment)	{attach_to<TGeomObjClass>(attachment, attachment.default_pass_on_behaviour());}

		inline void attach_to_vertices(IAttachment& attachment)	{attach_to<VertexBase>(attachment);}
		inline void attach_to_edges(IAttachment& attachment)	{attach_to<EdgeBase>(attachment);}
		inline void attach_to_faces(IAttachment& attachment)	{attach_to<Face>(attachment);}
		inline void attach_to_volumes(IAttachment& attachment)	{attach_to<Volume>(attachment);}

	////////////////////////////////////////////////////////////////////////
	//	attach with specified default value and default pass-on behaviour
		template <class TGeomObjClass, class TAttachment>
		void attach_to_dv(TAttachment& attachment, const typename TAttachment::ValueType& defaultValue);

		template <class TAttachment>
		inline void attach_to_vertices_dv(TAttachment& attachment, const typename TAttachment::ValueType& defaultValue)	{attach_to_dv<VertexBase>(attachment, defaultValue);}
		template <class TAttachment>
		inline void attach_to_edges_dv(TAttachment& attachment, const typename TAttachment::ValueType& defaultValue)	{attach_to_dv<EdgeBase>(attachment, defaultValue);}
		template <class TAttachment>
		inline void attach_to_faces_dv(TAttachment& attachment, const typename TAttachment::ValueType& defaultValue)	{attach_to_dv<Face>(attachment, defaultValue);}
		template <class TAttachment>
		inline void attach_to_volumes_dv(TAttachment& attachment, const typename TAttachment::ValueType& defaultValue)	{attach_to_dv<Volume>(attachment, defaultValue);}

	////////////////////////////////////////////////////////////////////////
	//	attach with specified default value and custom pass-on behaviour
		template <class TGeomObjClass, class TAttachment>
		void attach_to_dv(TAttachment& attachment, const typename TAttachment::ValueType& defaultValue, bool passOnValues);

		template <class TAttachment>
		inline void attach_to_vertices_dv(TAttachment& attachment, const typename TAttachment::ValueType& defaultValue, bool passOnValues)	{attach_to_dv<VertexBase>(attachment, defaultValue, passOnValues);}
		template <class TAttachment>
		inline void attach_to_edges_dv(TAttachment& attachment, const typename TAttachment::ValueType& defaultValue, bool passOnValues)		{attach_to_dv<EdgeBase>(attachment, defaultValue, passOnValues);}
		template <class TAttachment>
		inline void attach_to_faces_dv(TAttachment& attachment, const typename TAttachment::ValueType& defaultValue, bool passOnValues)		{attach_to_dv<Face>(attachment, defaultValue, passOnValues);}
		template <class TAttachment>
		inline void attach_to_volumes_dv(TAttachment& attachment, const typename TAttachment::ValueType& defaultValue, bool passOnValues)	{attach_to_dv<Volume>(attachment, defaultValue, passOnValues);}

	////////////////////////////////////////////////////////////////////////
	//	detach
		template <class TGeomObjClass>
		void detach_from(IAttachment& attachment);

		inline void detach_from_vertices(IAttachment& attachment)	{detach_from<VertexBase>(attachment);}
		inline void detach_from_edges(IAttachment& attachment)		{detach_from<EdgeBase>(attachment);}
		inline void detach_from_faces(IAttachment& attachment)		{detach_from<Face>(attachment);}
		inline void detach_from_volumes(IAttachment& attachment)	{detach_from<Volume>(attachment);}

		template <class TGeomObjClass>
		inline bool has_attachment(IAttachment& attachment)			{return get_attachment_pipe<TGeomObjClass>().has_attachment(attachment);}

		inline bool has_vertex_attachment(IAttachment& attachment)	{return has_attachment<VertexBase>(attachment);}
		inline bool has_edge_attachment(IAttachment& attachment)	{return has_attachment<EdgeBase>(attachment);}
		inline bool has_face_attachment(IAttachment& attachment)	{return has_attachment<Face>(attachment);}
		inline bool has_volume_attachment(IAttachment& attachment)	{return has_attachment<Volume>(attachment);}

	////////////////////////////////////////////////////////////////////////
	//	direct attachment access
		template <class TGeomObj>
		uint get_attachment_data_index(TGeomObj* pObj) const;

		template <class TGeomObj, class TAttachment>
		typename TAttachment::ContainerType*
		get_attachment_data_container(TAttachment& attachment);
		
	////////////////////////////////////////////////
	//	observers
		/**
		 * observerType may be any or-combination of constants enumerated in ObserverType.
		 */
		void register_observer(GridObserver* observer, uint observerType = OT_FULL_OBSERVER);
		void unregister_observer(GridObserver* observer);

/*
		template <class GeomObjClass>
		util::IAttachmentDataContainer* get_data_container(util::IAttachment& attachment);

		util::IAttachmentDataContainer* get_vertex_data_container(util::IAttachment& attachment)	{return get_data_container<VertexBase>(attachment);}
		util::IAttachmentDataContainer* get_edge_data_container(util::IAttachment& attachment)		{return get_data_container<EdgeBase>(attachment);}
		util::IAttachmentDataContainer* get_face_data_container(util::IAttachment& attachment)		{return get_data_container<Face>(attachment);}
		util::IAttachmentDataContainer* get_volume_data_container(util::IAttachment& attachment)	{return get_data_container<Volume>(attachment);}
*/

	////////////////////////////////////////////////
	//	pass_on_values
		void pass_on_values(VertexBase* objSrc, VertexBase* objDest);
		void pass_on_values(EdgeBase* objSrc, EdgeBase* objDest);
		void pass_on_values(Face* objSrc, Face* objDest);
		void pass_on_values(Volume* objSrc, Volume* objDest);

	//	subject to change!
		AssociatedEdgeIterator associated_edges_begin(VertexBase* vrt);///< DO NOT INVOKE! Subject to change.
		AssociatedEdgeIterator associated_edges_end(VertexBase* vrt);///< DO NOT INVOKE! Subject to change.
		AssociatedEdgeIterator associated_edges_begin(Face* face);///< DO NOT INVOKE! Subject to change.
		AssociatedEdgeIterator associated_edges_end(Face* face);///< DO NOT INVOKE! Subject to change.
		AssociatedEdgeIterator associated_edges_begin(Volume* vol);///< DO NOT INVOKE! Subject to change.
		AssociatedEdgeIterator associated_edges_end(Volume* vol);///< DO NOT INVOKE! Subject to change.

		AssociatedFaceIterator associated_faces_begin(VertexBase* vrt);///< DO NOT INVOKE! Subject to change.
		AssociatedFaceIterator associated_faces_end(VertexBase* vrt);///< DO NOT INVOKE! Subject to change.
		AssociatedFaceIterator associated_faces_begin(EdgeBase* edge);///< DO NOT INVOKE! Subject to change.
		AssociatedFaceIterator associated_faces_end(EdgeBase* edge);///< DO NOT INVOKE! Subject to change.
		AssociatedFaceIterator associated_faces_begin(Volume* vol);///< DO NOT INVOKE! Subject to change.
		AssociatedFaceIterator associated_faces_end(Volume* vol);///< DO NOT INVOKE! Subject to change.

		AssociatedVolumeIterator associated_volumes_begin(VertexBase* vrt);///< DO NOT INVOKE! Subject to change.
		AssociatedVolumeIterator associated_volumes_end(VertexBase* vrt);///< DO NOT INVOKE! Subject to change.
		AssociatedVolumeIterator associated_volumes_begin(EdgeBase* edge);///< DO NOT INVOKE! Subject to change.
		AssociatedVolumeIterator associated_volumes_end(EdgeBase* edge);///< DO NOT INVOKE! Subject to change.
		AssociatedVolumeIterator associated_volumes_begin(Face* face);///< DO NOT INVOKE! Subject to change.
		AssociatedVolumeIterator associated_volumes_end(Face* face);///< DO NOT INVOKE! Subject to change.

	////////////////////////////////////////////////
	//	advanced element manipulation
		inline void register_element(VertexBase* v, GeometricObject* pParent = NULL)	{register_vertex(v, pParent);}
		inline void unregister_element(VertexBase* v)									{unregister_vertex(v);}
		inline void register_element(EdgeBase* e, GeometricObject* pParent = NULL)		{register_edge(e, pParent);}
		inline void unregister_element(EdgeBase* e)										{unregister_edge(e);}
		inline void register_element(Face* f, GeometricObject* pParent = NULL)			{register_face(f, pParent);}
		inline void unregister_element(Face* f)											{unregister_face(f);}
		inline void register_element(Volume* v, GeometricObject* pParent = NULL)		{register_volume(v, pParent);}
		inline void unregister_element(Volume* v)										{unregister_volume(v);}

	///	registers the given element and replaces the old one. Calls pass_on_values.
	/// \{
		void register_and_replace_element(VertexBase* v, VertexBase* pReplaceMe);
		void register_and_replace_element(EdgeBase* e, EdgeBase* pReplaceMe);
		void register_and_replace_element(Face* f, Face* pReplaceMe);
		void register_and_replace_element(Volume* v, Volume* pReplaceMe);
	/// \}

	////////////////////////////////////////////////
	//	marks
	///	begin marking.
	/**	Call this method whenever you want to start a marking sequence.
	 *	On a call to this method all old marks are deleted.
	 *	When called for the first time, some preparations have to be taken,
	 *	which may consume some time. Successive calls however are very fast.*/
		void begin_marking();
	///	clears all marks
	/**	Calls are only valid between calls to Grid::begin_marking and Grid::end_marking.*/
		void clear_marks();

	///	marks the object. Calls are only valid between calls to Grid::begin_marking and Grid::end_marking.
	/**	Only pass objects that are contained by the grid.*/
		inline void mark(VertexBase* obj);
	///	marks the object. Calls are only valid between calls to Grid::begin_marking and Grid::end_marking.
	/**	Only pass objects that are contained by the grid.*/
		inline void mark(EdgeBase* obj);
	///	marks the object. Calls are only valid between calls to Grid::begin_marking and Grid::end_marking.
	/**	Only pass objects that are contained by the grid.*/
		inline void mark(Face* obj);
	///	marks the object. Calls are only valid between calls to Grid::begin_marking and Grid::end_marking.
	/**	Only pass objects that are contained by the grid.*/
		inline void mark(Volume* obj);

	///	unmarks the object. Calls are only valid between calls to Grid::begin_marking and Grid::end_marking.
	/**	Only pass objects that are contained by the grid.*/
		inline void unmark(VertexBase* obj);
	///	unmarks the object. Calls are only valid between calls to Grid::begin_marking and Grid::end_marking.
	/**	Only pass objects that are contained by the grid.*/
		inline void unmark(EdgeBase* obj);
	///	unmarks the object. Calls are only valid between calls to Grid::begin_marking and Grid::end_marking.
	/**	Only pass objects that are contained by the grid.*/
		inline void unmark(Face* obj);
	///	unmarks the object. Calls are only valid between calls to Grid::begin_marking and Grid::end_marking.
	/**	Only pass objects that are contained by the grid.*/
		inline void unmark(Volume* obj);

	///	returns true if the object is marked, false if not.
	/**	Only pass objects that are contained by the grid.*/
		inline bool is_marked(VertexBase* obj);
	///	returns true if the object is marked, false if not.
	/**	Only pass objects that are contained by the grid.*/
		inline bool is_marked(EdgeBase* obj);
	///	returns true if the object is marked, false if not.
	/**	Only pass objects that are contained by the grid.*/
		inline bool is_marked(Face* obj);
	///	returns true if the object is marked, false if not.
	/**	Only pass objects that are contained by the grid.*/
		inline bool is_marked(Volume* obj);

	///	ends a marking sequence. Call this method when you're done with marking.
		void end_marking();
		
	protected:
	///	copies the contents from the given grid to this grid.
	/**	Make sure that the grid on which this method is called is
	 *	empty before the method is called.*/
		void assign_grid(const Grid& grid);
		
	///	assigns a unique hash value to a Vertex.
	/**	overflow is not handled properly.
	 *	If sombody creates 2^32 elements, the uniquness can no longer be guaranteed.*/
		inline void assign_hash_value(VertexBase* vrt)	{vrt->m_hashValue = m_hashCounter++;}

		void register_vertex(VertexBase* v, GeometricObject* pParent = NULL);///< pDF specifies the element from which v derives its values
		void unregister_vertex(VertexBase* v);
		void register_edge(EdgeBase* e, GeometricObject* pParent = NULL);///< pDF specifies the element from which v derives its values
		void unregister_edge(EdgeBase* e);
		void register_face(Face* f, GeometricObject* pParent = NULL);///< pDF specifies the element from which v derives its values
		void unregister_face(Face* f);
		void register_volume(Volume* v, GeometricObject* pParent = NULL);///< pDF specifies the element from which v derives its values
		void unregister_volume(Volume* v);

		//void pass_on_values(util::AttachmentPipe& attachmentPipe, util::AttachmentPipeElement* elemFrom, util::AttachmentPipeElement* elemTo);
		void change_options(uint optsNew);
	/*
		void enable_vertex_options(uint opts);
		void enable_edge_options(uint opts);
		void enable_face_options(uint opts);
		void enable_volume_options(uint opts);
		void disable_vertex_options(uint opts);
		void disable_edge_options(uint opts);
		void disable_face_options(uint opts);
		void disable_volume_options(uint opts);
	*/
		void change_vertex_options(uint optsNew);
		void change_edge_options(uint optsNew);
		void change_face_options(uint optsNew);
		void change_volume_options(uint optsNew);

		void vertex_store_associated_edges(bool bStoreIt);
		void vertex_store_associated_faces(bool bStoreIt);
		void vertex_store_associated_volumes(bool bStoreIt);
		void edge_store_associated_faces(bool bStoreIt);
		void edge_store_associated_volumes(bool bStoreIt);
		void face_store_associated_edges(bool bStoreIt);
		void face_store_associated_volumes(bool bStoreIt);
		void face_autogenerate_edges(bool bAutogen);
		void volume_store_associated_edges(bool bStoreIt);
		void volume_store_associated_faces(bool bStoreIt);
		void volume_autogenerate_edges(bool bAutogen);
		void volume_autogenerate_faces(bool bAutogen);

		void pass_on_values(AttachmentPipe& attachmentPipe,
							GeometricObject* pSrc, GeometricObject* pDest);

		template <class TGeomObj>
		ug::AttachmentPipe<GeometricObject*, Grid>&
		get_attachment_pipe();

	//	some methods that simplify auto-enabling of grid options
		inline void autoenable_option(uint option, const char* caller, const char* optionName);

	//	neighbourhood access
		template <class TGeomObj>
		EdgeBase* find_edge_in_associated_edges(TGeomObj* obj,
												EdgeVertices& ev);
												
		template <class TGeomObj>
		Face* find_face_in_associated_faces(TGeomObj* obj,
											FaceVertices& fv);
												
		template <class TGeomObj>
		Volume* find_volume_in_associated_volumes(TGeomObj* obj,
												VolumeVertices& vv);
												
	//	marks
		void init_marks();
		void reset_marks();
		void remove_marks();
												
	protected:
		typedef ug::SectionContainer<GeometricObject*, std::list<GeometricObject*> > SectionContainer;

	///	This struct is used to hold GeometricObjects and their attachment pipes.
		struct ElementStorage
		{
			SectionContainer	m_sectionContainer;///	holds elements
			AttachmentPipe		m_attachmentPipe;///	holds the data of the stored elements.
		};

	protected:
	//	typedefs
		typedef std::vector<GridObserver*>	ObserverContainer;

		typedef Attachment<VertexContainer>	AVertexContainer;
		typedef Attachment<EdgeContainer>	AEdgeContainer;
		typedef Attachment<FaceContainer>	AFaceContainer;
		typedef Attachment<VolumeContainer>	AVolumeContainer;
		
		typedef Attachment<int>	AMark;

	protected:
		ElementStorage	m_elementStorage[NUM_GEOMETRIC_BASE_OBJECTS];
		uint			m_options;
		uint32			m_hashCounter;

	//	observer handling
		ObserverContainer	m_gridObservers;
		ObserverContainer	m_vertexObservers;
		ObserverContainer	m_edgeObservers;
		ObserverContainer	m_faceObservers;
		ObserverContainer	m_volumeObservers;

	//	interconnection managment
		AVertexContainer	m_aVertexContainer;
		AEdgeContainer		m_aEdgeContainer;
		AFaceContainer		m_aFaceContainer;
		AVolumeContainer	m_aVolumeContainer;

		AttachmentAccessor<VertexBase, AEdgeContainer>		m_aaEdgeContainerVERTEX;
		AttachmentAccessor<VertexBase, AFaceContainer>		m_aaFaceContainerVERTEX;
		AttachmentAccessor<VertexBase, AVolumeContainer>	m_aaVolumeContainerVERTEX;

		AttachmentAccessor<EdgeBase, AEdgeContainer>		m_aaEdgeContainerEDGE;
		AttachmentAccessor<EdgeBase, AFaceContainer>		m_aaFaceContainerEDGE;
		AttachmentAccessor<EdgeBase, AVolumeContainer>		m_aaVolumeContainerEDGE;

		AttachmentAccessor<Face, AEdgeContainer>		m_aaEdgeContainerFACE;
		AttachmentAccessor<Face, AFaceContainer>		m_aaFaceContainerFACE;
		AttachmentAccessor<Face, AVolumeContainer>		m_aaVolumeContainerFACE;

		AttachmentAccessor<Volume, AEdgeContainer>		m_aaEdgeContainerVOLUME;
		AttachmentAccessor<Volume, AFaceContainer>		m_aaFaceContainerVOLUME;
		AttachmentAccessor<Volume, AVolumeContainer>	m_aaVolumeContainerVOLUME;
		
	//	marks
		int m_currentMark;	// 0: marks inactive. -1: reset-marks (sets currentMark to 1)
		bool m_bMarking;
		AMark	m_aMark;
		VertexAttachmentAccessor<AMark>	m_aaMarkVRT;
		EdgeAttachmentAccessor<AMark>	m_aaMarkEDGE;
		FaceAttachmentAccessor<AMark>	m_aaMarkFACE;
		VolumeAttachmentAccessor<AMark>	m_aaMarkVOL;
};

///	@}
}//end of namespace

#include "grid_impl.hpp"

#endif
