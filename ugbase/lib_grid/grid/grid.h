//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y08 m10 d09

#ifndef __H__LIB_GRID__GRID__
#define __H__LIB_GRID__GRID__

#include <vector>
#include <list>
#include <memory>
#include <boost/function.hpp>
#include "common/util/message_hub.h"
#include "common/ug_config.h"
#include "grid_constants.h"
#include "geometric_base_objects.h"
#include "grid_observer.h"
#include "geometric_object_collection.h"
#include "element_storage.h"
#include "geometric_base_object_traits.h"

//	Define PROFILE_GRID to profile some often used gird-methods
//#define PROFILE_GRID
#ifdef PROFILE_GRID
	#include "common/profiler/profiler.h"
	#define GRID_PROFILE_FUNC()	PROFILE_FUNC()
	#define GRID_PROFILE(name)	PROFILE_BEGIN(name)
	#define GRID_PROFILE_END()	PROFILE_END()
#else
	#define GRID_PROFILE_FUNC()
	#define GRID_PROFILE(name)
	#define GRID_PROFILE_END()
#endif

namespace ug
{

//	predeclaration of the distributed grid manager, which handles parallelization
//	of grids. If you want to use it, you have to include
//	"lib_grid/parallelization/distributed_grid.h"
class DistributedGridManager;


/**
 * \brief Grid, MultiGrid and GeometricObjectCollection are contained in this group
 * \defgroup lib_grid_grid grid
 * \ingroup lib_grid
 * \{ */
 
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
 *
 * The Grid class features a message-hub, which can be used to distribute messages
 * to registered callbacks.
 */
class UG_API Grid
{
	public:
	///	The traits class holds some important types for each element-type
		template <class TElem>
		struct traits{
			typedef typename TElem::geometric_base_object			base_object;
			typedef ug::ElementStorage<base_object>					ElementStorage;
			typedef typename ElementStorage::AttachmentPipe			AttachmentPipe;
			typedef typename ElementStorage::AttachedElementList	AttachedElementList;
			typedef typename ElementStorage::SectionContainer		SectionContainer;

			typedef typename geometry_traits<TElem>::iterator		iterator;
			typedef typename geometry_traits<TElem>::const_iterator	const_iterator;

			typedef PointerConstArray<TElem*>						secure_container;

		//	CALLBACKS
		///	callback type for the elements base type.
			typedef boost::function<bool (base_object*)>			callback;

		///	convenience method which can be used as a callback that always returns true
			static bool cb_consider_all(base_object*)				{return true;}
		///	convenience method which can be used as a callback that always returns false
			static bool cb_consider_none(base_object*)				{return false;}
		};

	///	Convenience access to grid elements
	/** \{ */
		typedef traits<VertexBase>	vertex_traits;
		typedef traits<EdgeBase>	edge_traits;
		typedef traits<Face>		face_traits;
		typedef traits<Volume>		volume_traits;
	/** \} */

	///	Container to store associated vertices.
		typedef traits<VertexBase>::secure_container	SecureVertexContainer;
	///	Container to store associated edges.
		typedef traits<EdgeBase>::secure_container		SecureEdgeContainer;
	///	Container to store associated faces.
		typedef traits<Face>::secure_container			SecureFaceContainer;
	///	Container to store associated volumes.
		typedef traits<Volume>::secure_container		SecureVolumeContainer;

	///	the attachment-pipe used by Grid
		typedef ug::AttachmentPipe<VertexBase*, VertexElementStorage>	VertexAttachmentPipe;
		typedef ug::AttachmentPipe<EdgeBase*, EdgeElementStorage>		EdgeAttachmentPipe;
		typedef ug::AttachmentPipe<Face*, FaceElementStorage>			FaceAttachmentPipe;
		typedef ug::AttachmentPipe<Volume*, VolumeElementStorage>		VolumeAttachmentPipe;


	///	the generic attachment-accessor for access to grids attachment pipes.
		template <class TElem, class TAttachment>
		class AttachmentAccessor : public ug::AttachmentAccessor<typename TElem::geometric_base_object*,
																 TAttachment,
																 typename traits<TElem>::ElementStorage>
		{
			public:
				AttachmentAccessor();
				AttachmentAccessor(const AttachmentAccessor& aa);
				AttachmentAccessor(Grid& grid, TAttachment& a);
				AttachmentAccessor(Grid& grid, TAttachment& a, bool autoAttach);

				inline bool access(Grid& grid, TAttachment& a)
					{
						return ug::AttachmentAccessor<typename TElem::geometric_base_object*,
													  TAttachment,
													  typename traits<TElem>::ElementStorage>::
							access(grid.get_attachment_pipe<TElem>(), a);
					}
		};

	//	half-specialized AttachmentAccessors:
		template <class TAttachment>
		class VertexAttachmentAccessor : public AttachmentAccessor<VertexBase, TAttachment>
		{
			public:
				VertexAttachmentAccessor();
				VertexAttachmentAccessor(const VertexAttachmentAccessor& aa);
				VertexAttachmentAccessor(Grid& grid, TAttachment& a);
				VertexAttachmentAccessor(Grid& grid, TAttachment& a, bool autoAttach);
		};

		template <class TAttachment>
		class EdgeAttachmentAccessor : public AttachmentAccessor<EdgeBase, TAttachment>
		{
			public:
				EdgeAttachmentAccessor();
				EdgeAttachmentAccessor(const EdgeAttachmentAccessor& aa);
				EdgeAttachmentAccessor(Grid& grid, TAttachment& a);
				EdgeAttachmentAccessor(Grid& grid, TAttachment& a, bool autoAttach);
		};

		template <class TAttachment>
		class FaceAttachmentAccessor : public AttachmentAccessor<Face, TAttachment>
		{
			public:
				FaceAttachmentAccessor();
				FaceAttachmentAccessor(const FaceAttachmentAccessor& aa);
				FaceAttachmentAccessor(Grid& grid, TAttachment& a);
				FaceAttachmentAccessor(Grid& grid, TAttachment& a, bool autoAttach);
		};

		template <class TAttachment>
		class VolumeAttachmentAccessor : public AttachmentAccessor<Volume, TAttachment>
		{
			public:
				VolumeAttachmentAccessor();
				VolumeAttachmentAccessor(const VolumeAttachmentAccessor& aa);
				VolumeAttachmentAccessor(Grid& grid, TAttachment& a);
				VolumeAttachmentAccessor(Grid& grid, TAttachment& a, bool autoAttach);
		};

	///	Container used to store associated vertices
		typedef std::vector<VertexBase*>	VertexContainer;
	///	Container used to store associated edges
		typedef std::vector<EdgeBase*>		EdgeContainer;
	///	Container used to store associated faces
		typedef std::vector<Face*>			FaceContainer;
	///	Container used to store associated volumes
		typedef std::vector<Volume*>		VolumeContainer;

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

	////////////////////////////////////////////////
	//	parallelism
	///	tell the grid whether it will be used in a serial or in a parallel environment.
	/**	If parallelism is enabled, the grid will internally create a distributed grid
	 * manager, which will handle horizontal and vertical process interfaces. The
	 * manager can be queried through the method ug::Grid::distributed_grid_manager.
	 * If false is passed and a distributed grid manager already existed, it will
	 * be destroyed.
	 * parallelism is disabled by default.
	 * \note	set_parallel(true) may currently only be executed on ug::MultiGrid.
	 * 			There is currently no parallelization support for plain ug::Grid.
	 * \note	parallelism may only be activated if ug was compiled with PARALLEL=ON.*/
		void set_parallel(bool parallel);

	///	returns true if the grid is prepared for parallel computations.
	/**	If the method returns true, it is also clear, that a distributed grid manager
	 * exists in the grid. The manager can be queried through ug::Grid::distributed_grid_manager.*/
		inline bool is_parallel() const;

	///	returns a pointer to the associated distributed grid manager.
	/**	The method returns NULL, if no distributed grid manager for the given
	 * grid exists. This should be the case for serial environments or for
	 * serial grids in parallel environments.
	 * Use ug::Grid::set_parallel() to enable or disable parallelism in a grid.
	 * You have to include "lib_grid/parallelization/distributed_grid.h" to access
	 * methods of a distributed grid manager.
	 * \{ */
		inline DistributedGridManager* distributed_grid_manager();
		inline const DistributedGridManager* distributed_grid_manager() const;
	/** \} */

	////////////////////////////////////////////////
	//	clear
	///	clears the grids geometry and attachments
		void clear();
	///	clears the grids geometry. Registered attachments remain.
		void clear_geometry();
	///	clears the grids attachments. The geometry remains.
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

	///	Reserves memory for the creation of the given object type
	/**	Calls to this method are optional, but can improve runtime.
	 * Specify the total number of objects which the grid should be
	 * capable to hold (if more are required, the grid will automatically
	 * adjust sizes)
	 */
		template <class TGeomObj>
		void reserve(size_t num);

	////////////////////////////////////////////////
	//	element deletion
		void erase(GeometricObject* geomObj);
		void erase(VertexBase* vrt);
		void erase(EdgeBase* edge);
		void erase(Face* face);
		void erase(Volume* vol);

	/**\todo: This erase method can cause problems if used with multi-grids.*/
		template <class GeomObjIter>
		void erase(const GeomObjIter& iterBegin, const GeomObjIter& iterEnd);

		template <class TGeomObj>
		void clear();

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

	///	notifies the grid that two objects will be merged.
	/**	The grid forwards this notification to its GridObservers.
	 * The notification is not relevant for the grid itself.
	 * \{ */
		void objects_will_be_merged(VertexBase* target, VertexBase* elem1,
									VertexBase* elem2);
		void objects_will_be_merged(EdgeBase* target, EdgeBase* elem1,
									EdgeBase* elem2);
		void objects_will_be_merged(Face* target, Face* elem1,
									Face* elem2);
		void objects_will_be_merged(Volume* target, Volume* elem1,
									Volume* elem2);
	/**	\} */
	/*
	///	VrtPairIterator has to be an iterator with value-type std::pair<VertexBase*, VertexBase*>
		template <class VrtPairIter>
		void replace_vertices(VrtPairIter& iterBegin, VrtPairIter& iterEnd);
	*/

	////////////////////////////////////////////////
	//	geometric-object-collection
	///	returns the the GeometricObjectCollection of the grid:
		virtual GeometricObjectCollection get_geometric_objects();

	////////////////////////////////////////////////
	///	flips the orientation of an edge.
		void flip_orientation(EdgeBase* e);
		
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

	///	returns the first element of the given type.
	/**	Make sure that elements of the given type exist!
	 *	Behaviour is undefined, if not.*/
		template <class TGeomObj> TGeomObj* front();
		
	///	returns the last element of the given type.
	/**	Make sure that elements of the given type exist!
	 *	Behaviour is undefined, if not.*/
		template <class TGeomObj> TGeomObj* back();
		
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
	 *	should be enabled.
	 *	\note	If all edges of a face have to be processed, it may be faster
	 *			to get all edges in one single call.
	 *			This can be done using Grid::associated_elements.*/
		EdgeBase* get_edge(Face* f, int ind);
	///	If it exists, this method returns the i-th edge of the given Volume. If not NULL is returned.
	/**	To make sure that associated edges always exist, enable the grid-option
	 *	VOLOPT_AUTOGENERATE_EDGES.
	 *	For maximal performance, the option VOLOPT_STORE_ASSOCIATED_EDGES
	 *	should be enabled.
	 *	\note	If all edges of a volume have to be processed, it may be faster
	 *			to get all edges in one single call.
	 *			This can be done using Grid::associated_elements.*/
		EdgeBase* get_edge(Volume* v, int ind);
	///	returns the face that is described by fv.
	/**	Note that you may pass a FaceDescriptor to this method.*/
		Face* get_face(FaceVertices& fv);
	///	If it exists, this method returns the i-th face of the given Volume. If not NULL is returned.
	/**	To make sure that associated faces always exist, enable the grid-option
	 *	VOLOPT_AUTOGENERATE_FACES.
	 *	For maximal performance, the option VOLOPT_STORE_ASSOCIATED_FACES
	 *	should be enabled.
	 *	\note	If all faces of a volume have to be processed, it may be faster
	 *			to get all faces in one single call.
	 *			This can be done using Grid::associated_elements.*/
		Face* get_face(Volume* v, int ind);
	///	returns the volume that is described by ev.
	/**	Note that you may pass an VolumeDescriptor to this method.*/
		Volume* get_volume(VolumeVertices& vv);
		
	///	returns the element for the given vertices.
	/**	Note that you can either pass an element type (EdgeBase, Face, Volume)
	 * or a descriptor (EdgeDescriptor, FaceDescriptor, VolumeDescriptor).
	 * The method returns NULL, if the specified element does not exist.
	 * A special overload exists for VertexBase*, which simply returns the
	 * specified vertex. Useful for template programming...
	 * \{ */
		EdgeBase* get_element(EdgeVertices& ev)	{return get_edge(ev);}
		Face* get_element(FaceVertices& fv)		{return get_face(fv);}
		Volume* get_element(VolumeVertices& vv)	{return get_volume(vv);}
	/**	\} */

	////////////////////////////////////////////////
	//	access to the sides of an geometric object
	///	This method returns the i-th side of an EdgeBase, Face or Volume.
	/**	If obj has dimension d, then all associated elements of dimension d-1
	 *	are regarded as sides. (Face -> EdgeBase). Only derivates of Volume,
	 *	Face or EdgeBase may be queried for their sides. If you call this method
	 *	with VertexBase*, an assertion is triggered, since vertices do not have sides.
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
	 *
	 *	\note	If all sides of an element have to be processed, it may be faster
	 *			to get all sides in one single call.
	 *			This can be done using Grid::associated_elements.
	 */
		VertexBase::side* get_side(VertexBase* obj, size_t side);
		EdgeBase::side* get_side(EdgeBase* obj, size_t side);
		Face::side* get_side(Face* obj, size_t side);
		Volume::side* get_side(Volume* obj, size_t side);

	///	returns the geometric object on the opposing side of the given vertex regarding the given element.
	/**	\note Currently only implemented for Face and Volume.
	 * \{ */
		GeometricObject* get_opposing_object(VertexBase* vrt, Face* elem);
		GeometricObject* get_opposing_object(VertexBase* vrt, Volume* elem);
	/** \} */


	///	Puts all elements of type TAss which are contained in 'e' or which contain 'e' into elemsOut
	/**	One shouldn't depend on the order of elements in elemsOut.
	 * Use Grid::associated_elements_sorted if the order matters.
	 *
	 * \note	The returned container is only valid until changes to the queried
	 * 			grid are performed.
	 *
	 * \note	The following typedefs might ease your work with associated elements:
	 * 			AssociatedVertices, AssociatedEdges, AssociatedFaces and
	 * 			AssociatedVolumes.
	 * 			Those are typedefs for Grid::traits<VertexBase>::container, ...
	 *
	 * \note	Depending on the current grid options, this method may use Grid::mark.
	 * Valid arguments for TElem are VertexBase, EdgeBase, Face, Volume.
	 * \sa Grid::associated_elements_sorted
	 * \{ */
		template <class TElem>
		void associated_elements(traits<VertexBase>::secure_container& elemsOut, TElem* e);
		template <class TElem>
		void associated_elements(traits<EdgeBase>::secure_container& elemsOut, TElem* e);
		template <class TElem>
		void associated_elements(traits<Face>::secure_container& elemsOut, TElem* e);
		template <class TElem>
		void associated_elements(traits<Volume>::secure_container& elemsOut, TElem* e);
	/** \} */
	
	///	Puts all elements of type TAss which are contained in 'e' into elemsOut in the reference elements order
	/**	The order of elements in elemsOut is the same as the order in the reference element.
	 * Note that, depending on active grid options, this method may take more
	 * time than associated_elements, and should thus only be invoked, if the order
	 * of elements matters. Use Grid::associated_elements if the order does not matter.
	 *
	 * Valid arguments for TElem are VertexBase, EdgeBase, Face, Volume.
	 * Let TAss be the type of queried elements in elemsOut.
	 * Only associated elements of lower dimension can be sorted. The method
	 * thus behaves as follows:
	 * 	TAss::dim < TElem::dim	->	associated elements are written to 'elemsOut'
	 * 								with the same order as in the reference element
	 * 	TAss::dim == TElem::dim	->	onle 'e' itself will be written to 'elemsOut'
	 * 	TAss::dim > TElem::dim	->	elemsOut will be empty.
	 *
	 * \note	The returned container is only valid until changes to the queried
	 * 			grid are performed.
	 *
	 * \note	The following typedefs might ease your work with associated elements:
	 * 			AssociatedVertices, AssociatedEdges, AssociatedFaces and
	 * 			AssociatedVolumes.
	 * 			Those are typedefs for traits<VertexBase>::container, ...
	 *
	 * \note	Depending on the current grid options, this method may use Grid::mark.
	 * \sa Grid::associated_elements*/
		template <class TElem>
		void associated_elements_sorted(traits<VertexBase>::secure_container& elemsOut, TElem* e);
		template <class TElem>
		void associated_elements_sorted(traits<EdgeBase>::secure_container& elemsOut, TElem* e);
		template <class TElem>
		void associated_elements_sorted(traits<Face>::secure_container& elemsOut, TElem* e);
		template <class TElem>
		void associated_elements_sorted(traits<Volume>::secure_container& elemsOut, TElem* e);


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

		inline void attach_to_all(IAttachment& attachment, bool passOnValues);

	////////////////////////////////////////////////////////////////////////
	///	attach with default pass-on behaviour and unspecified default value.
		template <class TGeomObjClass>
		inline void attach_to(IAttachment& attachment)	{attach_to<TGeomObjClass>(attachment, attachment.default_pass_on_behaviour());}

		inline void attach_to_vertices(IAttachment& attachment)	{attach_to<VertexBase>(attachment);}
		inline void attach_to_edges(IAttachment& attachment)	{attach_to<EdgeBase>(attachment);}
		inline void attach_to_faces(IAttachment& attachment)	{attach_to<Face>(attachment);}
		inline void attach_to_volumes(IAttachment& attachment)	{attach_to<Volume>(attachment);}

	///	attaches to vertices, edges, faces and volumes at once.
		inline void attach_to_all(IAttachment& attachment);

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

	///	attaches to vertices, edges, faces and volumes at once.
		template <class TAttachment>
		inline void attach_to_all_dv(TAttachment& attachment, const typename TAttachment::ValueType& defaultValue);

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

		template <class TAttachment>
		inline void attach_to_all_dv(TAttachment& attachment, const typename TAttachment::ValueType& defaultValue, bool passOnValues);
	////////////////////////////////////////////////////////////////////////
	//	detach
		template <class TGeomObjClass>
		void detach_from(IAttachment& attachment);

		inline void detach_from_vertices(IAttachment& attachment)	{detach_from<VertexBase>(attachment);}
		inline void detach_from_edges(IAttachment& attachment)		{detach_from<EdgeBase>(attachment);}
		inline void detach_from_faces(IAttachment& attachment)		{detach_from<Face>(attachment);}
		inline void detach_from_volumes(IAttachment& attachment)	{detach_from<Volume>(attachment);}

		inline void detach_from_all(IAttachment& attachment);


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
		
	///	returns the attachment-pipe in which data associated with the given objects-types are stored.
	/**	This method is seldomly used, can however be useful. Use with care.
	 * If in doubt please use the methods featured by Grid instead of directly
	 * operating on the attachment pipe.
	 */
		template <class TGeomObj>
		typename traits<TGeomObj>::AttachmentPipe&
		get_attachment_pipe();

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
	/**	Only pass objects that are contained by the grid.
	 * \{ */
		inline void mark(GeometricObject* obj);
		inline void mark(VertexBase* obj);
		inline void mark(EdgeBase* obj);
		inline void mark(Face* obj);
		inline void mark(Volume* obj);
	/**	\} */

	///	marks all objects between begin and end
	/**	TIterator::value_type has to be either
	 *	VertexBase*, EdgeBase*, Face* or Volume*.*/
		template <class TIterator>
		void mark(TIterator begin, TIterator end);
		
	///	unmarks the object. Calls are only valid between calls to Grid::begin_marking and Grid::end_marking.
	/**	Only pass objects that are contained by the grid.
	 * \{ */
		inline void unmark(GeometricObject* obj);
		inline void unmark(VertexBase* obj);
		inline void unmark(EdgeBase* obj);
		inline void unmark(Face* obj);
		inline void unmark(Volume* obj);
	/** \} */

	///	unmarks all objects between begin and end
	/**	TIterator::value_type has to be either
	 *	VertexBase*, EdgeBase*, Face* or Volume*.*/
		template <class TIterator>
		void unmark(TIterator begin, TIterator end);
		
	///	returns true if the object is marked, false if not.
	/**	Only pass objects that are contained by the grid.
	 * \{ */
		inline bool is_marked(GeometricObject* obj);
		inline bool is_marked(VertexBase* obj);
		inline bool is_marked(EdgeBase* obj);
		inline bool is_marked(Face* obj);
		inline bool is_marked(Volume* obj);
	/** \} */

	///	ends a marking sequence. Call this method when you're done with marking.
		void end_marking();
		
	///	gives access to the grid's message-hub
		SPMessageHub message_hub()		{return m_messageHub;}

	///	a temporary testing method
		void test_attached_linked_lists();

	protected:
		typedef std::vector<GridObserver*>	ObserverContainer;

		typedef Attachment<VertexContainer>	AVertexContainer;
		typedef Attachment<EdgeContainer>	AEdgeContainer;
		typedef Attachment<FaceContainer>	AFaceContainer;
		typedef Attachment<VolumeContainer>	AVolumeContainer;

		typedef Attachment<int>	AMark;

	protected:
	///	unregisters all observers. Call this method in destructors of derived classes.
	/**	If the derived class is an observer itself and if you don't want it to be
	 * notified on grid-destruction, e.g., because you call this method in the
	 * destructor of your derived class, then pass a pointer to your class through
	 * the initiator parameter to this function.
	 * \param initiator:	The initiator won't be notified about grid destruction
	 */
		void notify_and_clear_observers_on_grid_destruction(GridObserver* initiator = NULL);

	///	returns the element storage for a given element type
		template <class TElem> inline
		typename traits<TElem>::ElementStorage&
		element_storage()
		{return ElementStorageSelector<typename geometry_traits<TElem>::geometric_base_object>::
				element_storage(m_vertexElementStorage, m_edgeElementStorage,
								m_faceElementStorage, m_volumeElementStorage);
		}

	///	returns the const element storage for a given element type
		template <class TElem> inline
		const typename traits<TElem>::ElementStorage&
		element_storage() const
		{return ElementStorageSelector<typename geometry_traits<TElem>::geometric_base_object>::
				element_storage(m_vertexElementStorage, m_edgeElementStorage,
								m_faceElementStorage, m_volumeElementStorage);
		}

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
		void register_edge(EdgeBase* e, GeometricObject* pParent = NULL,
						Face* createdByFace = NULL, Volume* createdByVol = NULL);///< pDF specifies the element from which v derives its values
		void unregister_edge(EdgeBase* e);
		void register_face(Face* f, GeometricObject* pParent = NULL,
						   Volume* createdByVol = NULL);///< pDF specifies the element from which v derives its values
		void unregister_face(Face* f);
		void register_volume(Volume* v, GeometricObject* pParent = NULL);///< pDF specifies the element from which v derives its values
		void unregister_volume(Volume* v);

		void change_options(uint optsNew);

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

		void volume_sort_associated_edge_container();

		template <class TAttachmentPipe, class TElem>
		void pass_on_values(TAttachmentPipe& attachmentPipe,
							TElem* pSrc, TElem* pDest);

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

	//	get associated elements
		void get_associated(SecureVertexContainer& vrts, VertexBase* v);
		void get_associated(SecureVertexContainer& vrts, EdgeBase* e);
		void get_associated(SecureVertexContainer& vrts, Face* f);
		void get_associated(SecureVertexContainer& vrts, Volume* v);

		void get_associated(SecureEdgeContainer& edges, VertexBase* v);
		void get_associated(SecureEdgeContainer& edges, EdgeBase* e);
		void get_associated(SecureEdgeContainer& edges, Face* f);
		void get_associated(SecureEdgeContainer& edges, Volume* v);

		void get_associated(SecureFaceContainer& faces, VertexBase* v);
		void get_associated(SecureFaceContainer& faces, EdgeBase* e);
		void get_associated(SecureFaceContainer& faces, Face* f);
		void get_associated(SecureFaceContainer& faces, Volume* v);

		void get_associated(SecureVolumeContainer& vols, VertexBase* v);
		void get_associated(SecureVolumeContainer& vols, EdgeBase* e);
		void get_associated(SecureVolumeContainer& vols, Face* f);
		void get_associated(SecureVolumeContainer& vols, Volume* v);
		
	/**	this method does not use possibly attached containers and can thus
	 * be used, when such containers are to be built.*/
		void get_associated_vols_raw(SecureVolumeContainer& vols, Face* f);

		void get_associated_sorted(SecureVertexContainer& vrts, EdgeBase* e) const;
		void get_associated_sorted(SecureVertexContainer& vrts, Face* f) const;
		void get_associated_sorted(SecureVertexContainer& vrts, Volume* v) const;

		void get_associated_sorted(SecureEdgeContainer& edges, VertexBase* v);
		void get_associated_sorted(SecureEdgeContainer& edges, Face* f);
		void get_associated_sorted(SecureEdgeContainer& edges, Volume* v);

		void get_associated_sorted(SecureFaceContainer& faces, VertexBase* v);
		void get_associated_sorted(SecureFaceContainer& faces, EdgeBase* e);
		void get_associated_sorted(SecureFaceContainer& faces, Volume* v);
		
		void get_associated_sorted(SecureVolumeContainer& vols, VertexBase* v);
		void get_associated_sorted(SecureVolumeContainer& vols, EdgeBase* e);
		void get_associated_sorted(SecureVolumeContainer& vols, Face* f);

		template <class TElem>
		void get_associated_sorted(typename traits<TElem>::secure_container& elems, TElem* e);


	///	helps in copying attachment pipes during assign_grid
	/**	Note that this method only copies attachments with m_userData==1.
	 * \todo	Copy behavior should be changed to all user-attachments.*/
		template <class TAttachmentPipe>
		void copy_user_attachments(const TAttachmentPipe& apSrc, TAttachmentPipe& apDest,
									std::vector<int>& srcDataIndices);

	//	marks
		void init_marks();
		void reset_marks();
		void remove_marks();

	protected:
	///	returns the iterator at which the given element lies in the section container
	/**	This method may only be called if the element has already been registered at the grid.
	 * \{
	 */
		inline traits<VertexBase>::SectionContainer::iterator
		get_iterator(VertexBase* o)
		{
			return m_vertexElementStorage.m_sectionContainer.
					get_container().get_iterator(o);
		}

		inline traits<EdgeBase>::SectionContainer::iterator
		get_iterator(EdgeBase* o)
		{
			return m_edgeElementStorage.m_sectionContainer.
					get_container().get_iterator(o);
		}

		inline traits<Face>::SectionContainer::iterator
		get_iterator(Face* o)
		{
			return m_faceElementStorage.m_sectionContainer.
					get_container().get_iterator(o);
		}

		inline traits<Volume>::SectionContainer::iterator
		get_iterator(Volume* o)
		{
			return m_volumeElementStorage.m_sectionContainer.
					get_container().get_iterator(o);
		}
	/**	\}	*/


	///	helper to clear_attachments
		template <class TElem>
		void clear_attachments();

	protected:
		VertexElementStorage	m_vertexElementStorage;
		EdgeElementStorage		m_edgeElementStorage;
		FaceElementStorage		m_faceElementStorage;
		VolumeElementStorage	m_volumeElementStorage;

		uint			m_options;
		uint32			m_hashCounter;

	//	observer handling
		ObserverContainer	m_gridObservers;
		ObserverContainer	m_vertexObservers;
		ObserverContainer	m_edgeObservers;
		ObserverContainer	m_faceObservers;
		ObserverContainer	m_volumeObservers;

	//	interconnection management
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

		SPMessageHub 							m_messageHub;
		std::auto_ptr<DistributedGridManager>	m_distGridMgr;
};

/** \} */
}//end of namespace

#include "grid_impl.hpp"

#endif
