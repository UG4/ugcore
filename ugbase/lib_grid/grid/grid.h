/*
 * Copyright (c) 2009-2015:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#ifndef __H__LIB_GRID__GRID__
#define __H__LIB_GRID__GRID__

#include <vector>
#include <list>
#include <memory>

#include <boost/function.hpp>

#include "common/util/message_hub.h"
#include "common/ug_config.h"
#include "common/util/pointer_const_array.h"
//#include "grid_constants.h"
#include "grid_base_objects.h"
#include "grid_observer.h"
#include "grid_object_collection.h"
#include "element_storage.h"
#include "grid_base_object_traits.h"

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

namespace ug {

//	predeclaration of the distributed grid manager, which handles parallelization
//	of grids. If you want to use it, you have to include
//	"lib_grid/parallelization/distributed_grid.h"
class DistributedGridManager;

//	predeclaration of the periodic boundary identifier, which handles the master/slave
//  identification of grid elements. If you want to use it, you have to include
//	"lib_grid/tools/periodic_boundary_identifier.h"
class PeriodicBoundaryManager;

/**
 * \brief Grid, MultiGrid and GridObjectCollection are contained in this group
 * \defgroup lib_grid_grid grid
 * \ingroup lib_grid
 * \{ */
 
////////////////////////////////////////////////////////////////////////////////////////////////
//	Grid
///	Manages the elements of a grid and their interconnection.
/**
 * The Grid class is the heart of libGrid. It can be used to create elements of
 * custom types, for example Vertices, Triangles or Tetrahedrons.
 * All elements have to be derived from either Vertex, Edge, Face or Volume
 * and have to specialize the template-class geometry_traits
 * (both defined in grid_base_objects.h).
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
 * and Volume (classes Vertex, Edge, Face, Volume).
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
 * base objects should be passed (Vertex, Edge, Face and Volume).
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
		template <typename TElem>
		struct traits{
			using base_object = typename TElem::grid_base_object;
			using ElementStorage = ElementStorage<base_object>; // ø using A = A<T>
			using AttachmentPipe = typename ElementStorage::AttachmentPipe;
			using AttachedElementList = typename ElementStorage::AttachedElementList;
			using SectionContainer = typename ElementStorage::SectionContainer;

			using iterator = typename geometry_traits<TElem>::iterator;
			using const_iterator = typename geometry_traits<TElem>::const_iterator;

			using secure_container = PointerConstArray<TElem*>;

		//	CALLBACKS
		///	callback type for the elements base type.
			using callback = boost::function<bool (base_object*)>;
		};

	///	Convenience access to grid elements
	/** \{ */
		using vertex_traits = traits<Vertex>;
		using edge_traits = traits<Edge>;
		using face_traits = traits<Face>;
		using volume_traits = traits<Volume>;
	/** \} */

	///	Container to store associated vertices.
		using SecureVertexContainer = traits<Vertex>::secure_container;
	///	Container to store associated edges.
		using SecureEdgeContainer = traits<Edge>::secure_container;
	///	Container to store associated faces.
		using SecureFaceContainer = traits<Face>::secure_container;
	///	Container to store associated volumes.
		using SecureVolumeContainer = traits<Volume>::secure_container;

	///	the attachment-pipe used by Grid
		using VertexAttachmentPipe = AttachmentPipe<Vertex*, VertexElementStorage>;
		using EdgeAttachmentPipe = AttachmentPipe<Edge*, EdgeElementStorage>;
		using FaceAttachmentPipe = AttachmentPipe<Face*, FaceElementStorage>;
		using VolumeAttachmentPipe = AttachmentPipe<Volume*, VolumeElementStorage>;


	///	the generic attachment-accessor for access to grids attachment pipes.
		template <typename TElem, typename TAttachment>
		class AttachmentAccessor : public ug::AttachmentAccessor<typename TElem::grid_base_object*,
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
						return ug::AttachmentAccessor<typename TElem::grid_base_object*,
													  TAttachment,
													  typename traits<TElem>::ElementStorage>::
							access(grid.get_attachment_pipe<TElem>(), a);
					}
		};

	//	half-specialized AttachmentAccessors:
		template <typename TAttachment>
		class VertexAttachmentAccessor : public AttachmentAccessor<Vertex, TAttachment>
		{
			public:
				VertexAttachmentAccessor();
				VertexAttachmentAccessor(const VertexAttachmentAccessor& aa);
				VertexAttachmentAccessor(Grid& grid, TAttachment& a);
				VertexAttachmentAccessor(Grid& grid, TAttachment& a, bool autoAttach);
		};

		template <typename TAttachment>
		class EdgeAttachmentAccessor : public AttachmentAccessor<Edge, TAttachment>
		{
			public:
				EdgeAttachmentAccessor();
				EdgeAttachmentAccessor(const EdgeAttachmentAccessor& aa);
				EdgeAttachmentAccessor(Grid& grid, TAttachment& a);
				EdgeAttachmentAccessor(Grid& grid, TAttachment& a, bool autoAttach);
		};

		template <typename TAttachment>
		class FaceAttachmentAccessor : public AttachmentAccessor<Face, TAttachment>
		{
			public:
				FaceAttachmentAccessor();
				FaceAttachmentAccessor(const FaceAttachmentAccessor& aa);
				FaceAttachmentAccessor(Grid& grid, TAttachment& a);
				FaceAttachmentAccessor(Grid& grid, TAttachment& a, bool autoAttach);
		};

		template <typename TAttachment>
		class VolumeAttachmentAccessor : public AttachmentAccessor<Volume, TAttachment>
		{
			public:
				VolumeAttachmentAccessor();
				VolumeAttachmentAccessor(const VolumeAttachmentAccessor& aa);
				VolumeAttachmentAccessor(Grid& grid, TAttachment& a);
				VolumeAttachmentAccessor(Grid& grid, TAttachment& a, bool autoAttach);
		};

	///	Container used to store associated vertices
		using VertexContainer = std::vector<Vertex*>;
	///	Container used to store associated edges
		using EdgeContainer = std::vector<Edge*>;
	///	Container used to store associated faces
		using FaceContainer = std::vector<Face*>;
	///	Container used to store associated volumes
		using VolumeContainer = std::vector<Volume*>;

	///	used to iterate over associated edges of vertices, faces and volumes
		using AssociatedEdgeIterator = EdgeContainer::iterator;
	///	used to iterate over associated faces of vertices, edges and volumes
		using AssociatedFaceIterator = FaceContainer::iterator;
	///	used to iterate over associated volumes of vertices, edges and faces
		using AssociatedVolumeIterator = VolumeContainer::iterator;

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
	/**	While all elements and the options are copied completely from the source-grid,
	 *	the attachments are only copied if their pass-on behaviour is set to true.*/
		Grid(const Grid& grid);

		virtual ~Grid();

	///	copies all elements and some attachments from the passed grid to this grid.
	/**	While all elements and the options are copied completely from the source-grid,
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
		[[nodiscard]] uint get_options() const;
		void enable_options(uint options);	///< see set_options for a description of valid parameters.
		void disable_options(uint options);	///< see set_options for a description of valid parameters.
		[[nodiscard]] bool option_is_enabled(uint option) const;///< see set_options for a description of valid parameters.

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
	/**	The method returns nullptr, if no distributed grid manager for the given
	 * grid exists. This should be the case for serial environments or for
	 * serial grids in parallel environments.
	 * Use ug::Grid::set_parallel() to enable or disable parallelism in a grid.
	 * You have to include "lib_grid/parallelization/distributed_grid.h" to access
	 * methods of a distributed grid manager.
	 * \{ */
		inline DistributedGridManager* distributed_grid_manager();
		[[nodiscard]] inline const DistributedGridManager* distributed_grid_manager() const;
	/** \} */

	////////////////////////////////////////////////
	//	periodic boundaries
	/// tell the grid whether it may contain periodic boundaries.
	/**
	 * If the grid may contain periodic boundaries, it instantiate a PeriodicBoundaryManager.*/
		void set_periodic_boundaries(bool);
	/// returns true, if grid has the possibility to handle periodic boundaries.
		[[nodiscard]] bool has_periodic_boundaries() const;

	///	returns a pointer to the associated periodic boundary manager.
	/**	The method returns nullptr, if no periodic boundary get_attachment_accessor for the given
	 * grid exists.
	 * Use ug::Grid::set_periodic_boundaries() to enable or disable periodic boundaries in a grid.
	 * You have to include "lib_grid/tools/periodic_boundary_manager.h" to access
	 * methods of the peridodic boundary manager.
	 * \{ */
		PeriodicBoundaryManager* periodic_boundary_manager();
		[[nodiscard]] const PeriodicBoundaryManager* periodic_boundary_manager() const;
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
	 * TGeomObj has to be a geometric object type as described in grid_base_objects.h.
	 * You may optionally specify a GridObject pParent.
	 * pParent may be used by observers to initialize the created object in a specific way.
	 */
		template <typename TGeomObj>
		typename geometry_traits<TGeomObj>::iterator
		create(GridObject* pParent = nullptr);

	///	create a custom element from a descriptor.
	/**
	 * TGeomObj has to be a geometric object type as described in grid_base_objects.h.
	 * You may optionally specify a GridObject pParent.
	 * pParent may be used by observers to initialize the created object in a specific way.
	 */
		template <typename TGeomObj>
		typename geometry_traits<TGeomObj>::iterator
		create(const typename geometry_traits<TGeomObj>::Descriptor& descriptor,
				GridObject* pParent = nullptr);

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
		template <typename TGeomObj>
		typename geometry_traits<TGeomObj>::iterator
		create_and_replace(typename geometry_traits<TGeomObj>::grid_base_object* pReplaceMe);

	///	this method creates a new vertex, which has the same type as pCloneMe.
		VertexIterator create_by_cloning(Vertex* pCloneMe, GridObject* pParent = nullptr);

	///	this method creates a new edge, which has the same type as pCloneMe.
		EdgeIterator create_by_cloning(Edge* pCloneMe, const IVertexGroup& ev, GridObject* pParent = nullptr);

	///	this method creates a new face, which has the same type as pCloneMe.
		FaceIterator create_by_cloning(Face* pCloneMe, const IVertexGroup& fv, GridObject* pParent = nullptr);

	///	this method creates a new volume, which has the same type as pCloneMe.
		VolumeIterator create_by_cloning(Volume* pCloneMe, const IVertexGroup& vv, GridObject* pParent = nullptr);

	///	Reserves memory for the creation of the given object type
	/**	Calls to this method are optional, but can improve runtime.
	 * Specify the total number of objects which the grid should be
	 * capable to hold (if more are required, the grid will automatically
	 * adjust sizes)
	 */
		template <typename TGeomObj>
		void reserve(size_t num);

	////////////////////////////////////////////////
	//	element deletion
		void erase(GridObject* geomObj);
		void erase(Vertex* vrt);
		void erase(Edge* edge);
		void erase(Face* face);
		void erase(Volume* vol);

	/**\todo: This erase method can cause problems if used with multi-grids.*/
		template <typename GeomObjIter>
		void erase(const GeomObjIter& iterBegin, const GeomObjIter& iterEnd);

		template <typename TGeomObj>
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
		bool replace_vertex(Vertex* vrtOld, Vertex* vrtNew);

	///	checks if replace_vertex would be a valid operation
	/**
	 * Checks if a call of replace_vertex with vertices vrtOld and vrtNew
	 * would lead to degenerate elements. If so false is returned.
	 * If not true is returned.
	 *
	 * requires options in GRIDOPT_VERTEXCENTRIC_INTERCONNECTION.
	 */
		bool replace_vertex_is_valid(Vertex* vrtOld, Vertex* vrtNew);

	///	notifies the grid that two objects will be merged.
	/**	The grid forwards this notification to its GridObservers.
	 * The notification is not relevant for the grid itself.
	 * \{ */
		void objects_will_be_merged(Vertex* target, Vertex* elem1,
									Vertex* elem2);
		void objects_will_be_merged(Edge* target, Edge* elem1,
									Edge* elem2);
		void objects_will_be_merged(Face* target, Face* elem1,
									Face* elem2);
		void objects_will_be_merged(Volume* target, Volume* elem1,
									Volume* elem2);
	/**	\} */
	/*
	///	VrtPairIterator has to be an iterator with value-type std::pair<Vertex*, Vertex*>
		template <typename VrtPairIter>
		void replace_vertices(VrtPairIter& iterBegin, VrtPairIter& iterEnd);
	*/

	////////////////////////////////////////////////
	//	geometric-object-collection
	///	returns the the GridObjectCollection of the grid:
		virtual GridObjectCollection get_grid_objects();

	////////////////////////////////////////////////
	///	flips the orientation of an edge.
		void flip_orientation(Edge* e); // ø why is this not a responsibility of the edge self?
		
	////////////////////////////////////////////////
	///	flips the orientation of a face.
		void flip_orientation(Face* f); // ø why is this not a responsibility of the face self?

	////////////////////////////////////////////////
	///	flips the orientation of a volume.
		void flip_orientation(Volume* vol); // ø why is this not a responsibility of the volume self?
		
	////////////////////////////////////////////////
	//	Iterators
		template <typename TGeomObj>
		typename geometry_traits<TGeomObj>::iterator
		begin();

		template <typename TGeomObj>
		typename geometry_traits<TGeomObj>::iterator
		end();

		inline VertexIterator vertices_begin() {return begin<Vertex>();}
		inline VertexIterator vertices_end() {return end<Vertex>();}
		inline EdgeIterator edges_begin() {return begin<Edge>();}
		inline EdgeIterator edges_end() {return end<Edge>();}
		inline FaceIterator faces_begin() {return begin<Face>();}
		inline FaceIterator faces_end() {return end<Face>();}
		inline VolumeIterator volumes_begin() {return begin<Volume>();}
		inline VolumeIterator volumes_end() {return end<Volume>();}
	
		template <typename TGeomObj>
		typename geometry_traits<TGeomObj>::const_iterator
		begin() const;

		template <typename TGeomObj>
		typename geometry_traits<TGeomObj>::const_iterator
		end() const;

	///	returns the first element of the given type.
	/**	Make sure that elements of the given type exist!
	 *	Behaviour is undefined, if not.*/
		template <typename TGeomObj> TGeomObj* front();
		
	///	returns the last element of the given type.
	/**	Make sure that elements of the given type exist!
	 *	Behaviour is undefined, if not.*/
		template <typename TGeomObj> TGeomObj* back();
		
	//	element manipulation
	/*
		void set_vertices(Edge* e, EdgeDesctiptor& ed);
		void set_vertices(Face* f, FaceDescriptor& fd);
		void set_vertices(Volume* v, VolumeDescriptor& vd);
	*/

	////////////////////////////////////////////////
	//	element numbers
		template <typename TGeomObj>
		[[nodiscard]] size_t num() const;
		[[nodiscard]] inline size_t num_vertices() const {return num<Vertex>();}
		[[nodiscard]] inline size_t num_edges() const {return num<Edge>();}
		[[nodiscard]] inline size_t num_faces()	const {return num<Face>();}
		[[nodiscard]] inline size_t num_volumes()const {return num<Volume>();}

		size_t vertex_fragmentation();	///< returns the number of unused vertex-data-entries.
		size_t edge_fragmentation();		///< returns the number of unused edge-data-entries.
		size_t face_fragmentation();		///< returns the number of unused face-data-entries.
		size_t volume_fragmentation();	///< returns the number of unused volume-data-entries.

	///	returns the size of the associated attachment containers.
	/**	valid types for TGeomObj are Vertex, Edge, Face, Volume.*/
		template <typename TGeomObj>
		[[nodiscard]] size_t attachment_container_size() const;

	////////////////////////////////////////////////
	//	connectivity-information
	///	returns the edge between v1 and v2, if it exists. Returns nullptr if not.
		Edge* get_edge(Vertex* v1, Vertex* v2);
	///	returns the edge that is described by ev.
	/**	Note that you may pass an EdgeDescriptor to this method.*/
		Edge* get_edge(const EdgeVertices& ev);
	///	If it exists, this method returns the i-th edge of the given Face. If not nullptr is returned.
	/**	To make sure that associated edges always exist, enable the grid-option
	 *	FACEOPT_AUTOGENERATE_EDGES.
	 *	For maximal performance, the option FACEOPT_STORE_ASSOCIATED_EDGES
	 *	should be enabled.
	 *	\note	If all edges of a face have to be processed, it may be faster
	 *			to get all edges in one single call.
	 *			This can be done using Grid::associated_elements.*/
		Edge* get_edge(Face* f, int ind);
	///	If it exists, this method returns the i-th edge of the given Volume. If not nullptr is returned.
	/**	To make sure that associated edges always exist, enable the grid-option
	 *	VOLOPT_AUTOGENERATE_EDGES.
	 *	For maximal performance, the option VOLOPT_STORE_ASSOCIATED_EDGES
	 *	should be enabled.
	 *	\note	If all edges of a volume have to be processed, it may be faster
	 *			to get all edges in one single call.
	 *			This can be done using Grid::associated_elements.*/
		Edge* get_edge(Volume* v, int ind);
	///	returns the face that is described by fv.
	/**	Note that you may pass a FaceDescriptor to this method.*/
		Face* get_face(const FaceVertices& fv);
	///	If it exists, this method returns the i-th face of the given Volume. If not nullptr is returned.
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
		Volume* get_volume(const VolumeVertices& vv);
		
	///	returns the element for the given vertices.
	/**	Note that you can either pass an element type (Edge, Face, Volume)
	 * or a descriptor (EdgeDescriptor, FaceDescriptor, VolumeDescriptor).
	 * The method returns nullptr, if the specified element does not exist.
	 * A special overload exists for Vertex*, which simply returns the
	 * specified vertex. Useful for template programming...
	 * \{ */
		Edge* get_element(const EdgeVertices& ev)		{return get_edge(ev);}
		Face* get_element(const FaceVertices& fv)		{return get_face(fv);}
		Volume* get_element(const VolumeVertices& vv)	{return get_volume(vv);}
	/**	\} */

	///	This overload is only useful to avoid compile issues in templated code
		Vertex* get_element(const VertexDescriptor& vd) const {return vd.vertex();}
		

	////////////////////////////////////////////////
	//	access to the sides of an geometric object
	///	This method returns the i-th side of an Edge, Face or Volume.
	/**	If obj has dimension d, then all associated elements of dimension d-1
	 *	are regarded as sides. (Face -> Edge). Only derivates of Volume,
	 *	Face or Edge may be queried for their sides. If you call this method
	 *	with Vertex*, an assertion is triggered, since vertices do not have sides.
	 *
	 *	It is not in all cases guaranteed that an object has sides.
	 *	If i.e. the FACEOPT_AUTOGENERATE_EDGES is not enabled in the grids options,
	 *	then it is not guaranteed that all side-edges of each face exist (in most
	 *	cases they wont exist). The method returns nullptr in this case.
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
		Vertex::side* get_side(Vertex* obj, size_t side);
		Edge::side* get_side(Edge* obj, size_t side);
		Face::side* get_side(Face* obj, size_t side);
		Volume::side* get_side(Volume* obj, size_t side);

	///	returns the geometric object on the opposing side of the given vertex regarding the given element.
	/**	\note Currently only implemented for Face and Volume.
	 * \{ */
		GridObject* get_opposing_object(Vertex* vrt, Face* elem);
		GridObject* get_opposing_object(Vertex* vrt, Volume* elem);
	/** \} */


	///	Puts all elements of type TAss which are contained in 'e' or which contain 'e' into elemsOut
	/**	One shouldn't depend on the order of elements in elemsOut.
	 * Use Grid::associated_elements_sorted if the order matters.
	 *
	 * \note	The returned container is only valid until changes to the queried
	 * 			grid are performed.
	 *
	 * \note	The following type definitions might ease your work with associated elements:
	 * 			AssociatedVertices, AssociatedEdges, AssociatedFaces and
	 * 			AssociatedVolumes.
	 * 			Those are type definitions for Grid::traits<Vertex>::container, ...
	 *
	 * \note	Depending on the current grid options, this method may use Grid::mark.
	 * Valid arguments for TElem are Vertex, Edge, Face, Volume.
	 * \sa Grid::associated_elements_sorted
	 * \{ */
		template <typename TElem>
		void associated_elements(traits<Vertex>::secure_container& elemsOut, TElem* e);
		template <typename TElem>
		void associated_elements(traits<Edge>::secure_container& elemsOut, TElem* e);
		template <typename TElem>
		void associated_elements(traits<Face>::secure_container& elemsOut, TElem* e);
		template <typename TElem>
		void associated_elements(traits<Volume>::secure_container& elemsOut, TElem* e);
	/** \} */
	
	///	Puts all elements of type TAss which are contained in 'e' into elemsOut in the reference elements order
	/**	The order of elements in elemsOut is the same as the order in the reference element.
	 * Note that, depending on active grid options, this method may take more
	 * time than associated_elements, and should thus only be invoked, if the order
	 * of elements matters. Use Grid::associated_elements if the order does not matter.
	 *
	 * Valid arguments for TElem are Vertex, Edge, Face, Volume.
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
	 * \note	The following type definitoins might ease your work with associated elements:
	 * 			AssociatedVertices, AssociatedEdges, AssociatedFaces and
	 * 			AssociatedVolumes.
	 * 			Those are type definitoins for traits<Vertex>::container, ...
	 *
	 * \note	Depending on the current grid options, this method may use Grid::mark.
	 * \sa Grid::associated_elements*/
		template <typename TElem>
		void associated_elements_sorted(traits<Vertex>::secure_container& elemsOut, TElem* e);
		template <typename TElem>
		void associated_elements_sorted(traits<Edge>::secure_container& elemsOut, TElem* e);
		template <typename TElem>
		void associated_elements_sorted(traits<Face>::secure_container& elemsOut, TElem* e);
		template <typename TElem>
		void associated_elements_sorted(traits<Volume>::secure_container& elemsOut, TElem* e);


	////////////////////////////////////////////////
	//	attachments
	////////////////////////////////////////////////////////////////////////
	///	attach with custom pass-on-behaviour and unspecified default value.
		template <typename TGeomObjClass>
		void attach_to(IAttachment& attachment, bool passOnValues);

		inline void attach_to_vertices(IAttachment& attachment, bool passOnValues)	{attach_to<Vertex>(attachment, passOnValues);}
		inline void attach_to_edges(IAttachment& attachment, bool passOnValues)		{attach_to<Edge>(attachment, passOnValues);}
		inline void attach_to_faces(IAttachment& attachment, bool passOnValues)		{attach_to<Face>(attachment, passOnValues);}
		inline void attach_to_volumes(IAttachment& attachment, bool passOnValues)		{attach_to<Volume>(attachment, passOnValues);}

		inline void attach_to_all(IAttachment& attachment, bool passOnValues);

	////////////////////////////////////////////////////////////////////////
	///	attach with default pass-on behaviour and unspecified default value.
		template <typename TGeomObjClass>
		inline void attach_to(IAttachment& attachment)	{attach_to<TGeomObjClass>(attachment, attachment.default_pass_on_behaviour());}

		inline void attach_to_vertices(IAttachment& attachment)	{attach_to<Vertex>(attachment);}
		inline void attach_to_edges(IAttachment& attachment)	{attach_to<Edge>(attachment);}
		inline void attach_to_faces(IAttachment& attachment)	{attach_to<Face>(attachment);}
		inline void attach_to_volumes(IAttachment& attachment)	{attach_to<Volume>(attachment);}

	///	attaches to vertices, edges, faces and volumes at once.
		inline void attach_to_all(IAttachment& attachment);

	////////////////////////////////////////////////////////////////////////
	//	attach with specified default value and default pass-on behaviour
		template <typename TGeomObjClass, typename TAttachment>
		void attach_to_dv(TAttachment& attachment, const typename TAttachment::ValueType& defaultValue);

		template <typename TAttachment>
		inline void attach_to_vertices_dv(TAttachment& attachment, const typename TAttachment::ValueType& defaultValue)	{attach_to_dv<Vertex>(attachment, defaultValue);}
		template <typename TAttachment>
		inline void attach_to_edges_dv(TAttachment& attachment, const typename TAttachment::ValueType& defaultValue)	{attach_to_dv<Edge>(attachment, defaultValue);}
		template <typename TAttachment>
		inline void attach_to_faces_dv(TAttachment& attachment, const typename TAttachment::ValueType& defaultValue)	{attach_to_dv<Face>(attachment, defaultValue);}
		template <typename TAttachment>
		inline void attach_to_volumes_dv(TAttachment& attachment, const typename TAttachment::ValueType& defaultValue)	{attach_to_dv<Volume>(attachment, defaultValue);}

	///	attaches to vertices, edges, faces and volumes at once.
		template <typename TAttachment>
		inline void attach_to_all_dv(TAttachment& attachment, const typename TAttachment::ValueType& defaultValue);

	////////////////////////////////////////////////////////////////////////
	//	attach with specified default value and custom pass-on behaviour
		template <typename TGeomObjClass, typename TAttachment>
		void attach_to_dv(TAttachment& attachment, const typename TAttachment::ValueType& defaultValue, bool passOnValues);

		template <typename TAttachment>
		inline void attach_to_vertices_dv(TAttachment& attachment, const typename TAttachment::ValueType& defaultValue, bool passOnValues)	{attach_to_dv<Vertex>(attachment, defaultValue, passOnValues);}
		template <typename TAttachment>
		inline void attach_to_edges_dv(TAttachment& attachment, const typename TAttachment::ValueType& defaultValue, bool passOnValues)		{attach_to_dv<Edge>(attachment, defaultValue, passOnValues);}
		template <typename TAttachment>
		inline void attach_to_faces_dv(TAttachment& attachment, const typename TAttachment::ValueType& defaultValue, bool passOnValues)		{attach_to_dv<Face>(attachment, defaultValue, passOnValues);}
		template <typename TAttachment>
		inline void attach_to_volumes_dv(TAttachment& attachment, const typename TAttachment::ValueType& defaultValue, bool passOnValues)	{attach_to_dv<Volume>(attachment, defaultValue, passOnValues);}

		template <typename TAttachment>
		inline void attach_to_all_dv(TAttachment& attachment, const typename TAttachment::ValueType& defaultValue, bool passOnValues);
	////////////////////////////////////////////////////////////////////////
	//	detach
		template <typename TGeomObjClass>
		void detach_from(IAttachment& attachment);

		inline void detach_from_vertices(IAttachment& attachment)	{detach_from<Vertex>(attachment);}
		inline void detach_from_edges(IAttachment& attachment)		{detach_from<Edge>(attachment);}
		inline void detach_from_faces(IAttachment& attachment)		{detach_from<Face>(attachment);}
		inline void detach_from_volumes(IAttachment& attachment)	{detach_from<Volume>(attachment);}

		inline void detach_from_all(IAttachment& attachment);


		template <typename TGeomObjClass>
		inline bool has_attachment(IAttachment& attachment)			{return get_attachment_pipe<TGeomObjClass>().has_attachment(attachment);}

		inline bool has_vertex_attachment(IAttachment& attachment)	{return has_attachment<Vertex>(attachment);}
		inline bool has_edge_attachment(IAttachment& attachment)	{return has_attachment<Edge>(attachment);}
		inline bool has_face_attachment(IAttachment& attachment)	{return has_attachment<Face>(attachment);}
		inline bool has_volume_attachment(IAttachment& attachment)	{return has_attachment<Volume>(attachment);}

	////////////////////////////////////////////////////////////////////////
	//	direct attachment access
		template <typename TGeomObj>
		uint get_attachment_data_index(TGeomObj* pObj) const;

		template <typename TGeomObj, typename TAttachment>
		typename TAttachment::ContainerType*
		get_attachment_data_container(TAttachment& attachment);
		
	///	returns the attachment-pipe in which data associated with the given objects-types are stored.
	/**	This method is seldomly used, can however be useful. Use with care.
	 * If in doubt please use the methods featured by Grid instead of directly
	 * operating on the attachment pipe.
	 */
		template <typename TGeomObj>
		typename traits<TGeomObj>::AttachmentPipe&
		get_attachment_pipe();

	////////////////////////////////////////////////
	//	observers
		/**
		 * observerType may be any or-combination of constants enumerated in ObserverType.
		 */
		void register_observer(GridObserver* observer, ObserverType observerType = ObserverType::OT_FULL_OBSERVER);
		void unregister_observer(GridObserver* observer);

/*
		template <typename GeomObjClass>
		util::IAttachmentDataContainer* get_data_container(util::IAttachment& attachment);

		util::IAttachmentDataContainer* get_vertex_data_container(util::IAttachment& attachment)	{return get_data_container<Vertex>(attachment);}
		util::IAttachmentDataContainer* get_edge_data_container(util::IAttachment& attachment)		{return get_data_container<Edge>(attachment);}
		util::IAttachmentDataContainer* get_face_data_container(util::IAttachment& attachment)		{return get_data_container<Face>(attachment);}
		util::IAttachmentDataContainer* get_volume_data_container(util::IAttachment& attachment)	{return get_data_container<Volume>(attachment);}
*/

	////////////////////////////////////////////////
	//	pass_on_values
		void pass_on_values(Vertex* objSrc, Vertex* objDest);
		void pass_on_values(Edge* objSrc, Edge* objDest);
		void pass_on_values(Face* objSrc, Face* objDest);
		void pass_on_values(Volume* objSrc, Volume* objDest);

	//	subject to change!
		AssociatedEdgeIterator associated_edges_begin(Vertex* vrt);///< DO NOT INVOKE! Subject to change.
		AssociatedEdgeIterator associated_edges_end(Vertex* vrt);///< DO NOT INVOKE! Subject to change.
		AssociatedEdgeIterator associated_edges_begin(Face* face);///< DO NOT INVOKE! Subject to change.
		AssociatedEdgeIterator associated_edges_end(Face* face);///< DO NOT INVOKE! Subject to change.
		AssociatedEdgeIterator associated_edges_begin(Volume* vol);///< DO NOT INVOKE! Subject to change.
		AssociatedEdgeIterator associated_edges_end(Volume* vol);///< DO NOT INVOKE! Subject to change.

		AssociatedFaceIterator associated_faces_begin(Vertex* vrt);///< DO NOT INVOKE! Subject to change.
		AssociatedFaceIterator associated_faces_end(Vertex* vrt);///< DO NOT INVOKE! Subject to change.
		AssociatedFaceIterator associated_faces_begin(Edge* edge);///< DO NOT INVOKE! Subject to change.
		AssociatedFaceIterator associated_faces_end(Edge* edge);///< DO NOT INVOKE! Subject to change.
		AssociatedFaceIterator associated_faces_begin(Volume* vol);///< DO NOT INVOKE! Subject to change.
		AssociatedFaceIterator associated_faces_end(Volume* vol);///< DO NOT INVOKE! Subject to change.

		AssociatedVolumeIterator associated_volumes_begin(Vertex* vrt);///< DO NOT INVOKE! Subject to change.
		AssociatedVolumeIterator associated_volumes_end(Vertex* vrt);///< DO NOT INVOKE! Subject to change.
		AssociatedVolumeIterator associated_volumes_begin(Edge* edge);///< DO NOT INVOKE! Subject to change.
		AssociatedVolumeIterator associated_volumes_end(Edge* edge);///< DO NOT INVOKE! Subject to change.
		AssociatedVolumeIterator associated_volumes_begin(Face* face);///< DO NOT INVOKE! Subject to change.
		AssociatedVolumeIterator associated_volumes_end(Face* face);///< DO NOT INVOKE! Subject to change.

	////////////////////////////////////////////////
	//	advanced element manipulation
		inline void register_element(Vertex* v, GridObject* pParent = nullptr)	{register_vertex(v, pParent);}
		inline void unregister_element(Vertex* v)									{unregister_vertex(v);}
		inline void register_element(Edge* e, GridObject* pParent = nullptr)		{register_edge(e, pParent);}
		inline void unregister_element(Edge* e)										{unregister_edge(e);}
		inline void register_element(Face* f, GridObject* pParent = nullptr)			{register_face(f, pParent);}
		inline void unregister_element(Face* f)											{unregister_face(f);}
		inline void register_element(Volume* v, GridObject* pParent = nullptr)		{register_volume(v, pParent);}
		inline void unregister_element(Volume* v)										{unregister_volume(v);}

	///	registers the given element and replaces the old one. Calls pass_on_values.
	/// \{
		void register_and_replace_element(Vertex* v, Vertex* pReplaceMe);
		void register_and_replace_element(Edge* e, Edge* pReplaceMe);
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
		inline void mark(GridObject* obj);
		inline void mark(Vertex* obj);
		inline void mark(Edge* obj);
		inline void mark(Face* obj);
		inline void mark(Volume* obj);
	/**	\} */

	///	marks all objects between begin and end
	/**	TIterator::value_type has to be either
	 *	Vertex*, Edge*, Face* or Volume*.*/
		template <typename TIterator>
		void mark(TIterator begin, TIterator end);
		
	///	unmarks the object. Calls are only valid between calls to Grid::begin_marking and Grid::end_marking.
	/**	Only pass objects that are contained by the grid.
	 * \{ */
		inline void unmark(GridObject* obj);
		inline void unmark(Vertex* obj);
		inline void unmark(Edge* obj);
		inline void unmark(Face* obj);
		inline void unmark(Volume* obj);
	/** \} */

	///	unmarks all objects between begin and end
	/**	TIterator::value_type has to be either
	 *	Vertex*, Edge*, Face* or Volume*.*/
		template <typename TIterator>
		void unmark(TIterator begin, TIterator end);
		
	///	returns true if the object is marked, false if not.
	/**	Only pass objects that are contained by the grid.
	 * \{ */
		inline bool is_marked(GridObject* obj) const;
		inline bool is_marked(Vertex* obj) const;
		inline bool is_marked(Edge* obj) const;
		inline bool is_marked(Face* obj) const;
		inline bool is_marked(Volume* obj) const;
	/** \} */

	///	ends a marking sequence. Call this method when you're done with marking.
		void end_marking();
		
	///	gives access to the grid's message-hub
		SPMessageHub message_hub()		{return m_messageHub;}

	///	a temporary testing method
		void test_attached_linked_lists();

	protected:
		using ObserverContainer = std::vector<GridObserver*>;

		using AVertexContainer = Attachment<VertexContainer>;
		using AEdgeContainer = Attachment<EdgeContainer>;
		using AFaceContainer = Attachment<FaceContainer>;
		using AVolumeContainer = Attachment<VolumeContainer>;

		using AMark = Attachment<int>;

	protected:
	///	unregisters all observers. Call this method in destructors of derived classes.
	/**	If the derived class is an observer itself and if you don't want it to be
	 * notified on grid-destruction, e.g., because you call this method in the
	 * destructor of your derived class, then pass a pointer to your class through
	 * the initiator parameter to this function.
	 * \param initiator:	The initiator won't be notified about grid destruction
	 */
		void notify_and_clear_observers_on_grid_destruction(GridObserver* initiator = nullptr);

	///	returns the element storage for a given element type
		template <typename TElem> inline
		typename traits<TElem>::ElementStorage&
		element_storage()
		{return ElementStorageSelector<typename geometry_traits<TElem>::grid_base_object>::
				element_storage(m_vertexElementStorage, m_edgeElementStorage,
								m_faceElementStorage, m_volumeElementStorage);
		}

	///	returns the const element storage for a given element type
		template <typename TElem> inline
		const typename traits<TElem>::ElementStorage&
		element_storage() const
		{return ElementStorageSelector<typename geometry_traits<TElem>::grid_base_object>::
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
		inline void assign_hash_value(Vertex* vrt)	{vrt->m_hashValue = m_hashCounter++;}

		void register_vertex(Vertex* v, GridObject* pParent = nullptr);///< pDF specifies the element from which v derives its values
		void unregister_vertex(Vertex* v);
		void register_edge(Edge* e, GridObject* pParent = nullptr,
						Face* createdByFace = nullptr, Volume* createdByVol = nullptr);///< pDF specifies the element from which v derives its values
		void unregister_edge(Edge* e);
		void register_face(Face* f, GridObject* pParent = nullptr,
						   Volume* createdByVol = nullptr);///< pDF specifies the element from which v derives its values
		void unregister_face(Face* f);
		void register_volume(Volume* v, GridObject* pParent = nullptr);///< pDF specifies the element from which v derives its values
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

		template <typename TAttachmentPipe, typename TElem>
		void pass_on_values(TAttachmentPipe& attachmentPipe,
							TElem* pSrc, TElem* pDest);

	//	some methods that simplify auto-enabling of grid options
		inline void autoenable_option(uint option, const char* caller, const char* optionName);

	//	neighbourhood access
		template <typename TGeomObj>
		Edge* find_edge_in_associated_edges(TGeomObj* obj,
											const EdgeVertices& ev);
												
		template <typename TGeomObj>
		Face* find_face_in_associated_faces(TGeomObj* obj,
											const FaceVertices& fv);
												
		template <typename TGeomObj>
		Volume* find_volume_in_associated_volumes(TGeomObj* obj,
												  const VolumeVertices& vv);

	//	get associated elements
		void get_associated(SecureVertexContainer& vrts, Edge* e);
		void get_associated(SecureVertexContainer& vrts, Face* f);
		void get_associated(SecureVertexContainer& vrts, Volume* v);

		void get_associated(SecureEdgeContainer& edges, Vertex* v);
		void get_associated(SecureEdgeContainer& edges, Face* f);
		void get_associated(SecureEdgeContainer& edges, Volume* v);

		void get_associated(SecureFaceContainer& faces, Vertex* v);
		void get_associated(SecureFaceContainer& faces, Edge* e);
		void get_associated(SecureFaceContainer& faces, Volume* v);

		void get_associated(SecureVolumeContainer& vols, Vertex* v);
		void get_associated(SecureVolumeContainer& vols, Edge* e);
		void get_associated(SecureVolumeContainer& vols, Face* f);

		template <typename TContainer>
		void get_associated(TContainer& container, GridObject* o);

		template <typename TElem>
		void get_associated(typename traits<typename TElem::grid_base_object>
							::secure_container& elems, TElem* e);
		
	/**	this method does not use possibly attached containers and can thus
	 * be used, when such containers are to be built.*/
		void get_associated_vols_raw(SecureVolumeContainer& vols, Face* f);

		void get_associated_sorted(SecureVertexContainer& vrts, Edge* e) const;
		void get_associated_sorted(SecureVertexContainer& vrts, Face* f) const;
		void get_associated_sorted(SecureVertexContainer& vrts, Volume* v) const;

		void get_associated_sorted(SecureEdgeContainer& edges, Vertex* v);
		void get_associated_sorted(SecureEdgeContainer& edges, Face* f);
		void get_associated_sorted(SecureEdgeContainer& edges, Volume* v);

		void get_associated_sorted(SecureFaceContainer& faces, Vertex* v);
		void get_associated_sorted(SecureFaceContainer& faces, Edge* e);
		void get_associated_sorted(SecureFaceContainer& faces, Volume* v);
		
		void get_associated_sorted(SecureVolumeContainer& vols, Vertex* v);
		void get_associated_sorted(SecureVolumeContainer& vols, Edge* e);
		void get_associated_sorted(SecureVolumeContainer& vols, Face* f);

		template <typename TElem>
		void get_associated_sorted(typename traits<typename TElem::grid_base_object>
								   ::secure_container& elems, TElem* e);


	///	helps in copying attachment pipes during assign_grid
	/**	Note that this method only copies attachments with m_userData==1.
	 * \todo	Copy behavior should be changed to all user-attachments.*/
		template <typename TAttachmentPipe>
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
		inline traits<Vertex>::SectionContainer::iterator
		get_iterator(Vertex* o)
		{
			return m_vertexElementStorage.m_sectionContainer.
					get_container().get_iterator(o);
		}

		inline traits<Edge>::SectionContainer::iterator
		get_iterator(Edge* o)
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
		template <typename TElem>
		void clear_attachments();

	protected:
		VertexElementStorage m_vertexElementStorage;
		EdgeElementStorage m_edgeElementStorage;
		FaceElementStorage m_faceElementStorage;
		VolumeElementStorage m_volumeElementStorage;

		uint m_options;
		uint32 m_hashCounter;

	//	observer handling
		ObserverContainer m_gridObservers;
		ObserverContainer m_vertexObservers;
		ObserverContainer m_edgeObservers;
		ObserverContainer m_faceObservers;
		ObserverContainer m_volumeObservers;

	//	interconnection management
		AVertexContainer m_aVertexContainer;
		AEdgeContainer m_aEdgeContainer;
		AFaceContainer m_aFaceContainer;
		AVolumeContainer m_aVolumeContainer;

		AttachmentAccessor<Vertex, AEdgeContainer> m_aaEdgeContainerVERTEX;
		AttachmentAccessor<Vertex, AFaceContainer> m_aaFaceContainerVERTEX;
		AttachmentAccessor<Vertex, AVolumeContainer> m_aaVolumeContainerVERTEX;

		AttachmentAccessor<Edge, AEdgeContainer> m_aaEdgeContainerEDGE;
		AttachmentAccessor<Edge, AFaceContainer> m_aaFaceContainerEDGE;
		AttachmentAccessor<Edge, AVolumeContainer> m_aaVolumeContainerEDGE;

		AttachmentAccessor<Face, AEdgeContainer> m_aaEdgeContainerFACE;
		AttachmentAccessor<Face, AFaceContainer> m_aaFaceContainerFACE;
		AttachmentAccessor<Face, AVolumeContainer> m_aaVolumeContainerFACE;

		AttachmentAccessor<Volume, AEdgeContainer> m_aaEdgeContainerVOLUME;
		AttachmentAccessor<Volume, AFaceContainer> m_aaFaceContainerVOLUME;
		AttachmentAccessor<Volume, AVolumeContainer> m_aaVolumeContainerVOLUME;
		
	//	marks
		int m_currentMark;	// 0: marks inactive. -1: reset-marks (sets currentMark to 1)
		bool m_bMarking;
		AMark m_aMark;
		VertexAttachmentAccessor<AMark>	m_aaMarkVRT;
		EdgeAttachmentAccessor<AMark> m_aaMarkEDGE;
		FaceAttachmentAccessor<AMark> m_aaMarkFACE;
		VolumeAttachmentAccessor<AMark>	m_aaMarkVOL;

		SPMessageHub m_messageHub;
		DistributedGridManager* m_distGridMgr;
		PeriodicBoundaryManager* m_periodicBndMgr;
};

/** \} */
}//end of namespace

#include "grid_impl.hpp"

#endif
