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

#ifndef __H__LIB_GRID__SERIALIZATION__
#define __H__LIB_GRID__SERIALIZATION__

#include <iostream>
#include "common/util/smart_pointer.h"
#include "common/util/binary_buffer.h"
#include "common/serialization.h"
#include "lib_grid/grid_objects/grid_objects.h"
#include "lib_grid/tools/subset_handler_interface.h"
#include "lib_grid/tools/selector_interface.h"
#include "lib_grid/multi_grid.h"
#include "lib_grid/common_attachments.h"
#include "lib_grid/algorithms/attachment_util.h"
#include "lib_grid/parallelization/grid_object_id.h"
// #include "lib_grid/refinement/projectors/refinement_projector.h"
// #include "lib_grid/refinement/projectors/projection_handler.h"

namespace ug
{

/**
 * Methods that perform grid-related serialization are grouped here.
 * \defgroup lib_grid_algorithms_serialization serialization
 * \ingroup lib_grid_algorithms
 * @{
 */

////////////////////////////////////////////////////////////////////////
//	Utilities

/**	\brief	Interface for handling serialization and deserialization of
 * 			data associated with geometric objects.
 *
 * The GeomObjDataSerializer allows to serialize data associated with
 * geometric objects. Before the data will be serialized, write_info is
 * called. Accordingly read_info is called before data is deserialized.
 *
 * Note that this class handles serialization and deserialization at once.
 *
 * Make sure to completely read all data written by the associated write calls.
 *
 * Note that the following type definition exist: VertexDataSerializer,
 * EdgeDataSerializer, FaceDataSerializer, VolumeDataSerializer.
 *
 * If one wants to serialize data of all objects in a grid, he should
 * take a look at GridDataSerializer.
 *
 * If you call read_info and/or read_data directly, make sure to also call
 * deserialization_done after deserialization has been performed for all
 * geometric objects.
 */
template <class TGeomObj>
class GeomObjDataSerializer
{
	public:
		virtual ~GeomObjDataSerializer()	{}

	///	can be used to write arbitrary info to the file.
	/**	Make sure to read everything you've written during read_data.
	 * Default implementation is empty.*/
		virtual void write_info(BinaryBuffer& out) const {};
	///	write data associated with the given object. Pure virtual.
		virtual void write_data(BinaryBuffer& out, TGeomObj* o) const = 0;

	///	Read the info written during write_info here. Default: empty implementation.
		virtual void read_info(BinaryBuffer& in)	{};
	///	read data associated with the given object. Pure virtual.
		virtual void read_data(BinaryBuffer& in, TGeomObj* o) = 0;

	///	this method is called after read_info has been called for all geometric objects.
		virtual void deserialization_starts()					{}

	///	this method will be called after read_info has been called for all geometric objects.
		virtual void deserialization_done()						{}
};

using VertexDataSerializer = GeomObjDataSerializer<Vertex>;
using EdgeDataSerializer = GeomObjDataSerializer<Edge>;
using FaceDataSerializer = GeomObjDataSerializer<Face>;
using VolumeDataSerializer = GeomObjDataSerializer<Volume>;

using SPVertexDataSerializer = SmartPtr<VertexDataSerializer>;
using SPEdgeDataSerializer = SmartPtr<EdgeDataSerializer>;
using SPFaceDataSerializer = SmartPtr<FaceDataSerializer>;
using SPVolumeDataSerializer = SmartPtr<VolumeDataSerializer>;

/**	\brief	Interface for handling serialization and deserialization of
 * 			data associated with all geometric objects in a grid.
 *
 * The GridDataSerializer allows to serialize data associated with
 * all geometric objects in a grid. Before the data will be serialized,
 * write_info is called. Accordingly read_info is called before data is
 * deserialized.
 *
 * Note that this class handles serialization and deserialization at once.
 *
 * Make sure to completely read all data written by the associated write calls.
 *
 * All methods have an empty implementation by default.
 *
 * If you call read_info and/or read_data directly, make sure to also call
 * deserialization_done after deserialization has been performed for all
 * geometric objects.
 */
class GridDataSerializer
{
	public:
		virtual ~GridDataSerializer() = default;

	///	can be used to write arbitrary info to the file.
	/**	Make sure to read everything you've written during read_data.
	 * Default implementation is empty.*/
		virtual void write_info(BinaryBuffer& out) const				{}

	///	Read the info written during write_info here. Default: empty implementation.
		virtual void read_info(BinaryBuffer& in)					{}

		virtual void write_data(BinaryBuffer& out, Vertex* o) const	{}
		virtual void write_data(BinaryBuffer& out, Edge* o) const	{}
		virtual void write_data(BinaryBuffer& out, Face* o) const		{}
		virtual void write_data(BinaryBuffer& out, Volume* o) const		{}

		virtual void read_data(BinaryBuffer& in, Vertex* o)	{}
		virtual void read_data(BinaryBuffer& in, Edge* o)	{}
		virtual void read_data(BinaryBuffer& in, Face* o)		{}
		virtual void read_data(BinaryBuffer& in, Volume* o)		{}

	///	this method is called after read_info has been called for all geometric objects.
		virtual void deserialization_starts()					{}

	///	this method is called after read_info has been called for all geometric objects.
		virtual void deserialization_done()						{}
};

using SPGridDataSerializer = SmartPtr<GridDataSerializer>;


///	Serialization of data associated with grid elements.
/**	Through the add-method callback-classes can be registered,
 * which will be called during serialization and deserialization to
 * write data associated with the given elements into a binary stream.
 *
 * Note when data for a given object-type not only registered
 * callback classes for the type are called, but also registered
 * instances of GridDataSerializer.
 *
 * Note that this class performs both serialization and deserialization.
 *
 * If you call the deserialize method directly, make sure to also call
 * deserialization_done after your last call to deserialize for a given
 * dezerialization.
 */
class GridDataSerializationHandler
{
	public:
		~GridDataSerializationHandler()	{}

	///	Adds a callback class for serialization and deserialization.
	/**	\{ */
		void add(SPVertexDataSerializer cb);
		void add(SPEdgeDataSerializer cb);
		void add(SPFaceDataSerializer cb);
		void add(SPVolumeDataSerializer cb);
		void add(SPGridDataSerializer cb);
	/**	\} */


	///	calls write_info on all registered serializers
		void write_infos(BinaryBuffer& out) const;

	/**	\{
	 * \brief Serializes data associated with the given object.*/
		inline void serialize(BinaryBuffer& out, Vertex* vrt) const;
		inline void serialize(BinaryBuffer& out, Edge* edge) const;
		inline void serialize(BinaryBuffer& out, Face* face) const;
		inline void serialize(BinaryBuffer& out, Volume* vol) const;
	/**	\} */

	///	Calls serialize on all elements between begin and end.
	/**	Make sure that TIterator::value_type is compatible with
	 * either Vertex*, Edge*, Face*, Volume*.*/
		template <class TIterator>
		void serialize(BinaryBuffer& out, TIterator begin, TIterator end) const;

	///	Calls serialize on all elements in the given geometric object collection
		void serialize(BinaryBuffer& out, GridObjectCollection goc) const;

	///	calls read_info on all registered serializers
		void read_infos(BinaryBuffer& in);

	/**	\{
	 * \brief Deserializes data associated with the given object.*/
		inline void deserialize(BinaryBuffer& in, Vertex* vrt);
		inline void deserialize(BinaryBuffer& in, Edge* edge);
		inline void deserialize(BinaryBuffer& in, Face* face);
		inline void deserialize(BinaryBuffer& in, Volume* vol);
	/**	\} */

	///	Calls deserialize on all elements between begin and end.
	/**	Make sure that TIterator::value_type is compatible with
	 * either Vertex*, Edge*, Face*, Volume*.*/
		template <class TIterator>
		void deserialize(BinaryBuffer& in, TIterator begin, TIterator end);

	///	Calls deserialize on all elements in the given geometric object collection
		void deserialize(BinaryBuffer& in, GridObjectCollection goc);

	///	this method will be called before read_infos is called for the first time
	///	in a deserialization run.
		void deserialization_starts();

	///	this method will be called after deserialize was called for the last time
	///	in a deserialization run.
		void deserialization_done();

	private:
	///	performs serialization on all given serializers.
		template<class TGeomObj, class TSerializers>
		void serialize(BinaryBuffer& out, TGeomObj* o,
					   TSerializers& serializers) const;

	///	performs deserialization on all given deserializers.
		template<class TGeomObj, class TDeserializers>
		void deserialize(BinaryBuffer& in, TGeomObj* o,
					   TDeserializers& deserializers);

		template<class TSerializers>
		void write_info(BinaryBuffer& out, TSerializers& serializers) const;

		template<class TSerializers>
		void read_info(BinaryBuffer& in, TSerializers& serializers);

		template<class TSerializers>
		void deserialization_starts(TSerializers& serializers);

		template<class TSerializers>
		void deserialization_done(TSerializers& serializers);

	private:
		std::vector<SPVertexDataSerializer>	m_vrtSerializers;
		std::vector<SPEdgeDataSerializer>	m_edgeSerializers;
		std::vector<SPFaceDataSerializer>	m_faceSerializers;
		std::vector<SPVolumeDataSerializer>	m_volSerializers;
		std::vector<SPGridDataSerializer>	m_gridSerializers;
};


////////////////////////////////////////////////////////////////////////
///	Serialization callback for grid attachments
/**	template class where TGeomObj should be one of the
 * following types: Vertex, Edge, Face, Volume.
 *
 * Note that the attachment is automatically attached, if not yet present.
 */
template <class TGeomObj, class TAttachment>
class GeomObjAttachmentSerializer :
	public GeomObjDataSerializer<TGeomObj>
{
	public:
		static SmartPtr<GeomObjDataSerializer<TGeomObj> >
		create(Grid& g, TAttachment a)
		{return SmartPtr<GeomObjDataSerializer<TGeomObj> >(new GeomObjAttachmentSerializer(g, a));}

		GeomObjAttachmentSerializer(Grid& g, TAttachment a) :
			m_aa(g, a, true)	{}

		virtual ~GeomObjAttachmentSerializer() {};

		virtual void write_data(BinaryBuffer& out, TGeomObj* o) const
		{Serialize(out, m_aa[o]);}

		virtual void read_data(BinaryBuffer& in, TGeomObj* o)
		{Deserialize(in, m_aa[o]);}

	private:
		Grid::AttachmentAccessor<TGeomObj, TAttachment>	m_aa;
};

class SubsetHandlerSerializer : public GridDataSerializer
{
	public:
		static SPGridDataSerializer create(ISubsetHandler& sh)
		{return SPGridDataSerializer(new SubsetHandlerSerializer(sh));}

		SubsetHandlerSerializer(ISubsetHandler& sh);

	///	writes subset-infos to the stream (subset names and colors)
		virtual void write_info(BinaryBuffer& out) const;

		///	Read the info written during write_info here. Default: empty implementation.
		virtual void read_info(BinaryBuffer& in);

		virtual void write_data(BinaryBuffer& out, Vertex* o) const;
		virtual void write_data(BinaryBuffer& out, Edge* o) const;
		virtual void write_data(BinaryBuffer& out, Face* o) const;
		virtual void write_data(BinaryBuffer& out, Volume* o) const;

		virtual void read_data(BinaryBuffer& in, Vertex* o);
		virtual void read_data(BinaryBuffer& in, Edge* o);
		virtual void read_data(BinaryBuffer& in, Face* o);
		virtual void read_data(BinaryBuffer& in, Volume* o);

	private:
		ISubsetHandler& m_sh;
};

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	GRID

////////////////////////////////////////////////////////////////////////
///	Writes a part of the grids elements to a binary-stream.
/**
 * The passed GridObjectCollection goc may only reference
 * elements of the given grid. It is important, that the goc
 * is complete - that means that all referenced vertices are
 * contained in the goc.
 *
 * If you pack several different parts of your grid, you should use
 * this method, since it is faster than calling SerializeGridElements
 * without the attachment.
 *
 * the integer attachment aInt is used during this method to store
 * an index in each vertex of the goc. The initial content of
 * the referenced attachment is ignored.
 *
 * The caller is responsible to attach the aIntVRT attachment to the 
 * to the vertices of the grid before calling this method.
 * The caller is also responsible to detach aIntVRT from the grids
 * vertices when it is no longer required.
 *
 * After termination the attachment holds the indices at which
 * the respcetive vertices are stored in the pack.
 */
bool SerializeGridElements(Grid& grid, GridObjectCollection goc,
						   AInt& aIntVRT, BinaryBuffer& out);

////////////////////////////////////////////////////////////////////////
///	Writes all grid elements into a binary-stream.
bool SerializeGridElements(Grid& grid, BinaryBuffer& out);

////////////////////////////////////////////////////////////////////////
///	Writes a part of the grids elements to a binary-stream.
/**
 * The passed GridObjectCollection goc may only reference
 * elements of the given grid. It is important, that the goc
 * is complete - that means that all referenced vertices are
 * contained in the goc.
 *
 * If you're planning to serialize multiple parts of one grid, you
 * should consider to use the full-featured serialization method.
 */
bool SerializeGridElements(Grid& grid, GridObjectCollection goc,
						   BinaryBuffer& out);

////////////////////////////////////////////////////////////////////////
///	Creates grid elements from a binary stream
/**	Old versions of SerializeGrid did not support grid-headers.
 * This is why you can specify via readGridHeader whether a
 * header should be read (default is true).
 */
bool DeserializeGridElements(Grid& grid, BinaryBuffer& in,
							bool readGridHeader = true);


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	MULTI-GRID

////////////////////////////////////////////////////////////////////////
//	SerializeMultiGridElements
///	writes a part of the elements of a MultiGrid to a binary stream.
/**
 * THIS METHOD USES Grid::mark.
 *
 * The passed GridObjectCollection goc may only
 * reference elements of the given grid. It is important, that the goc
 * is complete - that means that all referenced vertices are
 * contained in the goc.
 * The goc has also to be complete in regard to the multi-grid hierarchy.
 * This means that for each element of the goc, the parent has to be
 * also part of the goc.
 *
 * If you pack several different parts of your grid, you should use
 * this method, since it is faster than calling SerializeGridElements
 * without the attachment.
 *
 * the integer attachment aInt is used during this method to store
 * an index in each element of the goc. The initial content of
 * the referenced attachment is ignored.
 *
 * The caller is responsible to attach the aInt attachment to the
 * to the elements of the grid before calling this method.
 * The caller is also responsible to detach aInt from the grids
 * elements when it is no longer required.
 *
 * Optionally an accessor to global ids in mg can be specified through paaID.
 * Those ids have to be attached and correctly set before the method is called.
 * If you specify those ids, they are serialized together with the grid elements.
 * Make sure to pass a corresponding id-accessor on a call to deserialize.
 *
 * After termination the attachments hold the indices that were
 * assigned to the respective elements - starting from 0 for each
 * element type.
 *
 * \todo	add support for constrained/constraining faces
 * \todo	use ConstVertexArrays instead of virtual functions ...->vertex(...)
 */
bool SerializeMultiGridElements(MultiGrid& mg,
								GridObjectCollection goc,
								MultiElementAttachmentAccessor<AInt>&	aaInt,
								BinaryBuffer& out,
								MultiElementAttachmentAccessor<AGeomObjID>* paaID = nullptr);

////////////////////////////////////////////////////////////////////////
//	SerializeMultiGridElements
///	writes a part of the elements of a MultiGrid to a binary stream.
/**
 * The passed GridObjectCollection goc may only reference
 * elements of the given grid. It is important, that the goc
 * is complete - that means that all referenced vertices are
 * contained in the goc.
 * The goc has also to be complete in regard to the multi-grid hierarchy.
 * This means that for each element of the goc, the parent has to be
 * also part of the goc.
 *
 * If you're planning to serialize multiple parts of one grid, you
 * should consider to use the full-featured serialization method.
 */
bool SerializeMultiGridElements(MultiGrid& mg,
								GridObjectCollection goc,
								BinaryBuffer& out);

////////////////////////////////////////////////////////////////////////
//	SerializeMultiGridElements
///	writes the elements of a MultiGrid to a binary stream.
bool SerializeMultiGridElements(MultiGrid& mg,
								BinaryBuffer& out);
								
////////////////////////////////////////////////////////////////////////
///	Creates multi-grid elements from a binary stream
/**
 * If you pass a pointer to a std::vector using pvVrts, pvEdges,
 * pvFaces or pvVolumes, those vectors will contain the elements
 * of the grid in the order they were read.
 * Specifying those vectors does not lead to a performance loss.
 *
 * An accessor to global ids in the grids elements can optionally be specified
 * through paaID.
 * If it is specified, the existing grid is automatically merged with the new
 * elements based on the global ids. The accessor must only be specified if it
 * was also specified in SerializeMultiGridElements.
 *
 * \todo	add support for constrained/constraining faces
 */
bool DeserializeMultiGridElements(MultiGrid& mg, BinaryBuffer& in,
									std::vector<Vertex*>* pvVrts = nullptr,
									std::vector<Edge*>* pvEdges = nullptr,
									std::vector<Face*>* pvFaces = nullptr,
									std::vector<Volume*>* pvVols = nullptr,
									MultiElementAttachmentAccessor<AGeomObjID>* paaID = nullptr);



////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	ATTACHMENTS

////////////////////////////////////////////////////////////////////////
//	copies attached values to a binary stream.
/**
 * copies attached values of the grids elements of the given type
 * to the binary stream.
 *
 * Make sure that attachment is attached to the specified elements.
 */
template <class TElem, class TAttachment>
bool SerializeAttachment(Grid& grid, TAttachment& attachment,
						 BinaryBuffer& out);

////////////////////////////////////////////////////////////////////////
///	copies attached values to a binary stream.
/**
 * copies attached values of the elements between iterBegin and iterEnd
 * to the binary stream.
 *
 * Make sure that attachment is attached to the specified elements.
 */
template <class TElem, class TAttachment>
bool SerializeAttachment(Grid& grid, TAttachment& attachment,
						 typename geometry_traits<TElem>::iterator iterBegin,
						 typename geometry_traits<TElem>::iterator iterEnd,
						 BinaryBuffer& out);

////////////////////////////////////////////////////////////////////////
///	copies attached values from a binary stream
/**
 * copies values from the given binary stream to the given attachment of
 * elements between iterBegin and iterEnd.
 * If attachment was not attached to the grid, then it will be attached
 * automatically.
 */
template <class TElem, class TAttachment>
bool DeserializeAttachment(Grid& grid, TAttachment& attachment,
						 BinaryBuffer& in);

////////////////////////////////////////////////////////////////////////
///	copies attached values from a binary stream
/**
 * copies values from the given binary stream to the given attachment of
 * elements between iterBegin and iterEnd.
 * If attachment was not attached to the grid, then it will be attached
 * automatically.
 */
template <class TElem, class TAttachment>
bool DeserializeAttachment(Grid& grid, TAttachment& attachment,
						 typename geometry_traits<TElem>::iterator iterBegin,
						 typename geometry_traits<TElem>::iterator iterEnd,
						 BinaryBuffer& in);


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	SUBSET-HANDLER

////////////////////////////////////////////////////////////////////////
///	writes the subset-indices of all elements in the goc to a stream.
bool SerializeSubsetHandler(Grid& grid, ISubsetHandler& sh,
							GridObjectCollection goc,
							BinaryBuffer& out);
							
////////////////////////////////////////////////////////////////////////
///	writes the subset-indices of all elements in the grid to a stream.
bool SerializeSubsetHandler(Grid& grid, ISubsetHandler& sh,
							BinaryBuffer& out);

////////////////////////////////////////////////////////////////////////
///	assigns subset-indices to all elements in the goc from a stream.
/**	One has to be very careful that the given goc only contains
 * the elements that were passed to the serialization routine.
 * Problems could be caused by automatic element creation.
 * consider to set grid.set_option(GRIDOPT_NONE) before loading
 * the grid.
 *
 * readPropertyMap should always be true. It is only contained for backwards
 * compatibility with older binary files, which did not support property maps.
 */
bool DeserializeSubsetHandler(Grid& grid, ISubsetHandler& sh,
							GridObjectCollection goc,
							BinaryBuffer& in,
							bool readPropertyMap = true);

							
////////////////////////////////////////////////////////////////////////
///	assigns subset-indices to all elements in the grid from a stream.
/**	One has to be very careful that the given grid only contains
 * the elements that were passed to the serialization routine.
 * Problems could be caused by automatic element creation.
 * consider to set grid.set_option(GRIDOPT_NONE) before loading
 * the grid.
 *
 * readPropertyMap should always be true. It is only contained for backwards
 * compatibility with older binary files, which did not support property maps.
 */
bool DeserializeSubsetHandler(Grid& grid, ISubsetHandler& sh,
							BinaryBuffer& in,
							bool readPropertyMap = true);


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	SELECTOR

////////////////////////////////////////////////////////////////////////
///	writes the subset-indices of all elements in the goc to a stream.
bool SerializeSelector(Grid& grid, ISelector& sel,
					   GridObjectCollection goc,
					   BinaryBuffer& out);
							
////////////////////////////////////////////////////////////////////////
///	writes the subset-indices of all elements in the grid to a stream.
bool SerializeSelector(Grid& grid, ISelector& sel,
					   BinaryBuffer& out);

////////////////////////////////////////////////////////////////////////
///	assigns subset-indices to all elements in the goc from a stream.
/**	One has to be very careful that the given goc only contains
 * the elements that were passed to the serialization routine.
 * Problems could be caused by automatic element creation.
 * consider to set grid.set_option(GRIDOPT_NONE) before loading
 * the grid.
 *
 * readPropertyMap should always be true. It is only contained for backwards
 * compatibility with older binary files, which did not support property maps.
 */
bool DeserializeSelector(Grid& grid, ISelector& sel,
						 GridObjectCollection goc,
						 BinaryBuffer& in);

							
////////////////////////////////////////////////////////////////////////
///	assigns subset-indices to all elements in the grid from a stream.
/**	One has to be very careful that the given grid only contains
 * the elements that were passed to the serialization routine.
 * Problems could be caused by automatic element creation.
 * consider to set grid.set_option(GRIDOPT_NONE) before loading
 * the grid.
 *
 * readPropertyMap should always be true. It is only contained for backwards
 * compatibility with older binary files, which did not support property maps.
 */
bool DeserializeSelector(Grid& grid, ISelector& sel,
						 BinaryBuffer& in);


// void SerializeProjector(BinaryBuffer& out, RefinementProjector& proj);
// void SerializeProjectionHandler(BinaryBuffer& out, ProjectionHandler& ph);
// SPRefinementProjector DeserializeProjector(BinaryBuffer& in);
// void DeserializeProjectionHandler(BinaryBuffer& in, ProjectionHandler& ph);

/**@}*/ // end of doxygen defgroup command

}

////////////////////////////////
//	include implementation
#include "serialization_impl.hpp"

#endif
