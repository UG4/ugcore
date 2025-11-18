/*
 * Copyright (c) 2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG_global_attachments
#define __H__UG_global_attachments

#include "attachments/attachment_info_traits.h"
#include "attachments/attachment_io_traits.h"
#include "algorithms/serialization.h"
#include <algorithm>

namespace ug{

///	Global attachments are automatically read/written from/to files and are considered during redistribution
class GlobalAttachments {
	public:
		template <class TAttachment>
		static void declare_attachment (const std::string& name,
										bool passOnBehaviour = false)
		{
			const char* typeName = attachment_info_traits<TAttachment>::type_name();
			if(attachments()[name].attachment != nullptr) {
				UG_COND_THROW(
					dynamic_cast<TAttachment*> (attachments()[name].attachment) == nullptr,
						  "Attachment with name '" << name
						  << "' was already declared in GlobalAttachments with a different type. "
						  << "Old type: " << attachments()[name].type <<
						  ", new type: " << typeName);
				return;
			}

			int fi = static_cast<int> (attachment_names().size());
			attachment_names().push_back(name);
			attachments()[name] = AttachmentEntry(new TAttachment(passOnBehaviour), typeName, fi);
			functions<Vertex>().push_back(FunctionEntry<Vertex, TAttachment>());
			functions<Edge>().push_back(FunctionEntry<Edge, TAttachment>());
			functions<Face>().push_back(FunctionEntry<Face, TAttachment>());
			functions<Volume>().push_back(FunctionEntry<Volume, TAttachment>());
		}

		template <class TAttachment>
		static void undeclare_attachment(const std::string& name) {
			UG_COND_THROW(!is_declared(name), "Trying undeclaring a non-declared attachment.");
			AttachmentEntry& ae = attachment_entry(name);
			remove_function_entry<Volume>(ae);
			remove_function_entry<Face>(ae);
			remove_function_entry<Edge>(ae);
			remove_function_entry<Vertex>(ae);
			attachment_names().erase(
				std::remove(attachment_names().begin(), attachment_names().end(), name),
				attachment_names().end());
			attachments().erase(name);
			const char* typeName = attachment_info_traits<TAttachment>::type_name();
			attachment_types().erase(typeName);
		}

		static void declare_attachment (const std::string& name,
										const std::string& typeName,
										bool passOnBehaviour = false)
		{
			UG_COND_THROW(attachment_types()[typeName].declareFunc == 0,
						  "Unregistered attachment type used in "
						  << "GlobalAttachments::declare_attachment: '"
						  << typeName << "' during declaration of attachment '"
						  << name << "'.");
			attachment_types()[typeName].declareFunc(name, passOnBehaviour);
		}

		template <class TAttachment>
		static void register_attachment_type ()
		{
			std::string typeName = attachment_info_traits<TAttachment>::type_name();
			attachment_types()[typeName] = AttachmentType<TAttachment>();
		}

		static const std::vector<std::string>&
		declared_attachment_names ()
		{
			return attachment_names();
		}

		static
		bool is_declared(const std::string& name)
		{
			return attachments().find(name) != attachments().end();
		}

		static
		bool attachment_pass_on_behaviour(const std::string& name)
		{
			UG_COND_THROW(!is_declared(name), "Undeclared attachment queried: " << name);
			return attachments()[name].attachment->default_pass_on_behaviour();
		}

		static
		bool type_is_registered(const std::string& typeName)
		{
			return attachment_types().find(typeName) != attachment_types().end();
		}
		
		#ifdef UG_PARALLEL
		static void SynchronizeDeclaredGlobalAttachments(Grid& grid, int procId)
		{
			if (procId < 0) return; // this is not a parallel run
			
			//declare global attachments on all processors
			pcl::ProcessCommunicator procComm;
			std::vector<std::string> possible_attachment_names = declared_attachment_names();
			// only master proc loaded the grid
			if (procId == 0)
				procComm.broadcast<std::vector<std::string> >(possible_attachment_names, procId);
			else
				UG_THROW("There are more than one proc loading the grid"<<
						"please make sure all processes broadcast their GlobalAttachments");
						
			std::vector<byte_t> locDeclared(possible_attachment_names.size(), 0);
			std::vector<byte_t> globDeclared(possible_attachment_names.size(), 0);
			// record local info
			for(size_t i = 0; i < possible_attachment_names.size(); ++i){
				byte_t& b = locDeclared[i];
				if(is_declared(possible_attachment_names[i])){
					b |= 1;
					if(is_attached<Vertex>(grid, possible_attachment_names[i]))
						b |= 1<<1;
					if(is_attached<Edge>(grid, possible_attachment_names[i]))
						b |= 1<<2;
					if(is_attached<Face>(grid, possible_attachment_names[i]))
						b |= 1<<3;
					if(is_attached<Volume>(grid, possible_attachment_names[i]))
						b |= 1<<4;
				}
			}
			// sum up all the local to the global
			procComm.allreduce(locDeclared, globDeclared, PCL_RO_BOR);
			// update the local with the global
			for(size_t i = 0; i < possible_attachment_names.size(); ++i){
				byte_t& b = globDeclared[i];
				if(b & 1){
					if(!is_declared(possible_attachment_names[i]))
						declare_attachment(possible_attachment_names[i], "double", true);
					if(b & 1<<1)
						attach<Vertex>(grid, possible_attachment_names[i]);
					if(b & 1<<2)
						attach<Edge>(grid, possible_attachment_names[i]);
					if(b & 1<<3)
						attach<Face>(grid, possible_attachment_names[i]);
					if(b & 1<<4)
						attach<Volume>(grid, possible_attachment_names[i]);
				}
			}

		}
		#endif
		
		template <class TElem>
		static
		void attach(Grid& g, const std::string& name)
		{
			AttachmentEntry& ae = attachment_entry(name);
			IFunctionEntry& fe = function_entry<TElem>(ae);
			fe.attach(g, *ae.attachment);
		}

		template <class TElem>
		static
		bool is_attached(Grid& g, const std::string& name)
		{
			AttachmentEntry& ae = attachment_entry(name);
			return g.has_attachment<TElem>(*ae.attachment);
		}

		template <class TAttachment>
		static
		TAttachment attachment (const std::string& name)
		{
			AttachmentEntry& e = attachment_entry(name);
			auto* a = dynamic_cast<TAttachment*>(e.attachment);
			UG_COND_THROW(!a, "Attachment with invalid type queried. Given type "
						  "is " << e.type << ", queried type is " <<
						  attachment_info_traits<TAttachment>::type_name());
			return *a;
		}

		static
		const char* type_name (const std::string& name)
		{
			AttachmentEntry& e = attachment_entry(name);
			return e.type;
		}

		template <class TElem>
		static
		void read_attachment_values (std::istream& in,
									 Grid& grid,
									 const std::string& name)
		{
			AttachmentEntry& ae = attachment_entry(name);
			IFunctionEntry& fe = function_entry<TElem>(ae);
			fe.readFunc(in, grid, *ae.attachment);
		}


		template <class TElem>
		static
		void write_attachment_values (std::ostream& out,
									  Grid& grid,
									  const std::string& name)
		{
			AttachmentEntry& ae = attachment_entry(name);
			function_entry<TElem>(ae).writeFunc(out, grid, *ae.attachment);
		}

		template <class TElem>
		static
		void add_data_serializer(GridDataSerializationHandler& handler,
								 Grid& grid,
								 const std::string& name)
		{
			AttachmentEntry& ae = attachment_entry(name);
			IFunctionEntry& fe = function_entry<TElem>(ae);
			fe.addSerializer(handler, grid, *ae.attachment);
		}


	private:
	////////////////////////////////////////
	//	TYPES
		struct AttachmentEntry {
			AttachmentEntry () : attachment(nullptr), type(""), functionIndex(-1) {}
			AttachmentEntry (IAttachment* a, const char* t, int fi) :
				attachment(a), type(t), functionIndex(fi) {}
			IAttachment* 	attachment;
			const char*		type;
			int				functionIndex;
		};

		struct IFunctionEntry {
			IFunctionEntry() : readFunc(nullptr), writeFunc(nullptr), addSerializer(nullptr), attach(nullptr) {}
			void (*readFunc	)		(std::istream&, Grid&, IAttachment&);
			void (*writeFunc)		(std::ostream&, Grid&, IAttachment&);
			void (*addSerializer)	(GridDataSerializationHandler&, Grid&, IAttachment&);
			void (*attach)			(Grid&, IAttachment&);
		};

		template <class TElem, class TAttachment>
		struct FunctionEntry : public IFunctionEntry {
			FunctionEntry() {
				readFunc = 		&read_attachment_from_stream<TElem, TAttachment>;
				writeFunc = 	&write_attachment_to_stream<TElem, TAttachment>;
				addSerializer = &add_attachment_serializer<TElem, TAttachment>;
				attach = 		&cast_and_attach<TElem, TAttachment>;
			}
		};

		struct IAttachmentType {
			IAttachmentType() : declareFunc(nullptr)	{}
			void (*declareFunc)	(const std::string&, bool);
		};

		template <class TAttachment>
		struct AttachmentType : public IAttachmentType {
			AttachmentType() {
				declareFunc =	&declare_attachment<TAttachment>;
			}
		};


		using AttachmentMap = std::map<std::string, AttachmentEntry>;
		using AttachmentTypeMap = std::map<std::string, IAttachmentType>;
		using FunctionVec = std::vector<IFunctionEntry>;


	////////////////////////////////////////
	//	METHODS
		static
		GlobalAttachments& inst ()
		{
			static GlobalAttachments h;
			return h;
		}

		GlobalAttachments () = default;

		~GlobalAttachments ()
		{
			AttachmentMap& m = attachments();
			for(auto i = m.begin(); i != m.end(); ++i) {
				if(i->second.attachment)
					delete i->second.attachment;
			}
		}

		static
		std::vector<std::string>& attachment_names() {
			return inst().m_attachmentNames;
		}

		static
		AttachmentMap& attachments() {
			return inst().m_attachmentMap;
		}

		static
		AttachmentTypeMap& attachment_types() {
			static bool initialized = false;
			if(!initialized){
				initialized = true;
				register_standard_attachment_types();
			}
			return inst().m_attachmentTypeMap;
		}

		template <class TElem>
		static
		FunctionVec& functions() {
			return inst().m_functionVecs[TElem::BASE_OBJECT_ID];
		}

		static
		AttachmentEntry& attachment_entry(const std::string& name)
		{
			AttachmentEntry& e = attachments()[name];
			UG_COND_THROW(!e.attachment, "Undeclared attachment queried: " << name);
			return e;
		}

		template <class TElem>
		static
		IFunctionEntry& function_entry(const std::string& name)
		{
			return function_entry<TElem>(attachment_entry(name));
		}

		template <class TElem>
		static
		IFunctionEntry& function_entry(const AttachmentEntry& ae)
		{
			return functions<TElem>().at(ae.functionIndex);
		}

		template <class TElem>
		static
		void remove_function_entry(const AttachmentEntry& ae) {
			functions<TElem>().erase(functions<TElem>().begin() + ae.functionIndex);
		}

		template <class TElem, class TAttachment>
		static
		void read_attachment_from_stream (
				std::istream& in,
				Grid& grid,
				IAttachment& attachment)
		{
			auto& a = dynamic_cast<TAttachment&>(attachment);
			
			if(!grid.has_attachment<TElem>(a))
				grid.attach_to<TElem>(a);

			Grid::AttachmentAccessor<TElem, TAttachment> aaVal(grid, a);
			
			for(typename Grid::traits<TElem>::iterator iter = grid.begin<TElem>();
				iter != grid.end<TElem>(); ++iter)
			{
				attachment_io_traits<TAttachment>::read_value(in, aaVal[*iter]);
				UG_COND_THROW(!in, "Failed to read attachment entry.\n");
			}
		}


		template <class TElem, class TAttachment>
		static
		void write_attachment_to_stream (
				std::ostream& out,
				Grid& grid,
				IAttachment& attachment)
		{
			auto& a = dynamic_cast<TAttachment&>(attachment);
			
			if(!grid.has_attachment<TElem>(a))
				return;

			Grid::AttachmentAccessor<TElem, TAttachment> aaVal(grid, a);

			const typename Grid::traits<TElem>::iterator iterEnd = grid.end<TElem>();
			for(typename Grid::traits<TElem>::iterator iter = grid.begin<TElem>();
				iter != iterEnd;)
			{
				attachment_io_traits<TAttachment>::write_value(out, aaVal[*iter]);
				UG_COND_THROW(!out, "Failed to write attachment entry.\n");
				++iter;
				if(iter != iterEnd)
					out << " ";
			}
		}

		template <class TElem, class TAttachment>
		static
		void//SmartPtr<GeomObjDataSerializer<TElem> >
		add_attachment_serializer (
			GridDataSerializationHandler& handler,
			Grid& g,
			IAttachment& attachment)
		{
			auto& a = dynamic_cast<TAttachment&>(attachment);
			handler.add(GeomObjAttachmentSerializer<TElem, TAttachment>::
									create(g, a));
		}

		template <class TElem, class TAttachment>
		static
		void cast_and_attach (
				Grid& grid,
				IAttachment& attachment)
		{
			auto& a = dynamic_cast<TAttachment&>(attachment);
			
			if(!grid.has_attachment<TElem>(a))
				grid.attach_to<TElem>(a);
		}

		static
		void register_standard_attachment_types()
		{
		//todo:	explicit registration of common types in
		//		GlobalAttachments itself isn't really the best way to go
		//		(e.g. requires inclusion of 'ugmath_types.h').
			register_attachment_type<Attachment<bool> >();
			register_attachment_type<Attachment<char> >();
			register_attachment_type<Attachment<byte_t> >();
			register_attachment_type<Attachment<int> >();
			register_attachment_type<Attachment<uint> >();
			register_attachment_type<Attachment<float> >();
			register_attachment_type<Attachment<double> >();

			register_attachment_type<Attachment<vector1> >();
			register_attachment_type<Attachment<vector2> >();
			register_attachment_type<Attachment<vector3> >();
			register_attachment_type<Attachment<vector4> >();
		}
	////////////////////////////////////////
	//	VARIABLES
		std::vector<std::string>	m_attachmentNames;
		AttachmentMap				m_attachmentMap;
		AttachmentTypeMap			m_attachmentTypeMap;
		FunctionVec					m_functionVecs[4];
};

}//	end of namespace

#endif	//__H__UG_global_attachments
