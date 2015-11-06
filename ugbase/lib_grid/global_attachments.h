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

namespace ug{

///	Global attachments are automatically read/written from/to files and are considered during redistribution
class GlobalAttachments {
	public:
		template <class TAttachment>
		static void declare_attachment (const std::string& name,
										bool passOnBehaviour = false)
		{
			const char* typeName = attachment_info_traits<TAttachment>::type_name();
			if(attachments()[name].attachment != NULL) {
				UG_COND_THROW(
					dynamic_cast<TAttachment*> (attachments()[name].attachment) == NULL,
						  "Attachment with name '" << name
						  << "' was already declared in GlobalAttachments with a different type. "
						  << "Old type: " << attachments()[name].type <<
						  ", new type: " << typeName);
				return;
			}

			int fi = static_cast<int> (attachment_names().size());
			attachment_names().push_back(name);
			attachments()[name] =
				AttachmentEntry(new TAttachment(passOnBehaviour), typeName, fi);
			functions<Vertex>().push_back(FunctionEntry<Vertex, TAttachment>());
			functions<Edge>().push_back(FunctionEntry<Edge, TAttachment>());
			functions<Face>().push_back(FunctionEntry<Face, TAttachment>());
			functions<Volume>().push_back(FunctionEntry<Volume, TAttachment>());
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
			TAttachment* a = dynamic_cast<TAttachment*>(e.attachment);
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
			AttachmentEntry () : attachment(NULL), type(""), functionIndex(-1) {}
			AttachmentEntry (IAttachment* a, const char* t, int fi) :
				attachment(a), type(t), functionIndex(fi) {}
			IAttachment* 	attachment;
			const char*		type;
			int				functionIndex;
		};

		struct IFunctionEntry {
			IFunctionEntry() : readFunc(0), writeFunc(0), addSerializer(0), attach(0) {}
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

		typedef std::map<std::string, AttachmentEntry>	AttachmentMap;
		typedef std::vector<IFunctionEntry>				FunctionVec;


	////////////////////////////////////////
	//	METHODS
		static
		GlobalAttachments& inst ()
		{
			static GlobalAttachments h;
			return h;
		}

		GlobalAttachments ()
		{}

		~GlobalAttachments ()
		{
			AttachmentMap& m = attachments();
			for(AttachmentMap::iterator i = m.begin(); i != m.end(); ++i) {
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

		template <class TElem, class TAttachment>
		static
		void read_attachment_from_stream (
				std::istream& in,
				Grid& grid,
				IAttachment& attachment)
		{
			TAttachment& a = dynamic_cast<TAttachment&>(attachment);
			
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
			TAttachment& a = dynamic_cast<TAttachment&>(attachment);
			
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
			TAttachment& a = dynamic_cast<TAttachment&>(attachment);
			handler.add(GeomObjAttachmentSerializer<TElem, TAttachment>::
									create(g, a));
		}

		template <class TElem, class TAttachment>
		static
		void cast_and_attach (
				Grid& grid,
				IAttachment& attachment)
		{
			TAttachment& a = dynamic_cast<TAttachment&>(attachment);
			
			if(!grid.has_attachment<TElem>(a))
				grid.attach_to<TElem>(a);
		}

	////////////////////////////////////////
	//	VARIABLES
		std::vector<std::string>	m_attachmentNames;
		AttachmentMap				m_attachmentMap;
		FunctionVec					m_functionVecs[4];
};

}//	end of namespace

#endif	//__H__UG_global_attachments
