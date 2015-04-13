// created by Sebastian Reiter
// s.b.reiter@gmail.com

#ifndef __H__UG_attachment_io_handler
#define __H__UG_attachment_io_handler

#include "lib_grid/attachments/attachment_info_traits.h"
#include "lib_grid/attachments/attachment_io_traits.h"

namespace ug{


template <class TElem, class TAttachment>
void ReadAttachmentFromStream (
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
void WriteAttachmentToStream (
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



class AttachmentIOHandler {
	public:
		~AttachmentIOHandler ()
		{
			clear();
		}

	/**	\param	attachment	an attachment into which data shall be read or
	 *						from which data shall be written.
	 *						The attachment either has to exist longer than
	 *						the AttachmentIOHandler instance at which it is
	 *						registered, or an explicit call to
	 *						AttachmentIOHandler::unregister_attachment has to
	 *						be performed.
	 *
	 * \param name			The name of the attachment as it will appear in
	 *						a file
	 * \{ */
		template <class TElem, class TAttachment>
		void register_attachment (TAttachment& a, const std::string& name)
		{
			Entry e;
			e.attachment = new TAttachment(a);
			e.readFunc = &ReadAttachmentFromStream<TElem, TAttachment>;
			e.writeFunc = &WriteAttachmentToStream<TElem, TAttachment>;
			e.typeName = attachment_info_traits<TAttachment>::type_name();
			get_map<TElem>()[name] = e;
		}


		template <class TElem>
		void unregister_attachment (IAttachment& a)
		{
			EntryMap& m = get_map<TElem>();
			for(EntryMap::iterator iter = m.begin(); iter != m.end(); ++iter)
			{
				if((*iter->second.attachment) == a){
					delete iter->second.attachment;
					m.erase(iter);
					break;
				}
			}
		}


		void clear ()
		{
			for(int imap = 0; imap < NUM_GEOMETRIC_BASE_OBJECTS; ++imap){
				for(EntryMap::iterator iter = m_maps[imap].begin();
					iter != m_maps[imap].end(); ++iter)
				{
					delete iter->second.attachment;
				}
				m_maps[imap] = EntryMap();
			}
		}

		template <class TElem>
		bool attachment_is_registered(const std::string& name)
		{
			EntryMap& m = get_map<TElem>();
			return m.find(name) != m.end();
		}

		template <class TElem>
		size_t registered_attachments(std::vector<std::string>& namesOut)
		{
			namesOut.clear();
			EntryMap& m = get_map<TElem>();
			for(EntryMap::iterator iter = m.begin(); iter != m.end(); ++iter)
			{
				namesOut.push_back(iter->first);
			}
			return namesOut.size();
		}

		template <class TElem>
		const char* type_name (const std::string& name)
		{
			EntryMap& m = get_map<TElem>();
			EntryMap::iterator iter = m.find(name);
			if(iter != m.end()){
				return iter->second.typeName.c_str();
			}
			UG_THROW("Can't evaluate type-name of unregistered attachment "
					 << name << " for element-type: " << TElem::BASE_OBJECT_ID);
			return NULL;
		}

		template <class TElem>
		void read_attachment_values (std::istream& in,
									 Grid& grid,
									 const std::string& name)
		{
			EntryMap& m = get_map<TElem>();
			EntryMap::iterator iter = m.find(name);
			if(iter != m.end()){
				iter->second.readFunc(in, grid, *iter->second.attachment);
			}
		}


		template <class TElem>
		void write_attachment_values (std::ostream& out,
									  Grid& grid,
									  const std::string& name)
		{
			EntryMap& m = get_map<TElem>();
			EntryMap::iterator iter = m.find(name);
			if(iter != m.end()){
				iter->second.writeFunc(out, grid, *iter->second.attachment);
			}
		}

	private:
	////////////////////////
	//	TYPES
		struct Entry{
			IAttachment*	attachment;
			std::string		typeName;
			void (*			readFunc	)	(std::istream&, Grid&, IAttachment&);
			void (*			writeFunc	)	(std::ostream&, Grid&, IAttachment&);
		};

		typedef std::map<std::string, Entry> EntryMap;


	////////////////////////
	//	FUNCTIONS
		template <class TElem>
		EntryMap& get_map()		{return m_maps[TElem::BASE_OBJECT_ID];}


	////////////////////////
	//	VARIABLES
		EntryMap	m_maps[NUM_GEOMETRIC_BASE_OBJECTS];

};

}//	end of namespace

#endif	//__H__UG_attachment_io_handler
