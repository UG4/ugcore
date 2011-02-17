// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 09.02.2011 (m,d,y)

#ifndef __H__UG__compol_selection__
#define __H__UG__compol_selection__
#include "pcl/pcl_communication_structs.h"

namespace ug
{


template <class TLayout>
class ComPol_Selection : public pcl::ICommunicationPolicy<TLayout>
{
	public:
		typedef TLayout							Layout;
		typedef typename Layout::Type			GeomObj;
		typedef typename Layout::Element		Element;
		typedef typename Layout::Interface		Interface;
		typedef typename Interface::iterator	InterfaceIter;

	///	Construct the communication policy with a ug::Selector.
	/**	Through the parameters select and deselect one may specify whether
	 * a process selects and/or deselects elements based on the received
	 * selection status.*/
		ComPol_Selection(Selector& sel, bool select = true, bool deselect = true)
			 :	m_sel(sel), m_bSelectAllowed(select), m_bDeselectAllowed(deselect)
		{}

		virtual int
		get_required_buffer_size(Interface& interface)		{return interface.size() * sizeof(byte);}

	///	writes 1 for selected and 0 for unselected interface entries
		virtual bool
		collect(std::ostream& buff, Interface& interface)
		{
			byte zero = 0; byte one = 1;
		//	write the entry indices of marked elements.
			for(InterfaceIter iter = interface.begin();
				iter != interface.end(); ++iter)
			{
				Element elem = interface.get_element(iter);
				if(m_sel.is_selected(elem)){
					buff.write((char*)&one, sizeof(byte));
				}
				else{
					buff.write((char*)&zero, sizeof(byte));
				}
			}

			return true;
		}

	///	reads marks from the given stream
		virtual bool
		extract(std::istream& buff, Interface& interface)
		{
			byte val;
			for(InterfaceIter iter = interface.begin();
				iter != interface.end(); ++iter)
			{
				Element elem = interface.get_element(iter);
				buff.read((char*)&val, sizeof(byte));
				if(val && select_allowed())
					m_sel.select(elem);
				else if(!val && deselect_allowed())
					m_sel.deselect(elem);
			}
			return true;
		}

	protected:
		inline bool select_allowed()	{return m_bSelectAllowed;}
		inline bool deselect_allowed()	{return m_bDeselectAllowed;}

		Selector& m_sel;
		bool m_bSelectAllowed;
		bool m_bDeselectAllowed;
};

}//	end of namespace

#endif
