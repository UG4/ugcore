#ifndef __H__UG__compol_selection__
#define __H__UG__compol_selection__
#include "pcl/pcl_communication_structs.h"

namespace ug
{


template <class TLayout>
class ComPol_Selection : public pcl::ICommunicationPolicy<TLayout>
{
	public:
		typedef TLayout								Layout;
		typedef typename Layout::Type				GeomObj;
		typedef typename Layout::Element			Element;
		typedef typename Layout::Interface			Interface;
		typedef typename Interface::const_iterator	InterfaceIter;

	///	Construct the communication policy with a ug::Selector.
	/**	Through the parameters select and deselect one may specify whether
	 * a process selects and/or deselects elements based on the received
	 * selection status.*/
		ComPol_Selection(ISelector& sel, bool select = true, bool deselect = true)
			 :	m_sel(sel), m_bSelectAllowed(select), m_bDeselectAllowed(deselect)
		{}

		virtual int
		get_required_buffer_size(const Interface& interface)
		{return interface.size() * sizeof(byte);}

	///	writes 1 for selected and 0 for unselected interface entries
		virtual bool
		collect(ug::BinaryBuffer& buff, const Interface& interface)
		{
		//	write the entry indices of marked elements.
			for(InterfaceIter iter = interface.begin();
				iter != interface.end(); ++iter)
			{
				Element elem = interface.get_element(iter);
				byte refMark = m_sel.get_selection_status(elem);
				buff.write((char*)&refMark, sizeof(byte));
			}

			return true;
		}

	///	reads marks from the given stream
		virtual bool
		extract(ug::BinaryBuffer& buff, const Interface& interface)
		{
			byte val;
			for(InterfaceIter iter = interface.begin();
				iter != interface.end(); ++iter)
			{
				Element elem = interface.get_element(iter);
				buff.read((char*)&val, sizeof(byte));
				if(val && select_allowed())
					m_sel.select(elem, val);
				else if(!val && deselect_allowed())
					m_sel.deselect(elem);
			}
			return true;
		}

	protected:
		inline bool select_allowed()	{return m_bSelectAllowed;}
		inline bool deselect_allowed()	{return m_bDeselectAllowed;}

		ISelector& m_sel;
		bool m_bSelectAllowed;
		bool m_bDeselectAllowed;
};



////////////////////////////////////////////////////////////////////////////////
/**	Enables only the part of the selection-state for which the specified stateBits
 * hold 1.
 */
template <class TLayout>
class ComPol_EnableSelectionStateBits : public pcl::ICommunicationPolicy<TLayout>
{
	public:
		typedef TLayout								Layout;
		typedef typename Layout::Type				GeomObj;
		typedef typename Layout::Element			Element;
		typedef typename Layout::Interface			Interface;
		typedef typename Interface::const_iterator	InterfaceIter;

		ComPol_EnableSelectionStateBits(ISelector& sel, byte stateBits)
			 :	m_sel(sel), m_stateBits(stateBits)
		{}

		virtual int
		get_required_buffer_size(const Interface& interface)
		{return interface.size() * sizeof(byte);}

	///	writes writes the selection states of the interface entries
		virtual bool
		collect(ug::BinaryBuffer& buff, const Interface& interface)
		{
		//	write the entry indices of marked elements.
			for(InterfaceIter iter = interface.begin();
				iter != interface.end(); ++iter)
			{
				Element elem = interface.get_element(iter);
				byte refMark = m_sel.get_selection_status(elem);
				buff.write((char*)&refMark, sizeof(byte));
			}

			return true;
		}

	///	reads marks from the given stream
		virtual bool
		extract(ug::BinaryBuffer& buff, const Interface& interface)
		{
			byte val;
			for(InterfaceIter iter = interface.begin();
				iter != interface.end(); ++iter)
			{
				Element elem = interface.get_element(iter);
				buff.read((char*)&val, sizeof(byte));
				m_sel.select(elem, m_sel.get_selection_status(elem)
								   | (val & m_stateBits));
			}
			return true;
		}

		ISelector& m_sel;
		byte m_stateBits;
};

}//	end of namespace

#endif
