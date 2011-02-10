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

		ComPol_Selection(Selector& sel)
			 :	m_sel(sel)
		{}

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
				if(m_sel.is_selected(elem))
					buff.write((char*)&one, sizeof(byte));
				else
					buff.write((char*)&zero, sizeof(byte));
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
				if(val)
					m_sel.select(elem);
				else
					m_sel.deselect(elem);
			}

			return true;
		}

	protected:
		Selector& m_sel;
};

}//	end of namespace

#endif
