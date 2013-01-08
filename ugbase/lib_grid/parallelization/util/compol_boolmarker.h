// created by Andreas Vogel

#ifndef __H__UG__compol_boolmarker__
#define __H__UG__compol_boolmarker__

#include "pcl/pcl_communication_structs.h"
#include "lib_grid/tools/bool_marker.h"

namespace ug
{

///	adds marking at extracting side
template <class TLayout>
class ComPol_BoolMarker_AddMarks : public pcl::ICommunicationPolicy<TLayout>
{
	public:
		typedef TLayout							Layout;
		typedef typename Layout::Type			GeomObj;
		typedef typename Layout::Element		Element;
		typedef typename Layout::Interface		Interface;
		typedef typename Interface::iterator	InterfaceIter;

	///	Construct the communication policy with a ug::BoolMarker.
		ComPol_BoolMarker_AddMarks(BoolMarker& marker)
			 :	m_rMarker(marker)
		{}

		virtual int
		get_required_buffer_size(Interface& interface)		{return interface.size() * sizeof(byte);}

	///	writes 1 for marked and 0 for unmarked interface entries
		virtual bool
		collect(ug::BinaryBuffer& buff, Interface& interface)
		{
		//	write the entry indices of marked elements.
			for(InterfaceIter iter = interface.begin();
				iter != interface.end(); ++iter)
			{
				Element elem = interface.get_element(iter);
				bool isMarked = m_rMarker.is_marked(elem);
				byte refMark = (isMarked ? 1 : 0);
				buff.write((char*)&refMark, sizeof(byte));
			}

			return true;
		}

	///	reads marks from the given stream
		virtual bool
		extract(ug::BinaryBuffer& buff, Interface& interface)
		{
			byte val;
			for(InterfaceIter iter = interface.begin();
				iter != interface.end(); ++iter)
			{
				Element elem = interface.get_element(iter);
				buff.read((char*)&val, sizeof(byte));
				if(val)	m_rMarker.mark(elem);
			}
			return true;
		}

	protected:
		BoolMarker& m_rMarker;
};

///	removes marks at extracting side, if no mark was received
template <class TLayout>
class ComPol_BoolMarker_RemoveMarks : public pcl::ICommunicationPolicy<TLayout>
{
	public:
		typedef TLayout							Layout;
		typedef typename Layout::Type			GeomObj;
		typedef typename Layout::Element		Element;
		typedef typename Layout::Interface		Interface;
		typedef typename Interface::iterator	InterfaceIter;

	///	Construct the communication policy with a ug::BoolMarker.
		ComPol_BoolMarker_RemoveMarks(BoolMarker& marker)
			 :	m_rMarker(marker)
		{}

		virtual int
		get_required_buffer_size(Interface& interface)		{return interface.size() * sizeof(byte);}

	///	writes 1 for marked and 0 for unmarked interface entries
		virtual bool
		collect(ug::BinaryBuffer& buff, Interface& interface)
		{
		//	write the entry indices of marked elements.
			for(InterfaceIter iter = interface.begin();
				iter != interface.end(); ++iter)
			{
				Element elem = interface.get_element(iter);
				bool isMarked = m_rMarker.is_marked(elem);
				byte refMark = (isMarked ? 1 : 0);
				buff.write((char*)&refMark, sizeof(byte));
			}

			return true;
		}

	///	reads marks from the given stream
		virtual bool
		extract(ug::BinaryBuffer& buff, Interface& interface)
		{
			byte val;
			for(InterfaceIter iter = interface.begin();
				iter != interface.end(); ++iter)
			{
				Element elem = interface.get_element(iter);
				buff.read((char*)&val, sizeof(byte));
				if(!val)	m_rMarker.unmark(elem);
			}
			return true;
		}

	protected:
		BoolMarker& m_rMarker;
};

}//	end of namespace

#endif
