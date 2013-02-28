// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 25.7.2011 (m,d,y)

#ifndef __H__UG__compol_subset__
#define __H__UG__compol_subset__
#include "pcl/pcl_communication_structs.h"

namespace ug
{


template <class TLayout>
class ComPol_Subset : public pcl::ICommunicationPolicy<TLayout>
{
	public:
		typedef TLayout								Layout;
		typedef typename Layout::Type				GeomObj;
		typedef typename Layout::Element			Element;
		typedef typename Layout::Interface			Interface;
		typedef typename Interface::const_iterator	InterfaceIter;

	///	Construct the communication policy with a ug::SubsetHandler.
	/**	If overwrite is specified, the subset index in the target handler is
	 * simply overwritten (default is false).*/
		ComPol_Subset(ISubsetHandler& sel, bool overwrite = false)
			 :	m_sh(sel), m_overwriteEnabled(overwrite)
		{}

		virtual int
		get_required_buffer_size(const Interface& interface)
		{
			return interface.size() * sizeof(int);
		}

	///	writes 1 for selected and 0 for unassigned interface entries
		virtual bool
		collect(ug::BinaryBuffer& buff, const Interface& interface)
		{
		//	write the entry indices of marked elements.
			for(InterfaceIter iter = interface.begin();
				iter != interface.end(); ++iter)
			{
				Element elem = interface.get_element(iter);
				int si = m_sh.get_subset_index(elem);
				buff.write((char*)&si, sizeof(int));
			}

			return true;
		}

	///	reads marks from the given stream
		virtual bool
		extract(ug::BinaryBuffer& buff, const Interface& interface)
		{
			int nsi;
			bool retVal = true;
			for(InterfaceIter iter = interface.begin();
				iter != interface.end(); ++iter)
			{
				Element elem = interface.get_element(iter);
				buff.read((char*)&nsi, sizeof(int));
				if(m_overwriteEnabled){
					m_sh.assign_subset(elem, nsi);
				}
				else{
					if(m_sh.get_subset_index(elem) == -1){
						m_sh.assign_subset(elem, nsi);
					}
					else if(m_sh.get_subset_index(elem) != nsi){
					//	if the subset indices do not match, we have a problem here.
						retVal = false;
					}
				}
			}
			return retVal;
		}

	protected:
		ISubsetHandler&	m_sh;
		bool m_overwriteEnabled;
};

}//	end of namespace

#endif
