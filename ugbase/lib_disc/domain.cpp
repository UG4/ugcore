/*
 * Copyright (c) 2013:  G-CSC, Goethe University Frankfurt
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

#include "domain.h"
#include "common/util/table.h"

#ifdef UG_PARALLEL
	#include "lib_grid/refinement/projectors/projectors.h"
	#include "common/boost_serialization_routines.h"
	#include "common/util/archivar.h"
	#include "common/util/factory.h"
	#include <boost/archive/text_oarchive.hpp>
	#include <boost/archive/text_iarchive.hpp>
#endif


using namespace std;

namespace ug{

std::string DomainInfo::
to_string() const
{
	if((m_numElems.size() != m_numLocalElems.size())
	   || (m_numElems.size() != m_numLocalGhosts.size()))
	{
		UG_THROW("elem-arrays have to have the same number of entries!");
	}

	StringStreamTable t;

	t(0, 0) << "lvl";
	t(0, 1) << "#total-elems";
	t(0, 2) << "#local-elems";
	t(0, 3) << "(% of total)";
	t(0, 4) << "#local-ghosts";
	t(0, 5) << "#min-local-elems";
	t(0, 6) << "#max-local-elems";

	for(size_t i = 0; i < m_numElems.size(); ++i){
		int r = i+1;
		t(r, 0) << i;
		t(r, 1) << m_numElems[i];
		t(r, 2) << m_numLocalElems[i];
		if(m_numElems[i] > 0)
			t(r, 3) << (float)m_numLocalElems[i] / (float)m_numElems[i];
		else
			t(r, 3) << "-";
		t(r, 4) << m_numLocalGhosts[i];
		t(r, 5) << m_minNumLocalElems[i];
		t(r, 6) << m_maxNumLocalElems[i];
	}

	return t.to_string();
}


#ifdef UG_PARALLEL
namespace detail{

SPRefinementProjector
BroadcastRefinementProjector(
		int rootProc,
		pcl::ProcessCommunicator& procCom,
		SmartPtr<ISubsetHandler> subsetHandler,
		SPIGeometry3d geometry,
		SPRefinementProjector projector)
{
	BinaryBuffer 	buf;
	const int		magicNumber	= 3243578;
	const bool 		isRoot		= (pcl::ProcRank() == rootProc);

	static Factory<RefinementProjector, ProjectorTypes>	projFac;

	if(isRoot){
		Archivar<boost::archive::text_oarchive,
				RefinementProjector,
				ProjectorTypes>
			archivar;

	//	if the specified projector is a projection handler, we'll perform a
	//	special operation.
		ProjectionHandler* ph = NULL;
		int projectorType = -1;// -1: none, 0: normal projector, 1: projection handler
		if(projector.valid()){
			ph = dynamic_cast<ProjectionHandler*>(projector.get());
			if(ph)
				projectorType = 1;
			else
				projectorType = 0;
		}

		Serialize(buf, projectorType);
		if(ph){
			size_t numProjectors = ph->num_projectors();
			Serialize(buf, numProjectors);
			for(size_t iproj = 0; iproj < numProjectors; ++iproj){
				SPRefinementProjector	proj		= ph->projector(iproj);
				const string&			projName 	= projFac.class_name(*proj);
				Serialize(buf, projName);

				stringstream ss;
				boost::archive::text_oarchive ar(ss, boost::archive::no_header);
				archivar.archive(ar, *proj);
				Serialize(buf, ss.str());
			}
		}
		else if(projector.valid()){
			RefinementProjector&	proj		= *projector;
			const string&			projName 	= projFac.class_name(proj);
			Serialize(buf, projName);

			stringstream ss;
			boost::archive::text_oarchive ar(ss, boost::archive::no_header);
			archivar.archive(ar, proj);
			Serialize(buf, ss.str());
		}

		Serialize(buf, magicNumber);
	}

	procCom.broadcast(buf, rootProc);

	if(!isRoot){
		Archivar<boost::archive::text_iarchive,
				RefinementProjector,
				ProjectorTypes>
			archivar;

		int projectorType;
		Deserialize(buf, projectorType);
		if(projectorType == 1){
			ProjectionHandler* ph = new ProjectionHandler(geometry, subsetHandler);
			SPProjectionHandler projHandler = make_sp(ph);

			size_t numProjectors;
			Deserialize(buf, numProjectors);

			for(size_t iproj = 0; iproj < numProjectors; ++iproj){
				string name;
				Deserialize(buf, name);
				SPRefinementProjector proj = projFac.create(name);

				string data;
				Deserialize(buf, data);
				stringstream ss(data, ios_base::in);
				boost::archive::text_iarchive ar(ss, boost::archive::no_header);
				archivar.archive(ar, *proj);

				ph->set_projector(iproj, proj);
			}

			projector = projHandler;
		}
		else if(projectorType == 0){
			string name;
			Deserialize(buf, name);
			SPRefinementProjector proj = projFac.create(name);

			string data;
			Deserialize(buf, data);
			stringstream ss(data, ios_base::in);
			boost::archive::text_iarchive ar(ss, boost::archive::no_header);
			archivar.archive(ar, *proj);

			projector = proj;
		}
		else if(projectorType == -1){
			projector = SPNULL;
		}
		else{
			UG_THROW("Invalid projector type in 'BroadcastRefinementProjector': "
					 << projectorType);
		}

		int tmp;
		Deserialize(buf, tmp);
		UG_COND_THROW(tmp != magicNumber, "Magic number mismatch in "
					  "'BroadcastRefinementProjector'. Received "
					  << tmp << ", but expected " << magicNumber);
	}

	return projector;
}
}// end of namespace detail

#endif
}//	end of namespace
