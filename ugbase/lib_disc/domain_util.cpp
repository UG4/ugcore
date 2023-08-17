/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel
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

#include "domain_util.h"
#include "domain_traits.h"
#include "common/util/string_util.h"
#include "common/util/file_util.h"
#include "lib_grid/file_io/file_io.h"
#include "lib_grid/file_io/file_io_ugx.h"
#include "lib_grid/algorithms/geom_obj_util/misc_util.h"
#include "lib_grid/refinement/projectors/projection_handler.h"
#include "common/profiler/profiler.h"

using namespace std;

namespace ug{

template <typename TDomain>
void LoadDomain(TDomain& domain, const char* filename)
{
	LoadDomain(domain, filename, 0);
}

// Use procId = -2 (i.e., 4294967294 for 32bit int) to achieve loading on all procs.
template <typename TDomain>
void LoadDomain(TDomain& domain, const char* filename, int procId)
{
	PROFILE_FUNC_GROUP("grid");
	size_t num_ph = 0;
	//domain.create_additional_subset_handler("markSH");
	vector<string> additionalSHNames = domain.additional_subset_handler_names();
	SPProjectionHandler ph = make_sp(new ProjectionHandler(domain.geometry3d(), domain.subset_handler()));

	if(additionalSHNames.size()>0){
		SmartPtr<ISubsetHandler> sh;
		vector<SmartPtr<ISubsetHandler>> ash;
		for(size_t i_name = 0; i_name < additionalSHNames.size(); ++i_name){
			SmartPtr<ISubsetHandler> sh = domain.additional_subset_handler(additionalSHNames[i_name]);
			ash[i_name]= sh;
		}

		if(!LoadGridFromFile(*domain.grid(), ph, num_ph, *domain.subset_handler(), additionalSHNames, ash,
								filename, domain.position_attachment(), procId))
		{
			UG_THROW("LoadDomain: Could not load file: "<<filename);
		}
	}
	else if(!LoadGridFromFile(*domain.grid(), *domain.subset_handler(),
								filename, domain.position_attachment(), procId))
	{
		UG_THROW("LoadDomain: Could not load file: "<<filename);
	}
	
	//	declare global attachments on all processors
	#ifdef UG_PARALLEL
		pcl::ProcessCommunicator procComm;
		const vector<string> possible_attachment_names = {"PORO", "PERM", "PERMX", "PERMY", "PERMZ"};
		vector<byte> locDeclared(possible_attachment_names.size(), 0);
		vector<byte> globDeclared(possible_attachment_names.size(), 0);
		// record local info
		for(size_t i = 0; i < possible_attachment_names.size(); ++i){
			byte& b = locDeclared[i];
			if(GlobalAttachments::is_declared(possible_attachment_names[i])){
				b |= 1;
				if(GlobalAttachments::is_attached<Vertex>(*domain.grid(), possible_attachment_names[i]))
					b |= 1<<1;
				if(GlobalAttachments::is_attached<Edge>(*domain.grid(), possible_attachment_names[i]))
					b |= 1<<2;
				if(GlobalAttachments::is_attached<Face>(*domain.grid(), possible_attachment_names[i]))
					b |= 1<<3;
				if(GlobalAttachments::is_attached<Volume>(*domain.grid(), possible_attachment_names[i]))
					b |= 1<<4;
			}
		}
		// sum up all the local to the global
		procComm.allreduce(locDeclared, globDeclared, PCL_RO_BOR);
		// update the local with the global
		for(size_t i = 0; i < possible_attachment_names.size(); ++i){
			byte& b = globDeclared[i];
			if(b & 1){
				if(!GlobalAttachments::is_declared(possible_attachment_names[i]))
					GlobalAttachments::declare_attachment(possible_attachment_names[i], "double", true);
				if(b & 1<<1)
					GlobalAttachments::attach<Vertex>(*domain.grid(), possible_attachment_names[i]);
				if(b & 1<<2)
					GlobalAttachments::attach<Edge>(*domain.grid(), possible_attachment_names[i]);	
				if(b & 1<<3)
					GlobalAttachments::attach<Face>(*domain.grid(), possible_attachment_names[i]);	
				if(b & 1<<4)
					GlobalAttachments::attach<Volume>(*domain.grid(), possible_attachment_names[i]);	
			}
		}
	#endif
	
	if(num_ph > 0)
	{
		domain.set_refinement_projector(ph);
	}
}


template <typename TDomain>
void SaveDomain(TDomain& domain, const char* filename)
{
	PROFILE_FUNC_GROUP("grid");
	if(GetFilenameExtension(string(filename)) == string("ugx")){
		GridWriterUGX ugxWriter;
		ugxWriter.add_grid(*domain.grid(), "defGrid", domain.position_attachment());
		ugxWriter.add_subset_handler(*domain.subset_handler(), "defSH", 0);

		vector<string> additionalSHNames = domain.additional_subset_handler_names();
		for(size_t i_name = 0; i_name < additionalSHNames.size(); ++i_name){
			const char* shName = additionalSHNames[i_name].c_str();
			ugxWriter.add_subset_handler(*domain.additional_subset_handler(shName), shName, 0);
		}

		if(!ugxWriter.write_to_file(filename)){
			UG_THROW("Couldn't save domain to the specified file: " << filename);
		}
	}
	else if(!SaveGridToFile(*domain.grid(), *domain.subset_handler(),
						  filename, domain.position_attachment()))
		UG_THROW("SaveDomain: Could not save to file: "<<filename);
}


template <typename TDomain>
number MaxElementDiameter(TDomain& domain, int level)
{
	typedef typename domain_traits<TDomain::dim>::grid_base_object TElem;
	return  MaxElementDiameter(*domain.grid(), domain.position_accessor(),
	                           domain.grid()->template begin<TElem>(level),
	                           domain.grid()->template end<TElem>(level));
}

template <typename TDomain>
number MinElementDiameter(TDomain& domain, int level)
{
	typedef typename domain_traits<TDomain::dim>::grid_base_object TElem;
	return  MinElementDiameter(*domain.grid(), domain.position_accessor(),
	                           domain.grid()->template begin<TElem>(level),
	                           domain.grid()->template end<TElem>(level));
}

template void LoadDomain<Domain1d>(Domain1d& domain, const char* filename);
template void LoadDomain<Domain2d>(Domain2d& domain, const char* filename);
template void LoadDomain<Domain3d>(Domain3d& domain, const char* filename);

template void LoadDomain<Domain1d>(Domain1d& domain, const char* filename, int procId);
template void LoadDomain<Domain2d>(Domain2d& domain, const char* filename, int procId);
template void LoadDomain<Domain3d>(Domain3d& domain, const char* filename, int procId);

template void SaveDomain<Domain1d>(Domain1d& domain, const char* filename);
template void SaveDomain<Domain2d>(Domain2d& domain, const char* filename);
template void SaveDomain<Domain3d>(Domain3d& domain, const char* filename);

template number MaxElementDiameter<Domain1d>(Domain1d& domain, int level);
template number MaxElementDiameter<Domain2d>(Domain2d& domain, int level);
template number MaxElementDiameter<Domain3d>(Domain3d& domain, int level);

template number MinElementDiameter<Domain1d>(Domain1d& domain, int level);
template number MinElementDiameter<Domain2d>(Domain2d& domain, int level);
template number MinElementDiameter<Domain3d>(Domain3d& domain, int level);

} // end namespace ug

