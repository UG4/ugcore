/*
 * domain_util.cpp
 *
 *  Created on: 05.03.2012
 *      Author: andreasvogel
 */

#include "domain_util.h"
#include "domain_traits.h"
#include "common/util/string_util.h"
#include "lib_grid/file_io/file_io.h"
#include "lib_grid/algorithms/geom_obj_util/misc_util.h"

using namespace std;

namespace ug{

template <typename TDomain>
void LoadDomain(TDomain& domain, const char* filename)
{
	LoadDomain(domain, filename, 0);
}


template <typename TDomain>
void LoadDomain(TDomain& domain, const char* filename, int procId)
{
	if(GetFilenameExtension(string(filename)) == string("ugx")){
		GridReaderUGX ugxReader;
		if(!ugxReader.parse_file(filename)){
			UG_THROW("ERROR in LoadDomain: File not found: " << filename);
		}

		if(ugxReader.num_grids() < 1){
			UG_THROW("ERROR in LoadGridFromUGX: File contains no grid.");
		}

		ugxReader.grid(*domain.grid(), 0, domain.position_attachment());

		if(ugxReader.num_subset_handlers(0) > 0)
			ugxReader.subset_handler(*domain.subset_handler(), 0, 0);

		vector<string> additionalSHNames = domain.additional_subset_handler_names();
		for(size_t i_name = 0; i_name < additionalSHNames.size(); ++i_name){
			string shName = additionalSHNames[i_name];
			for(size_t i_sh = 0; i_sh < ugxReader.num_subset_handlers(0); ++i_sh){
				if(shName == ugxReader.get_subset_handler_name(0, i_sh)){
					ugxReader.subset_handler(*domain.additional_subset_handler(shName), i_sh, 0);
				}
			}
		}
	}
	else if(!LoadGridFromFile(*domain.grid(), *domain.subset_handler(),
						 filename, domain.position_attachment(), procId))
	{
		UG_THROW("LoadDomain: Could not load file: "<<filename);
	}
}


template <typename TDomain>
void SaveDomain(TDomain& domain, const char* filename)
{
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
	typedef typename domain_traits<TDomain::dim>::geometric_base_object TElem;
	return  MaxElementDiameter(*domain.grid(), domain.position_accessor(),
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

} // end namespace ug

