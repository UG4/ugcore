/*
 * domain_util.cpp
 *
 *  Created on: 05.03.2012
 *      Author: andreasvogel
 */

#include "domain_util.h"
#include "lib_grid/file_io/file_io.h"

namespace ug{

template <typename TDomain>
void LoadDomain(TDomain& domain, const char* filename)
{
	LoadDomain(domain, filename, 0);
}


template <typename TDomain>
void LoadDomain(TDomain& domain, const char* filename, int procId)
{
#ifdef UG_PARALLEL
	if((procId >= 0 ) && (pcl::GetProcRank() != procId))
		return;
#endif

	if(!LoadGridFromFile(*domain.grid(), *domain.subset_handler(),
						 filename, domain.position_attachment()))
		UG_THROW("LoadDomain: Could not load file: "<<filename);

	domain.update_local_subset_dim_property();
}


template <typename TDomain>
void SaveDomain(TDomain& domain, const char* filename)
{
	if(!SaveGridToFile(*domain.grid(), *domain.subset_handler(),
						  filename, domain.position_attachment()))
		UG_THROW("SaveDomain: Could not save to file: "<<filename);
}

// writes domain to *.ugx file
template <typename TDomain>
void WriteDomainToUGX(const char* filename, const TDomain& domain)
{
	// filename
	std::string strName = filename;

	// check filename
	if(strName.find(" ") != std::string::npos)
		UG_THROW("Filename must not include spaces. Cannot write domain.");

	// check if filename has already ending (if not add it)
	if(strName.find(".ugx") == std::string::npos)
	{
		if(strName.find(".") != std::string::npos)
		{
			UG_THROW("Filename must not include dots. Cannot write domain.");
		}
		else
		{
			strName = strName + ".ugx";
		}
	}

	TDomain* pDomain = const_cast<TDomain*>(&domain);

	// save grid
	if(!SaveGridToUGX(*pDomain->grid(), *pDomain->subset_handler(), strName.c_str()))
	{
		UG_THROW("WriteDomainToUGX: Cannot save grid.");
	}
}

template void WriteDomainToUGX<Domain1d>(const char* filename, const Domain1d& domain);
template void WriteDomainToUGX<Domain2d>(const char* filename, const Domain2d& domain);
template void WriteDomainToUGX<Domain3d>(const char* filename, const Domain3d& domain);

template void LoadDomain<Domain1d>(Domain1d& domain, const char* filename);
template void LoadDomain<Domain2d>(Domain2d& domain, const char* filename);
template void LoadDomain<Domain3d>(Domain3d& domain, const char* filename);

template void LoadDomain<Domain1d>(Domain1d& domain, const char* filename, int procId);
template void LoadDomain<Domain2d>(Domain2d& domain, const char* filename, int procId);
template void LoadDomain<Domain3d>(Domain3d& domain, const char* filename, int procId);

template void SaveDomain<Domain1d>(Domain1d& domain, const char* filename);
template void SaveDomain<Domain2d>(Domain2d& domain, const char* filename);
template void SaveDomain<Domain3d>(Domain3d& domain, const char* filename);

} // end namespace ug

