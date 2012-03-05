/*
 * domain_util.cpp
 *
 *  Created on: 05.03.2012
 *      Author: andreasvogel
 */

#include "domain_util.h"

#include "lib_grid/lib_grid.h"

namespace ug{

template <typename TDomain>
bool LoadDomain(TDomain& domain, const char* filename)
{
	return LoadDomain(domain, filename, 0);
}


template <typename TDomain>
bool LoadDomain(TDomain& domain, const char* filename, int procId)
{
#ifdef UG_PARALLEL
	if((procId >= 0 ) && (pcl::GetProcRank() != procId))
		return true;
#endif

	if(!LoadGridFromFile(*domain.grid(), *domain.subset_handler(),
						 filename, domain.position_attachment()))
	{
		return false;
	}

	domain.update_local_subset_dim_property();

	return true;
}


template <typename TDomain>
bool SaveDomain(TDomain& domain, const char* filename)
{
	return SaveGridToFile(*domain.grid(), *domain.subset_handler(),
						  filename, domain.position_attachment());
}

// writes domain to *.ugx file
template <typename TDomain>
bool WriteDomainToUGX(const char* filename, const TDomain& domain)
{
	// filename
	std::string strName = filename;

	// check filename
	if(strName.find(" ") != std::string::npos)
		{UG_LOG("Filename must not include spaces. Cannot write domain."); return false;}

	// check if filename has already ending (if not add it)
	if(strName.find(".ugx") == std::string::npos)
	{
		if(strName.find(".") != std::string::npos)
		{
			UG_LOG("Filename must not include dots. Cannot write domain.");
			return false;
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
		UG_LOG("WriteDomainToUGX: Cannot save grid.\n");
		return false;
	}

	return true;
}

template bool WriteDomainToUGX<Domain1d>(const char* filename, const Domain1d& domain);
template bool WriteDomainToUGX<Domain2d>(const char* filename, const Domain2d& domain);
template bool WriteDomainToUGX<Domain3d>(const char* filename, const Domain3d& domain);

template bool LoadDomain<Domain1d>(Domain1d& domain, const char* filename);
template bool LoadDomain<Domain2d>(Domain2d& domain, const char* filename);
template bool LoadDomain<Domain3d>(Domain3d& domain, const char* filename);

template bool LoadDomain<Domain1d>(Domain1d& domain, const char* filename, int procId);
template bool LoadDomain<Domain2d>(Domain2d& domain, const char* filename, int procId);
template bool LoadDomain<Domain3d>(Domain3d& domain, const char* filename, int procId);

template bool SaveDomain<Domain1d>(Domain1d& domain, const char* filename);
template bool SaveDomain<Domain2d>(Domain2d& domain, const char* filename);
template bool SaveDomain<Domain3d>(Domain3d& domain, const char* filename);

} // end namespace ug

