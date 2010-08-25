//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y08 m11 d13

#include <string>
#include "file_io.h"

using namespace std;

namespace ug
{

//	this method performs the actual loading.
static bool LoadGrid(Grid& grid, const char* filename,
					AVector3& aPos,
					ISubsetHandler* pSH)
{
	string strName = filename;
	bool bAutoassignFaces = false;
	bool bSuccess = false;
	if(strName.find(".ugx") != string::npos){
		if(pSH){
			bSuccess = LoadGridFromUGX(grid, *pSH, filename, aPos);
		}
		else {
			SubsetHandler emptySH(grid);
			bSuccess = LoadGridFromUGX(grid, emptySH, filename, aPos);
		}

	}
	else if(strName.find(".txt") != string::npos)
	{
		bAutoassignFaces = true;
		bSuccess = LoadGridFromTXT(grid, filename, aPos);
	}
	else if(strName.find(".obj") != string::npos)
		bSuccess = LoadGridFromOBJ(grid, filename, aPos, NULL, pSH);
	else if(strName.find(".lgb") != string::npos)
	{
		SubsetHandler* shGrid = dynamic_cast<SubsetHandler*>(pSH);
		int numSHs = 0;
		if(shGrid)
			numSHs = 1;

		bSuccess = LoadGridFromLGB(grid, filename, &shGrid, numSHs, aPos);
	}
	else if(strName.find(".net") != string::npos)
		bSuccess = LoadGridFromART(grid, filename, pSH, aPos);
	else if(strName.find(".art") != string::npos)
		bSuccess = LoadGridFromART(grid, filename, pSH, aPos);
	else if(strName.find(".dat") != string::npos)
		bSuccess = LoadGridFromART(grid, filename, pSH, aPos);
	else if(strName.find(".lgm") != string::npos)
		bSuccess = ImportGridFromLGM(grid, filename, aPos, pSH);
	else if(strName.find(".ng") != string::npos)
		bSuccess = ImportGridFromNG(grid, filename, aPos, pSH);
	else if(strName.find(".dump") != string::npos)
	{
		bSuccess = LoadGridFromDUMP(grid, filename, pSH, aPos);
	}
	else if(strName.find(".ele") != string::npos)
		return LoadGridFromELE(grid, filename, pSH, aPos);

	if(bAutoassignFaces && pSH)
		pSH->assign_subset(grid.faces_begin(), grid.faces_end(), 0);

	return bSuccess;
}

//	this method performs the actual save.
static bool SaveGrid(Grid& grid, const char* filename,
					AVector3& aPos,
					SubsetHandler* pSH)
{
	string strName = filename;
	if(strName.find(".ugx") != string::npos){
		SubsetHandler* shGrid = dynamic_cast<SubsetHandler*>(pSH);
		if(shGrid)
			return SaveGridToUGX(grid, *shGrid, filename, aPos);
		else {
			SubsetHandler emptySH(grid);
			return SaveGridToUGX(grid, emptySH, filename, aPos);
		}
	}
	else if(strName.find(".txt") != string::npos)
		return SaveGridToTXT(grid, filename, aPos);
	else if(strName.find(".obj") != string::npos)
		return SaveGridToOBJ(grid, filename, aPos, NULL, pSH);
	else if(strName.find(".lgb") != string::npos)
	{
		int numSHs = 0;
		if(pSH)
			numSHs = 1;

		return SaveGridToLGB(grid, filename, &pSH, numSHs, aPos);
	}
	else if(strName.find(".ele") != string::npos)
		return SaveGridToELE(grid, filename, pSH, aPos);
	else if(strName.find(".net") != string::npos)
		return SaveGridToART(grid, filename, pSH, aPos);
	else if(strName.find(".art") != string::npos)
		return SaveGridToART(grid, filename, pSH, aPos);
	else if(strName.find(".ncdf") != string::npos)
		return SaveGridToNCDF(grid, filename, pSH, aPos);

	return false;
}


////////////////////////////////////////////////////////////////////////
bool LoadGridFromFile(Grid& grid, const char* filename, AVector3& aPos)
{
	return LoadGrid(grid, filename, aPos, NULL);
}

////////////////////////////////////////////////////////////////////////
bool SaveGridToFile(Grid& grid, const char* filename, AVector3& aPos)
{
	return SaveGrid(grid, filename, aPos, NULL);
}

////////////////////////////////////////////////////////////////////////
bool LoadGridFromFile(Grid& grid, const char* filename,
						ISubsetHandler& sh, AVector3& aPos)
{
	return LoadGrid(grid, filename, aPos, &sh);
}

////////////////////////////////////////////////////////////////////////
bool SaveGridToFile(Grid& grid, const char* filename,
					SubsetHandler& sh, AVector3& aPos)
{
	return SaveGrid(grid, filename, aPos, &sh);
}

}//	end of namespace
