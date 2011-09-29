//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y08 m11 d13

#include <string>
#include "file_io.h"
#include "lib_grid/algorithms/attachment_util.h"

using namespace std;

namespace ug
{

////////////////////////////////////////////////////////////////////////////////
//	this method performs the actual loading.
static bool LoadGrid3d_IMPL(Grid& grid, ISubsetHandler* pSH,
					   const char* filename, AVector3& aPos)
{
	string strName = filename;
	bool bAutoassignFaces = false;
	bool bSuccess = false;
	if(strName.find(".txt") != string::npos)
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
	else if(strName.find(".msh") != string::npos)
		bSuccess = LoadGridFromMSH(grid, filename, pSH, aPos);

	if(bAutoassignFaces && pSH)
		pSH->assign_subset(grid.faces_begin(), grid.faces_end(), 0);

	return bSuccess;
}


static bool LoadGrid3d(Grid& grid, ISubsetHandler* psh,
					   const char* filename, APosition1& aPos)
{
	APosition aPosTMP;
	grid.attach_to_vertices(aPosTMP);
	if(LoadGrid3d_IMPL(grid, psh, filename, aPosTMP)){
	//	convert the position data from 3d to the required dimension.
		ConvertMathVectorAttachmentValues<VertexBase>(grid, aPosTMP, aPos);
		grid.detach_from_vertices(aPosTMP);
		return true;
	}
	grid.detach_from_vertices(aPosTMP);
	return false;
}

static bool LoadGrid3d(Grid& grid, ISubsetHandler* psh,
					   const char* filename, APosition2& aPos)
{
	APosition aPosTMP;
	grid.attach_to_vertices(aPosTMP);
	if(LoadGrid3d_IMPL(grid, psh, filename, aPosTMP)){
	//	convert the position data from 3d to the required dimension.
		ConvertMathVectorAttachmentValues<VertexBase>(grid, aPosTMP, aPos);
		grid.detach_from_vertices(aPosTMP);
		return true;
	}
	grid.detach_from_vertices(aPosTMP);
	return false;
}

static bool LoadGrid3d(Grid& grid, ISubsetHandler* psh,
					   const char* filename, APosition3& aPos)
{
	return LoadGrid3d_IMPL(grid, psh, filename, aPos);
}

////////////////////////////////////////////////////////////////////////////////
///	This method calls specific load routines or delegates loading to LoadGrid3d
template <class TAPos>
static bool LoadGrid(Grid& grid, ISubsetHandler* psh,
					 const char* filename, TAPos& aPos)
{
//	first all load methods, which do accept template position types are
//	handled. Then all those which only work with 3d position types are processed.
	string strName = filename;
	if(strName.find(".ugx") != string::npos){
		if(psh)
			return LoadGridFromUGX(grid, *psh, filename, aPos);
		else{
		//	we have to create a temporary subset handler, since
			SubsetHandler shTmp(grid);
			return LoadGridFromUGX(grid, shTmp, filename, aPos);
		}
	}
	else{
	//	now we'll handle those methods, which only support 3d position types.
	//	Use a temporary position attachment for dimension 1 and 2.
		return LoadGrid3d(grid, psh, filename, aPos);
	}
}

////////////////////////////////////////////////////////////////////////////////
template <class TAPos>
bool LoadGridFromFile(Grid& grid, ISubsetHandler& sh,
					  const char* filename, TAPos& aPos)
{
	return LoadGrid(grid, &sh, filename, aPos);
}

template <class TAPos>
bool LoadGridFromFile(Grid& grid, const char* filename, TAPos& aPos)
{
	return LoadGrid(grid, NULL, filename, aPos);
}

bool LoadGridFromFile(Grid& grid, ISubsetHandler& sh, const char* filename)
{
	return LoadGrid(grid, &sh, filename, aPosition);
}

bool LoadGridFromFile(Grid& grid, const char* filename)
{
	return LoadGrid(grid, NULL, filename, aPosition);
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//	this method performs the actual save.
static bool SaveGrid3d_IMPL(Grid& grid, SubsetHandler* pSH,
							const char* filename, AVector3& aPos)
{
	string strName = filename;
	if(strName.find(".txt") != string::npos)
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

////////////////////////////////////////////////////////////////////////////////
static bool SaveGrid3d(Grid& grid, SubsetHandler* psh,
					   const char* filename, APosition1& aPos)
{
	APosition aPosTMP;
	grid.attach_to_vertices(aPosTMP);
//	convert the position data from the given dimension to 3d
	ConvertMathVectorAttachmentValues<VertexBase>(grid, aPos, aPosTMP);
	if(SaveGrid3d_IMPL(grid, psh, filename, aPosTMP)){
		grid.detach_from_vertices(aPosTMP);
		return true;
	}
	grid.detach_from_vertices(aPosTMP);
	return false;
}

////////////////////////////////////////////////////////////////////////////////
static bool SaveGrid3d(Grid& grid, SubsetHandler* psh,
					   const char* filename, APosition2& aPos)
{
	APosition aPosTMP;
	grid.attach_to_vertices(aPosTMP);
//	convert the position data from the given dimension to 3d
	ConvertMathVectorAttachmentValues<VertexBase>(grid, aPos, aPosTMP);
	if(SaveGrid3d_IMPL(grid, psh, filename, aPosTMP)){
		grid.detach_from_vertices(aPosTMP);
		return true;
	}
	grid.detach_from_vertices(aPosTMP);
	return false;
}

////////////////////////////////////////////////////////////////////////////////
static bool SaveGrid3d(Grid& grid, SubsetHandler* psh,
					   const char* filename, APosition3& aPos)
{
	return SaveGrid3d_IMPL(grid, psh, filename, aPos);
}

////////////////////////////////////////////////////////////////////////////////
template <class TAPos>
static bool SaveGrid(Grid& grid, SubsetHandler* psh,
					 const char* filename, TAPos& aPos)
{
	string strName = filename;
	if(strName.find(".ugx") != string::npos){
		if(psh)
			return SaveGridToUGX(grid, *psh, filename, aPos);
		else {
			SubsetHandler shTmp(grid);
			return SaveGridToUGX(grid, shTmp, filename, aPos);
		}
	}
	else
		return SaveGrid3d(grid, psh, filename, aPos);
}


////////////////////////////////////////////////////////////////////////////////
template <class TAPos>
bool SaveGridToFile(Grid& grid, SubsetHandler& sh,
					const char* filename, TAPos& aPos)
{
	return SaveGrid3d(grid, &sh, filename, aPos);
}

template <class TAPos>
bool SaveGridToFile(Grid& grid, const char* filename, TAPos& aPos)
{
	return SaveGrid3d(grid, NULL, filename, aPos);
}

bool SaveGridToFile(Grid& grid, SubsetHandler& sh, const char* filename)
{
	return SaveGrid3d(grid, &sh, filename, aPosition);
}

bool SaveGridToFile(Grid& grid, const char* filename)
{
	return SaveGrid3d(grid, NULL, filename, aPosition);
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//	explicit template instantiation
template bool LoadGridFromFile(Grid&, ISubsetHandler&, const char*, AVector1&);
template bool LoadGridFromFile(Grid&, ISubsetHandler&, const char*, AVector2&);
template bool LoadGridFromFile(Grid&, ISubsetHandler&, const char*, AVector3&);

template bool LoadGridFromFile(Grid&, const char*, AVector1&);
template bool LoadGridFromFile(Grid&, const char*, AVector2&);
template bool LoadGridFromFile(Grid&, const char*, AVector3&);

template bool SaveGridToFile(Grid&, SubsetHandler&, const char*, AVector1&);
template bool SaveGridToFile(Grid&, SubsetHandler&, const char*, AVector2&);
template bool SaveGridToFile(Grid&, SubsetHandler&, const char*, AVector3&);

template bool SaveGridToFile(Grid&, const char*, AVector1&);
template bool SaveGridToFile(Grid&, const char*, AVector2&);
template bool SaveGridToFile(Grid&, const char*, AVector3&);

}//	end of namespace
