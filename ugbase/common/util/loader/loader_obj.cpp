/*
 * Copyright (c) 2009-2015:  G-CSC, Goethe University Frankfurt
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

#include <vector>
#include <fstream>
#include <cstring>
#include <string>
#include "loader_obj.h"
#include "loader_util.h"

using namespace std;

namespace ug
{
////////////////////////////////////////////////////////////////////////////////////////////////
//	some string functions
//TODO: IMPROVE THIS FUNCTION!!! MOVE IT TO A COMMON UTIL-FILE
typedef vector<string>		ParameterList;
static void split_parameters(ParameterList* pParamList, const char* pParamString)
{
	string param(pParamString);

	while(param.size())
	{
	//	delete empty spaces at the front of the string
		while(param.find(" ", 0) == 0)
		{
			param.erase(0, 1);
		}
	//	find the next empty space
		int iSpace = param.find(" ", 0);
		if(iSpace == -1)
		{
		//	there is no -> the rest of the string is a parameter
			if(param.size() > 0)
			{
			/*
				SParameter nParam;
				nParam.pFirst = pParamString + PosIndex;
				nParam.size = param.size();
			*/
				string nParam = param;
				pParamList->push_back(nParam);
				break;
			}
		}
		else
		{
		//	extract the parameter from param and push back to ParamList
		/*
			SParameter nParam;
			nParam.pFirst = pParameterString + PosIndex;
			nParam.size = iSpace;
		*/
			string nParam;
			nParam.assign(param, 0, iSpace);
			pParamList->push_back(nParam);
			param.erase(0, iSpace + 1);
		}
	}
}

//TODO: IMPROVE THIS FUNCTION!!! MOVE IT TO A COMMON UTIL-FILE
string replace_chars(string& str, char cToR, char cNew)
{
	char* nStr = new char[str.length() + 1];
	memcpy(nStr, str.c_str(), str.length());
	nStr[str.length()] = 0;
	char* tPtr = nStr;
	while(*tPtr)
	{
		if(*tPtr == cToR)
			*tPtr = cNew;
		tPtr++;
	}
	string Ret = nStr;
	delete[] nStr;
	return Ret;
}

//TODO: IMPROVE THIS FUNCTION!!! MOVE IT TO A COMMON UTIL-FILE
string extract_path(const string& filename)
{
	string::size_type idx1 = filename.rfind("/");
	string::size_type idx2 = filename.rfind("\\");

	if(idx2 < idx1)
		idx1 = idx2;
	return filename.substr(0, idx1+1);
}

/*
**
*/
////////////////////////////////////////////////////////////////////////////////////////////////
//	LoaderObj implementation
LoaderObj::~LoaderObj()
{
	clear();
}

void LoaderObj::clear()
{
	for(LoaderObj::ObjectIterator iter = m_vObjects.begin(); iter != m_vObjects.end(); iter++)
		delete (*iter);
	m_vObjects.clear();
}

bool LoaderObj::load_file(const char* strFilename, bool convertQuadsToTris)
{
	clear();

	ifstream in(strFilename);
	if(!in)
		return false;

	string strObj = "o";
	string strPos = "v";
	string strNorm = "vn";
	string strTex1 = "vt";
	string strFace = "f";
	string strMtrl = "usemtl";
	string strMtrlLib = "mtllib";

	bool bGotPosition = false;
	bool bGotTexture = false;
	bool bGotNormal = false;

	Object* pActiveObject = NULL;

	string strCommand;
	char BUFFER[256];
	BUFFER[255] = 0;
	for(uint i = 0; i < 256; ++i)
		BUFFER[i] = 0;
	ParameterList	lstParams;
	ParameterList::iterator PIter;

	while(!in.eof())
	{
		in.getline(BUFFER, 256);
		if(BUFFER[255] != 0)
			return false;

//	files that come from a windows environment always have '\r' at the end of each line.
//	the following code removes those '\r'. The implementation is far from optimal!!
		/*
		string tStr(BUFFER);
		tStr = replace_chars(tStr, '\r', '\0');

		lstParams.clear();

		split_parameters(&lstParams, tStr.c_str());
		*/
		split_parameters(lstParams, BUFFER, " \r\t");

		if(lstParams.empty())
			continue;

		strCommand = *lstParams.begin();
		lstParams.erase(lstParams.begin());

		if(strCommand == strObj)
		{
			if(pActiveObject)
				m_vObjects.push_back(pActiveObject);

			pActiveObject = new LoaderObj::Object;
			pActiveObject->m_strName = *lstParams.begin();
			pActiveObject->m_iMaterialIndex = -1;
		}
		else if(strCommand == strMtrlLib)
		{//	read materials
			string mtlLib = extract_path(string(strFilename));
			mtlLib.append(*lstParams.begin());

			string newMaterial = "newmtl";
			string colorDiffuse = "Kd";
			string textureDiffuse = "map_Kd";
			string colorAlpha = "d";

			ifstream mtlIn(mtlLib.c_str());
			if(!mtlIn)
				continue;

			Material* pActMaterial = NULL;
			while(!mtlIn.eof())
			{
				mtlIn.getline(BUFFER, 256);
				if(BUFFER[255] != 0)
					return false;

				lstParams.clear();
				split_parameters(&lstParams, BUFFER);
				if(lstParams.empty())
					continue;

				strCommand = *lstParams.begin();
				lstParams.erase(lstParams.begin());

				if(strCommand == newMaterial)
				{
					m_vMaterials.push_back(Material());
					pActMaterial = &m_vMaterials[m_vMaterials.size() - 1];
					if(!lstParams.empty())
						pActMaterial->m_strName = *lstParams.begin();
				}
				else if(pActMaterial)
				{
					if(strCommand == colorDiffuse)
					{
						if(lstParams.size() != 3)
							continue;
						PIter = lstParams.begin();
						pActMaterial->m_vDiffuse.x() = atof((*PIter).c_str());
						PIter++;
						pActMaterial->m_vDiffuse.y() = atof((*PIter).c_str());
						PIter++;
						pActMaterial->m_vDiffuse.z() = atof((*PIter).c_str());
					}
					else if(strCommand == colorAlpha)
					{
						if(lstParams.size() != 1)
							continue;
						PIter = lstParams.begin();
						pActMaterial->m_fAlpha = pActMaterial->m_vDiffuse.w() = atof((*PIter).c_str());
					}
					else if(strCommand == textureDiffuse)
					{
						if(lstParams.size() != 1)
							continue;
						pActMaterial->m_strTextureDiffuse = *lstParams.begin();
					}
				}
			}
		}
		else if(strCommand == strPos)
		{
			bGotPosition = true;
			if(lstParams.size() != 3)
				continue;
			vector3 v;
			PIter = lstParams.begin();
			v.x() = atof((*PIter).c_str());
			PIter++;
			v.y() = atof((*PIter).c_str());
			PIter++;
			v.z() = atof((*PIter).c_str());
			m_vPoints.push_back(v);
		}
		else if(strCommand == strNorm)
		{
			bGotNormal = true;
		}
		else if(strCommand == strTex1)
		{
			bGotTexture = true;
			if(lstParams.size() < 2)
				continue;
			vector2 v;
			PIter = lstParams.begin();
			v.x() = atof((*PIter).c_str());
			PIter++;
			v.y() = -atof((*PIter).c_str());
			m_vTexCoords.push_back(v);
		}
		else if(strCommand == strFace){
			if(!pActiveObject){
				pActiveObject = new LoaderObj::Object;
				pActiveObject->m_strName = "default";
				pActiveObject->m_iMaterialIndex = -1;
			}

			PIter = lstParams.begin();

			if(lstParams.size() == 2)
			{
				for(; PIter != lstParams.end(); PIter++)
				{
					string tStr = replace_chars((*PIter), '/', ' ');
					ParameterList	lstIndexParams;
					split_parameters(&lstIndexParams, tStr.c_str());
					uint iCount = 0;
					if(bGotPosition && (iCount < lstIndexParams.size()))
					{
						pActiveObject->m_vEdgeList.push_back(atoi(lstIndexParams[iCount].c_str()) - 1);
						iCount++;
					}
				}
			}
			else if(lstParams.size() == 3)
			{
				for(; PIter != lstParams.end(); PIter++)
				{
					string tStr = replace_chars((*PIter), '/', ' ');
					ParameterList	lstIndexParams;
					split_parameters(&lstIndexParams, tStr.c_str());

					uint iCount = 0;
					if(bGotPosition)
					{
						pActiveObject->m_vTriangleList.push_back(atoi(lstIndexParams[iCount].c_str()) - 1);
						iCount++;
					}
					if(bGotTexture && iCount < lstIndexParams.size())
					{
						pActiveObject->m_vTriangleListTex.push_back(atoi(lstIndexParams[iCount].c_str()) - 1);
						iCount++;
					}
					if(bGotNormal && iCount < lstIndexParams.size())
					{
						//pActiveObject->m_vTriangleList.push_back(atoi(lstIndexParams[iCount].c_str()) - 1);
						iCount++;
					}


				}
			}
			else if(lstParams.size() == 4)
			{
				unsigned int tInd[4];
				unsigned int tIndTex[4];

				int counter = 0;
				for(; PIter != lstParams.end(); PIter++)
				{
					string tStr = replace_chars((*PIter), '/', ' ');
					ParameterList	lstIndexParams;
					split_parameters(&lstIndexParams, tStr.c_str());
					tInd[counter] = atoi((*lstIndexParams.begin()).c_str());
					uint iCount = 0;
					if(bGotPosition)
					{
						tInd[counter] = atoi(lstIndexParams[iCount].c_str());
						iCount++;
					}
					if(bGotTexture && iCount < lstIndexParams.size())
					{
						tIndTex[counter] = atoi(lstIndexParams[iCount].c_str());
						iCount++;
					}
					if(bGotNormal && iCount < lstIndexParams.size())
					{
						//iIndNorm
						iCount++;
					}
					counter++;
				}
				if(convertQuadsToTris)
				{
					if(bGotPosition)
					{
						pActiveObject->m_vTriangleList.push_back(tInd[0] - 1);
						pActiveObject->m_vTriangleList.push_back(tInd[1] - 1);
						pActiveObject->m_vTriangleList.push_back(tInd[2] - 1);
						pActiveObject->m_vTriangleList.push_back(tInd[0] - 1);
						pActiveObject->m_vTriangleList.push_back(tInd[2] - 1);
						pActiveObject->m_vTriangleList.push_back(tInd[3] - 1);
					}
					if(bGotTexture)
					{
						pActiveObject->m_vTriangleListTex.push_back(tIndTex[0] - 1);
						pActiveObject->m_vTriangleListTex.push_back(tIndTex[1] - 1);
						pActiveObject->m_vTriangleListTex.push_back(tIndTex[2] - 1);
						pActiveObject->m_vTriangleListTex.push_back(tIndTex[0] - 1);
						pActiveObject->m_vTriangleListTex.push_back(tIndTex[2] - 1);
						pActiveObject->m_vTriangleListTex.push_back(tIndTex[3] - 1);
					}
				}
				else
				{
					if(bGotPosition)
					{
						pActiveObject->m_vQuadList.push_back(tInd[0] - 1);
						pActiveObject->m_vQuadList.push_back(tInd[1] - 1);
						pActiveObject->m_vQuadList.push_back(tInd[2] - 1);
						pActiveObject->m_vQuadList.push_back(tInd[3] - 1);
					}
					if(bGotTexture)
					{
						pActiveObject->m_vQuadListTex.push_back(tIndTex[0] - 1);
						pActiveObject->m_vQuadListTex.push_back(tIndTex[1] - 1);
						pActiveObject->m_vQuadListTex.push_back(tIndTex[2] - 1);
						pActiveObject->m_vQuadListTex.push_back(tIndTex[3] - 1);
					}
				}

			}
		}
		else if(strCommand == strMtrl)
		{
			if(!pActiveObject){
				pActiveObject = new LoaderObj::Object;
				pActiveObject->m_strName = "default";
				pActiveObject->m_iMaterialIndex = -1;
			}

			if(lstParams.empty())
			{
				pActiveObject->m_strMaterialName = *lstParams.begin();
				pActiveObject->m_iMaterialIndex = get_material_index_by_name(pActiveObject->m_strMaterialName.c_str());
			}
		}
	}

	if(pActiveObject)
		m_vObjects.push_back(pActiveObject);

	return true;
}

int LoaderObj::get_material_index_by_name(const char* name) const
{
	int counter = 0;
	for(MaterialVector::const_iterator iter = m_vMaterials.begin(); iter != m_vMaterials.end(); iter++)
	{
		if((*iter).m_strName.compare(name) == 0)
			return counter;
		counter++;
	}
	return -1;
}

}
