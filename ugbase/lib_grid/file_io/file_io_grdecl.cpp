/*
 * Copyright (c) 2009-2015:  G-CSC, Goethe University Frankfurt
 * Author: Shuai Lu
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

#include "file_io_grdecl.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "../lg_base.h"
#include "common/util/string_util.h"
#include "lib_grid/global_attachments.h"

using namespace std;

namespace ug {

	struct xy
	{
	   double x, y;
	};
	
	struct xyz
	{
	   double x, y, z;
	};

	struct ab
	{
		size_t a1, a2, a3, a4, b1, b2, b3, b4;
	};


void GetDim(string str, vector<size_t>& dim)
{
    string word = "";
    for (auto x : str) 
    {
        if (x == ' ')
        {
			if (word!="" && dim.size()<3)
			{
				size_t a;
				stringstream ss;
				ss << word;
				ss >> a;
				dim.push_back(a);
			}  
            word = "";
        }
        else {
            word = word + x;
        }
    }
	if (word!="" && dim.size()<3)
		{
			size_t a;
			stringstream ss;
			ss << word;
			ss >> a;
			dim.push_back(a);
		}
}

void GetCoord(string str, vector<double>& coord)
{
    string word = "";
    for (auto x : str) 
    {
        if (x == ' ')
        {
			if (word!="")
			{
				string snum = "";
				size_t num=1;
				for (auto y : word)
				{
					if (y=='*')
					{
						stringstream ss;
						ss << snum;
						ss >> num;
						snum="";
					}
					else
						snum=snum+y;
				}
				
				double value;
				stringstream ss;
				ss << snum;
				ss >> value;
				for (size_t i=0; i<num; i++)
				{
					coord.push_back(value);
				}
			}  
            word = "";
        }
        else {
            word = word + x;
        }
    }
	if (word!="" && word.compare(0,1,"/")!=0 && word!="\r")
	{
		string snum = "";
		size_t num=1;
		for (auto y : word)
		{
			if (y=='*')
			{
				stringstream ss;
				ss << snum;
				ss >> num;
				snum="";
			}
			else
				snum=snum+y;
		}
				
		double value;
		stringstream ss;
		ss << snum;
		ss >> value;
		for (size_t i=0; i<num; i++)
		{
			coord.push_back(value);
		}
	}
}
	
void GetZcorn(string str, vector<double>& zcorn)
{
    string word = "";
    for (auto x : str) 
    {
        if (x == ' ')
        {
			if (word!="")
			{
				string snum = "";
				size_t num=1;
				for (auto y : word)
				{
					if (y=='*')
					{
						stringstream ss;
						ss << snum;
						ss >> num;
						snum="";
					}
					else
						snum=snum+y;
				}
				
				double value;
				stringstream ss;
				ss << snum;
				ss >> value;
				for (size_t i=0; i<num; i++)
				{
					zcorn.push_back(value);
				}
			}  
            word = "";
        }
        else {
            word = word + x;
        }
    }
	if (word!="" && word.compare(0,1,"/")!=0)
	{
		string snum = "";
		size_t num=1;
		for (auto y : word)
		{
			if (y=='*')
			{
				stringstream ss;
				ss << snum;
				ss >> num;
				snum="";
			}
			else
				snum=snum+y;
		}
				
		double value;
		stringstream ss;
		ss << snum;
		ss >> value;
		for (size_t i=0; i<num; i++)
		{
			zcorn.push_back(value);
		}
	}
}

void GetAct(string str, vector<size_t>& act)
{
    string word = "";
    for (auto x : str) 
    {
        if (x == ' ')
        {
			if (word!="")
			{
				string snum = "";
				size_t num=1;
				for (auto y : word)
				{
					if (y=='*')
					{
						stringstream ss;
						ss << snum;
						ss >> num;
						snum="";
					}
					else
						snum=snum+y;
				}
				
				size_t value;
				stringstream ss;
				ss << snum;
				ss >> value;
				for (size_t i=0; i<num; i++)
				{
					act.push_back(value);
				}
			}  
            word = "";
        }
        else {
            word = word + x;
        }
    }
	if (word!="" && word.compare(0,1,"/")!=0)
	{
		string snum = "";
		size_t num=1;
		for (auto y : word)
		{
			if (y=='*')
			{
				stringstream ss;
				ss << snum;
				ss >> num;
				snum="";
			}
			else
				snum=snum+y;
		}
				
		size_t value;
		stringstream ss;
		ss << snum;
		ss >> value;
		for (size_t i=0; i<num; i++)
		{
			act.push_back(value);
		}
	}
}

void GetProperty(string str, vector<double>& prop)
{
    string word = "";
    for (auto x : str) 
    {
        if (x == ' ')
        {
			if (word!="")
			{
				string snum = "";
				size_t num=1;
				for (auto y : word)
				{
					if (y=='*')
					{
						stringstream ss;
						ss << snum;
						ss >> num;
						snum="";
					}
					else
						snum=snum+y;
				}
				
				double value;
				stringstream ss;
				ss << snum;
				ss >> value;
				for (size_t i=0; i<num; i++)
				{
					prop.push_back(value);
				}
			}  
            word = "";
        }
        else {
            word = word + x;
        }
    }
	if (word!="" && word.compare(0,1,"/")!=0)
	{
		string snum = "";
		size_t num=1;
		for (auto y : word)
		{
			if (y=='*')
			{
				stringstream ss;
				ss << snum;
				ss >> num;
				snum="";
			}
			else
				snum=snum+y;
		}
				
		double value;
		stringstream ss;
		ss << snum;
		ss >> value;
		for (size_t i=0; i<num; i++)
		{
			prop.push_back(value);
		}
	}
}


bool AttachAct(Grid& grid, const char* filename, string name)
{
	vector<size_t> actnum;
	bool d = false;
	string buf;
	//string name2 = name +' ';
	ifstream ifs(filename);
	if(!ifs)
		return false;

	while (getline(ifs, buf))
	{
		if (d==true)
			{
				if (buf.compare(2,1,"/")==0)
					d=false;
				else if (buf.compare(buf.size()-2,1,"/")==0)
				{
					GetAct(buf, actnum);
					d=false;
				}
				else
					GetAct(buf, actnum);
			}
		if (buf.compare(0,name.length(),name)==0)
			d=true;
	}
	
	if (!GlobalAttachments::is_declared(name))
	{
		try {GlobalAttachments::declare_attachment<ANumber>(name);}
		catch (ug::UGError& err)
		{
			UG_LOGN(err.get_msg());
			return false;
		}
	}
	ANumber aVols = GlobalAttachments::attachment<ANumber>(name);
	if(!grid.has_volume_attachment(aVols))
		grid.attach_to_volumes(aVols);

	Grid::VolumeAttachmentAccessor<ANumber> aaVols(grid, aVols);
	
	size_t i = 0;
	for(VolumeIterator iter = grid.volumes_begin(); iter != grid.volumes_end(); ++iter, i++)
		{
			aaVols[*iter]=actnum[i];
//			ofs<<i<<' '<<aaVols[*iter]<<endl;
		}
	return true;
}


bool AttachProperty(Grid& grid, const char* filename, string name)
{
	vector<double> prop;
	bool d = false;
	bool e = false;
	string buf;
	string name2 = name +' ';
	ifstream ifs(filename);
	if(!ifs)
		return false;

	while (getline(ifs, buf))
	{
		if (d==true)
			{
				if (buf.compare(2,1,"/")==0)
					d=false;
				else if (buf.compare(buf.size()-2,1,"/")==0)
				{
					GetProperty(buf, prop);
					d=false;
				}
				else
					GetProperty(buf, prop);
			}
		if ((buf.compare(0,11,"-- Property")==0)&&(e==true))
			d=true;
		if (buf.compare(0,name2.length(),name2)==0)
			e=true;
	}
	
	if (!GlobalAttachments::is_declared(name))
	{
		try {GlobalAttachments::declare_attachment<ANumber>(name);}
		catch (ug::UGError& err)
		{
			UG_LOGN(err.get_msg());
			return false;
		}
	}
	ANumber aVols = GlobalAttachments::attachment<ANumber>(name);
	if(!grid.has_volume_attachment(aVols))
		grid.attach_to_volumes(aVols);

	Grid::VolumeAttachmentAccessor<ANumber> aaVols(grid, aVols);
	
//	ofstream ofs;
//	ofs.open ("Volumes.txt");
	
	size_t i = 0;
	for(VolumeIterator iter = grid.volumes_begin(); iter != grid.volumes_end(); ++iter, i++)
		{
			aaVols[*iter]=prop[i];
//			ofs<<i<<' '<<aaVols[*iter]<<endl;
		}
	return true;
}

bool LoadGridFromGRDECL(Grid& grid, const char* filename, AVector3& aPos)
{
	
	string buf;
	vector<size_t> dim;
	vector<double> coord;
	vector<double> zcorn;
	vector<size_t> act;
	vector<string> prop;
	vector<string> prop_name;
	
	bool a=false;
	bool b=false;
	bool c=false;
	bool d=false;
	
	ifstream ifs(filename);
	if(!ifs)
		return false;
	
	while (getline(ifs, buf))
		{
			//Get the dimension of the cells
			if (a==true)
			{
				GetDim(buf, dim);
				a=false;
			}
			if (buf.compare(0,8,"SPECGRID")==0)
				a=true;
			
			//Get the coordinate list of the top and bot surfaces
			if (b==true)
				{
					if (buf.compare(2,1,"/")==0)
						b=false;
					else if (buf.compare(buf.size()-2,1,"/")==0)
					{
						GetCoord(buf, coord);
						b=false;
					}
					else
						GetCoord(buf, coord);
				}
			if (buf.compare(0,6,"COORD ")==0 || buf.compare(0,7,"COORD\r")==0)
				b=true;
				
			//Get the depth list
			if (c==true)
				{
					if (buf.compare(2,1,"/")==0)
						c=false;
					else if (buf.compare(buf.size()-2,1,"/")==0)
					{
						GetCoord(buf, zcorn);
						c=false;
					}
					else
						GetCoord(buf, zcorn);
				}
			if (buf.compare(0,5,"ZCORN")==0)
				c=true;
				
			//Get the ACTNUM list
			if (d==true)
				{
					if (buf.compare(2,1,"/")==0)
						d=false;
					else if (buf.compare(buf.size()-2,1,"/")==0)
					{
						GetAct(buf, act);
						d=false;
					}
					else
						GetAct(buf, act);
				}
			if (buf.compare(0,6,"ACTNUM")==0)
				d=true;
	
			//Get property name list
			prop.push_back(buf);
			if (buf.compare(0,11,"-- Property")==0)
			{
				prop.pop_back();
				string name = "";
				for (auto x : prop.back()) 
				{
					if (x != ' ')
						name = name + x;
					else
						break;
				}
				prop_name.push_back(name);
			}
			
		}

		
	vector<xy> top;
	vector<xy> bot;
	// loop over pillar
	for (size_t j=0; j<dim[1]+1; j++)
	{
		for (size_t i=0; i<dim[0]+1; i++)
		{
			//build pillar coords
			top.push_back({coord[(j*(dim[0]+1)+i)*6], coord[(j*(dim[0]+1)+i)*6+1]});
			bot.push_back({coord[(j*(dim[0]+1)+i)*6+3], coord[(j*(dim[0]+1)+i)*6+4]});
			//check pillar is vertical
			if ( coord[(j*(dim[0]+1)+i)*6] != coord[(j*(dim[0]+1)+i)*6+3]|| coord[(j*(dim[0]+1)+i)*6+1] != coord[(j*(dim[0]+1)+i)*6+4] )
				UG_THROW ("LoadGridFromGRDECL: pillar is not vertical");
		}
	}
	
	vector<vector<double> > zcorn_ij((dim[0]+1)*(dim[1]+1), vector<double>(1,0));
	
	for (size_t k=0; k<dim[2]; k++)
	{
		for (size_t j=0; j<dim[1]; j++)
		{
			for (size_t i=0; i<dim[0]; i++)
			{
				size_t ij=j*(dim[0]+1)+i;
				size_t ijk4=k*dim[1]*dim[0]*8+j*dim[0]*4+i*2;
				
				zcorn_ij[ij].push_back(zcorn[ijk4]);
				zcorn_ij[ij+1].push_back(zcorn[ijk4+1]);
				zcorn_ij[ij+dim[0]+2].push_back(zcorn[ijk4+dim[0]*2+1]);
				zcorn_ij[ij+dim[0]+1].push_back(zcorn[ijk4+dim[0]*2]);
				
				zcorn_ij[ij].push_back(zcorn[ijk4+dim[1]*dim[0]*4]);
				zcorn_ij[ij+1].push_back(zcorn[ijk4+dim[1]*dim[0]*4+1]);
				zcorn_ij[ij+dim[0]+2].push_back(zcorn[ijk4+dim[1]*dim[0]*4+dim[0]*2+1]);
				zcorn_ij[ij+dim[0]+1].push_back(zcorn[ijk4+dim[1]*dim[0]*4+dim[0]*2]);			
			}
		}
	}
	
	vector<vector<size_t> > Index_zcorn((dim[0]+1)*(dim[1]+1), vector<size_t>(1,0));
	vector<vector<size_t> > Index_zcorn_Global((dim[0]+1)*(dim[1]+1), vector<size_t>(1,0));
	vector<vector<double> > New_zcorn((dim[0]+1)*(dim[1]+1), vector<double>(1,0));
	size_t numVrts=0;
	
	vector<xyz> coord_list;

	for (size_t j=0; j<dim[1]+1; j++)
	{
		for (size_t i=0; i<dim[0]+1; i++)
		{	
			size_t ij=j*(dim[0]+1)+i;
			
			//remove redundant coords
			New_zcorn[ij].push_back(zcorn_ij[ij][1]);
			for (size_t t=2; t<zcorn_ij[ij].size(); t++)
			{
				bool New_coord=true;
				for (size_t s =1; s<New_zcorn[ij].size(); s++)
				{
					if (zcorn_ij[ij][t]==New_zcorn[ij][s])
					{
						New_coord = false;
						Index_zcorn[ij].push_back(s-1);
						break;
					}
				}
				
				if (New_coord==true)
				{
					Index_zcorn[ij].push_back(New_zcorn[ij].size()-1);
					New_zcorn[ij].push_back(zcorn_ij[ij][t]);
				}
			}
			
			for (size_t k=1; k<New_zcorn[ij].size(); k++)
			{
				coord_list.push_back({top[ij].x, top[ij].y, -New_zcorn[ij][k]});
				Index_zcorn_Global[ij].push_back(numVrts);
				numVrts++;
			}
			
			reverse(New_zcorn[ij].begin(), New_zcorn[ij].end());
			New_zcorn[ij].pop_back();
			reverse(zcorn_ij[ij].begin(), zcorn_ij[ij].end());
			zcorn_ij[ij].pop_back();
			reverse(Index_zcorn[ij].begin(), Index_zcorn[ij].end());
			
			reverse(Index_zcorn_Global[ij].begin(), Index_zcorn_Global[ij].end());
			Index_zcorn_Global[ij].pop_back();
			reverse(Index_zcorn_Global[ij].begin(), Index_zcorn_Global[ij].end());
		}
		
	}
	
	vector<ab> ele_list;
	
	for (size_t k=0; k<dim[2]; k++)
	{
		for (size_t j=0; j<dim[1]; j++)
		{
			for (size_t i=0; i<dim[0]; i++)
			{
				size_t ij=j*(dim[0]+1)+i;
				
				size_t a1=Index_zcorn_Global[ij][Index_zcorn[ij].back()];
				Index_zcorn[ij].pop_back();
				size_t a2=Index_zcorn_Global[ij+1][Index_zcorn[ij+1].back()];
				Index_zcorn[ij+1].pop_back();
				size_t a3=Index_zcorn_Global[ij+dim[0]+2][Index_zcorn[ij+dim[0]+2].back()];
				Index_zcorn[ij+dim[0]+2].pop_back();
				size_t a4=Index_zcorn_Global[ij+dim[0]+1][Index_zcorn[ij+dim[0]+1].back()];
				Index_zcorn[ij+dim[0]+1].pop_back();
				
				size_t b1=Index_zcorn_Global[ij][Index_zcorn[ij].back()];
				Index_zcorn[ij].pop_back();
				size_t b2=Index_zcorn_Global[ij+1][Index_zcorn[ij+1].back()];
				Index_zcorn[ij+1].pop_back();
				size_t b3=Index_zcorn_Global[ij+dim[0]+2][Index_zcorn[ij+dim[0]+2].back()];
				Index_zcorn[ij+dim[0]+2].pop_back();
				size_t b4=Index_zcorn_Global[ij+dim[0]+1][Index_zcorn[ij+dim[0]+1].back()];
				Index_zcorn[ij+dim[0]+1].pop_back();
				
				ele_list.push_back({a1, a2, a3, a4, b1, b2, b3, b4});
			}
		}
	}
		
	size_t numElems = ele_list.size();

//	create points
//	store pointers to the vertices on the fly in a vector.
	vector<Vertex*>	vVrts;
	vector<size_t> vVrtIds;
	vVrts.resize(numVrts); vVrtIds.resize(numVrts);

	for(size_t i = 0; i < numVrts; ++i)
		vVrts[i] = *grid.create<RegularVertex>();

	if(!grid.has_vertex_attachment(aPos))
		grid.attach_to_vertices(aPos);
	Grid::VertexAttachmentAccessor<APosition> aaPos(grid, aPos);

//	read the points
	{
		size_t i = 0;
		for(VertexIterator iter = grid.vertices_begin(); iter != grid.vertices_end(); ++iter, ++i)
		{
			vVrtIds[i]=i;
			aaPos[*iter].x()=coord_list[i].x;
			aaPos[*iter].y()=coord_list[i].y;
			aaPos[*iter].z()=coord_list[i].z;
		}

	}

//	read the hexahedrons
	{
		for(size_t i = 0; i < numElems; ++i)
		{
			grid.create<Hexahedron>
				(HexahedronDescriptor
					(vVrts[ele_list[i].a1], vVrts[ele_list[i].a2], vVrts[ele_list[i].a3], vVrts[ele_list[i].a4], vVrts[ele_list[i].b1], vVrts[ele_list[i].b2], vVrts[ele_list[i].b3], vVrts[ele_list[i].b4]));
		}
	}

	AttachAct(grid, filename, "ACTNUM");	

//	read the volumes(properties)
	
	for (size_t i = 0; i<prop_name.size(); i++)
		{
			AttachProperty(grid, filename, prop_name[i]);
		}


	return true;
}

}//	end of namespace