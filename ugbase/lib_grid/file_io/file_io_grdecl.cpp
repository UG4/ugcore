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
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "../lg_base.h"
#include "common/util/string_util.h"
#include "lib_grid/global_attachments.h"

#include "file_io_grdecl.h"

using namespace std;

namespace ug
{

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
					GetAct(buf, actnum);
					d=false;
				}
				else
					GetAct(buf, actnum);
			}
		if (buf.compare(0,name2.length(),name2)==0)
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
	
//	ofstream ofs;
//	ofs.open ("Actnum.txt");
	
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
			if (buf.compare(0,6,"COORD ")==0)
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
			if (buf.compare(0,6,"ZCORN ")==0)
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
			if (buf.compare(0,7,"ACTNUM ")==0)
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
		
/*	ofstream ofs;
	ofs.open ("prop_name.txt");	
	for (auto i = 0; i<zcorn.size(); i++)
		{
			ofs<<i<<' '<<zcorn[i]<<endl;
			//ofs<<i<<' '<<buf<<endl;
		}*/
		
		
		
	vector<xyz> top;
	vector<xyz> bot;
	for (size_t j=0; j<dim[1]+1; j++)
	{
		for (size_t i=0; i<dim[0]+1; i++)
		{
			top.push_back({coord[(j*(dim[0]+1)+i)*6], coord[(j*(dim[0]+1)+i)*6+1], coord[(j*(dim[0]+1)+i)*6+2]});
			bot.push_back({coord[(j*(dim[0]+1)+i)*6+3], coord[(j*(dim[0]+1)+i)*6+4], coord[(j*(dim[0]+1)+i)*6+5]});
		}
	}
	
	vector<xyz> coord_list;
	vector<ab> ele_list;
	
	for (size_t k=0; k<dim[2]; k++)
	{
		for (size_t j=0; j<dim[1]; j++)
		{
			for (size_t i=0; i<dim[0]; i++)
			{
				
				size_t ij=j*(dim[0]+1)+i;
				size_t ijk4=k*dim[1]*dim[0]*8+j*dim[0]*4+i*2;
				size_t ijk8=(k*dim[1]*dim[0]+j*dim[0]+i)*8;
				
				ele_list.push_back({ijk8+1, ijk8+2, ijk8+3, ijk8+4, ijk8+5, ijk8+6, ijk8+7, ijk8+8});
				
				coord_list.push_back({top[ij].x, top[ij].y, zcorn[ijk4]});
				coord_list.push_back({top[ij+1].x, top[ij+1].y, zcorn[ijk4+1]});
				coord_list.push_back({top[ij+dim[0]+2].x, top[ij+dim[0]+2].y, zcorn[ijk4+dim[0]*2+1]});
				coord_list.push_back({top[ij+dim[0]+1].x, top[ij+dim[0]+1].y, zcorn[ijk4+dim[0]*2]});
				
				coord_list.push_back({top[ij].x, top[ij].y, zcorn[ijk4+dim[1]*dim[0]*4]});
				coord_list.push_back({top[ij+1].x, top[ij+1].y, zcorn[ijk4+dim[1]*dim[0]*4+1]});
				coord_list.push_back({top[ij+dim[0]+2].x, top[ij+dim[0]+2].y, zcorn[ijk4+dim[1]*dim[0]*4+dim[0]*2+1]});
				coord_list.push_back({top[ij+dim[0]+1].x, top[ij+dim[0]+1].y, zcorn[ijk4+dim[1]*dim[0]*4+dim[0]*2]});
			}
		}
	}
	
	//remove redundant coords
	vector<xyz> New_coord_list;
	vector<ab> New_ele_list;
	vector<size_t> Index_coord;
	New_coord_list.push_back(coord_list[0]);
	Index_coord.push_back(0);

	for (size_t i=1; i<coord_list.size(); i++)
	{
		bool New_coord=true;
		for (size_t j=0; j<New_coord_list.size(); j++)
		{
			if ((coord_list[i].x==New_coord_list[j].x) && (coord_list[i].y==New_coord_list[j].y) && (coord_list[i].z==New_coord_list[j].z))
				{
					New_coord=false;
					Index_coord.push_back(j);
					break;
				}
		}
		if (New_coord==true)
			{
				Index_coord.push_back(New_coord_list.size());
				New_coord_list.push_back(coord_list[i]);
			}
	}
	
	ofstream ofs;
	ofs.open ("coord_list.txt");	
	for (size_t i = 0; i<New_coord_list.size(); i++)
		{
			ofs<<i<<' '<<New_coord_list[i].x<<' '<<New_coord_list[i].y<<' '<<New_coord_list[i].z<<endl;
		}
	
	
	
	size_t numVrts, numElems;

	numVrts=New_coord_list.size();
	numElems=ele_list.size();

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
			aaPos[*iter].x()=New_coord_list[i].x;
			aaPos[*iter].y()=New_coord_list[i].y;
			aaPos[*iter].z()=New_coord_list[i].z;
		}

	}



//	read the hexahedrons
	{
		for(size_t i = 0; i < numElems; ++i)
		{
			size_t Index, vrt_id_1, vrt_id_2, vrt_id_3, vrt_id_4, vrt_id_5, vrt_id_6, vrt_id_7, vrt_id_8;
			vector<size_t>::iterator i1, i2, i3, i4, i5, i6, i7, i8;
			
			Index=i;
			vrt_id_1=Index_coord[ele_list[i].a1-1];
			vrt_id_2=Index_coord[ele_list[i].a2-1];
			vrt_id_3=Index_coord[ele_list[i].a3-1];
			vrt_id_4=Index_coord[ele_list[i].a4-1];
			vrt_id_5=Index_coord[ele_list[i].b1-1];
			vrt_id_6=Index_coord[ele_list[i].b2-1];
			vrt_id_7=Index_coord[ele_list[i].b3-1];
			vrt_id_8=Index_coord[ele_list[i].b4-1];
			
			i1 = find (vVrtIds.begin(), vVrtIds.end(), vrt_id_1);
			i2 = find (vVrtIds.begin(), vVrtIds.end(), vrt_id_2);
			i3 = find (vVrtIds.begin(), vVrtIds.end(), vrt_id_3);
			i4 = find (vVrtIds.begin(), vVrtIds.end(), vrt_id_4);
			i5 = find (vVrtIds.begin(), vVrtIds.end(), vrt_id_5);
			i6 = find (vVrtIds.begin(), vVrtIds.end(), vrt_id_6);
			i7 = find (vVrtIds.begin(), vVrtIds.end(), vrt_id_7);
			i8 = find (vVrtIds.begin(), vVrtIds.end(), vrt_id_8);
			
	if (i1 == vVrtIds.end() || i2 == vVrtIds.end() || i3 == vVrtIds.end() || i4 == vVrtIds.end() || i5 == vVrtIds.end() || i6 == vVrtIds.end() || i7 == vVrtIds.end() || i8 == vVrtIds.end() )
				UG_THROW ("LoadGridFromGRDECL: Wrong index of a vertex of element " << Index << "in '" << filename << "'");
			grid.create<Hexahedron>
				(HexahedronDescriptor
					(vVrts[i1 - vVrtIds.begin()], vVrts[i2 - vVrtIds.begin()], vVrts[i3 - vVrtIds.begin()], vVrts[i4 - vVrtIds.begin()], vVrts[i5 - vVrtIds.begin()], vVrts[i6 - vVrtIds.begin()], vVrts[i7 - vVrtIds.begin()], vVrts[i8 - vVrtIds.begin()]));
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