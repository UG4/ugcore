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

#include <fstream>
#include "lib_grid/lg_base.h"
#include "lib_grid/algorithms/serialization.h"
#include "file_io_lgb.h"

#include <sstream>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include "common/boost_serialization_routines.h"
#include "common/util/archivar.h"
#include "common/util/factory.h"
#include "lib_grid/refinement/projectors/projection_handler.h"
#include "lib_grid/refinement/projectors/projectors.h"

using namespace std;

namespace ug
{

enum LGBConstants
{
	LGBC_NONE = 0,
	LGBC_POS2D = 1,
	LGBC_POS3D = 1 << 1,
	LGBC_SUBSET_HANDLER = 1 << 2,
	LGBC_SELECTOR = 1 << 3,
	LGBC_PROJECTION_HANDLER = 1 << 4
};

void SerializeProjector(BinaryBuffer& out, RefinementProjector& proj)
{
	static Factory<RefinementProjector, ProjectorTypes>	projFac;
	static Archivar<boost::archive::text_oarchive, RefinementProjector, ProjectorTypes>	archivar;
	
	const string& projName = projFac.class_name(proj);

	Serialize(out, projName);

	stringstream ss;
	boost::archive::text_oarchive ar(ss, boost::archive::no_header);
	archivar.archive(ar, proj);

	Serialize(out, ss.str());
}


void SerializeProjectionHandler(BinaryBuffer& out, ProjectionHandler& ph)
{
	const int magicNumber = 978523;
	out.write((char*)&magicNumber, sizeof(int));

	if(ph.default_projector().valid()){
		byte_t b = 1;
		out.write((char*)&b, sizeof(b));
		SerializeProjector(out, *ph.default_projector());
	}
	else{
		byte_t b = 0;
		out.write((char*)&b, sizeof(b));
	}

	int numProjectors = (int)ph.num_projectors();
	out.write((char*)&numProjectors, sizeof(int));

	for(int i = -1; i < numProjectors; ++i){
		if(!ph.projector(i).valid()){
			const int invInd = -2;
			out.write((char*)& invInd, sizeof(int));
			continue;
		}

		RefinementProjector& proj= *ph.projector(i);

		out.write((char*)& i, sizeof(int));
		SerializeProjector(out, proj);
	}
	out.write((char*)&magicNumber, sizeof(int));
}


SPRefinementProjector DeserializeProjector(BinaryBuffer& in)
{
	static Archivar<boost::archive::text_iarchive,
				RefinementProjector,
				ProjectorTypes>
			archivar;

	static Factory<RefinementProjector, ProjectorTypes>	projFac;

	std::string name;
	Deserialize(in, name);

	if(name.empty())
		return SPRefinementProjector();
	
	SPRefinementProjector proj = projFac.create(name);

	std::string data;
	Deserialize(in, data);
	std::stringstream ss(data, std::ios_base::in);
	boost::archive::text_iarchive ar(ss, boost::archive::no_header);
	archivar.archive(ar, *proj);
	return proj;
}


void DeserializeProjectionHandler(BinaryBuffer& in, ProjectionHandler& ph)
{
	const int magicNumber = 978523;
	int tmpMagicNumber;
	in.read((char*)&tmpMagicNumber, sizeof(int));
	UG_COND_THROW(tmpMagicNumber != magicNumber,
	              "Magic number mismatch in DeserializeProjectionHandler (1)!");

	ph.clear();
	
	byte_t b;
	in.read((char*)&b, sizeof(b));
	if(b){
		ph.set_default_projector(DeserializeProjector(in));
	}
	else
		ph.set_default_projector(SPRefinementProjector());

	int numProjectors;
	in.read((char*)&numProjectors, sizeof(int));

	for(int i = -1; i < numProjectors; ++i){
		int index;
		in.read((char*)& index, sizeof(int));
		if(index == -2){
			continue;
		}

		ph.set_projector(index, DeserializeProjector(in));
	}

	in.read((char*)&tmpMagicNumber, sizeof(int));
	UG_COND_THROW(tmpMagicNumber != magicNumber,
	              "Magic number mismatch in DeserializeProjectionHandler (2)!");
}




bool SaveGridToLGB(Grid& grid, const char* filename,
				   ISubsetHandler** ppSH, int numSHs,
				   ProjectionHandler* pPH,
				   APosition aPos)
{
	return SaveGridToLGB (grid, filename,ppSH, numSHs, nullptr, 0, pPH, aPos);
}


bool SaveGridToLGB(Grid& grid, const char* filename,
				   ISubsetHandler** ppSH, int numSHs,
				   ISelector** ppSel, int numSels,
				   ProjectionHandler* pPH,
				   APosition aPos)
{
//	make sure that aPos is attached to the grid
	if(!grid.has_vertex_attachment(aPos))
		return false;

	BinaryBuffer tbuf;

//	write the header
	byte_t endianess = 1;
	byte_t intSize = (byte_t)sizeof(int);
	byte_t numberSize = (byte_t)sizeof(number);
	int versionNumber = 4;
	
	tbuf.write((char*)&endianess, sizeof(byte_t));
	tbuf.write((char*)&intSize, sizeof(byte_t));
	tbuf.write((char*)&numberSize, sizeof(byte_t));
	tbuf.write((char*)&versionNumber, sizeof(int));

//	the options
	uint opts = LGBC_POS3D;
	if(numSHs > 0)
		opts |= LGBC_SUBSET_HANDLER;
	if(numSels > 0)
		opts |= LGBC_SELECTOR;
	if(pPH)
		opts |= LGBC_PROJECTION_HANDLER;

//	write the options
	tbuf.write((char*)&opts, sizeof(uint));

//	serialize the grid
	SerializeGridElements(grid, tbuf);

//	serialize the position-attachment
	SerializeAttachment<Vertex>(grid, aPos, tbuf);

//	Serialize the subset-handler
	if(numSHs > 0){
	//	write the number of subset handlers which shall be serialized
		tbuf.write((char*)&numSHs, sizeof(int));
		for(int i = 0; i< numSHs; ++i)
			SerializeSubsetHandler(grid, *ppSH[i], tbuf);
	}

//	Serialize the selectors
	if(numSels > 0){
		tbuf.write((char*)&numSels, sizeof(int));
		for(int i = 0; i< numSels; ++i)
			SerializeSelector(grid, *ppSel[i], tbuf);
	}

//	Serializie the projection handler
	if(pPH){
		ProjectionHandler& ph = *pPH;

		//	find corresponding subset handler
		{
			UG_COND_THROW(numSHs < 1, "If a ProjectionHandler is specified, the corresponding SubsetHandler also has to be specified!");

			int i = 0;
			for (; i < numSHs; ++i)
				if (ppSH[i] == ph.subset_handler())
					break;

			UG_COND_THROW(i == numSHs, "ERROR in 'SaveGridToLGB': "
				"No matching SubsetHandler could be found.\n"
				"Please make sure to add the associated SubsetHandler before adding a ProjectionHandler");

		//	write index of corresponding subset handler
			tbuf.write((char*)&i, sizeof(int));
		}

		SerializeProjectionHandler(tbuf, ph);
	}

//	write a magic-number that allows us to check during read
//	whether everything went ok.
	int magicNumber = 3478384;
	tbuf.write((char*)&magicNumber, sizeof(int));

//	open the outstream
	ofstream out(filename, ios::binary);
	if(!out) return false;

	out.write(tbuf.buffer(), tbuf.write_pos());

//	done
	out.close();
	return true;
}



bool LoadGridFromLGB(Grid& grid, const char* filename,
				   ISubsetHandler** ppSH, int numSHs,
				   ProjectionHandler* pPH,
				   APosition aPos)
{
	return LoadGridFromLGB (grid, filename, ppSH, numSHs, nullptr, 0, pPH, aPos);
}


bool LoadGridFromLGB(Grid& grid, const char* filename,
				   ISubsetHandler** ppSH, int numSHs,
				   ISelector** ppSel, int numSels,
				   ProjectionHandler* pPH,
				   APosition aPos)
{
	grid.clear_geometry();

//	if aPos is not yet attached to the grid, we'll do it now.
	if(!grid.has_vertex_attachment(aPos))
		grid.attach_to_vertices(aPos);

//	open the outstream
	ifstream in(filename, ios::binary);
	if(!in){
		UG_LOG("ERROR in LoadGridFromLGB: couldn't open file: " << filename << endl);
		return false;
	}

//	read the whole file into a binary buffer
	BinaryBuffer tbuf;
	in.seekg(0, ios::end);
	streampos fileLen = in.tellg();
	in.seekg(0, ios::beg);

	tbuf.reserve(fileLen);
	in.read(tbuf.buffer(), fileLen);
	tbuf.set_write_pos(fileLen);

	in.close();

//	read the header
	byte_t endianess;
	byte_t intSize;
	byte_t numberSize;
	int versionNumber;

	tbuf.read((char*)&endianess, sizeof(byte_t));
	tbuf.read((char*)&intSize, sizeof(byte_t));
	tbuf.read((char*)&numberSize, sizeof(byte_t));
	tbuf.read((char*)&versionNumber, sizeof(int));

//	check whether the values are ok
	if(endianess != 1)
	{
		LOG("ERROR in LoadGridFromLGB: wrong endianess\n");
		return false;
	}
	if(intSize != sizeof(int))
	{
		LOG("ERROR in LoadGridFromLGB: bad integer-size\n");
		return false;
	}
	if(numberSize != sizeof(number))
	{
		LOG("ERROR in LoadGridFromLGB: bad number-size\n");
		return false;
	}
	if((versionNumber < 2) || (versionNumber > 4))
	{
		LOG("ERROR in LoadGridFromLGB: bad file-version: " << versionNumber << ". Expected 2 or 3.\n");
		return false;
	}

//	from version 3 on grids write a small header
	bool readGridHeader = (versionNumber >= 3);

//	read the options
	uint opts;
	tbuf.read((char*)&opts, sizeof(uint));

//	to avoid problems with autgenerated elements we'll deactivate
//	all options and reactivate them later on
	uint gridOptions = grid.get_options();
	grid.set_options(GRIDOPT_NONE);

//	deserialize the grid
	DeserializeGridElements(grid, tbuf, readGridHeader);

//	deserialize the position-attachment
	DeserializeAttachment<Vertex>(grid, aPos, tbuf);

//	Serialize the subset-handler
	if((opts & LGBC_SUBSET_HANDLER) == LGBC_SUBSET_HANDLER)
	{
	//	read number of subset handlers
		int numSrcSHs;
		tbuf.read((char*)&numSrcSHs, sizeof(int));
		
	//	starting from version 4, subset-infos contain a property-map
		bool readPropertyMap = (versionNumber >= 4);

		int i;
		for(i = 0; i < min(numSrcSHs, numSHs); ++i){
			DeserializeSubsetHandler(grid, *ppSH[i], tbuf, readPropertyMap);
		}
		
	//	read the rest
		for(; i < numSrcSHs; ++i)
		{
			SubsetHandler sh(grid);
			DeserializeSubsetHandler(grid, sh, tbuf, readPropertyMap);
		}
	}

	if((opts & LGBC_SELECTOR) == LGBC_SELECTOR)
	{
	//	read number of subset handlers
		int numSrcSels;
		tbuf.read((char*)&numSrcSels, sizeof(int));
		
		int i;
		for(i = 0; i < min(numSrcSels, numSels); ++i){
			DeserializeSelector(grid, *ppSel[i], tbuf);
		}
		
	//	read the rest
		for(; i < numSrcSels; ++i)
		{
			Selector sel(grid);
			DeserializeSelector(grid, sel, tbuf);
		}
	}


	if((opts & LGBC_PROJECTION_HANDLER) == LGBC_PROJECTION_HANDLER)
	{
	//	read number of subset handlers
		int shIndex;
		tbuf.read((char*)&shIndex, sizeof(int));
		
		if(pPH)
			DeserializeProjectionHandler(tbuf, *pPH);
		else{
			ProjectionHandler ph;
			DeserializeProjectionHandler(tbuf, ph);
		}
	}

//	reactivate the grid-options
	grid.set_options(gridOptions);

//	check the magic-number
	int magicNumber;
	tbuf.read((char*)&magicNumber, sizeof(int));

	if(magicNumber != 3478384){
		LOG("ERROR in LoadGridFromLGB: Bad magic number at end of file: " << magicNumber << endl);
		return false;
	}
//	done
	return true;
}

}//	end of namespace
