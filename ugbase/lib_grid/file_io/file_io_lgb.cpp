//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m11 d11

#include <fstream>
#include "lib_grid/lib_grid.h"
#include "file_io_lgb.h"

using namespace std;

namespace ug
{

enum LGBConstants
{
	LGBC_NONE = 0,
	LGBC_POS2D = 1,
	LGBC_POS3D = 1 << 1,
	LGBC_SUBSET_HANDLER = 1 << 2
};

bool SaveGridToLGB(Grid& grid, const char* filename,
				   SubsetHandler* pSH, APosition aPos)
{
//	make sure that aPos is attached to the grid
	if(!grid.has_vertex_attachment(aPos))
		return false;

//	open the outstream
	ofstream out(filename, ios::binary);
	if(!out) return false;

//	write the header
	byte endianess = 1;
	byte intSize = (byte)sizeof(int);
	byte numberSize = (byte)sizeof(number);
	int versionNumber = 1;
	
	out.write((char*)&endianess, sizeof(byte));
	out.write((char*)&intSize, sizeof(byte));
	out.write((char*)&numberSize, sizeof(byte));
	out.write((char*)&versionNumber, sizeof(int));

//	the options
	uint opts = LGBC_POS3D;
	if(pSH)
		opts |= LGBC_SUBSET_HANDLER;

//	write the options
	out.write((char*)&opts, sizeof(uint));

//	serialize the grid
	SerializeGridElements(grid, out);

//	serialize the position-attachment
	SerializeAttachment<VertexBase>(grid, aPos, out);

//	Serialize the subset-handler
	if(pSH)
		SerializeSubsetHandler(grid, *pSH, out);

//	write a magic-number that allows us to check during read
//	whether everything went ok.
	int magicNumber = 3478384;
	out.write((char*)&magicNumber, sizeof(int));

//	done
	out.close();
	return true;
}



bool LoadGridFromLGB(Grid& grid, const char* filename,
				   SubsetHandler* pSH, APosition aPos)
{
//	if aPos is not yet attached to the grid, we'll do it now.
	if(!grid.has_vertex_attachment(aPos))
		grid.attach_to_vertices(aPos);

//	open the outstream
	ifstream in(filename, ios::binary);
	if(!in) return false;

//	read the header
	byte endianess;
	byte intSize;
	byte numberSize;
	int versionNumber;

	in.read((char*)&endianess, sizeof(byte));
	in.read((char*)&intSize, sizeof(byte));
	in.read((char*)&numberSize, sizeof(byte));
	in.read((char*)&versionNumber, sizeof(int));

//	check whether the values are ok
	if(endianess != 1)
	{
		LOG("ERROR during LoadGridFromLGB: wrong endianess\n");
		return false;
	}
	if(intSize != sizeof(int))
	{
		LOG("ERROR during LoadGridFromLGB: bad integer-size\n");
		return false;
	}
	if(numberSize != sizeof(number))
	{
		LOG("ERROR during LoadGridFromLGB: bad number-size\n");
		return false;
	}
	if(versionNumber != 1)
	{
		LOG("ERROR during LoadGridFromLGB: bad file-version: " << versionNumber << ". Expected 1.\n");
		return false;
	}

//	read the options
	uint opts;
	in.read((char*)&opts, sizeof(uint));

//	to avoid problems with autgenerated elements we'll deactivate
//	all options and reactivate them later on
	uint gridOptions = grid.get_options();
	grid.set_options(GRIDOPT_NONE);

//	serialize the grid
	DeserializeGridElements(grid, in);

//	serialize the position-attachment
	DeserializeAttachment<VertexBase>(grid, aPos, in);

//	Serialize the subset-handler
	if((opts & LGBC_SUBSET_HANDLER) == LGBC_SUBSET_HANDLER)
	{
		if(pSH)
			DeserializeSubsetHandler(grid, *pSH, in);
		else
		{
			SubsetHandler sh(grid);
			DeserializeSubsetHandler(grid, sh, in);
		}
	}

//	reactivate the grid-options
	grid.set_options(gridOptions);

//	check the magic-number
	int magicNumber;
	in.read((char*)&magicNumber, sizeof(int));

	if(magicNumber != 3478384)
		return false;
//	done
	in.close();
	return true;
}

}//	end of namespace
