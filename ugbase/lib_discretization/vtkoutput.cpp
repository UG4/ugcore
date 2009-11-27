/*
 * vtkoutput.cpp
 *
 *  Created on: 06.07.2009
 *      Author: andreasvogel
 */


#define __SWAPBYTES__            /* if using LittleEndian */

#include "vtkoutput.h"
#include <iostream>
#include <cstring>

namespace ug{


	enum{OSIZE=32,            /* must be a multiple of 16 */
		 BSIZE=(3*OSIZE/4)};  /* will be just right then  */

static struct {
    char buffer[BSIZE];
	char output[OSIZE];
    int front;
    int size;
} BStream;


/*BEGIN: Helper Functions */
static void EncodeTriplet(char *_in, char *out, int n)
{
	static char digits[] =
		"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
	unsigned char *in = (unsigned char *)_in;

	out[0] = digits[in[0] >> 2];
	out[1] = digits[((in[0] & 0x03) << 4) | (in[1] >> 4)];
	out[2] = n > 1? digits[((in[1] & 0x0F) << 2) | (in[2] >> 6)] : '=';
	out[3] = n > 2? digits[in[2] & 0x3F] : '=';
}

static void BStreamWrite(FILE* File, void *item)
{
	int i;
	char *p;

	memcpy(BStream.buffer + BStream.front, item, BStream.size);
	BStream.front += BStream.size;
	if (BStream.front == BSIZE) {
		p = BStream.output;
		for (i = 0; i < BSIZE; i += 3) {
			EncodeTriplet(BStream.buffer + i, p, 3);
			p += 4;
		}
		fwrite(BStream.output, 1, OSIZE, File);
		BStream.front = 0;
	}
}

static void BStreamFlush(FILE* File)
{
	int i, r, to;
	char *p;

	if (BStream.front != 0) {
		p = BStream.output;
		r = BStream.front % 3;
		to = BStream.front - r;
		for (i = 0; i < to; i += 3) {
			EncodeTriplet(BStream.buffer + i, p, 3);
			p += 4;
		}
		if (r) {
			memset(BStream.buffer + BStream.front, 0, 3-r);
			EncodeTriplet(BStream.buffer + to, p, r);
			p += 4;
		}
		fwrite(BStream.output, 1, p - BStream.output, File);
		BStream.front = 0;
	}
}
 /* END: Helper Functions */

bool VTKOutput::InitNumbers(SubsetHandler& sh, uint subsetIndex)
{
	// TODO: Test Version for 2D only!
	int n = 0;

	Numbers.noVertices = sh.num_elements<Vertex>(subsetIndex);
	Numbers.noElements = sh.num_elements<Triangle>(subsetIndex) + sh.num_elements<Quadrilateral>(subsetIndex); //TODO: This is 2D only!
	Numbers.noConnections = 3*sh.num_elements<Triangle>(subsetIndex) +
							4*sh.num_elements<Quadrilateral>(subsetIndex); //TODO: May be wrong for hybrid grids, 2D only!

	return true;
}

bool VTKOutput::print(NumericalSolution& u, const char* filename, double Time)
{
	SubsetHandler* sh = u.get_pattern()->get_assigned_subset();
	Grid* grid = sh->get_assigned_grid();

	// attach help indices
	grid->attach_to_vertices(m_aDOFIndex);
	m_aaDOFIndexVRT.access(*grid, m_aDOFIndex);

	// open stream
	FILE* File = fopen(filename, "w");
	if(File == NULL)
	{
		std::cout << "ERROR (in VTKOutput::print(...)): Can not open Output File" << std::endl;
		return false;
	}

	BStream.front = 0;


	// Write to File
	if(write_prolog(File, Time) != true)
	{
		std::cout << "ERROR (in VTKOutput::print(...)): Can not write Prolog" << std::endl;
		fclose(File);return false;
	}

	// loop subsets
	for(int subsetIndex = 0; subsetIndex < sh->num_subsets(); ++subsetIndex)
	{
		if(write_subset(*sh, subsetIndex, u, File)!= true)
		{
			std::cout << "ERROR (in VTKOutput::print(...)): Can not write Subset" << std::endl;
			fclose(File);return false;
		}
	}

	if(write_epilog(File) != true)
	{
		std::cout << "ERROR (in VTKOutput::print(...)): Can not write Epilog" << std::endl;
		fclose(File);return false;
	}

	// close stream
	fclose(File);

	// detach help indices
	grid->detach_from_vertices(m_aDOFIndex);
	m_aaDOFIndexVRT.invalidate();

	return true;
}

bool VTKOutput::print_subset(NumericalSolution& u, int subsetIndex, const char* filename, double Time)
{
	SubsetHandler* sh = u.get_pattern()->get_assigned_subset();
	Grid* grid = sh->get_assigned_grid();

	// attach help indices
	grid->attach_to_vertices(m_aDOFIndex);
	m_aaDOFIndexVRT.access(*grid, m_aDOFIndex);

	// open stream
	FILE* File = fopen(filename, "w");
	if(File == NULL)
	{
		std::cout << "ERROR (in VTKOutput::print(...)): Can not open Output File" << std::endl;
		return false;
	}

	BStream.front = 0;


	// Write to File
	if(write_prolog(File, Time) != true)
	{
		std::cout << "ERROR (in VTKOutput::print(...)): Can not write Prolog" << std::endl;
		fclose(File);return false;
	}

	// Write Subset
	if(write_subset(*sh, subsetIndex, u, File) != true)
	{
		std::cout << "ERROR (in VTKOutput::print(...)): Can not write Subset" << std::endl;
		fclose(File);return false;
	}

	if(write_epilog(File) != true)
	{
		std::cout << "ERROR (in VTKOutput::print(...)): Can not write Epilog" << std::endl;
		fclose(File);return false;
	}

	// close stream
	fclose(File);

	// detach help indices
	grid->detach_from_vertices(m_aDOFIndex);
	m_aaDOFIndexVRT.invalidate();

	return true;
}


bool VTKOutput::write_subset(SubsetHandler& sh, int subsetIndex, NumericalSolution& u, FILE* File)
{
	// Read sizes
	InitNumbers(sh, subsetIndex);

	if(write_piece_prolog(File) != true)
	{
		std::cout << "ERROR (in VTKOutput::print(...)): Can not write Prolog" << std::endl;
		fclose(File);return false;
	}

	if(write_points(File, sh, subsetIndex) != true)
	{
		std::cout << "ERROR (in VTKOutput::print(...)): Can not write Points" << std::endl;
		fclose(File);return false;
	}
	if(write_elements(File, sh, subsetIndex) != true)
	{
		std::cout << "ERROR (in VTKOutput::print(...)): Can not write Elements" << std::endl;
		fclose(File);return false;
	}

	fprintf(File, "      <PointData>\n");
	for(int comp = 0; comp < u.get_pattern()->num_comp(); ++comp)
	{
		if(u.get_pattern()->comp_def_in_subset(comp, subsetIndex) == false) continue;

		if(write_scalar(File, u, comp, sh, subsetIndex) != true)
		{
			std::cout << "ERROR (in VTKOutput::print(...)): Can not write Scalar Values" << std::endl;
			fclose(File);return false;
		}
	}
	fprintf(File, "      </PointData>\n");

	if(write_piece_epilog(File) != true)
	{
		std::cout << "ERROR (in VTKOutput::print(...)): Can not write Epilog" << std::endl;
		fclose(File);return false;
	}

	return true;
}

bool VTKOutput::write_points(FILE* File, SubsetHandler& sh, uint subsetIndex)
{
	float co;
	int n;

	fprintf(File, "      <Points>\n");
	fprintf(File, "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"binary\">\n");
	BStream.size = sizeof(int);
	n = 3*sizeof(float)*Numbers.noVertices;
	BStreamWrite(File, &n);
	BStreamFlush(File);
	BStream.size = sizeof(float);

	VertexIterator iterBegin, iterEnd, iter;
	iterBegin = sh.begin<Vertex>(subsetIndex);
	iterEnd = sh.end<Vertex>(subsetIndex);

	Grid* grid = sh.get_assigned_grid();
	Grid::VertexAttachmentAccessor<APosition> aaPos(*grid, aPosition);

	assert(m_aaDOFIndexVRT.valid());

	// write points and remember numbering
	vector3 Pos;
	n = 0;
	for(iter = iterBegin; iter != iterEnd; iter++)
	{
		Vertex *v = *iter;
		m_aaDOFIndexVRT[v] = n++;
		Pos = aaPos[v];

		co = Pos[0];
		BStreamWrite(File, &co);
		co = Pos[1];
		BStreamWrite(File, &co);
		co = Pos[2];
		BStreamWrite(File, &co);
	}

	BStreamFlush(File);
	fprintf(File, "\n        </DataArray>\n");
	fprintf(File, "      </Points>\n");

	return true;
}

template <class TElem>
bool VTKOutput::write_elements_connectivity(FILE* File, SubsetHandler& sh, uint subsetIndex)
{
	int id;
	typename geometry_traits<TElem>::Descriptor TDesc;
	typename geometry_traits<TElem>::iterator iterBegin, iterEnd, iter;

	iterBegin = sh.begin<TElem>(subsetIndex);
	iterEnd = sh.end<TElem>(subsetIndex);

	assert(m_aaDOFIndexVRT.valid());

	for(iter = iterBegin; iter != iterEnd; iter++)
	{
		TElem *t = *iter;
		for(int i=0; i< TDesc.num_vertices(); i++)
		{
			VertexBase* vert = t->vertex(i);
			id = m_aaDOFIndexVRT[vert];
			BStreamWrite(File, &id);
		}
	}

	return true;
}

template <class TElem>
bool VTKOutput::write_elements_offsets(FILE* File, SubsetHandler& sh, uint subsetIndex, int& n)
{
	typename geometry_traits<TElem>::iterator iterBegin, iterEnd, iter;
	typename geometry_traits<TElem>::Descriptor TDesc;

	iterBegin = sh.begin<TElem>(subsetIndex);
	iterEnd = sh.end<TElem>(subsetIndex);
	for(iter = iterBegin; iter != iterEnd; iter++)
	{
		n += TDesc.num_vertices();
		BStreamWrite(File, &n);
	}
	BStreamFlush(File);

	return true;
}

template <class TElem>
bool VTKOutput::write_elements_types(FILE* File, SubsetHandler& sh, uint subsetIndex)
{
	char type;
	typename geometry_traits<TElem>::iterator iterBegin, iterEnd, iter;
	typename geometry_traits<TElem>::Descriptor TDesc;

	iterBegin = sh.begin<TElem>(subsetIndex);
	iterEnd = sh.end<TElem>(subsetIndex);

	BStream.size = sizeof(char);
	for(iter = iterBegin; iter != iterEnd; iter++)
	{
		if(TDesc.num_vertices() == 3)
		{
			type = 5;
		}
		else if(TDesc.num_vertices() == 4)
		{
			type = 9;
		}
	    BStreamWrite(File, &type);
	}
	BStreamFlush(File);
	return true;
}


bool VTKOutput::write_elements(FILE* File, SubsetHandler& sh, uint subsetIndex)
{
	int n;

	fprintf(File, "      <Cells>\n");

	/*** connectivity ***/
	fprintf(File, "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"binary\">\n");
	BStream.size = sizeof(int);
	n = sizeof(int)*Numbers.noConnections;
	BStreamWrite(File, &n);
	BStreamFlush(File);

	if(write_elements_connectivity<Triangle>(File, sh, subsetIndex) != true) return false;
	if(write_elements_connectivity<Quadrilateral>(File, sh, subsetIndex) != true) return false;

	BStreamFlush(File);
	fprintf(File, "\n        </DataArray>\n");

	/*** offsets ***/
	fprintf(File, "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"binary\">\n");
	n = sizeof(int)*Numbers.noElements;
	BStreamWrite(File, &n);
	BStreamFlush(File);
	n = 0;
	if(write_elements_offsets<Triangle>(File, sh, subsetIndex, n) == false) return false;
	if(write_elements_offsets<Quadrilateral>(File, sh, subsetIndex, n) == false) return false;
	fprintf(File, "\n        </DataArray>\n");

	/*** types ***/
	fprintf(File, "        <DataArray type=\"Int8\" Name=\"types\" format=\"binary\">\n");
	BStreamWrite(File, &Numbers.noElements);
	BStreamFlush(File);

	if(write_elements_types<Triangle>(File, sh, subsetIndex) == false) return false;
	if(write_elements_types<Quadrilateral>(File, sh, subsetIndex) == false) return false;

	fprintf(File, "\n        </DataArray>\n");

	fprintf(File, "      </Cells>\n");

	return true;
}

bool VTKOutput::write_scalar(FILE* File, NumericalSolution& u, int comp, SubsetHandler& sh, uint subsetIndex)
{
	double val[3];
	int id[3];
	float valf;
	int n;

	fprintf(File, "        <DataArray type=\"Float32\" Name=\"%s\" "
			"NumberOfComponents=\"%d\" format=\"binary\">\n", u.get_name(comp).c_str(), 1);

	BStream.size = sizeof(int);
	n = sizeof(float)*Numbers.noVertices;
	BStreamWrite(File, &n);
	BStreamFlush(File);

	BStream.size = sizeof(float);

	VertexIterator iterBegin, iterEnd, iter;
	iterBegin = sh.begin<Vertex>(subsetIndex);
	iterEnd = sh.end<Vertex>(subsetIndex);

	for(iter = iterBegin; iter != iterEnd; iter++)
	{
		Vertex *v = *iter;
		id[0] = u.get_pattern()->get_index(v, comp);
		if(u.GridVector()->get_values(1, id, val)!=true) return false;
		valf = (float)val[0];
		BStreamWrite(File, &valf);
	}
	BStreamFlush(File);
	fprintf(File, "\n        </DataArray>\n");

	return true;
}


bool VTKOutput::write_prolog(FILE* File, double Time)
{
	fprintf(File, "<?xml version=\"1.0\"?>\n");
	fprintf(File, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" "
			"byte_order=\"%s\">\n",
#ifdef __SWAPBYTES__
			"LittleEndian"
#else
			"BigEndian"
#endif
	);
	fprintf(File, "  <Time timestep=\"%g\"/>\n", Time);
	fprintf(File, "  <UnstructuredGrid>\n");

	return true;
}

bool VTKOutput::write_piece_prolog(FILE* File)
{
	fprintf(File, "    <Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n",
			Numbers.noVertices, Numbers.noElements);

	return true;
}

bool VTKOutput::write_piece_epilog(FILE* File)
{
	fprintf(File, "    </Piece>\n");
	return true;
}

bool VTKOutput::write_epilog(FILE* File)
{
	fprintf(File, "  </UnstructuredGrid>\n");
	fprintf(File, "</VTKFile>\n");

	return true;

}

}
