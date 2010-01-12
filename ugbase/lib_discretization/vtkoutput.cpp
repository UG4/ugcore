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
template <int d>
bool VTKOutput<d>::print(NumericalSolution<d>& u, int level, const char* filename, double Time)
{

	ISubsetHandler* ish = u.get_pattern()->get_assigned_subset();
	Grid* grid = ish->get_assigned_grid();

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
	for(int subsetIndex = 0; subsetIndex < ish->num_subsets(); ++subsetIndex)
	{
		if(write_subset(*ish, subsetIndex, u, level, File)!= true)
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

template <int d>
bool VTKOutput<d>::print_subset(NumericalSolution<d>& u, int level, int subsetIndex, const char* filename, double Time)
{

	ISubsetHandler* ish = u.get_pattern()->get_assigned_subset();
	Grid* grid = ish->get_assigned_grid();

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
	if(write_subset(*ish, subsetIndex, u, level, File) != true)
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


template <int d>
bool VTKOutput<d>::write_subset(ISubsetHandler& ish, int subsetIndex, NumericalSolution<d>& u, int level, FILE* File)
{
	// Read sizes
	if(init_subset(ish, subsetIndex, level) != true)
	{
		std::cout << "ERROR (in VTKOutput::print(...)): Can not init subset" << std::endl;
		fclose(File);return false;
	}

	if(write_piece_prolog(File) != true)
	{
		std::cout << "ERROR (in VTKOutput::print(...)): Can not write Prolog" << std::endl;
		fclose(File);return false;
	}

	if(write_points(File, ish, subsetIndex, u.get_domain()->get_position_attachment()) != true)
	{
		std::cout << "ERROR (in VTKOutput::print(...)): Can not write Points" << std::endl;
		fclose(File);return false;
	}
	if(write_elements(File, ish, subsetIndex) != true)
	{
		std::cout << "ERROR (in VTKOutput::print(...)): Can not write Elements" << std::endl;
		fclose(File);return false;
	}

	fprintf(File, "      <PointData>\n");
	for(int comp = 0; comp < u.get_pattern()->num_comp(); ++comp)
	{
		if(u.get_pattern()->comp_def_in_subset(comp, subsetIndex) == false) continue;

		if(write_scalar(File, u, comp, ish, subsetIndex) != true)
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


template <int d>
bool VTKOutput<d>::init_subset(ISubsetHandler& ish, uint subsetIndex, int level)
{
	// TODO: Test Version for 2D only!
	Grid* grid = ish.get_assigned_grid();
	MultiGrid* mg = dynamic_cast<MultiGrid*>(grid);

	if(mg == NULL)
	{
		if(level != 0)
		{
			std::cout << "ERROR (in VTKOutput::init_subset(...)): level > 0 but not a MultiGrid.\n" << std::endl;
			return false;
		}
		_iterBeginVRT = grid->begin<Vertex>();
		_iterEndVRT = grid->end<Vertex>();

		_iterBeginTRIANGLE = grid->begin<Triangle>();
		_iterEndTRIANGLE = grid->end<Triangle>();

		_iterBeginQUADRILATERAL = grid->begin<Quadrilateral>();
		_iterEndQUADRILATERAL = grid->end<Quadrilateral>();
		_level = 0;
	}
	else
	{
		_iterBeginVRT = mg->begin<Vertex>(level);
		_iterEndVRT = mg->end<Vertex>(level);

		_iterBeginTRIANGLE = mg->begin<Triangle>(level);
		_iterEndTRIANGLE = mg->end<Triangle>(level);

		_iterBeginQUADRILATERAL = mg->begin<Quadrilateral>(level);
		_iterEndQUADRILATERAL = mg->end<Quadrilateral>(level);
		_level = level;
	}

	Numbers.noVertices = 0;
	for( _iterVRT = _iterBeginVRT; _iterVRT != _iterEndVRT; ++_iterVRT)
	{
		if(subsetIndex == ish.get_subset_index(*_iterVRT)) Numbers.noVertices++;
	}
	Numbers.noElements = 0;
	Numbers.noConnections = 0;
	for( _iterTRIANGLE = _iterBeginTRIANGLE; _iterTRIANGLE != _iterEndTRIANGLE; ++_iterTRIANGLE)
	{
		if(subsetIndex == ish.get_subset_index(*_iterTRIANGLE))
		{
			Numbers.noElements++;
			Numbers.noConnections += 3;
		}
	}
	for( _iterQUADRILATERAL = _iterBeginQUADRILATERAL; _iterQUADRILATERAL != _iterEndQUADRILATERAL; ++_iterQUADRILATERAL)
	{
		if(subsetIndex == ish.get_subset_index(*_iterQUADRILATERAL))
		{
			Numbers.noElements++;
			Numbers.noConnections += 4;
		}
	}



	return true;
}


template <int d>
bool VTKOutput<d>::write_points(FILE* File, ISubsetHandler& ish, uint subsetIndex, typename Domain<d>::position_attachment_type* aPos)
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

	VertexIterator iter;

	Grid* grid = ish.get_assigned_grid();
	Grid::VertexAttachmentAccessor<typename Domain<d>::position_attachment_type > aaPos(*grid, *aPos);

	assert(m_aaDOFIndexVRT.valid());

	// write points and remember numbering
	MathVector<d> Pos;
	n = 0;
	for(iter = _iterBeginVRT; iter != _iterEndVRT; iter++)
	{
		if(subsetIndex != ish.get_subset_index(*iter)) continue;

		Vertex *v = *iter;
		m_aaDOFIndexVRT[v] = n++;
		Pos = aaPos[v];

		for(uint i = 0; i < d; ++i)
		{
			co = Pos[i];
			BStreamWrite(File, &co);
		}
		for(uint i = d; i < 3; ++i)
		{
			co = 0.0;
			BStreamWrite(File, &co);
		}
	}

	BStreamFlush(File);
	fprintf(File, "\n        </DataArray>\n");
	fprintf(File, "      </Points>\n");

	return true;
}

template <int d>
template <class TElem>
bool VTKOutput<d>::write_elements_connectivity(FILE* File,
											typename geometry_traits<TElem>::iterator iterBegin,
											typename geometry_traits<TElem>::iterator iterEnd,
											ISubsetHandler& ish, uint subsetIndex)
{
	int id;
	typename geometry_traits<TElem>::Descriptor TDesc;
	typename geometry_traits<TElem>::iterator /*iterBegin, iterEnd,*/ iter;

	//iterBegin = sh.begin<TElem>(subsetIndex);
	//iterEnd = sh.end<TElem>(subsetIndex);

	assert(m_aaDOFIndexVRT.valid());

	for(iter = iterBegin; iter != iterEnd; iter++)
	{
		if(subsetIndex != ish.get_subset_index(*iter)) continue;

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

template <int d>
template <class TElem>
bool VTKOutput<d>::write_elements_offsets(	FILE* File,
										typename geometry_traits<TElem>::iterator iterBegin,
										typename geometry_traits<TElem>::iterator iterEnd,
										ISubsetHandler& ish, uint subsetIndex, int& n)
{
	typename geometry_traits<TElem>::iterator /*iterBegin, iterEnd,*/ iter;
	typename geometry_traits<TElem>::Descriptor TDesc;

	//iterBegin = sh.begin<TElem>(subsetIndex);
	//iterEnd = sh.end<TElem>(subsetIndex);
	for(iter = iterBegin; iter != iterEnd; iter++)
	{
		if(subsetIndex != ish.get_subset_index(*iter)) continue;

		n += TDesc.num_vertices();
		BStreamWrite(File, &n);
	}
	BStreamFlush(File);

	return true;
}

template <int d>
template <class TElem>
bool VTKOutput<d>::write_elements_types(FILE* File,
									typename geometry_traits<TElem>::iterator iterBegin,
									typename geometry_traits<TElem>::iterator iterEnd,
									ISubsetHandler& ish, uint subsetIndex)
{
	char type;
	typename geometry_traits<TElem>::iterator iter;
	typename geometry_traits<TElem>::Descriptor TDesc;

	BStream.size = sizeof(char);
	for(iter = iterBegin; iter != iterEnd; iter++)
	{
		if(subsetIndex != ish.get_subset_index(*iter)) continue;

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


template <int d>
bool VTKOutput<d>::write_elements(FILE* File, ISubsetHandler& ish, uint subsetIndex)
{
	int n;

	fprintf(File, "      <Cells>\n");

	/*** connectivity ***/
	fprintf(File, "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"binary\">\n");
	BStream.size = sizeof(int);
	n = sizeof(int)*Numbers.noConnections;
	BStreamWrite(File, &n);
	BStreamFlush(File);

	if(write_elements_connectivity<Triangle>(File, _iterBeginTRIANGLE, _iterEndTRIANGLE,  ish, subsetIndex) != true) return false;
	if(write_elements_connectivity<Quadrilateral>(File, _iterBeginQUADRILATERAL, _iterEndQUADRILATERAL, ish, subsetIndex) != true) return false;

	BStreamFlush(File);
	fprintf(File, "\n        </DataArray>\n");

	/*** offsets ***/
	fprintf(File, "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"binary\">\n");
	n = sizeof(int)*Numbers.noElements;
	BStreamWrite(File, &n);
	BStreamFlush(File);
	n = 0;
	if(write_elements_offsets<Triangle>(File, _iterBeginTRIANGLE, _iterEndTRIANGLE,  ish, subsetIndex,  n) == false) return false;
	if(write_elements_offsets<Quadrilateral>(File, _iterBeginQUADRILATERAL, _iterEndQUADRILATERAL,ish, subsetIndex, n) == false) return false;
	fprintf(File, "\n        </DataArray>\n");

	/*** types ***/
	fprintf(File, "        <DataArray type=\"Int8\" Name=\"types\" format=\"binary\">\n");
	BStreamWrite(File, &Numbers.noElements);
	BStreamFlush(File);

	if(write_elements_types<Triangle>(File, _iterBeginTRIANGLE, _iterEndTRIANGLE, ish, subsetIndex) == false) return false;
	if(write_elements_types<Quadrilateral>(File,  _iterBeginQUADRILATERAL, _iterEndQUADRILATERAL, ish, subsetIndex) == false) return false;

	fprintf(File, "\n        </DataArray>\n");

	fprintf(File, "      </Cells>\n");

	return true;
}

template <int d>
bool VTKOutput<d>::write_scalar(FILE* File, NumericalSolution<d>& u, int comp, ISubsetHandler& ish, uint subsetIndex)
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

	VertexIterator iter;

	for(iter = _iterBeginVRT; iter != _iterEndVRT; iter++)
	{
		if(subsetIndex != ish.get_subset_index(*iter)) continue;

		Vertex *v = *iter;
		id[0] = u.get_pattern()->get_index(v, comp);
		if(u.GridVector(_level)->get_values(1, id, val)!=true) return false;
		valf = (float)val[0];
		BStreamWrite(File, &valf);
	}
	BStreamFlush(File);
	fprintf(File, "\n        </DataArray>\n");

	return true;
}


template <int d>
bool VTKOutput<d>::write_prolog(FILE* File, double Time)
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

template <int d>
bool VTKOutput<d>::write_piece_prolog(FILE* File)
{
	fprintf(File, "    <Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n",
			Numbers.noVertices, Numbers.noElements);

	return true;
}

template <int d>
bool VTKOutput<d>::write_piece_epilog(FILE* File)
{
	fprintf(File, "    </Piece>\n");
	return true;
}

template <int d>
bool VTKOutput<d>::write_epilog(FILE* File)
{
	fprintf(File, "  </UnstructuredGrid>\n");
	fprintf(File, "</VTKFile>\n");

	return true;

}

// force code generation
template class VTKOutput<1>;
template class VTKOutput<2>;
template class VTKOutput<3>;


}


