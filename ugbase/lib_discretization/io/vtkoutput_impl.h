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
template <typename TDiscreteFunction>
bool
VTKOutput<TDiscreteFunction>::
print(discrete_function_type& u, const char* filename, double Time)
{
	Grid* grid = dynamic_cast<Grid*>(&u.get_domain().get_grid());

	// attach help indices
	grid->attach_to_vertices(m_aDOFIndex);
	m_aaDOFIndexVRT.access(*grid, m_aDOFIndex);

	// open stream
	FILE* File = fopen(filename, "w");
	if(File == NULL)
	{
		UG_LOG("ERROR (in VTKOutput::print(...)): Can not open Output File" << std::endl);
		return false;
	}

	BStream.front = 0;

	// Write to File
	if(write_prolog(File, Time) != true)
	{
		UG_LOG("ERROR (in VTKOutput::print(...)): Can not write Prolog" << std::endl);
		fclose(File);
		return false;
	}

	// loop subsets
	for(int subsetIndex = 0; subsetIndex < u.num_subsets(); ++subsetIndex)
	{
		if(write_subset(File, u, subsetIndex)!= true)
		{
			UG_LOG("ERROR (in VTKOutput::print(...)): Can not write Subset" << std::endl);
			fclose(File);
			return false;
		}
	}

	if(write_epilog(File) != true)
	{
		UG_LOG("ERROR (in VTKOutput::print(...)): Can not write Epilog" << std::endl);
		fclose(File);
		return false;
	}

	// close stream
	fclose(File);

	// detach help indices
	grid->detach_from_vertices(m_aDOFIndex);
	m_aaDOFIndexVRT.invalidate();

	return true;
}

template <typename TDiscreteFunction>
bool
VTKOutput<TDiscreteFunction>::
print_subset(discrete_function_type& u, int subsetIndex, const char* filename, double Time)
{
	Grid* grid = dynamic_cast<Grid*>(&u.get_domain().get_grid());

	// attach help indices
	grid->attach_to_vertices(m_aDOFIndex);
	m_aaDOFIndexVRT.access(*grid, m_aDOFIndex);

	// open stream
	FILE* File = fopen(filename, "w");
	if(File == NULL)
	{
		UG_LOG("ERROR (in VTKOutput::print(...)): Can not open Output File" << std::endl);
		return false;
	}

	BStream.front = 0;

	// Write to File
	if(write_prolog(File, Time) != true)
	{
		UG_LOG("ERROR (in VTKOutput::print(...)): Can not write Prolog" << std::endl);
		fclose(File);
		return false;
	}

	// Write Subset
	if(write_subset(u, subsetIndex, File) != true)
	{
		UG_LOG("ERROR (in VTKOutput::print(...)): Can not write Subset" << std::endl);
		fclose(File);
		return false;
	}

	if(write_epilog(File) != true)
	{
		UG_LOG("ERROR (in VTKOutput::print(...)): Can not write Epilog" << std::endl);
		fclose(File);
		return false;
	}

	// close stream
	fclose(File);

	// detach help indices
	grid->detach_from_vertices(m_aDOFIndex);
	m_aaDOFIndexVRT.invalidate();

	return true;
}


template <typename TDiscreteFunction>
bool
VTKOutput<TDiscreteFunction>::
write_subset(FILE* File, discrete_function_type& u, int subsetIndex)
{
	// Read sizes
	UG_DLOG(LIB_DISC_OUTPUT, 1, " Initiallizing subset.\n");
	if(init_subset(u, subsetIndex) != true)
	{
		UG_LOG("ERROR (in VTKOutput::print(...)): Can not init subset" << std::endl);
		return false;
	}

	UG_DLOG(LIB_DISC_OUTPUT, 1, " Writing prolog.\n");
	if(write_piece_prolog(File) != true)
	{
		UG_LOG("ERROR (in VTKOutput::print(...)): Can not write Prolog" << std::endl);
		return false;
	}

	UG_DLOG(LIB_DISC_OUTPUT, 1, " Writing points.\n");
	if(write_points(File, u, subsetIndex) != true)
	{
		UG_LOG("ERROR (in VTKOutput::print(...)): Can not write Points" << std::endl);
		return false;
	}

	UG_DLOG(LIB_DISC_OUTPUT, 1, " Writing elements.\n");
	if(write_elements(File, u, subsetIndex) != true)
	{
		UG_LOG("ERROR (in VTKOutput::print(...)): Can not write Elements" << std::endl);
		return false;
	}

	UG_DLOG(LIB_DISC_OUTPUT, 1, " Writing nodal values of " << u.num_fct() << " solutions.\n");
	fprintf(File, "      <PointData>\n");
	for(uint fct = 0; fct < u.num_fct(); ++fct)
	{
		if(u.fct_is_def_in_subset(fct, subsetIndex) == false) continue;

		UG_DLOG(LIB_DISC_OUTPUT, 1, "  - Writing nodal values of solutions ' " << u.get_name(fct) <<"'.\n");
		if(write_scalar(File, u, fct, subsetIndex) != true)
		{
			UG_LOG("ERROR (in VTKOutput::print(...)): Can not write Scalar Values" << std::endl);
			return false;
		}
	}
	fprintf(File, "      </PointData>\n");

	if(write_piece_epilog(File) != true)
	{
		UG_LOG("ERROR (in VTKOutput::print(...)): Can not write Epilog" << std::endl);
		return false;
	}

	return true;
}

template <typename TDiscreteFunction>
template <typename TElem>
bool
VTKOutput<TDiscreteFunction>::
count_elem_conn(int& num_elem, int& num_connections,
				typename geometry_traits<TElem>::iterator iterBegin,
				typename geometry_traits<TElem>::iterator iterEnd)
{
	num_elem = num_connections = 0;
	for( ; iterBegin != iterEnd; ++iterBegin)
	{
			num_elem++;
			num_connections += reference_element_traits<TElem>::num_corners;
	}
	return true;
}

template <typename TDiscreteFunction>
bool
VTKOutput<TDiscreteFunction>::
init_subset(discrete_function_type& u, int subsetIndex)
{
	UG_DLOG(LIB_DISC_OUTPUT, 2, "\n ---- Init Numbers ----\n");

	typename geometry_traits<VertexBase>::iterator iterVRT, iterEndVRT;
	iterEndVRT = u.template end<Vertex>(subsetIndex);

	// count vertices on subset
	Numbers.noVertices = 0;
	for( iterVRT = u.template begin<Vertex>(subsetIndex); iterVRT != iterEndVRT; ++iterVRT)
	{
		Numbers.noVertices++;
	}
	UG_DLOG(LIB_DISC_OUTPUT, 2, "Number of Vertices: " << Numbers.noVertices << "\n");

	// count elements on subset
	int num_elem = 0, num_connections = 0;
	Numbers.noElements = 0;
	Numbers.noConnections = 0;
	if(count_elem_conn<Triangle>(num_elem, num_connections, u.template begin<Triangle>(subsetIndex), u.template end<Triangle>(subsetIndex)) == true)
	{
		Numbers.noElements += num_elem;
		Numbers.noConnections += num_connections;
	}
	else return false;
	if(count_elem_conn<Quadrilateral>(num_elem, num_connections, u.template begin<Quadrilateral>(subsetIndex), u.template end<Quadrilateral>(subsetIndex)) == true)
	{
		Numbers.noElements += num_elem;
		Numbers.noConnections += num_connections;
	}
	else return false;
	UG_DLOG(LIB_DISC_OUTPUT, 2, "Number of Elements: " << Numbers.noElements << "\n");
	UG_DLOG(LIB_DISC_OUTPUT, 2, "Number of Connections: " << Numbers.noConnections << "\n");

	UG_DLOG(LIB_DISC_OUTPUT, 2, " ---- End ----\n");
	return true;
}


template <typename TDiscreteFunction>
bool
VTKOutput<TDiscreteFunction>::
write_points(FILE* File, discrete_function_type& u, int subsetIndex)
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

	typedef typename discrete_function_type::domain_type domain_type;
	typename domain_type::position_accessor_type& aaPos = u.get_domain().get_position_accessor();

	assert(m_aaDOFIndexVRT.valid());

	// write points and remember numbering
	typename domain_type::position_type Pos;
	n = 0;
	typename geometry_traits<VertexBase>::iterator iterVRT, iterEndVRT;
	iterEndVRT = u.template end<Vertex>(subsetIndex);

	UG_DLOG(LIB_DISC_OUTPUT, 3, "\n ---- Start: Writing Vertices to file ----\n");

	for(iterVRT = u.template begin<Vertex>(subsetIndex); iterVRT != iterEndVRT; ++iterVRT)
	{
		VertexBase *v = *iterVRT;
		m_aaDOFIndexVRT[v] = n++;
		Pos = aaPos[v];

		UG_DLOG(LIB_DISC_OUTPUT, 3, "Writing Vertex Nr. " << n-1 << " with position " << ".\n");

		for(int i = 0; i < domain_type::dim; ++i)
		{
			co = Pos[i];
			BStreamWrite(File, &co);
		}
		for(int i = domain_type::dim; i < 3; ++i)
		{
			co = 0.0;
			BStreamWrite(File, &co);
		}
	}

	UG_DLOG(LIB_DISC_OUTPUT, 3, " ---- " << n << " Vertices (Nr. 0 - " << (n-1) << ") written to file ----\n");

	BStreamFlush(File);
	fprintf(File, "\n        </DataArray>\n");
	fprintf(File, "      </Points>\n");

	return true;
}

template <typename TDiscreteFunction>
template <class TElem>
bool
VTKOutput<TDiscreteFunction>::
write_elements_connectivity(FILE* File,
							typename geometry_traits<TElem>::iterator iterBegin,
							typename geometry_traits<TElem>::iterator iterEnd)
{
	int id;
	typename geometry_traits<TElem>::Descriptor TDesc;

	assert(m_aaDOFIndexVRT.valid());

	UG_DLOG(LIB_DISC_OUTPUT, 3, "\n ---- Start: Writing Connectivity (i.e. Vertices of each cell) to file ----");

	for( ; iterBegin != iterEnd; iterBegin++)
	{
		TElem *t = *iterBegin;
		UG_DLOG(LIB_DISC_OUTPUT, 3, "\nWriting Vertices of Finite Element: ");

		for(uint i=0; i< TDesc.num_vertices(); i++)
		{
			VertexBase* vert = t->vertex(i);
			id = m_aaDOFIndexVRT[vert];

			UG_DLOG(LIB_DISC_OUTPUT, 3, id << " ");

			BStreamWrite(File, &id);
		}
	}

	UG_DLOG(LIB_DISC_OUTPUT, 3, "\n ---- End ----\n");

	return true;
}

template <typename TDiscreteFunction>
template <class TElem>
bool
VTKOutput<TDiscreteFunction>::
write_elements_offsets(	FILE* File,
										typename geometry_traits<TElem>::iterator iterBegin,
										typename geometry_traits<TElem>::iterator iterEnd,
										int& n)
{
	typename geometry_traits<TElem>::Descriptor TDesc;

	UG_DLOG(LIB_DISC_OUTPUT, 2, "\n ---- Start: Writing element offsets to file ----\n");

	for( ; iterBegin != iterEnd; ++iterBegin)
	{
		n += TDesc.num_vertices();
		UG_DLOG(LIB_DISC_OUTPUT, 3, n << " ");
		BStreamWrite(File, &n);
	}
	BStreamFlush(File);

	UG_DLOG(LIB_DISC_OUTPUT, 2, "\n ---- End ----\n");

	return true;
}

template <typename TDiscreteFunction>
template <class TElem>
bool
VTKOutput<TDiscreteFunction>::
write_elements_types(FILE* File,
					typename geometry_traits<TElem>::iterator iterBegin,
					typename geometry_traits<TElem>::iterator iterEnd)
{
	char type;
	typename geometry_traits<TElem>::Descriptor TDesc;

	// TODO: This is 2D only
	if(TDesc.num_vertices() == 3)
	{
		type = (char) 5;
	}
	else if(TDesc.num_vertices() == 4)
	{
		type = (char) 9;
	}
	else UG_ASSERT(0, "Element Type not known.\n");

	UG_DLOG(LIB_DISC_OUTPUT, 3, "\n ---- Start: Writing element types (VTK element types) to file ----\n");
	BStream.size = sizeof(char);
	for( ; iterBegin != iterEnd; ++iterBegin)
	{
		UG_DLOG(LIB_DISC_OUTPUT, 3, (int) type << ",");
	    BStreamWrite(File, &type);
	}
	BStreamFlush(File);

	UG_DLOG(LIB_DISC_OUTPUT, 3, "\n ---- End ----\n");

	return true;
}


template <typename TDiscreteFunction>
bool VTKOutput<TDiscreteFunction>::
write_elements(FILE* File, discrete_function_type& u, int subsetIndex)
{
	int n;

	fprintf(File, "      <Cells>\n");

	/*** connectivity ***/
	fprintf(File, "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"binary\">\n");
	BStream.size = sizeof(int);
	n = sizeof(int)*Numbers.noConnections;
	BStreamWrite(File, &n);
	BStreamFlush(File);

	if(write_elements_connectivity<Triangle>(File, u.template begin<Triangle>(subsetIndex), u.template end<Triangle>(subsetIndex)) != true) return false;
	if(write_elements_connectivity<Quadrilateral>(File,  u.template begin<Quadrilateral>(subsetIndex), u.template end<Quadrilateral>(subsetIndex)) != true) return false;

	BStreamFlush(File);
	fprintf(File, "\n        </DataArray>\n");

	/*** offsets ***/
	fprintf(File, "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"binary\">\n");
	n = sizeof(int)*Numbers.noElements;
	BStreamWrite(File, &n);
	BStreamFlush(File);
	n = 0;
	if(write_elements_offsets<Triangle>(File, u.template begin<Triangle>(subsetIndex), u.template end<Triangle>(subsetIndex),  n) == false) return false;
	if(write_elements_offsets<Quadrilateral>(File, u.template begin<Quadrilateral>(subsetIndex), u.template end<Quadrilateral>(subsetIndex), n) == false) return false;
	fprintf(File, "\n        </DataArray>\n");

	/*** types ***/
	fprintf(File, "        <DataArray type=\"Int8\" Name=\"types\" format=\"binary\">\n");
	BStreamWrite(File, &Numbers.noElements);
	BStreamFlush(File);

	if(write_elements_types<Triangle>(File, u.template begin<Triangle>(subsetIndex), u.template end<Triangle>(subsetIndex)) == false) return false;
	if(write_elements_types<Quadrilateral>(File, u.template begin<Quadrilateral>(subsetIndex), u.template end<Quadrilateral>(subsetIndex)) == false) return false;

	fprintf(File, "\n        </DataArray>\n");

	fprintf(File, "      </Cells>\n");

	return true;
}

template <typename TDiscreteFunction>
bool
VTKOutput<TDiscreteFunction>::
write_scalar(FILE* File, discrete_function_type& u, uint fct, int subsetIndex)
{
	float valf;
	int n;

	fprintf(File, "        <DataArray type=\"Float32\" Name=\"%s\" "
			"NumberOfComponents=\"%d\" format=\"binary\">\n", u.get_name(fct).c_str(), 1);

	BStream.size = sizeof(int);
	n = sizeof(float)*Numbers.noVertices;
	BStreamWrite(File, &n);
	BStreamFlush(File);

	BStream.size = sizeof(float);

	typename geometry_traits<VertexBase>::iterator iterVRT, iterEndVRT;
	iterEndVRT = u.template end<Vertex>(subsetIndex);

	UG_DLOG(LIB_DISC_OUTPUT, 3, "\n ---- Start: Writing nodal values to file for function " << u.get_name(fct) << " ----\n");

	for(iterVRT = u.template begin<Vertex>(subsetIndex); iterVRT != iterEndVRT; ++iterVRT)
	{
		typename discrete_function_type::local_vector_type val(1);
		VertexBase *v = *iterVRT;
		u.get_dof_values_of_geom_obj(v, fct, val);
		valf = (float) val[0];

		UG_DLOG(LIB_DISC_OUTPUT, 3, "Writing value: " << valf << "\n");
		BStreamWrite(File, &valf);
	}

	UG_DLOG(LIB_DISC_OUTPUT, 3, "\n ---- End  ----\n");

	BStreamFlush(File);
	fprintf(File, "\n        </DataArray>\n");

	return true;
}


template <typename TDiscreteFunction>
bool
VTKOutput<TDiscreteFunction>::
write_prolog(FILE* File, double Time)
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

template <typename TDiscreteFunction>
bool
VTKOutput<TDiscreteFunction>::
write_piece_prolog(FILE* File)
{
	fprintf(File, "    <Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n",
			Numbers.noVertices, Numbers.noElements);

	return true;
}

template <typename TDiscreteFunction>
bool
VTKOutput<TDiscreteFunction>::
write_piece_epilog(FILE* File)
{
	fprintf(File, "    </Piece>\n");
	return true;
}

template <typename TDiscreteFunction>
bool
VTKOutput<TDiscreteFunction>::
write_epilog(FILE* File)
{
	fprintf(File, "  </UnstructuredGrid>\n");
	fprintf(File, "</VTKFile>\n");

	return true;

}

}


