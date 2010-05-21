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
	m_grid = dynamic_cast<Grid*>(&u.get_domain().get_grid());

	// attach help indices
	m_grid->attach_to_vertices(m_aDOFIndex);
	m_aaDOFIndexVRT.access(*m_grid, m_aDOFIndex);

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
		int dim = 1;
		if( (u.template num<Triangle>(subsetIndex)) > 0 || (u.template num<Quadrilateral>(subsetIndex)) > 0) dim = 2;

		if(write_subset(File, u, subsetIndex, dim)!= true)
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
	m_grid->detach_from_vertices(m_aDOFIndex);
	m_aaDOFIndexVRT.invalidate();

	return true;
}

template <typename TDiscreteFunction>
bool
VTKOutput<TDiscreteFunction>::
print_subset(discrete_function_type& u, int subsetIndex, const char* filename, double Time)
{
	m_grid = dynamic_cast<Grid*>(&u.get_domain().get_grid());

	// attach help indices
	m_grid->attach_to_vertices(m_aDOFIndex);
	m_aaDOFIndexVRT.access(*m_grid, m_aDOFIndex);

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
	int dim = 1;
	if(u.template num<Triangle>(subsetIndex) > 0 || u.template num<Quadrilateral>(subsetIndex) > 0) dim = 2;
	if(write_subset(u, subsetIndex, File, dim) != true)
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
	m_grid->detach_from_vertices(m_aDOFIndex);
	m_aaDOFIndexVRT.invalidate();

	return true;
}


template <typename TDiscreteFunction>
bool
VTKOutput<TDiscreteFunction>::
write_subset(FILE* File, discrete_function_type& u, int subsetIndex, int dim)
{
	// Read sizes
	UG_DLOG(LIB_DISC_OUTPUT, 1, " Initiallizing subset.\n");
	if(init_subset(u, subsetIndex, dim) != true)
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
	if(write_points(File, u, subsetIndex, dim) != true)
	{
		UG_LOG("ERROR (in VTKOutput::print(...)): Can not write Points" << std::endl);
		return false;
	}

	UG_DLOG(LIB_DISC_OUTPUT, 1, " Writing elements.\n");
	if(write_elements(File, u, subsetIndex, dim) != true)
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
		if(write_scalar(File, u, fct, subsetIndex, dim) != true)
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
bool
VTKOutput<TDiscreteFunction>::
init_subset(discrete_function_type& u, int subsetIndex, int dim)
{
	UG_DLOG(LIB_DISC_OUTPUT, 2, "\n ---- Init Numbers ----\n");

	// count elements on subset
	Numbers.noElements = 0; Numbers.noConnections = 0;	Numbers.noVertices = 0;

	if(dim == 1)
	{
		count_elem_conn<Edge>(u, subsetIndex, u.template begin<Edge>(subsetIndex), u.template end<Edge>(subsetIndex));
	}
	if(dim == 2)
	{
		count_elem_conn<Triangle>(u, subsetIndex, u.template begin<Triangle>(subsetIndex), u.template end<Triangle>(subsetIndex));
		count_elem_conn<Quadrilateral>(u, subsetIndex, u.template begin<Quadrilateral>(subsetIndex), u.template end<Quadrilateral>(subsetIndex));
	}

	UG_DLOG(LIB_DISC_OUTPUT, 2, "Number of Vertices: " << Numbers.noVertices << "\n");
	UG_DLOG(LIB_DISC_OUTPUT, 2, "Number of Elements: " << Numbers.noElements << "\n");
	UG_DLOG(LIB_DISC_OUTPUT, 2, "Number of Connections: " << Numbers.noConnections << "\n");
	UG_DLOG(LIB_DISC_OUTPUT, 2, " ---- End ----\n");
	return true;
}


template <typename TDiscreteFunction>
bool
VTKOutput<TDiscreteFunction>::
write_points(FILE* File, discrete_function_type& u, int subsetIndex, int dim)
{

	fprintf(File, "      <Points>\n");
	fprintf(File, "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"binary\">\n");
	BStream.size = sizeof(int);
	int n = 3*sizeof(float)*Numbers.noVertices;
	BStreamWrite(File, &n);
	BStreamFlush(File);
	BStream.size = sizeof(float);

	UG_DLOG(LIB_DISC_OUTPUT, 3, "\n ---- Start: Writing Vertices to file ----\n");

	n = 0;
	switch(dim)
	{
		case 1:
			write_points_elementwise<Edge>(File, u, u.template begin<Edge>(subsetIndex), u.template end<Edge>(subsetIndex), n);
			break;
		case 2:
			write_points_elementwise<Triangle>(File, u, u.template begin<Triangle>(subsetIndex), u.template end<Triangle>(subsetIndex), n);
			write_points_elementwise<Quadrilateral>(File, u, u.template begin<Quadrilateral>(subsetIndex), u.template end<Quadrilateral>(subsetIndex), n);
			break;
	}

	fprintf(File, "\n        </DataArray>\n");
	fprintf(File, "      </Points>\n");

	return true;
}


template <typename TDiscreteFunction>
bool VTKOutput<TDiscreteFunction>::
write_elements(FILE* File, discrete_function_type& u, int subsetIndex, int dim)
{
	int n;

	fprintf(File, "      <Cells>\n");

	/*** connectivity ***/
	fprintf(File, "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"binary\">\n");
	BStream.size = sizeof(int);
	n = sizeof(int)*Numbers.noConnections;
	BStreamWrite(File, &n);
	BStreamFlush(File);

	switch(dim)
	{
		case 1:
			if(write_elements_connectivity<Edge>(File, u.template begin<Edge>(subsetIndex), u.template end<Edge>(subsetIndex)) != true) return false;
			break;
		case 2:
			if(write_elements_connectivity<Triangle>(File, u.template begin<Triangle>(subsetIndex), u.template end<Triangle>(subsetIndex)) != true) return false;
			if(write_elements_connectivity<Quadrilateral>(File,  u.template begin<Quadrilateral>(subsetIndex), u.template end<Quadrilateral>(subsetIndex)) != true) return false;
			break;
	}
	BStreamFlush(File);
	fprintf(File, "\n        </DataArray>\n");


	/*** offsets ***/
	fprintf(File, "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"binary\">\n");
	n = sizeof(int)*Numbers.noElements;
	BStreamWrite(File, &n);
	BStreamFlush(File);
	n = 0;
	switch(dim)
	{
		case 1:
			if(write_elements_offsets<Edge>(File, u.template begin<Edge>(subsetIndex), u.template end<Edge>(subsetIndex),  n) == false) return false;
			break;
		case 2:
			if(write_elements_offsets<Triangle>(File, u.template begin<Triangle>(subsetIndex), u.template end<Triangle>(subsetIndex),  n) == false) return false;
			if(write_elements_offsets<Quadrilateral>(File, u.template begin<Quadrilateral>(subsetIndex), u.template end<Quadrilateral>(subsetIndex), n) == false) return false;
			break;
	}
	BStreamFlush(File);
	fprintf(File, "\n        </DataArray>\n");


	/*** types ***/
	fprintf(File, "        <DataArray type=\"Int8\" Name=\"types\" format=\"binary\">\n");
	BStreamWrite(File, &Numbers.noElements);
	BStreamFlush(File);
	switch(dim)
	{
		case 1:
			if(write_elements_types<Edge>(File, u.template begin<Edge>(subsetIndex), u.template end<Edge>(subsetIndex)) == false) return false;
			break;
		case 2:
			if(write_elements_types<Triangle>(File, u.template begin<Triangle>(subsetIndex), u.template end<Triangle>(subsetIndex)) == false) return false;
			if(write_elements_types<Quadrilateral>(File, u.template begin<Quadrilateral>(subsetIndex), u.template end<Quadrilateral>(subsetIndex)) == false) return false;
			break;
	}
	fprintf(File, "\n        </DataArray>\n");
	fprintf(File, "      </Cells>\n");

	return true;
}


template <typename TDiscreteFunction>
bool
VTKOutput<TDiscreteFunction>::
write_scalar(FILE* File, discrete_function_type& u, uint fct, int subsetIndex, int dim)
{
	int n;

	fprintf(File, "        <DataArray type=\"Float32\" Name=\"%s\" "
			"NumberOfComponents=\"%d\" format=\"binary\">\n", u.get_name(fct).c_str(), 1);

	BStream.size = sizeof(int);
	n = sizeof(float)*Numbers.noVertices;
	BStreamWrite(File, &n);
	BStreamFlush(File);

	switch(dim)
	{
		case 1:
			if(write_scalar_elementwise<Edge>(File, u, fct, u.template begin<Edge>(subsetIndex), u.template end<Edge>(subsetIndex)) == false) return false;
			break;
		case 2:
			if(write_scalar_elementwise<Triangle>(File, u, fct, u.template begin<Triangle>(subsetIndex), u.template end<Triangle>(subsetIndex)) == false) return false;
			if(write_scalar_elementwise<Quadrilateral>(File, u, fct, u.template begin<Quadrilateral>(subsetIndex), u.template end<Quadrilateral>(subsetIndex)) == false) return false;
			break;
	}

	BStreamFlush(File);
	fprintf(File, "\n        </DataArray>\n");

	return true;
};


template <typename TDiscreteFunction>
template <typename TElem>
void
VTKOutput<TDiscreteFunction>::
count_elem_conn(discrete_function_type& u, int subsetIndex,
						typename geometry_traits<TElem>::iterator iterBegin,
						typename geometry_traits<TElem>::iterator iterEnd)
{
	Numbers.noElements += u.template num<TElem>(subsetIndex);
	Numbers.noConnections += u.template num<TElem>(subsetIndex) * reference_element_traits<TElem>::num_corners;

	m_grid->begin_marking();
	for( ; iterBegin != iterEnd; ++iterBegin)
	{
		TElem *elem = *iterBegin;
		for(uint i = 0; i < (uint)reference_element_traits<TElem>::num_corners; ++i)
		{
			VertexBase* v = elem->vertex(i);
			if(m_grid->is_marked(v)) continue;

			Numbers.noVertices++;
			m_grid->mark(v);
		}
	}
	m_grid->end_marking();
};


template <typename TDiscreteFunction>
template <typename TElem>
bool
VTKOutput<TDiscreteFunction>::
write_points_elementwise(FILE* File, discrete_function_type& u,
							typename geometry_traits<TElem>::iterator iterBegin,
							typename geometry_traits<TElem>::iterator iterEnd, int& n)
{
	typedef typename discrete_function_type::domain_type domain_type;
	typename domain_type::position_accessor_type& aaPos = u.get_domain().get_position_accessor();

	// write points and remember numbering
	assert(m_aaDOFIndexVRT.valid());
	typename domain_type::position_type Pos;
	float co;
	m_grid->begin_marking();
	for( ; iterBegin != iterEnd; ++iterBegin)
	{
		TElem *elem = *iterBegin;
		for(uint i = 0; i < (uint) reference_element_traits<TElem>::num_corners; ++i)
		{
			VertexBase* v = elem->vertex(i);
			if(m_grid->is_marked(v)) continue;

			m_grid->mark(v);

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
	}
	UG_DLOG(LIB_DISC_OUTPUT, 3, " ---- " << n << " Vertices (Nr. 0 - " << (n-1) << ") written to file ----\n");
	BStreamFlush(File);
	m_grid->end_marking();

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
	assert(m_aaDOFIndexVRT.valid());

	UG_DLOG(LIB_DISC_OUTPUT, 3, "\n ---- Start: Writing Connectivity (i.e. Vertices of each cell) to file ----");

	for( ; iterBegin != iterEnd; iterBegin++)
	{
		TElem *t = *iterBegin;
		UG_DLOG(LIB_DISC_OUTPUT, 3, "\nWriting Vertices of Finite Element: ");

		for(uint i=0; i< (uint) reference_element_traits<TElem>::num_corners; i++)
		{
			VertexBase* vert = t->vertex(i);
			int id = m_aaDOFIndexVRT[vert];

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
	UG_DLOG(LIB_DISC_OUTPUT, 2, "\n ---- Start: Writing element offsets to file ----\n");

	for( ; iterBegin != iterEnd; ++iterBegin)
	{
		n += reference_element_traits<TElem>::num_corners;
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
	static const ReferenceElementType reftype = reference_element_traits<TElem>::REFERENCE_ELEMENT_TYPE;
	char type;

	// TODO: This is 1D and 2D only
	switch(reftype)
	{
	case RET_EDGE: type = (char) 3; break;
	case RET_TRIANGLE: type = (char) 5; break;
	case RET_QUADRILATERAL: type = (char) 9; break;
	default:UG_ASSERT(0, "Element Type not known.");
	}

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
template <class TElem>
bool
VTKOutput<TDiscreteFunction>::
write_scalar_elementwise(FILE* File,
						discrete_function_type& u, uint fct,
						typename geometry_traits<TElem>::iterator iterBegin,
						typename geometry_traits<TElem>::iterator iterEnd)
{
	float valf;
	BStream.size = sizeof(float);

	UG_DLOG(LIB_DISC_OUTPUT, 3, "\n ---- Start: Writing nodal values to file for function " << u.get_name(fct) << " ----\n");

	m_grid->begin_marking();
	for( ; iterBegin != iterEnd; ++iterBegin)
	{
		TElem *elem = *iterBegin;
		for(uint i = 0; i < (uint) reference_element_traits<TElem>::num_corners; ++i)
		{
			VertexBase* v = elem->vertex(i);
			if(m_grid->is_marked(v)) continue;

			m_grid->mark(v);

			typename discrete_function_type::local_vector_type val(1);
			u.get_dof_values_of_geom_obj(v, fct, val);
			valf = (float) val[0];

			UG_DLOG(LIB_DISC_OUTPUT, 3, "Writing value: " << valf << "\n");
			BStreamWrite(File, &valf);
		}
	}
	m_grid->end_marking();

	UG_DLOG(LIB_DISC_OUTPUT, 3, "\n ---- End  ----\n");
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


