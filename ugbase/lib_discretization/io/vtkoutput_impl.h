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

template <typename TDiscreteFunction>
bool
VTKOutput<TDiscreteFunction>::
begin_timeseries(const char* filename, discrete_function_type& u)
{
	m_u = &u;
	strcpy(m_seriesname, filename);
	m_timestep.clear();
	return true;
}

template <typename TDiscreteFunction>
bool
VTKOutput<TDiscreteFunction>::
end_timeseries(const char* filename, discrete_function_type& u)
{
	// not same name
	if(strcmp(m_seriesname, filename) != 0) return false;
	if(m_u != &u) return false;

#ifdef UG_PARALLEL
	if(write_time_pvd(u, filename) != true)
	{
		UG_LOG("ERROR (in VTKOutput::print(...)): Can not write pvd - file. \n");
		return false;
	}
#endif

	return true;
}


/* END: Helper Functions */
template <typename TDiscreteFunction>
bool
VTKOutput<TDiscreteFunction>::
print(const char* filename, discrete_function_type& u, size_t step, number time)
{
	// loop subsets
	for(int si = 0; si < u.num_subsets(); ++si)
	{
		if(print_subset(filename, u, si, step, time)!= true)
		{
			UG_LOG("ERROR (in VTKOutput::print(...)): Can not write Subset "<< si << ".\n");
			return false;
		}
	}

	if(!write_pvd(u, filename, step, time))
		{UG_LOG("Cannot write pvd file.\n"); return false;}

	return true;
}

template <typename TDiscreteFunction>
bool
VTKOutput<TDiscreteFunction>::
print_subset(const char* filename, discrete_function_type& u, int si, size_t step, number time)
{
	m_grid = dynamic_cast<Grid*>(&u.get_approximation_space().get_domain().get_grid());

	// attach help indices
	m_grid->attach_to_vertices(m_aDOFIndex);
	m_aaDOFIndexVRT.access(*m_grid, m_aDOFIndex);

	// open stream
	char name[NAMESIZE];

#ifdef UG_PARALLEL
	int rank = pcl::GetProcRank();
#else
	int rank = 0;
#endif
	if(vtu_filename(name, filename, rank, si, step) != true) return false;

	FILE* File = fopen(name, "w");
	if(File == NULL)
	{
		UG_LOG("ERROR (in VTKOutput::print(...)): Can not open Output File" << std::endl);
		return false;
	}

	BStream.front = 0;

	// Write to File
	if(write_prolog(File, time) != true)
	{
		UG_LOG("ERROR (in VTKOutput::print(...)): Can not write Prolog" << std::endl);
		fclose(File);
		return false;
	}

	// Write Subset
	int dim = 1;
	if(u.template num<Triangle>(si) > 0 || u.template num<Quadrilateral>(si) > 0) dim = 2;
	if(u.template num<Tetrahedron>(si) > 0 || u.template num<Pyramid>(si) > 0 ||
		u.template num<Prism>(si) > 0	|| u.template num<Hexahedron>(si) > 0) dim = 3;
	if(write_subset(File, u, si, dim) != true)
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


#ifdef UG_PARALLEL
	if(write_pvtu(u, filename, si, step, time) != true)
	{
		UG_LOG("ERROR (in VTKOutput::print(...)): Can not write pvtu - file" << std::endl);
		fclose(File);
		return false;
	}

	if((m_u == &u) && (strcmp(filename, m_seriesname) == 0))
		if(m_timestep.empty() || (m_timestep.back() != time))
			m_timestep.push_back(time);
#endif

	return true;
}


template <typename TDiscreteFunction>
bool
VTKOutput<TDiscreteFunction>::
write_subset(FILE* File, discrete_function_type& u, int si, int dim)
{
	// Read sizes
	UG_DLOG(LIB_DISC_OUTPUT, 1, " Initiallizing subset.\n");
	if(init_subset(u, si, dim) != true)
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
	if(write_points(File, u, si, dim) != true)
	{
		UG_LOG("ERROR (in VTKOutput::print(...)): Can not write Points" << std::endl);
		return false;
	}

	UG_DLOG(LIB_DISC_OUTPUT, 1, " Writing elements.\n");
	if(write_elements(File, u, si, dim) != true)
	{
		UG_LOG("ERROR (in VTKOutput::print(...)): Can not write Elements" << std::endl);
		return false;
	}

	UG_DLOG(LIB_DISC_OUTPUT, 1, " Writing nodal values of " << u.num_fct() << " solutions.\n");
	fprintf(File, "      <PointData>\n");
	for(uint fct = 0; fct < u.num_fct(); ++fct)
	{
		if(u.is_def_in_subset(fct, si) == false) continue;

		UG_DLOG(LIB_DISC_OUTPUT, 1, "  - Writing nodal values of solutions ' " << u.name(fct) <<"'.\n");
		if(write_scalar(File, u, fct, si, dim) != true)
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
init_subset(discrete_function_type& u, int si, int dim)
{
	UG_DLOG(LIB_DISC_OUTPUT, 2, "\n ---- Init Numbers ----\n");

	// count elements on subset
	Numbers.noElements = 0; Numbers.noConnections = 0;	Numbers.noVertices = 0;

	if(dim == 1)
	{
		count_elem_conn<Edge>(u, si, u.template begin<Edge>(si), u.template end<Edge>(si));
	}
	if(dim == 2)
	{
		count_elem_conn<Triangle>(u, si, u.template begin<Triangle>(si), u.template end<Triangle>(si));
		count_elem_conn<Quadrilateral>(u, si, u.template begin<Quadrilateral>(si), u.template end<Quadrilateral>(si));
	}
	if(dim == 3)
	{
		count_elem_conn<Tetrahedron>(u, si, u.template begin<Tetrahedron>(si), u.template end<Tetrahedron>(si));
		count_elem_conn<Pyramid>(u, si, u.template begin<Pyramid>(si), u.template end<Pyramid>(si));
		count_elem_conn<Prism>(u, si, u.template begin<Prism>(si), u.template end<Prism>(si));
		count_elem_conn<Hexahedron>(u, si, u.template begin<Hexahedron>(si), u.template end<Hexahedron>(si));
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
write_points(FILE* File, discrete_function_type& u, int si, int dim)
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
			write_points_elementwise<Edge>(File, u, u.template begin<Edge>(si), u.template end<Edge>(si), n);
			break;
		case 2:
			write_points_elementwise<Triangle>(File, u, u.template begin<Triangle>(si), u.template end<Triangle>(si), n);
			write_points_elementwise<Quadrilateral>(File, u, u.template begin<Quadrilateral>(si), u.template end<Quadrilateral>(si), n);
			break;
		case 3:
			write_points_elementwise<Tetrahedron>(File, u, u.template begin<Tetrahedron>(si), u.template end<Tetrahedron>(si), n);
			write_points_elementwise<Pyramid>(File, u, u.template begin<Pyramid>(si), u.template end<Pyramid>(si), n);
			write_points_elementwise<Prism>(File, u, u.template begin<Prism>(si), u.template end<Prism>(si), n);
			write_points_elementwise<Hexahedron>(File, u, u.template begin<Hexahedron>(si), u.template end<Hexahedron>(si), n);
			break;
	}

	fprintf(File, "\n        </DataArray>\n");
	fprintf(File, "      </Points>\n");

	return true;
}


template <typename TDiscreteFunction>
bool VTKOutput<TDiscreteFunction>::
write_elements(FILE* File, discrete_function_type& u, int si, int dim)
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
			if(write_elements_connectivity<Edge>(File, u.template begin<Edge>(si), u.template end<Edge>(si)) != true) return false;
			break;
		case 2:
			if(write_elements_connectivity<Triangle>(File, u.template begin<Triangle>(si), u.template end<Triangle>(si)) != true) return false;
			if(write_elements_connectivity<Quadrilateral>(File,  u.template begin<Quadrilateral>(si), u.template end<Quadrilateral>(si)) != true) return false;
			break;
		case 3:
			if(write_elements_connectivity<Tetrahedron>(File, u.template begin<Tetrahedron>(si), u.template end<Tetrahedron>(si)) != true) return false;
			if(write_elements_connectivity<Pyramid>(File, u.template begin<Pyramid>(si), u.template end<Pyramid>(si)) != true) return false;
			if(write_elements_connectivity<Prism>(File, u.template begin<Prism>(si), u.template end<Prism>(si)) != true) return false;
			if(write_elements_connectivity<Hexahedron>(File, u.template begin<Hexahedron>(si), u.template end<Hexahedron>(si)) != true) return false;
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
			if(write_elements_offsets<Edge>(File, u.template begin<Edge>(si), u.template end<Edge>(si),  n) == false) return false;
			break;
		case 2:
			if(write_elements_offsets<Triangle>(File, u.template begin<Triangle>(si), u.template end<Triangle>(si),  n) == false) return false;
			if(write_elements_offsets<Quadrilateral>(File, u.template begin<Quadrilateral>(si), u.template end<Quadrilateral>(si), n) == false) return false;
			break;
		case 3:
			if(write_elements_offsets<Tetrahedron>(File, u.template begin<Tetrahedron>(si), u.template end<Tetrahedron>(si),  n) == false) return false;
			if(write_elements_offsets<Pyramid>(File, u.template begin<Pyramid>(si), u.template end<Pyramid>(si),  n) == false) return false;
			if(write_elements_offsets<Prism>(File, u.template begin<Prism>(si), u.template end<Prism>(si),  n) == false) return false;
			if(write_elements_offsets<Hexahedron>(File, u.template begin<Hexahedron>(si), u.template end<Hexahedron>(si),  n) == false) return false;
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
			if(write_elements_types<Edge>(File, u.template begin<Edge>(si), u.template end<Edge>(si)) == false) return false;
			break;
		case 2:
			if(write_elements_types<Triangle>(File, u.template begin<Triangle>(si), u.template end<Triangle>(si)) == false) return false;
			if(write_elements_types<Quadrilateral>(File, u.template begin<Quadrilateral>(si), u.template end<Quadrilateral>(si)) == false) return false;
			break;
		case 3:
			if(write_elements_types<Tetrahedron>(File, u.template begin<Tetrahedron>(si), u.template end<Tetrahedron>(si)) == false) return false;
			if(write_elements_types<Pyramid>(File, u.template begin<Pyramid>(si), u.template end<Pyramid>(si)) == false) return false;
			if(write_elements_types<Prism>(File, u.template begin<Prism>(si), u.template end<Prism>(si)) == false) return false;
			if(write_elements_types<Hexahedron>(File, u.template begin<Hexahedron>(si), u.template end<Hexahedron>(si)) == false) return false;
			break;
	}
	fprintf(File, "\n        </DataArray>\n");
	fprintf(File, "      </Cells>\n");

	return true;
}


template <typename TDiscreteFunction>
bool
VTKOutput<TDiscreteFunction>::
write_scalar(FILE* File, discrete_function_type& u, uint fct, int si, int dim)
{
	int n;

	fprintf(File, "        <DataArray type=\"Float32\" Name=\"%s\" "
			"NumberOfComponents=\"%d\" format=\"binary\">\n", u.name(fct).c_str(), 1);

	BStream.size = sizeof(int);
	n = sizeof(float)*Numbers.noVertices;
	BStreamWrite(File, &n);
	BStreamFlush(File);

	switch(dim)
	{
		case 1:
			if(write_scalar_elementwise<Edge>(File, u, fct, u.template begin<Edge>(si), u.template end<Edge>(si), si) == false) return false;
			break;
		case 2:
			if(write_scalar_elementwise<Triangle>(File, u, fct, u.template begin<Triangle>(si), u.template end<Triangle>(si), si) == false) return false;
			if(write_scalar_elementwise<Quadrilateral>(File, u, fct, u.template begin<Quadrilateral>(si), u.template end<Quadrilateral>(si), si) == false) return false;
			break;
		case 3:
			if(write_scalar_elementwise<Tetrahedron>(File, u, fct, u.template begin<Tetrahedron>(si), u.template end<Tetrahedron>(si), si) == false) return false;
			if(write_scalar_elementwise<Pyramid>(File, u, fct, u.template begin<Pyramid>(si), u.template end<Pyramid>(si), si) == false) return false;
			if(write_scalar_elementwise<Prism>(File, u, fct, u.template begin<Prism>(si), u.template end<Prism>(si), si) == false) return false;
			if(write_scalar_elementwise<Hexahedron>(File, u, fct, u.template begin<Hexahedron>(si), u.template end<Hexahedron>(si), si) == false) return false;
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
count_elem_conn(discrete_function_type& u, int si,
						typename geometry_traits<TElem>::iterator iterBegin,
						typename geometry_traits<TElem>::iterator iterEnd)
{
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;
	Numbers.noElements += u.template num<TElem>(si);
	Numbers.noConnections += u.template num<TElem>(si) * ref_elem_type::num_corners;

	m_grid->begin_marking();
	for( ; iterBegin != iterEnd; ++iterBegin)
	{
		TElem *elem = *iterBegin;
		for(uint i = 0; i < (uint) ref_elem_type::num_corners; ++i)
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
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;
	typedef typename discrete_function_type::domain_type domain_type;
	typename domain_type::position_accessor_type& aaPos = u.get_approximation_space().get_domain().get_position_accessor();

	// write points and remember numbering
	assert(m_aaDOFIndexVRT.valid());
	typename domain_type::position_type Pos;
	float co;
	m_grid->begin_marking();
	for( ; iterBegin != iterEnd; ++iterBegin)
	{
		TElem *elem = *iterBegin;
		for(uint i = 0; i < (uint) ref_elem_type::num_corners; ++i)
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
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;
	static const ReferenceObjectID refID = ref_elem_type::REFERENCE_OBJECT_ID;
	assert(m_aaDOFIndexVRT.valid());

	UG_DLOG(LIB_DISC_OUTPUT, 3, "\n ---- Start: Writing Connectivity (i.e. Vertices of each cell) to file ----");

	for( ; iterBegin != iterEnd; iterBegin++)
	{
		TElem *t = *iterBegin;
		UG_DLOG(LIB_DISC_OUTPUT, 3, "\nWriting Vertices of Finite Element: ");

		int id = 0;
		if(refID != ROID_PRISM)
		{
			for(uint i=0; i< (uint) ref_elem_type::num_corners; i++)
			{
				VertexBase* vert = t->vertex(i);
				id = m_aaDOFIndexVRT[vert];
				UG_DLOG(LIB_DISC_OUTPUT, 3, id << " ");
				BStreamWrite(File, &id);
			}
		}
		else
		{
			id = m_aaDOFIndexVRT[t->vertex(0)]; BStreamWrite(File, &id);
			id = m_aaDOFIndexVRT[t->vertex(2)]; BStreamWrite(File, &id);
			id = m_aaDOFIndexVRT[t->vertex(1)]; BStreamWrite(File, &id);
			id = m_aaDOFIndexVRT[t->vertex(3)]; BStreamWrite(File, &id);
			id = m_aaDOFIndexVRT[t->vertex(5)]; BStreamWrite(File, &id);
			id = m_aaDOFIndexVRT[t->vertex(4)]; BStreamWrite(File, &id);
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
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;
	UG_DLOG(LIB_DISC_OUTPUT, 2, "\n ---- Start: Writing element offsets to file ----\n");

	for( ; iterBegin != iterEnd; ++iterBegin)
	{
		n += ref_elem_type::num_corners;
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
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;
	static const ReferenceObjectID refID = ref_elem_type::REFERENCE_OBJECT_ID;

	char type;

	// TODO: This is 1D and 2D only
	switch(refID)
	{
	case ROID_EDGE: type = (char) 3; break;
	case ROID_TRIANGLE: type = (char) 5; break;
	case ROID_QUADRILATERAL: type = (char) 9; break;
	case ROID_TETRAHEDRON: type = (char) 10; break;
	case ROID_PYRAMID: type = (char) 14; break;
	case ROID_PRISM: type = (char) 13; break;
	case ROID_HEXAHEDRON: type = (char) 12; break;
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
						typename geometry_traits<TElem>::iterator iterEnd,
						int si)
{
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;
	float valf;
	BStream.size = sizeof(float);

	UG_DLOG(LIB_DISC_OUTPUT, 3, "\n ---- Start: Writing nodal values to file for function " << u.name(fct) << " ----\n");

	typename TDiscreteFunction::multi_index_vector_type multInd;
	typename TDiscreteFunction::vector_type& u_vec = u.get_vector();

	m_grid->begin_marking();
	for( ; iterBegin != iterEnd; ++iterBegin)
	{
		TElem *elem = *iterBegin;
		for(uint i = 0; i < (uint) ref_elem_type::num_corners; ++i)
		{
			VertexBase* v = elem->vertex(i);
			if(m_grid->is_marked(v)) continue;

			m_grid->mark(v);

			if(u.get_inner_multi_indices(v, fct, multInd) != 1)
				return false;

			const size_t index = multInd[0][0];
			const size_t alpha = multInd[0][1];

			valf = (float) BlockRef(u_vec[index], alpha);

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
write_prolog(FILE* File, double time)
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
	fprintf(File, "  <time timestep=\"%g\"/>\n", time);
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


template <typename TDiscreteFunction>
bool
VTKOutput<TDiscreteFunction>::
is_valid_filename(const char *nameIn)
{
	// filename: name_p0000_s0000_t0000.vtu
	// i.e. name given by user plus 22 char
	if(strlen(nameIn) + 22 > NAMESIZE)
	{
		UG_LOG("Filename to long. Cannot print to file.\n");
		return false;
	}

	// no dots allowed in file name
	const char* p = strrchr(nameIn, '.');
	if (p != NULL)
	{
		UG_LOG("Filename must not contain '.'. Cannot print to file.\n");
		return false;
	}

	return true;
}

template <typename TDiscreteFunction>
bool
VTKOutput<TDiscreteFunction>::
vtu_filename(char *nameOut, const char *nameIn, int rank, int si, size_t step)
{
	char ext[NAMESIZE];

	if(!is_valid_filename(nameIn)) return false;

	strcpy(nameOut, nameIn);

	// process index
#ifdef UG_PARALLEL
	sprintf(ext, "_p%04d", rank);
	strcat(nameOut, ext);
#endif

	// subset index
	sprintf(ext, "_s%04d", si);
	strcat(nameOut, ext);

	// time index
	sprintf(ext, "_t%04d", (int)step);
	strcat(nameOut, ext);

	// add file extension
	strcat(nameOut, ".vtu");
	return true;

}

template <typename TDiscreteFunction>
bool
VTKOutput<TDiscreteFunction>::
pvtu_filename(char *nameOut, const char *nameIn, int si, size_t step)
{
	char ext[NAMESIZE];

	if(!is_valid_filename(nameIn)) return false;

	strcpy(nameOut, nameIn);

	// subset index
	sprintf(ext, "_s%04d", si);
	strcat(nameOut, ext);

	// time index
	sprintf(ext, "_t%04d", (int)step);
	strcat(nameOut, ext);

	// add file extension
	strcat(nameOut, ".pvtu");
	return true;

}

template <typename TDiscreteFunction>
bool
VTKOutput<TDiscreteFunction>::
pvd_filename(char *nameOut, const char *nameIn)
{
	if(!is_valid_filename(nameIn)) return false;

	strcpy(nameOut, nameIn);

	// add file extension
	strcat(nameOut, ".pvd");
	return true;

}

template <typename TDiscreteFunction>
bool
VTKOutput<TDiscreteFunction>::
pvd_time_filename(char *nameOut, const char *nameIn, size_t timestep)
{
	char ext[NAMESIZE];
	if(!is_valid_filename(nameIn)) return false;

	strcpy(nameOut, nameIn);

	// timestep index
	sprintf(ext, "_t%04d", (int) timestep);
	strcat(nameOut, ext);

	// add file extension
	strcat(nameOut, ".pvd");
	return true;

}

#ifdef UG_PARALLEL
template <typename TDiscreteFunction>
bool
VTKOutput<TDiscreteFunction>::
write_pvtu(discrete_function_type& u, const char* filename, int si, size_t step, number time)
{
	FILE* file;
	char sname[NAMESIZE], pname[NAMESIZE];

	if (pcl::IsOutputProc()) {
		if(pvtu_filename(sname, filename, si, step) != true) return false;

		file = fopen(sname, "w");
		if (file == NULL)
		{
			UG_LOG("Cannot print to file.\n");
			return false;
		}

		// Write to file
		fprintf(file, "<?xml version=\"1.0\"?>\n");
		fprintf(file, "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\">\n");
		fprintf(file, "  <time timestep=\"%g\"/>\n", time);
		fprintf(file, "  <PUnstructuredGrid GhostLevel=\"0\">\n");
		fprintf(file, "    <PPoints>\n");
		fprintf(file, "      <PDataArray type=\"Float32\" NumberOfComponents=\"3\"/>\n");
		fprintf(file, "    </PPoints>\n");

		// Node Data
		fprintf(file, "    <PPointData>\n");
		for(uint fct = 0; fct < u.num_fct(); ++fct)
		{
			if(u.is_def_in_subset(fct, si) == false) continue;

			fprintf(file, "      <PDataArray type=\"Float32\" Name=\"%s\" "
					"NumberOfComponents=\"%d\"/>\n",
					u.name(fct).c_str(), 1);
		}
		fprintf(file, "    </PPointData>\n");

		// Element Data (currently not supported)

		// include files from all procs
		for (int i = 0; i < pcl::GetNumProcesses(); i++) {
			vtu_filename(pname, filename, i, si, step);
			fprintf(file, "    <Piece Source=\"%s\"/>\n", pname);
		}
		fprintf(file, "  </PUnstructuredGrid>\n");
		fprintf(file, "</VTKFile>\n");
		fclose(file);
	}

	return true;
}

template <typename TDiscreteFunction>
bool
VTKOutput<TDiscreteFunction>::
write_time_pvd(discrete_function_type& u, const char* filename)
{
	FILE* file;
	char sname[NAMESIZE], pname[NAMESIZE], procname[NAMESIZE];

	if (pcl::IsOutputProc()) {
		if(pvd_filename(sname, filename) != true) return false;

		file = fopen(sname, "w");
		if (file == NULL)
		{
			UG_LOG("Cannot print to file.\n");
			return false;
		}

		// Write to file
		fprintf(file, "<?xml version=\"1.0\"?>\n");
		fprintf(file, "<VTKFile type=\"Collection\" version=\"0.1\">\n");
		fprintf(file, "  <Collection>\n");

		// include files from all procs
		for(size_t t = 0; t < m_timestep.size(); ++t)
		{
			for(int si = 0; si < u.num_subsets(); ++si)
			{
				pvtu_filename(pname, filename, si, t);
				fprintf(file, "  <DataSet timestep=\"%g\" part=\"%d\" file=\"%s\"/>\n", m_timestep[t], si, pname);
			}
		}

		fprintf(file, "  </Collection>\n");
		fprintf(file, "</VTKFile>\n");
		fclose(file);
	}

	if (pcl::IsOutputProc()) {
		strcpy(sname, filename);
		strcat(sname, "_processwise");

		if(pvd_filename(procname, sname) != true) return false;

		file = fopen(procname, "w");
		if (file == NULL)
		{
			UG_LOG("Cannot print to file.\n");
			return false;
		}

		// Write to file
		fprintf(file, "<?xml version=\"1.0\"?>\n");
		fprintf(file, "<VTKFile type=\"Collection\" version=\"0.1\">\n");
		fprintf(file, "  <Collection>\n");

		// include files from all procs
		for(size_t t = 0; t < m_timestep.size(); ++t)
		{
			for (int i = 0; i < pcl::GetNumProcesses(); i++)
			{
				for(int si = 0; si < u.num_subsets(); ++si)
				{
					vtu_filename(pname, filename, i, si, t);
					fprintf(file, "  <DataSet timestep=\"%g\" part=\"%d\" file=\"%s\"/>\n", m_timestep[t], i, pname);
				}
			}
		}

		fprintf(file, "  </Collection>\n");
		fprintf(file, "</VTKFile>\n");
		fclose(file);
	}

	return true;
}

#endif


template <typename TDiscreteFunction>
bool
VTKOutput<TDiscreteFunction>::
write_pvd(discrete_function_type& u, const char* filename, size_t timestep, number time)
{
	FILE* file;
	char sname[NAMESIZE], pname[NAMESIZE], procname[NAMESIZE];

#ifdef UG_PARALLEL
	if (pcl::IsOutputProc()) {
#endif
		if(pvd_time_filename(sname, filename, timestep) != true) return false;

		file = fopen(sname, "w");
		if (file == NULL)
		{
			UG_LOG("Cannot print to file.\n");
			return false;
		}

		// Write to file
		fprintf(file, "<?xml version=\"1.0\"?>\n");
		fprintf(file, "<VTKFile type=\"Collection\" version=\"0.1\">\n");
		fprintf(file, "  <Collection>\n");

		// include files from all procs
		for(int si = 0; si < u.num_subsets(); ++si)
		{
			pvtu_filename(pname, filename, si, timestep);
			fprintf(file, "  <DataSet timestep=\"%g\" part=\"%d\" file=\"%s\"/>\n", time, si, pname);
		}

		fprintf(file, "  </Collection>\n");
		fprintf(file, "</VTKFile>\n");
		fclose(file);
#ifdef UG_PARALLEL
	}
#endif

#ifdef UG_PARALLEL
	if (pcl::IsOutputProc()) {
		strcpy(sname, filename);
		strcat(sname, "_processwise");

		if(!pvd_time_filename(procname, sname, timestep)) return false;

		file = fopen(procname, "w");
		if (file == NULL)
		{
			UG_LOG("Cannot print to file.\n");
			return false;
		}

		// Write to file
		fprintf(file, "<?xml version=\"1.0\"?>\n");
		fprintf(file, "<VTKFile type=\"Collection\" version=\"0.1\">\n");
		fprintf(file, "  <Collection>\n");

		// include files from all procs
		for (int i = 0; i < pcl::GetNumProcesses(); i++)
		{
			for(int si = 0; si < u.num_subsets(); ++si)
			{
				vtu_filename(pname, filename, i, si, timestep);
				fprintf(file, "  <DataSet timestep=\"%g\" part=\"%d\" file=\"%s\"/>\n", time, i, pname);
			}
		}

		fprintf(file, "  </Collection>\n");
		fprintf(file, "</VTKFile>\n");
		fclose(file);
	}
#endif

	return true;
}

}


