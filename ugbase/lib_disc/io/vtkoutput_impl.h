/*
 * vtkoutput.cpp
 *
 *  Created on: 06.07.2009
 *      Author: andreasvogel
 */


#define __SWAPBYTES__            /* if using LittleEndian */

//other libraries
#include <cstdio>
#include <iostream>
#include <cstring>
#include <string>
#include <algorithm>

// ug4 libraries
#include "common/log.h"
#include "common/util/string_util.h"
#include "lib_disc/common/function_group.h"
#include "lib_disc/common/groups_util.h"

// own interface
#include "vtkoutput.h"

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

inline void BStreamWrite(FILE* File, void *item)
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

inline void BStreamFlush(FILE* File)
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
print(const char* filename, function_type& u, int step, number time, bool makeConsistent)
{
#ifdef UG_PARALLEL
	if(makeConsistent)
	{
		if(!u.change_storage_type(PST_CONSISTENT))
		{
			UG_LOG("ERROR in 'VTK::print': "
				   "Cannot change storage type to consistent.\n");
			return false;
		}
	}
#endif

//	check functions
	bool bEverywhere = true;
	for(size_t fct = 0; fct < u.num_fct(); ++fct)
	{
	//	check if function is defined everywhere
		if(!u.is_def_everywhere(fct))
			bEverywhere = false;
	}

//	in case that all functions are defined everywhere, we write the grid as
//	a whole. If not, we must write each subset separately and group the files
//	later using a *.pvd file.
	if(bEverywhere)
	{
	//	si == -1 indicates whole grid
		int si = -1;

	//	write whole grid to a single file
		if(!print_subset(filename, u, si, step, time, makeConsistent))
		{
			UG_LOG("ERROR in 'VTK::print': Can not write grid.\n");
			return false;
		}
	}
	else
	{
		// 	loop subsets
		for(int si = 0; si < u.num_subsets(); ++si)
		{
		//	write each subset to a single file
			if(!print_subset(filename, u, si, step, time, makeConsistent))
			{
				UG_LOG("ERROR in 'VTK::print': Can not write Subset "<< si << ".\n");
				return false;
			}
		}

		//	write grouping pvd file
		if(!write_subset_pvd(u, filename, step, time))
		{
			UG_LOG("ERROR in 'VTK::print': Can not write pvd file.\n");
			return false;
		}
	}

//	we're done
	return true;
}

template <typename TDiscreteFunction>
bool
VTKOutput<TDiscreteFunction>::
print_subset(const char* filename, function_type& u, int si, int step, number time, bool makeConsistent)
{
#ifdef UG_PARALLEL
	if(makeConsistent)
	{
		if(!u.change_storage_type(PST_CONSISTENT))
		{
			UG_LOG("ERROR in 'VTK::print_subset': "
				   "Cannot change storage type to consistent.\n");
			return false;
		}
	}
#endif

//	get the grid associated to the solution
	m_pGrid = dynamic_cast<Grid*>(&u.get_approximation_space().get_domain().get_grid());

//	check grid
	if(!m_pGrid)
	{
		UG_LOG("ERROR in 'VTK::print_subset': Cannot find underlying Grid.\n");
		return false;
	}

// 	attach help indices
	m_pGrid->attach_to_vertices(m_aDOFIndex);
	m_aaDOFIndexVRT.access(*m_pGrid, m_aDOFIndex);

//	get rank of process
	int rank = 0;
#ifdef UG_PARALLEL
	rank = pcl::GetProcRank();
#endif

//	get name for *.vtu file
	std::string name;
	if(!vtu_filename(name, filename, rank, si, u.num_subsets()-1, step)) return false;

//	open the file
	FILE* File = fopen(name.c_str(), "w");
	if(File == NULL)
	{
		UG_LOG("ERROR in 'VTK::print_subset': Can not open Output File.\n");
		return false;
	}

//	reset stream
	BStream.front = 0;

//	bool if time point should be written to *.vtu file
//	in parallel we must not (!) write it to the *.vtu file, but to the *.pvtu
	bool bTimeDep = (step >= 0);
#ifdef UG_PARALLEL
	if(pcl::GetNumProcesses() > 1) bTimeDep = false;
#endif

//	header
	fprintf(File, "<?xml version=\"1.0\"?>\n");
	fprintf(File, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" "
			"byte_order=\"%s\">\n",
#ifdef __SWAPBYTES__
			"LittleEndian"
#else
			"BigEndian"
#endif
	);

//	writing time point
	if(bTimeDep)
		fprintf(File, "  <Time timestep=\"%g\"/>\n", time);

//	opening the grid
	fprintf(File, "  <UnstructuredGrid>\n");

// 	get dimension of grid-piece
	int dim = -1;
	if(si >= 0) dim = DimensionOfSubset(u.get_domain(), si);
	else dim = DimensionOfDomain(u.get_domain());

//	write piece of grid
	if(dim >= 0){
		if(!write_piece(File, u, si, dim))
		{
			UG_LOG("ERROR in 'VTK::print_subset': Can not write Subset.\n");
			fclose(File);
			return false;
		}
	}
	else
	{
	//	if dim < 0, some is wrong with grid, except no element is in the grid
		if( ((si < 0) && m_pGrid->num<VertexBase>() != 0) ||
			((si >=0) && u.get_domain().get_subset_handler().template num<VertexBase>(si) != 0))
		{
			UG_LOG("ERROR in 'VTK::print_subset': Dimension of grid/subset not"
					" detected correctly although grid objects present.\n");
			fclose(File);
			return false;
		}

	//	write that no elements are in the grid
		fprintf(File, "    <Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n",
														0, 0);
		if(!write_points(File, u, si, dim, 0) ||
			!write_cells(File, u, si, dim, 0, 0))
		{
			UG_LOG("ERROR in 'VTK::print_subset': Can not write empty Points/Cells.\n");
			return false;
		}
		fprintf(File, "    </Piece>\n");

	}

//	write closing xml tags
	fprintf(File, "  </UnstructuredGrid>\n");
	fprintf(File, "</VTKFile>\n");

//	close stream
	fclose(File);

// 	detach help indices
	m_pGrid->detach_from_vertices(m_aDOFIndex);
	m_aaDOFIndexVRT.invalidate();

#ifdef UG_PARALLEL
//	write grouping *.pvtu file in parallel case
	if(!write_pvtu(u, filename, si, step, time))
	{
		UG_LOG("ERROR in 'VTK::print_subset': Can not write pvtu - file.\n");
		fclose(File);
		return false;
	}
#endif

//	remember time step
	if(step >= 0)
	{
	//	get vector of time points for the name
		std::vector<number>& vTimestep = m_mTimestep[filename];

	//	resize the vector
		vTimestep.resize(step+1);

	//	write time point
		vTimestep[step] = time;
	}

//	we're done
	return true;
}


template <typename TDiscreteFunction>
bool
VTKOutput<TDiscreteFunction>::
write_piece(FILE* File, function_type& u, int si, int dim)
{
//	counters
	int numVert = 0, numElem = 0, numConn = 0;

// 	Count needed sizes for vertices, elements and connections
	if(!count_piece_sizes(u, si, dim, numVert, numElem, numConn))
	{
		UG_LOG("ERROR in 'VTK::write_subset': Can not count piece sizes.\n");
		return false;
	}

//	write the beginning of the piece, indicating the number of vertices
//	and the number of elements for this piece of the grid.
	fprintf(File, "    <Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n",
													numVert, numElem);

//	write vertices of this piece
	if(!write_points(File, u, si, dim, numVert))
	{
		UG_LOG("ERROR in 'VTK::write_subset': Can not write Points.\n");
		return false;
	}

//	write elements of this piece
	if(!write_cells(File, u, si, dim, numElem, numConn))
	{
		UG_LOG("ERROR in 'VTK::write_subset': Can not write Elements.\n");
		return false;
	}

//	write opening tag to indicate point data
	fprintf(File, "      <PointData>\n");

//	get function pattern
	const FunctionPattern& fctPatt = u.get_dof_distribution().get_function_pattern();

//	add all components if 'selectAll' chosen
	if(m_bSelectAll)
		for(size_t fct = 0; fct < u.num_fct(); ++fct)
			select_nodal_scalar(u.name(fct).c_str(), u.name(fct).c_str());

//	loop all selected symbolic names
	for(size_t sym = 0; sym < m_vSymbFct.size(); ++sym)
	{
	//	create functio group
		FunctionGroup fctGrp;

	//	get symb function
		const std::string& symbNames = m_vSymbFct[sym].first;
		const std::string& vtkName = m_vSymbFct[sym].second;

	//	extract names
		if(!ConvertStringToFunctionGroup(fctGrp, fctPatt, symbNames.c_str()))
			return false;

	//	check that all functions are contained in subset
		bool bContained = true;
		for(size_t i = 0; i < fctGrp.num_fct(); ++i)
		{
		//	get function
			const size_t fct = fctGrp[i];

		//	check
			if(!u.is_def_in_subset(fct, si)) bContained = false;
		}
		if(!bContained) continue;

	//	write scalar value of this function
		if(!write_nodal_values(File, u, fctGrp, vtkName, si, dim, numVert))
		{
			UG_LOG("ERROR in 'VTK::write_subset': Can not write Scalar Values.\n");
			return false;
		}
	}

//	write closing tag
	fprintf(File, "      </PointData>\n");

//	write closing tag
	fprintf(File, "    </Piece>\n");

//	we're done
	return true;
}

template <typename TDiscreteFunction>
bool
VTKOutput<TDiscreteFunction>::
count_piece_sizes(function_type& u, int si, int dim,
            int& numVert, int& numElem, int& numConn)
{
//	debug output
	UG_DLOG(LIB_DISC_OUTPUT, 2, "\n ---- Init Numbers ----\n");

//	switch dimension
	switch(dim)
	{
		case 0: if(si>=0) numVert = u.template num<VertexBase>(si);
				else {
					numVert = 0;
					for(si = 0; si < u.num_subsets(); ++si)
						numVert += u.template num<VertexBase>(si);
				}
				break;
		case 1: count_sizes<Edge>(u, si, numVert, numElem, numConn);
				break;
		case 2: count_sizes<Triangle>(u, si, numVert, numElem, numConn);
				count_sizes<Quadrilateral>(u, si, numVert, numElem, numConn);
				break;
		case 3: count_sizes<Tetrahedron>(u, si, numVert, numElem, numConn);
				count_sizes<Pyramid>(u, si, numVert, numElem, numConn);
				count_sizes<Prism>(u, si, numVert, numElem, numConn);
				count_sizes<Hexahedron>(u, si, numVert, numElem, numConn);
				break;
		default: UG_LOG("ERROR in 'VTK::init_subset': Dimension " << dim <<
		                " is not supported.\n");
				return false;
	}

//	debug output
	UG_DLOG(LIB_DISC_OUTPUT, 2, "Number of Vertices: " << numVert << "\n");
	UG_DLOG(LIB_DISC_OUTPUT, 2, "Number of Elements: " << numElem << "\n");
	UG_DLOG(LIB_DISC_OUTPUT, 2, "Number of Connections: " << numConn << "\n");
	UG_DLOG(LIB_DISC_OUTPUT, 2, " ---- End ----\n");

//	we're done
	return true;
}

template <typename TDiscreteFunction>
bool
VTKOutput<TDiscreteFunction>::
write_points(FILE* File, function_type& u, int si, int dim,
             int numVert)
{
//	write starting xml tag for points
	fprintf(File, "      <Points>\n");
	fprintf(File, "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"binary\">\n");
	BStream.size = sizeof(int);
	int n = 3*sizeof(float) * numVert;
	BStreamWrite(File, &n);
	BStreamFlush(File);
	BStream.size = sizeof(float);

//	reset counter for vertices
	n = 0;

//	switch dimension
	if(numVert > 0){
		switch(dim)
		{
			case 0: write_points_elementwise<VertexBase>(File, u, si, n);
					break;
			case 1: write_points_elementwise<Edge>(File, u, si, n);
					break;
			case 2: write_points_elementwise<Triangle>(File, u, si, n);
					write_points_elementwise<Quadrilateral>(File, u, si, n);
					break;
			case 3:	write_points_elementwise<Tetrahedron>(File, u, si, n);
					write_points_elementwise<Pyramid>(File, u, si, n);
					write_points_elementwise<Prism>(File, u, si, n);
					write_points_elementwise<Hexahedron>(File, u, si, n);
					break;
			default: UG_LOG("ERROR in 'VTK::write_points': Dimension " << dim <<
							" is not supported.\n");
					return false;
		}
	}

//	write closing tags
	fprintf(File, "\n        </DataArray>\n");
	fprintf(File, "      </Points>\n");

//	everything fine
	return true;
}

template <typename TDiscreteFunction>
bool VTKOutput<TDiscreteFunction>::
write_cells(FILE* File, function_type& u, int si, int dim,
               int numElem, int numConn)
{
//	write opening tag to indicate that elements will be written
	fprintf(File, "      <Cells>\n");

	///////////////////////////
	// connectivity
	///////////////////////////
//	write opening tag to indicate that connections will be written
	fprintf(File, "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"binary\">\n");
	BStream.size = sizeof(int);
	int n = sizeof(int) * numConn;
	BStreamWrite(File, &n);
	BStreamFlush(File);

//	switch dimension
	if(numConn > 0){
		switch(dim)
		{
			case 0: break; // no elements -> nothing to do
			case 1: write_cell_connectivity<Edge>(File, u, si);
					break;
			case 2: write_cell_connectivity<Triangle>(File, u, si);
					write_cell_connectivity<Quadrilateral>(File, u, si);
					break;
			case 3: write_cell_connectivity<Tetrahedron>(File, u, si);
					write_cell_connectivity<Pyramid>(File, u, si);
					write_cell_connectivity<Prism>(File, u, si);
					write_cell_connectivity<Hexahedron>(File, u, si);
					break;
			default: UG_LOG("ERROR in 'VTK::write_elements': Dimension " << dim <<
							" is not supported.\n");
					return false;
		}
	}
	BStreamFlush(File);

//	write closing tag
	fprintf(File, "\n        </DataArray>\n");

	///////////////////////////
	// offsets
	///////////////////////////
//	write opening tag indicating that offsets are going to be written
	fprintf(File, "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"binary\">\n");
	n = sizeof(int) * numElem;
	BStreamWrite(File, &n);
	BStreamFlush(File);

	n = 0;
//	switch dimension
	if(numElem > 0){
		switch(dim)
		{
			case 0: break; // no elements -> nothing to do
			case 1: write_cell_offsets<Edge>(File, u, si, n);
					break;
			case 2: write_cell_offsets<Triangle>(File, u, si, n);
					write_cell_offsets<Quadrilateral>(File, u, si, n);
					break;
			case 3: write_cell_offsets<Tetrahedron>(File, u, si, n);
					write_cell_offsets<Pyramid>(File, u, si, n);
					write_cell_offsets<Prism>(File, u, si, n);
					write_cell_offsets<Hexahedron>(File, u, si, n);
					break;
			default: UG_LOG("ERROR in 'VTK::write_elements': Dimension " << dim <<
							" is not supported.\n");
					return false;
		}
	}
	BStreamFlush(File);
	fprintf(File, "\n        </DataArray>\n");


	///////////////////////////
	// types of elements
	///////////////////////////
//	write opening tag to indicate that types will be written
	fprintf(File, "        <DataArray type=\"Int8\" Name=\"types\" format=\"binary\">\n");
	BStreamWrite(File, &numElem);
	BStreamFlush(File);

//	switch dimension
	if(numElem > 0) {
		switch(dim)
		{
			case 0: break; // no elements -> nothing to do
			case 1: write_cell_types<Edge>(File, u, si);
					break;
			case 2: write_cell_types<Triangle>(File, u, si);
					write_cell_types<Quadrilateral>(File, u, si);
					break;
			case 3: write_cell_types<Tetrahedron>(File, u, si);
					write_cell_types<Pyramid>(File, u, si);
					write_cell_types<Prism>(File, u, si);
					write_cell_types<Hexahedron>(File, u, si);
					break;
			default: UG_LOG("ERROR in 'VTK::write_elements': Dimension " << dim <<
							" is not supported.\n");
					return false;
		}
	}

//	write closing tag
	fprintf(File, "\n        </DataArray>\n");
	fprintf(File, "      </Cells>\n");

//	we're done
	return true;
}


template <typename TDiscreteFunction>
bool
VTKOutput<TDiscreteFunction>::
write_nodal_values(FILE* File, function_type& u,
                   const FunctionGroup& vFct, const std::string& name,
                   int si, int dim, int numVert)
{
//	write opening tag
	fprintf(File, "        <DataArray type=\"Float32\" Name=\"%s\" "
					"NumberOfComponents=\"%d\" format=\"binary\">\n",
					name.c_str(), (vFct.num_fct() == 1 ? 1 : 3));
	BStream.size = sizeof(int);
	int n = sizeof(float) * numVert * (vFct.num_fct() == 1 ? 1 : 3);
	BStreamWrite(File, &n);
	BStreamFlush(File);

//	switch dimension
	switch(dim)
	{
		case 0:
			if(!write_nodal_values_elementwise<VertexBase>(File, u, vFct, si)) return false;
			break;
		case 1:
			if(!write_nodal_values_elementwise<Edge>(File, u, vFct, si)) return false;
			break;
		case 2:
			if(!write_nodal_values_elementwise<Triangle>(File, u, vFct, si)) return false;
			if(!write_nodal_values_elementwise<Quadrilateral>(File, u, vFct, si)) return false;
			break;
		case 3:
			if(!write_nodal_values_elementwise<Tetrahedron>(File, u, vFct, si)) return false;
			if(!write_nodal_values_elementwise<Pyramid>(File, u, vFct, si)) return false;
			if(!write_nodal_values_elementwise<Prism>(File, u, vFct, si)) return false;
			if(!write_nodal_values_elementwise<Hexahedron>(File, u, vFct, si)) return false;
			break;
		default: UG_LOG("ERROR in 'VTK::write_scalar': Dimension " << dim <<
		                " is not supported.\n");
				return false;
	}

//	flush stream
	BStreamFlush(File);

//	write closing tag
	fprintf(File, "\n        </DataArray>\n");

//	everything fine
	return true;
};


template <typename TDiscreteFunction>
template <typename TElem>
void
VTKOutput<TDiscreteFunction>::
count_sizes(function_type& u, int si,
                int& numVert, int& numElem, int& numConn)
{
//	get reference element
	typedef typename reference_element_traits<TElem>::reference_element_type
																ref_elem_type;

//	iterator for the elements
	typename geometry_traits<TElem>::const_iterator iterBegin, iterEnd;

//	reset all marks
	m_pGrid->begin_marking();

//	loop all subsets for whole domain, only the given subset if si >= 0
	int sistart = si, siend = si+1;
	if(si < 0) {sistart = 0; siend = u.num_subsets();}
	for(int si = sistart; si < siend; ++si)
	{
	//	get iterators
		iterBegin = u.template begin<TElem>(si);
		iterEnd = u.template end<TElem>(si);

	//	count number of elements and number of connections
		numElem += u.template num<TElem>(si);
		numConn += u.template num<TElem>(si) * ref_elem_type::num_corners;

	//	loop over elements of this subset
		for( ; iterBegin != iterEnd; ++iterBegin)
		{
		//	get the element
			TElem *elem = *iterBegin;

		//	loop vertices of the element
			for(size_t i = 0; i < (size_t) ref_elem_type::num_corners; ++i)
			{
			//	get vertex of the element
				VertexBase* v = elem->vertex(i);

			//	if this vertex has already been counted, skip it
				if(m_pGrid->is_marked(v)) continue;

			// count vertex and mark it
				++numVert;
				m_pGrid->mark(v);
			}
		}
	}

//	signal end of marking
	m_pGrid->end_marking();
};


template <typename TDiscreteFunction>
template <typename TElem>
void
VTKOutput<TDiscreteFunction>::
write_points_elementwise(FILE* File, function_type& u, int si, int& n)
{
//	get reference element
	typedef typename reference_element_traits<TElem>::reference_element_type
																ref_elem_type;

//	get domain type
	typedef typename function_type::domain_type domain_type;

//	get position attachment
	typename domain_type::position_accessor_type& aaPos =
				u.get_approximation_space().get_domain().get_position_accessor();

// 	write points and remember numbering
	UG_ASSERT(m_aaDOFIndexVRT.valid(), "Missing attachment");

//	position vector
	typename domain_type::position_type Pos;

//	corner counter
	float co;

//	get iterators
	typename geometry_traits<TElem>::const_iterator iterBegin, iterEnd;

//	start marking of vertices
	m_pGrid->begin_marking();

//	loop all subsets for whole domain, only the given subset if si >= 0
	int sistart = si, siend = si+1;
	if(si < 0) {sistart = 0; siend = u.num_subsets();}
	for(int si = sistart; si < siend; ++si)
	{
		iterBegin = u.template begin<TElem>(si);
		iterEnd = u.template end<TElem>(si);

	//	loop all elements of the subset
		for( ; iterBegin != iterEnd; ++iterBegin)
		{
		//	get the element
			TElem *elem = *iterBegin;

		//	loop vertices of the element
			for(size_t i = 0; i < (size_t) ref_elem_type::num_corners; ++i)
			{
			//	get vertex of element
				VertexBase* v = GetVertex(elem, i);

			//	if vertex has already be handled, skip it
				if(m_pGrid->is_marked(v)) continue;

			//	mark the vertex as processed
				m_pGrid->mark(v);

			//	number vertex
				m_aaDOFIndexVRT[v] = n++;

			//	get position of vertex
				Pos = aaPos[v];

			//	write position to stream
				for(int i = 0; i < domain_type::dim; ++i)
				{
					co = Pos[i];
					BStreamWrite(File, &co);
				}

			//	fill with missing zeros (if dim < 3)
				for(int i = domain_type::dim; i < 3; ++i)
				{
					co = 0.0;
					BStreamWrite(File, &co);
				}
			}
		}
	}

	if(n > 0){
		UG_DLOG(LIB_DISC_OUTPUT, 3, " ---- " << n <<
		        " Vertices (Nr. 0 - " << (n-1) << ") written to file ----\n");}
	else{
		UG_DLOG(LIB_DISC_OUTPUT, 3, " ---- " << n <<
		        " Vertices written to file ----\n");}

//	flush the stream
	BStreamFlush(File);

//	signal end of marking the grid
	m_pGrid->end_marking();
}


template <typename TDiscreteFunction>
template <class TElem>
void
VTKOutput<TDiscreteFunction>::
write_cell_connectivity(FILE* File, function_type& u, int si)
{
//	get reference element type
	typedef typename reference_element_traits<TElem>::reference_element_type
																ref_elem_type;

//	get reference element id
	static const ReferenceObjectID refID = ref_elem_type::REFERENCE_OBJECT_ID;

//	check id accessor
	UG_ASSERT(m_aaDOFIndexVRT.valid(), "ID access invalid");

//	get iterators
	typename geometry_traits<TElem>::const_iterator iterBegin, iterEnd;

//	loop all subsets for whole domain, only the given subset if si >= 0
	int sistart = si, siend = si+1;
	if(si < 0) {sistart = 0; siend = u.num_subsets();}
	for(int si = sistart; si < siend; ++si)
	{
	//	get iterators
		iterBegin = u.template begin<TElem>(si);
		iterEnd = u.template end<TElem>(si);

	//	loop all elements
		for( ; iterBegin != iterEnd; iterBegin++)
		{
		//	get element
			TElem* elem = *iterBegin;

		//	write ids of the element
			if(refID != ROID_PRISM)
			{
				for(size_t i=0; i< (size_t) ref_elem_type::num_corners; i++)
				{
					VertexBase* vert = elem->vertex(i);
					int id = m_aaDOFIndexVRT[vert];
					BStreamWrite(File, &id);
				}
			}
			else
			{
				int id = m_aaDOFIndexVRT[elem->vertex(0)]; BStreamWrite(File, &id);
				id = m_aaDOFIndexVRT[elem->vertex(2)]; BStreamWrite(File, &id);
				id = m_aaDOFIndexVRT[elem->vertex(1)]; BStreamWrite(File, &id);
				id = m_aaDOFIndexVRT[elem->vertex(3)]; BStreamWrite(File, &id);
				id = m_aaDOFIndexVRT[elem->vertex(5)]; BStreamWrite(File, &id);
				id = m_aaDOFIndexVRT[elem->vertex(4)]; BStreamWrite(File, &id);
			}
		}
	}
}

template <typename TDiscreteFunction>
template <class TElem>
void
VTKOutput<TDiscreteFunction>::
write_cell_offsets(	FILE* File, function_type& u, int si, int& n)
{
//	get reference element
	typedef typename reference_element_traits<TElem>::reference_element_type
																ref_elem_type;

//	get iterators
	typename geometry_traits<TElem>::const_iterator iterBegin, iterEnd;

//	loop all subsets for whole domain, only the given subset if si >= 0
	int sistart = si, siend = si+1;
	if(si < 0) {sistart = 0; siend = u.num_subsets();}
	for(int si = sistart; si < siend; ++si)
	{
	//	get iterators
		iterBegin = u.template begin<TElem>(si);
		iterEnd = u.template end<TElem>(si);

	//	loop all elements
		for( ; iterBegin != iterEnd; ++iterBegin)
		{
		//	increase counter of vertices
			n += ref_elem_type::num_corners;

		//	write offset
			BStreamWrite(File, &n);
		}
	}

//	flush stream
	BStreamFlush(File);
}

template <typename TDiscreteFunction>
template <class TElem>
void
VTKOutput<TDiscreteFunction>::
write_cell_types(FILE* File, function_type& u, int si)
{
//	get reference element
	typedef typename reference_element_traits<TElem>::reference_element_type
																ref_elem_type;
//	get object type
	static const ReferenceObjectID refID = ref_elem_type::REFERENCE_OBJECT_ID;

//	type
	char type;

//	get type, based on reference element type
	switch(refID)
	{
		case ROID_EDGE: type = (char) 3; break;
		case ROID_TRIANGLE: type = (char) 5; break;
		case ROID_QUADRILATERAL: type = (char) 9; break;
		case ROID_TETRAHEDRON: type = (char) 10; break;
		case ROID_PYRAMID: type = (char) 14; break;
		case ROID_PRISM: type = (char) 13; break;
		case ROID_HEXAHEDRON: type = (char) 12; break;
		default: throw(UGFatalError("Element Type not known."));
	}

	BStream.size = sizeof(char);

//	get iterators
	typename geometry_traits<TElem>::const_iterator iterBegin, iterEnd;

//	loop all subsets for whole domain, only the given subset if si >= 0
	int sistart = si, siend = si+1;
	if(si < 0) {sistart = 0; siend = u.num_subsets();}
	for(int si = sistart; si < siend; ++si)
	{
	//	get iterators
		iterBegin = u.template begin<TElem>(si);
		iterEnd = u.template end<TElem>(si);

	//	loop all elements, write type for each element to stream
		for( ; iterBegin != iterEnd; ++iterBegin)
			BStreamWrite(File, &type);
	}

	BStreamFlush(File);
}


template <typename TDiscreteFunction>
template <class TElem>
bool
VTKOutput<TDiscreteFunction>::
write_nodal_values_elementwise(FILE* File, function_type& u,
                               const FunctionGroup& vFct, int si)
{
//	get reference element
	typedef typename reference_element_traits<TElem>::reference_element_type
																ref_elem_type;

	float valf;
	BStream.size = sizeof(float);

//	index vector
	typename TDiscreteFunction::multi_index_vector_type vMultInd;

//	start marking of grid
	m_pGrid->begin_marking();

//	get iterators
	typename geometry_traits<TElem>::const_iterator iterBegin, iterEnd;

//	loop all subsets for whole domain, only the given subset if si >= 0
	int sistart = si, siend = si+1;
	if(si < 0) {sistart = 0; siend = u.num_subsets();}
	for(int si = sistart; si < siend; ++si)
	{
	//	get iterators
		iterBegin = u.template begin<TElem>(si);
		iterEnd = u.template end<TElem>(si);

	//	loop all elements
		for( ; iterBegin != iterEnd; ++iterBegin)
		{
		//	get element
			TElem *elem = *iterBegin;

		//	loop vertices of element
			for(size_t co = 0; co < (size_t) ref_elem_type::num_corners; ++co)
			{
			//	get vertex of element
				VertexBase* v = GetVertex(elem, co);

			//	if vertex has been handled before, skip
				if(m_pGrid->is_marked(v)) continue;

			//	mark as used
				m_pGrid->mark(v);

			//	loop all compenents
				for(size_t i = 0; i < vFct.num_fct(); ++i)
				{
				//	get multi index of vertex for the function
					if(u.inner_multi_indices(v, vFct[i], vMultInd) != 1)
					{
						UG_LOG("ERROR in 'VTK:write_nodal_values_elementwise': "
								"The function component "<<vFct[i]<<" has "<<
								vMultInd.size()<<" DoFs in  a vertex. To write a "
								"component to vtk, exactly one DoF must be "
								"given in a vertex.\n");
						return false;
					}

				//	get index and subindex
					const size_t index = vMultInd[0][0];
					const size_t alpha = vMultInd[0][1];

				//	read value from vector
					valf = (float) BlockRef(u[index], alpha);

				//	flush stream
					BStreamWrite(File, &valf);
				}

			//	fill with zeros up to 3d if vector type
				if(vFct.num_fct() != 1)
					for(size_t i = vFct.num_fct(); i < 3; ++i)
					{
						valf = (float) 0.0;
						BStreamWrite(File, &valf);
					}
			}
		}
	}

//	end marking
	m_pGrid->end_marking();

//	everything fine
	return true;
}

template <typename TDiscreteFunction>
bool
VTKOutput<TDiscreteFunction>::
is_valid_filename(std::string& nameIn)
{
// 	search for dots, they are not allowed in file name
	nameIn = nameIn.substr(0, nameIn.find_first_of('.'));

//	everything ok
	return true;
}


template <typename TDiscreteFunction>
bool
VTKOutput<TDiscreteFunction>::
vtu_filename(std::string& nameOut, std::string nameIn, int rank,
             int si, int maxSi, int step)
{
	if(!is_valid_filename(nameIn)) return false;

//	copy name
	nameOut = nameIn;

#ifdef UG_PARALLEL
// 	process index
	if(pcl::GetNumProcesses() > 1)
		AppendCounterToString(nameOut, "_p", rank, pcl::GetNumProcesses() - 1);
#endif

// 	subset index
	if(si >= 0)
		AppendCounterToString(nameOut, "_s", si, maxSi);

// 	time index
	if(step >= 0)
		AppendCounterToString(nameOut, "_t", (int) step);

// 	add file extension
	nameOut.append(".vtu");

//	we're done
	return true;
}

template <typename TDiscreteFunction>
bool
VTKOutput<TDiscreteFunction>::
pvtu_filename(std::string& nameOut, std::string nameIn,
              int si, int maxSi, int step)
{
//	check name
	if(!is_valid_filename(nameIn)) return false;

//	copy name
	nameOut = nameIn;

// 	subset index
	if(si >= 0)
		AppendCounterToString(nameOut, "_s", si, maxSi);

// 	time index
	if(step >= 0)
		AppendCounterToString(nameOut, "_t", (int) step);

// 	add file extension
	nameOut.append(".pvtu");

//	we're done
	return true;
}

template <typename TDiscreteFunction>
bool
VTKOutput<TDiscreteFunction>::
pvd_filename(std::string& nameOut, std::string nameIn)
{
//	check name
	if(!is_valid_filename(nameIn)) return false;

//	copy name
	nameOut = nameIn;

// 	add file extension
	nameOut.append(".pvd");

//	we're done
	return true;
}

template <typename TDiscreteFunction>
bool
VTKOutput<TDiscreteFunction>::
pvd_time_filename(std::string& nameOut, std::string nameIn, int step)
{
//	check name
	if(!is_valid_filename(nameIn)) return false;

//	copy name
	nameOut = nameIn;

// 	time index
	if(step >= 0)
		AppendCounterToString(nameOut, "_t", (int) step);

// 	add file extension
	nameOut.append(".pvd");

//	we're done
	return true;
}

template <typename TDiscreteFunction>
bool
VTKOutput<TDiscreteFunction>::
write_pvtu(function_type& u, const std::string& filename,
           int si, int step, number time)
{
#ifdef UG_PARALLEL
//	File pointer
	FILE* file;

//	file name
	std::string name;

//	get and check number of procs (only for numProcs > 1 we write the pvtu)
	int numProcs = pcl::GetNumProcesses();
	if(numProcs == 1) return true;

//	check if this proc is output proc
	bool isOutputProc = pcl::IsOutputProc();

//	max subset
	int maxSi = u.num_subsets() - 1;

//	only the master process writes this file
	if (isOutputProc)
	{
	//	get name for *.pvtu file
		if(!pvtu_filename(name, filename, si, maxSi, step)) return false;

	//	open file
		file = fopen(name.c_str(), "w");
		if (file == NULL)
		{
			UG_LOG("Cannot print to file.\n");
			return false;
		}

	//	Write to file
		fprintf(file, "<?xml version=\"1.0\"?>\n");
		fprintf(file, "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\">\n");
		fprintf(file, "  <Time timestep=\"%g\"/>\n", time);
		fprintf(file, "  <PUnstructuredGrid GhostLevel=\"0\">\n");
		fprintf(file, "    <PPoints>\n");
		fprintf(file, "      <PDataArray type=\"Float32\" NumberOfComponents=\"3\"/>\n");
		fprintf(file, "    </PPoints>\n");

	// 	Node Data
		fprintf(file, "    <PPointData>\n");
		for(size_t fct = 0; fct < u.num_fct(); ++fct)
		{
		//	skip functions not defined in the subset
			if(!u.is_def_in_subset(fct, si)) continue;

			fprintf(file, "      <PDataArray type=\"Float32\" Name=\"%s\" "
						  "NumberOfComponents=\"%d\"/>\n",
						  u.name(fct).c_str(), 1);
		}
		fprintf(file, "    </PPointData>\n");

		// Element Data (currently not supported)

	// 	include files from all procs
		for (int i = 0; i < numProcs; i++) {
			vtu_filename(name, filename, i, si, maxSi, step);
			name = FilenameWithoutPath(name);
			fprintf(file, "    <Piece Source=\"%s\"/>\n", name.c_str());
		}

	//	write closing tags
		fprintf(file, "  </PUnstructuredGrid>\n");
		fprintf(file, "</VTKFile>\n");
		fclose(file);
	}
#endif

//	we're done
	return true;
}

template <typename TDiscreteFunction>
bool
VTKOutput<TDiscreteFunction>::
write_time_pvd(const char* filename, function_type& u)
{
//	File
	FILE* file;

//	filename
	std::string name;

// 	get some numbers
	bool isOutputProc = true;
	int numProcs = 1;
#ifdef UG_PARALLEL
	isOutputProc = pcl::IsOutputProc();
	numProcs = pcl::GetNumProcesses();
#endif

//	get time steps
	std::vector<number>& vTimestep = m_mTimestep[filename];

	if (isOutputProc)
	{
	//	get file name
		if(!pvd_filename(name, filename)) return false;

	//	open file
		file = fopen(name.c_str(), "w");
		if (file == NULL)
		{
			UG_LOG("Cannot print to file.\n");
			return false;
		}

	// 	Write to file
		fprintf(file, "<?xml version=\"1.0\"?>\n");
		fprintf(file, "<VTKFile type=\"Collection\" version=\"0.1\">\n");
		fprintf(file, "  <Collection>\n");

	//	check functions
		bool bEverywhere = true;
		for(size_t fct = 0; fct < u.num_fct(); ++fct)
		{
		//	check if function is defined everywhere
			if(!u.is_def_everywhere(fct))
				bEverywhere = false;
		}

	// 	include files from all procs
		if(bEverywhere)
		{
			for(int step = 0; step < (int)vTimestep.size(); ++step)
			{
				vtu_filename(name, filename, 0, -1, 0, step);
#ifdef UG_PARALLEL
				if(numProcs > 1) pvtu_filename(name, filename, -1, 0, step);
#endif
				name = FilenameWithoutPath(name);
				fprintf(file, "  <DataSet timestep=\"%g\" part=\"%d\" file=\"%s\"/>\n",
				        		vTimestep[step], 0, name.c_str());
			}
		}
		else
		{
			for(int step = 0; step < (int)vTimestep.size(); ++step)
				for(int si = 0; si < u.num_subsets(); ++si)
				{
					vtu_filename(name, filename, 0, si, u.num_subsets()-1, step);
#ifdef UG_PARALLEL
					if(numProcs > 1) pvtu_filename(name, filename, si, u.num_subsets()-1, step);
#endif
					name = FilenameWithoutPath(name);
					fprintf(file, "  <DataSet timestep=\"%g\" part=\"%d\" file=\"%s\"/>\n",
					        	vTimestep[step], si, name.c_str());
				}
			}

	//	close file
		fprintf(file, "  </Collection>\n");
		fprintf(file, "</VTKFile>\n");
		fclose(file);
	}

	if (isOutputProc && numProcs > 1)
	{
	//	adjust filename
		std::string procName = filename;
		procName.append("_processwise");
		if(!pvd_filename(name, procName)) return false;

	//	open file
		file = fopen(name.c_str(), "w");
		if (file == NULL)
		{
			UG_LOG("Cannot print to file.\n");
			return false;
		}

	// 	Write to file
		fprintf(file, "<?xml version=\"1.0\"?>\n");
		fprintf(file, "<VTKFile type=\"Collection\" version=\"0.1\">\n");
		fprintf(file, "  <Collection>\n");

	// 	include files from all procs
		for(int step = 0; step < (int)vTimestep.size(); ++step)
			for (int rank = 0; rank < numProcs; rank++)
				for(int si = 0; si < u.num_subsets(); ++si)
				{
					vtu_filename(name, filename, rank, si, u.num_subsets()-1, step);
#ifdef UG_PARALLEL
					if(numProcs > 1) pvtu_filename(name, filename, si, u.num_subsets()-1, step);
#endif
					name = FilenameWithoutPath(name);
					fprintf(file, "  <DataSet timestep=\"%g\" part=\"%d\" file=\"%s\"/>\n",
					        vTimestep[step], rank, name.c_str());
				}

	//	end file
		fprintf(file, "  </Collection>\n");
		fprintf(file, "</VTKFile>\n");
		fclose(file);
	}

//	we're done
	return true;
}


template <typename TDiscreteFunction>
bool
VTKOutput<TDiscreteFunction>::
write_subset_pvd(function_type& u, const std::string& filename, int step, number time)
{
//	file pointer
	FILE* file;

//	to store name of file
	std::string name;

//	get rank, outproc bool and number of processes
	bool isOutputProc = true;
	int rank = 0;
	int numProcs = 1;

#ifdef UG_PARALLEL
	isOutputProc = pcl::IsOutputProc();
	rank = pcl::GetProcRank();
	numProcs = pcl::GetNumProcesses();
#endif

//	only output proc writes this file
	if (isOutputProc)
	{
	//	get file name
		if(step >= 0){
			if(!pvd_time_filename(name, filename, step)) return false;
		}
		else{
			if(!pvd_filename(name, filename)) return false;
		}

	//	open file
		file = fopen(name.c_str(), "w");
		if (file == NULL)
		{
			UG_LOG("Cannot print to file.\n");
			return false;
		}

	// 	Write beginning of file
		fprintf(file, "<?xml version=\"1.0\"?>\n");
		fprintf(file, "<VTKFile type=\"Collection\" version=\"0.1\">\n");
		fprintf(file, "  <Collection>\n");

	// 	Include files for all subsets
		for(int si = 0; si < u.num_subsets(); ++si)
		{
			vtu_filename(name, filename, rank, si, u.num_subsets()-1, step);
#ifdef UG_PARALLEL
			if(numProcs > 1) pvtu_filename(name, filename, si, u.num_subsets()-1, step);
#endif
			name = FilenameWithoutPath(name);
			fprintf(file, "  <DataSet timestep=\"%g\" part=\"%d\" file=\"%s\"/>\n",
			        		time, si, name.c_str());
		}

	//	write closing tag
		fprintf(file, "  </Collection>\n");
		fprintf(file, "</VTKFile>\n");
		fclose(file);
	}

	if (isOutputProc && numProcs > 1)
	{
		std::string procName(filename);
		procName.append("_processwise");

	//	get file name
		if(step >= 0){
			if(!pvd_time_filename(name, filename, step)) return false;
		}
		else{
			if(!pvd_filename(name, filename)) return false;
		}

	//	open File
		file = fopen(name.c_str(), "w");
		if (file == NULL)
		{
			UG_LOG("Cannot print to file.\n");
			return false;
		}

	// 	Write to file
		fprintf(file, "<?xml version=\"1.0\"?>\n");
		fprintf(file, "<VTKFile type=\"Collection\" version=\"0.1\">\n");
		fprintf(file, "  <Collection>\n");

	// 	include files from all procs
		for (int r = 0; r < numProcs; r++)
			for(int si = 0; si < u.num_subsets(); ++si)
			{
				vtu_filename(name, filename, rank, si, u.num_subsets()-1, step);
#ifdef UG_PARALLEL
				if(numProcs > 1) pvtu_filename(name, filename, si, u.num_subsets()-1, step);
#endif
				name = FilenameWithoutPath(name);
				fprintf(file, "  <DataSet timestep=\"%g\" part=\"%d\" file=\"%s\"/>\n",
				        	time, r, name.c_str());
			}

	//	end file
		fprintf(file, "  </Collection>\n");
		fprintf(file, "</VTKFile>\n");
		fclose(file);
	}

//	we're done
	return true;
}

template <typename TDiscreteFunction>
bool
VTKOutput<TDiscreteFunction>::
select_nodal_scalar(const char* fctName, const char* name)
{
	std::vector<std::string> tokens;
	std::string fctString(fctName);
	TokenizeString(fctString, tokens, ',');
	if(tokens.size() != 1)
	{
		UG_LOG("ERROR in 'VTK:select_nodal_scalar': In order to select"
				" a nodal scalar for output to vtk,"
				" exactly one function components must be chosen.\n");
		return false;
	}

//	skip already selected
	for(size_t i = 0; i < tokens.size(); ++i)
	{
		RemoveWhitespaceFromString(tokens[i]);
		for(size_t j = 0; j < m_vSymbFct.size(); ++j)
		{
		//	skip function selected twice, if and only if same name given
			if(m_vSymbFct[j].first.size() == tokens[i].size() &&
				m_vSymbFct[j].first == tokens[i])
			{
				if(m_vSymbFct[j].second == name &&
					m_vSymbFct[j].second.size() == strlen(name)){
					return true;
				}
				else{
					UG_LOG("ERROR in 'VTK:select_nodal_scalar': Selecting component "
							<< tokens[i] << " again, but with different name " <<
							name << " instead of already scheduled "
							<< m_vSymbFct[j].second << "\n");
					return false;
				}
			}

		//	check if name is not in use
			if(m_vSymbFct[j].second == name &&
				m_vSymbFct[j].second.size() == strlen(name)){
				UG_LOG("ERROR in 'VTK:select_nodal_scalar': Selecting component "
						<< tokens[i] << ", but with already used name " <<
						name << ". This is not allowed, use different name.\n");
				return false;
			}
		}
	}

	m_vSymbFct.push_back(std::pair<std::string, std::string>(fctName, name));
	return true;
}

template <typename TDiscreteFunction>
bool
VTKOutput<TDiscreteFunction>::
select_nodal_vector(const char* fctNames, const char* name)
{
	std::vector<std::string> tokens;
	std::string fctString(fctNames);
	TokenizeString(fctString, tokens, ',');
	if(tokens.size() != (size_t)function_type::dim)
	{
		UG_LOG("ERROR in 'VTK:select_nodal_vector': In order to select"
				" a nodal vector for output to vtk,"
				" #dim function components must be chosen.\n");
		return false;
	}

	for(size_t j = 0; j < m_vSymbFct.size(); ++j)
	{
	//	check if name is not in use
		if(m_vSymbFct[j].second == name &&
			m_vSymbFct[j].second.size() == strlen(name)){
			UG_LOG("ERROR in 'VTK:select_nodal_vector': Using name " << name <<
			       " that is already used by other data is not allowed.\n");
			return false;
		}
	}

	m_vSymbFct.push_back(std::pair<std::string, std::string>(fctNames, name));
	return true;
}


}


