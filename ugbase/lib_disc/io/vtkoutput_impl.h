/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel
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

//other libraries
#include <cstdio>
#include <iostream>
#include <cstring>
#include <string>
#include <algorithm>

// ug4 libraries
#include "common/log.h"
#include "common/util/string_util.h"
#include "common/util/endian_detection.h"
#include "common/profiler/profiler.h"
#include "common/util/provider.h"
#include "lib_disc/common/function_group.h"
#include "lib_disc/common/groups_util.h"
#include "lib_disc/common/multi_index.h"
#include "lib_disc/reference_element/reference_element.h"
#include "lib_disc/local_finite_element/local_finite_element_provider.h"
#include "lib_disc/spatial_disc/disc_util/fv_util.h"
#include "lib_grid/algorithms/debug_util.h"

#ifdef UG_PARALLEL
#include "pcl/pcl_base.h"
#include "lib_algebra/parallelization/parallel_storage_type.h"
#endif

// own interface
#include "vtkoutput.h"

namespace ug{

/**
 * Writes a float in the VTU-compatible ascii format:
 */
template <int TDim>
VTKFileWriter& VTKOutput<TDim>::
write_asc_float(VTKFileWriter& File, float data)
{
	if (std::abs (data) < std::numeric_limits<float>::min ()) // a protection against the denormalized floats
		File << '0';
	else
		File << data;
	return File;
}


/* Functions that write data into the VTU file depending on the type.
 * Note that the data may be either scalars, 3-dimensional vectors or 3x3 tensors.
 * The binary values should be represented in the Float32 format, the ascii ones must
 * be representable in this format.
 */

template <int TDim>
void VTKOutput<TDim>::
write_item_to_file(VTKFileWriter& File, float data)
{
	if(m_bBinary)
		File << (float) data;
	else
		write_asc_float (File, data) << ' ';
}

// fill position data up with zeros if dim < 3.
template <int TDim>
void VTKOutput<TDim>::
write_item_to_file(VTKFileWriter& File, const ug::MathVector<1>& data)
{
	if(m_bBinary)
		File << (float) data[0] << (float) 0.f << (float) 0.f;
	else
		write_asc_float (File, data[0]) << " 0 0 ";
}

template <int TDim>
void VTKOutput<TDim>::
write_item_to_file(VTKFileWriter& File, const ug::MathVector<2>& data)
{
	if(m_bBinary)
		File << (float) data[0] << (float) data[1] << (float) 0.f;
	else
	{
		write_asc_float (File, data[0]) << ' ';
		write_asc_float (File, data[1]) << " 0 ";
	}
}

template <int TDim>
void VTKOutput<TDim>::
write_item_to_file(VTKFileWriter& File, const ug::MathVector<3>& data)
{
	if(m_bBinary)
		File << (float) data[0] << (float) data[1] << (float) data[2];
	else
	{
		write_asc_float (File, data[0]) << ' ';
		write_asc_float (File, data[1]) << ' ';
		write_asc_float (File, data[2]) << ' ';
	}
}

template <int TDim>
void VTKOutput<TDim>::
write_item_to_file(VTKFileWriter& File, const ug::MathMatrix<1,1>& data)
{
	if(m_bBinary)
		File << (float) data(0,0) << (float) 0.f << (float) 0.f << (float) 0.f << (float) 0.f << (float) 0.f << (float) 0.f << (float) 0.f << (float) 0.f;
	else
		write_asc_float (File, data(0,0)) << " 0 0 0 0 0 0 0 0 ";
}

template <int TDim>
void VTKOutput<TDim>::
write_item_to_file(VTKFileWriter& File, const ug::MathMatrix<2,2>& data)
{
	if(m_bBinary)
	{
		File << (float) data(0,0) << (float) data(0,1) << (float) 0.f;
		File << (float) data(1,0) << (float) data(1,1) << (float) 0.f << (float) 0.f << (float) 0.f << (float) 0.f;
	}
	else
	{
		write_asc_float (File, data(0,0)) << ' '; write_asc_float (File, data(0,1)) << " 0 ";
		write_asc_float (File, data(1,0)) << ' '; write_asc_float (File, data(1,1)) << " 0 0 0 0 ";
	}
}

template <int TDim>
void VTKOutput<TDim>::
write_item_to_file(VTKFileWriter& File, const ug::MathMatrix<3,3>& data)
{
	if(m_bBinary)
	{
		File << (float) data(0,0) << (float) data(0,1) << (float) data(0,2);
		File << (float) data(1,0) << (float) data(1,1) << (float) data(1,2);
		File << (float) data(2,0) << (float) data(2,1) << (float) data(2,2);
	}
	else
	{
		write_asc_float (File, data(0,0)) << ' '; write_asc_float (File, data(0,1)) << ' '; write_asc_float (File, data(0,2)) << ' ';
		write_asc_float (File, data(1,0)) << ' '; write_asc_float (File, data(1,1)) << ' '; write_asc_float (File, data(1,2)) << ' ';
		write_asc_float (File, data(2,0)) << ' '; write_asc_float (File, data(2,1)) << ' '; write_asc_float (File, data(2,2)) << ' ';
	}
}

/**
 * VTU output for the data based on the given grid function (u), for all the subsets
 * where this function is defined.
 */
template <int TDim>
template <typename TFunction>
void VTKOutput<TDim>::
print(const char* filename, TFunction& u, int step, number time, bool makeConsistent)
{
	PROFILE_FUNC();
#ifdef UG_PARALLEL
	if(makeConsistent)
		if(!u.change_storage_type(PST_CONSISTENT))
			UG_THROW("VTK::print: Cannot change storage type to consistent.");
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
	//	write whole grid to a single file
		try
		{
			SubsetGroup all_ss (u.domain()->subset_handler ());
			all_ss.add_all ();
			print_subsets(filename, u, all_ss, step, time, makeConsistent);
		}
		UG_CATCH_THROW("VTK::print: Failed to write grid.");
	}
	else
	{
		// 	loop subsets
		for(int si = 0; si < u.num_subsets(); ++si)
		{
		//	write each subset to a single file
			try{
				print_subset(filename, u, si, step, time, makeConsistent);
			}
			UG_CATCH_THROW("VTK::print: Failed to write Subset "<< si << ".");
		}

		//	write grouping pvd file
		try{
			write_subset_pvd(u.num_subsets(), filename, step, time);
		}
		UG_CATCH_THROW("VTK::print: Failed to write pvd-file.");
	}
	#ifdef UG_PARALLEL
		PCL_DEBUG_BARRIER_ALL();
	#endif
}


/**
 * VTU output for the data based on the given grid function (u), for only one given subset
 * si. The file name is composed with the number of the subset.
 */
template <int TDim>
template <typename TFunction>
void VTKOutput<TDim>::
print_subset(const char* filename, TFunction& u, int si, int step, number time, bool makeConsistent)
{
	PROFILE_FUNC();

#ifdef UG_PARALLEL
	if(makeConsistent)
		if(!u.change_storage_type(PST_CONSISTENT))
			UG_THROW("VTK::print_subset: Cannot change storage type to consistent.");
#endif

//	get the grid associated to the solution
	Grid& grid = *u.domain()->grid();
	SubsetGroup ssg (u.domain()->subset_handler ());
	ssg.add (si);

// 	attach help indices
	typedef ug::Attachment<int> AVrtIndex;
	AVrtIndex aVrtIndex;
	Grid::VertexAttachmentAccessor<AVrtIndex> aaVrtIndex;
	grid.attach_to_vertices(aVrtIndex);
	aaVrtIndex.access(grid, aVrtIndex);

//	get rank of process
	int rank = 0;
#ifdef UG_PARALLEL
	rank = pcl::ProcRank();
#endif

//	get name for *.vtu file
	std::string name;
	try{
		vtu_filename(name, filename, rank, si, u.num_subsets()-1, step);
	}
	UG_CATCH_THROW("VTK::print_subset: Failed to write vtu-file.");


//	open the file
	try
	{
		VTKFileWriter File(name.c_str());

	//	bool if time point should be written to *.vtu file
	//	in parallel we must not (!) write it to the *.vtu file, but to the *.pvtu
		bool bTimeDep = (step >= 0);
#ifdef UG_PARALLEL
		if(pcl::NumProcs() > 1) bTimeDep = false;
#endif

	//	header
		File << VTKFileWriter::normal;
		File << "<?xml version=\"1.0\"?>\n";

		write_comment(File);

		File << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"";
		if(IsLittleEndian()) File << "LittleEndian";
		else File << "BigEndian";
		File << "\">\n";

	//	writing time point
		if(bTimeDep)
		{
			File << "  <Time timestep=\""<<time<<"\"/>\n";
		}

	//	opening the grid
		File << "  <UnstructuredGrid>\n";

	// 	get dimension of grid-piece
		int dim = DimensionOfSubset(*u.domain()->subset_handler(), si);

	//	write piece of grid
		if(dim >= 0)
		{
			try{
				write_grid_solution_piece(File, aaVrtIndex, grid, u, time, ssg, dim);
			}
			UG_CATCH_THROW("VTK::print_subset: Can not write Subset: "<<si);
		}
		else
		{
		//	if dim < 0, something is wrong with grid, unless there are no elements in the grid
			if((si >=0) && u.domain()->subset_handler()->template num<Vertex>(si) != 0)
			{
				UG_THROW("VTK::print_subset: Dimension of grid/subset not"
						" detected correctly although grid objects present.");
			}

			write_empty_grid_piece(File);
		}

	//	write closing xml tags
		File << VTKFileWriter::normal;
		File << "  </UnstructuredGrid>\n";
		File << "</VTKFile>\n";

	// 	detach help indices
		grid.detach_from_vertices(aVrtIndex);

#ifdef UG_PARALLEL
	//	write grouping *.pvtu file in parallel case
		try{
			write_pvtu(u, filename, si, step, time);
		}
		UG_CATCH_THROW("VTK::print_subset: Failed to write pvtu-file.");
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

	}
	UG_CATCH_THROW("VTK::print_subset: Failed to open file: '" << filename << "' for output");
	#ifdef UG_PARALLEL
		PCL_DEBUG_BARRIER_ALL();
	#endif
}


/**
 * VTU output for the data based on the given grid function (u), only for the subsets from
 * a given subset group.
 */
template <int TDim>
template <typename TFunction>
void VTKOutput<TDim>::
print_subsets(const char* filename, TFunction& u, const SubsetGroup& ssGrp, int step, number time, bool makeConsistent)
{
	PROFILE_FUNC();

#ifdef UG_PARALLEL
	if(makeConsistent)
		if(!u.change_storage_type(PST_CONSISTENT))
			UG_THROW("VTK::print_subset: Cannot change storage type to consistent.");
#endif

//	get the grid associated to the solution
	Grid& grid = *u.domain()->grid();
	
// 	attach help indices
	typedef ug::Attachment<int> AVrtIndex;
	AVrtIndex aVrtIndex;
	Grid::VertexAttachmentAccessor<AVrtIndex> aaVrtIndex;
	grid.attach_to_vertices(aVrtIndex);
	aaVrtIndex.access(grid, aVrtIndex);

//	get rank of process
	int rank = 0;
#ifdef UG_PARALLEL
	rank = pcl::ProcRank();
#endif

//	get name for *.vtu file
	std::string name;
	try{
		vtu_filename(name, filename, rank, -1, u.num_subsets()-1, step); // "si == -1" because we do not want any subset prefixes!
	}
	UG_CATCH_THROW("VTK::print_subsets: Failed to write vtu-file.");


//	open the file
	try
	{
		VTKFileWriter File(name.c_str());

	//	bool if time point should be written to *.vtu file
	//	in parallel we must not (!) write it to the *.vtu file, but to the *.pvtu
		bool bTimeDep = (step >= 0);
#ifdef UG_PARALLEL
		if(pcl::NumProcs() > 1) bTimeDep = false;
#endif

	//	header
		File << VTKFileWriter::normal;
		File << "<?xml version=\"1.0\"?>\n";

		write_comment(File);

		File << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"";
		if(IsLittleEndian()) File << "LittleEndian";
		else File << "BigEndian";
		File << "\">\n";

	//	writing time point
		if(bTimeDep)
		{
			File << "  <Time timestep=\""<<time<<"\"/>\n";
		}

	//	opening the grid
		File << "  <UnstructuredGrid>\n";

	// 	get dimension of grid-piece: the highest dimension of the specified subsets
		int dim = -1;
		for(size_t i = 0; i < ssGrp.size(); i++)
		{
			int ssDim = DimensionOfSubset(*u.domain()->subset_handler(), ssGrp[i]);
			if(dim < ssDim)
				dim = ssDim;
		}

	//	write piece of grid
		if(dim >= 0)
		{
			try{
				write_grid_solution_piece(File, aaVrtIndex, grid, u, time, ssGrp, dim);
			}
			UG_CATCH_THROW("VTK::print_subsets: Can not write the subsets");
		}
		else
		{
			UG_THROW("VTK::print_subsets: Dimension of grid/subset not"
					" detected correctly although grid objects present.");
		}

	//	write closing xml tags
		File << VTKFileWriter::normal;
		File << "  </UnstructuredGrid>\n";
		File << "</VTKFile>\n";

	// 	detach help indices
		grid.detach_from_vertices(aVrtIndex);

#ifdef UG_PARALLEL
	//	write grouping *.pvtu file in parallel case
		try{
			write_pvtu(u, filename, -1, step, time); // "-1" because we do not want any subset prefixes!
		}
		UG_CATCH_THROW("VTK::print_subsets: Failed to write pvtu-file.");
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

	}
	UG_CATCH_THROW("VTK::print_subset: Failed open file: '" << filename << "' for output");
	#ifdef UG_PARALLEL
		PCL_DEBUG_BARRIER_ALL();
	#endif
}


/**
 * VTU output for the data based on the given grid function (u), only for the subsets from
 * a given list of subset names.
 */
template <int TDim>
template <typename TFunction>
void VTKOutput<TDim>::
print_subsets(const char* filename, TFunction& u, const char* ssNames, int step, number time, bool makeConsistent)
{
	std::vector<std::string> v_ss_names;
	SubsetGroup ssGrp (u.domain()->subset_handler());
	
	TokenizeTrimString (ssNames, v_ss_names);
	ssGrp.add (v_ss_names);
	print_subsets(filename, u, ssGrp, step, time, makeConsistent);
}


////////////////////////////////////////////////////////////////////////////////
//	writing pieces
////////////////////////////////////////////////////////////////////////////////


/**
 * This function composes the vtu sections with the coordinates of the grid points and
 * the connectivities (cells), a part of the unstructured grid description in the vtu file.
 */
template <int TDim>
template <typename T>
void VTKOutput<TDim>::
write_points_cells_piece(VTKFileWriter& File,
                         Grid::VertexAttachmentAccessor<Attachment<int> >& aaVrtIndex,
                         const Grid::VertexAttachmentAccessor<Attachment<MathVector<TDim> > >& aaPos,
                         Grid& grid, const T& iterContainer, const SubsetGroup& ssGrp, const int dim,
                         const int numVert, const int numElem, const int numConn)
{
//	write vertices of this piece
	try{
		write_points<T>(File, aaVrtIndex, aaPos, grid, iterContainer, ssGrp, dim, numVert);
	}
	UG_CATCH_THROW("VTK::write_points_cells_piece: Failed to write vertices (Points).");

//	write elements of this piece
	try{
		write_cells(File, aaVrtIndex, grid, iterContainer, ssGrp, dim, numElem, numConn);
	}
	UG_CATCH_THROW("VTK::write_points_cells_piece: Failed to write grid elements (Cells).");
}

/**
 * Composes only the vtu-piece of the description of the unstructured grid.
 */
template <int TDim>
template <typename T>
void VTKOutput<TDim>::
write_grid_piece(VTKFileWriter& File,
                 Grid::VertexAttachmentAccessor<Attachment<int> >& aaVrtIndex,
                 const Grid::VertexAttachmentAccessor<Attachment<MathVector<TDim> > >& aaPos,
                 Grid& grid, const T& iterContainer, SubsetGroup& ssGrp, int dim)
{
//	counters
	int numVert = 0, numElem = 0, numConn = 0;

// 	Count needed sizes for vertices, elements and connections
	try{
		count_piece_sizes(grid, iterContainer, ssGrp, dim, numVert, numElem, numConn);
	}
	UG_CATCH_THROW("VTK::write_piece: Failed to count piece sizes.");

//	write the beginning of the piece, indicating the number of vertices
//	and the number of elements for this piece of the grid.
	File << VTKFileWriter::normal;
	File << "    <Piece NumberOfPoints=\""<<numVert<<
	"\" NumberOfCells=\""<<numElem<<"\">\n";

//	write grid
	write_points_cells_piece<T>
	(File, aaVrtIndex, aaPos, grid, iterContainer, ssGrp, dim, numVert, numElem, numConn);

//	write closing tag
	File << VTKFileWriter::normal;
	File << "    </Piece>\n";
}

/**
 * This function writes a piece of the grid to the vtk file. First the
 * geometric data is written (points, cells, connectivity). Then the data
 * associated to the points or cells is written.
 */
template <int TDim>
template <typename TFunction>
void VTKOutput<TDim>::
write_grid_solution_piece(VTKFileWriter& File,
                          Grid::VertexAttachmentAccessor<Attachment<int> >& aaVrtIndex,
                          Grid& grid, TFunction& u, number time,
                          const SubsetGroup& ssGrp, const int dim)
{
//	counters
	int numVert = 0, numElem = 0, numConn = 0;

// 	Count needed sizes for vertices, elements and connections
	try{
		count_piece_sizes(grid, u, ssGrp, dim, numVert, numElem, numConn);
	}
	UG_CATCH_THROW("VTK::write_piece: Failed to count piece sizes.");

//	write the beginning of the piece, indicating the number of vertices
//	and the number of elements for this piece of the grid.
	File << VTKFileWriter::normal;
	File << "    <Piece NumberOfPoints=\""<<numVert<<
	"\" NumberOfCells=\""<<numElem<<"\">\n";

//	write grid
	write_points_cells_piece<TFunction>
	(File, aaVrtIndex, u.domain()->position_accessor(), grid, u, ssGrp, dim, numVert, numElem, numConn);

//	add all components if 'selectAll' chosen
	if(m_bSelectAll){
		for(size_t fct = 0; fct < u.num_fct(); ++fct){
			if(!vtk_name_used(u.name(fct).c_str())){
				if(LocalFiniteElementProvider::continuous(u.local_finite_element_id(fct))){
					select_nodal(u.name(fct).c_str(), u.name(fct).c_str());
				}else{
					select_element(u.name(fct).c_str(), u.name(fct).c_str());
				}
			}
		}
	}

//	write nodal data
	write_nodal_values_piece(File, u, time, grid, ssGrp, dim, numVert);

//	write cell data
	write_cell_values_piece(File, u, time, grid, ssGrp, dim, numElem);

//	write closing tag
	File << VTKFileWriter::normal;
	File << "    </Piece>\n";
}

////////////////////////////////////////////////////////////////////////////////
// Sizes
////////////////////////////////////////////////////////////////////////////////

template <int TDim>
template <typename TElem, typename T>
void VTKOutput<TDim>::
count_sizes(Grid& grid, const T& iterContainer, const int si,
            int& numVert, int& numElem, int& numConn)
{
//	get reference element
	typedef typename reference_element_traits<TElem>::reference_element_type
																ref_elem_type;

//	number of corners of element
	static const int numCo = ref_elem_type::numCorners;

//	get iterators
	typedef typename IteratorProvider<T>::template traits<TElem>::const_iterator const_iterator;
	const_iterator iterBegin = IteratorProvider<T>::template begin<TElem>(iterContainer, si);
	const_iterator iterEnd = IteratorProvider<T>::template end<TElem>(iterContainer, si);

//	loop elements
	for( ; iterBegin != iterEnd; ++iterBegin)
	{
	//	get the element
		TElem *elem = *iterBegin;

	//	count number of elements and number of connections;
	//	handle octahedrons separately by splitting into a top and bottom pyramid
		if(ref_elem_type::REFERENCE_OBJECT_ID != ROID_OCTAHEDRON)
		{
			++numElem;
			numConn += numCo;
		}
		else
		{
		// 	counting top and bottom pyramid
			numElem += 2;
			numConn += 10;
		}

	//	loop vertices of the element
		for(int i = 0; i < numCo; ++i)
		{
		//	get vertex of the element
			Vertex* v = GetVertex(elem,i);

		//	if this vertex has already been counted, skip it
			if(grid.is_marked(v)) continue;

		// count vertex and mark it
			++numVert;
			grid.mark(v);
		}
	}
};

/**
 * Counts VTK grid elements (points, full-dim. elements and connections) in a given
 * piece of the ug grid.
 */
template <int TDim>
template <typename T>
void VTKOutput<TDim>::
count_piece_sizes(Grid& grid, const T& iterContainer, const SubsetGroup& ssGrp, const int dim,
                  int& numVert, int& numElem, int& numConn)
{
	numVert = numElem = numConn = 0;

//	reset all marks (we use vertex marking to avoid counting vertices twice)
	grid.begin_marking();

	for(size_t i = 0; i < ssGrp.size(); i++)
	{
		const int si = ssGrp[i];
		switch(dim) // switch dimension
		{
			case 0: count_sizes<Vertex, T>(grid, iterContainer, si, numVert, numElem, numConn); break;
			case 1: count_sizes<RegularEdge, T>(grid, iterContainer, si, numVert, numElem, numConn);
					count_sizes<ConstrainingEdge, T>(grid, iterContainer, si, numVert, numElem, numConn); break;
			case 2: count_sizes<Triangle, T>(grid, iterContainer, si, numVert, numElem, numConn);
					count_sizes<Quadrilateral, T>(grid, iterContainer, si, numVert, numElem, numConn);
					count_sizes<ConstrainingTriangle, T>(grid, iterContainer, si, numVert, numElem, numConn);
					count_sizes<ConstrainingQuadrilateral, T>(grid, iterContainer, si, numVert, numElem, numConn); break;
			case 3: count_sizes<Tetrahedron, T>(grid, iterContainer, si, numVert, numElem, numConn);
					count_sizes<Pyramid, T>(grid, iterContainer, si, numVert, numElem, numConn);
					count_sizes<Prism, T>(grid, iterContainer, si, numVert, numElem, numConn);
					count_sizes<Octahedron, T>(grid, iterContainer, si, numVert, numElem, numConn);
					count_sizes<Hexahedron, T>(grid, iterContainer, si, numVert, numElem, numConn); break;
			default: UG_THROW("VTK::count_piece_sizes: Dimension " << dim <<
									" is not supported.");
		}
	}

//	signal end of marking
	grid.end_marking();
}

////////////////////////////////////////////////////////////////////////////////
// Points
////////////////////////////////////////////////////////////////////////////////

/**
 * This method loops all elements of a subset and writes the vertex
 * coordinates to a file. All vertices of the element are written reguardless
 * if the vertex itself is contained in the subset or another. Written
 * vertices are marked and not written again.
 */
template <int TDim>
template <typename TElem, typename T>
void VTKOutput<TDim>::
write_points_elementwise(VTKFileWriter& File,
                         Grid::VertexAttachmentAccessor<Attachment<int> >& aaVrtIndex,
                         const Grid::VertexAttachmentAccessor<Attachment<MathVector<TDim> > >& aaPos,
                         Grid& grid, const T& iterContainer, const int si, int& n)
{
//	get reference element
	typedef typename reference_element_traits<TElem>::reference_element_type
																ref_elem_type;
//	get iterators
	typedef typename IteratorProvider<T>::template traits<TElem>::const_iterator const_iterator;
	const_iterator iterBegin = IteratorProvider<T>::template begin<TElem>(iterContainer, si);
	const_iterator iterEnd = IteratorProvider<T>::template end<TElem>(iterContainer, si);
	if(m_bBinary)
		File << VTKFileWriter::base64_binary;
	else
		File << VTKFileWriter::normal;

//	loop all elements of the subset
	for( ; iterBegin != iterEnd; ++iterBegin)
	{
	//	get the element
		TElem *elem = *iterBegin;

	//	loop vertices of the element
		for(size_t i = 0; i < (size_t) ref_elem_type::numCorners; ++i)
		{
		//	get vertex of element
			Vertex* v = GetVertex(elem, i);

		//	if vertex has already be handled, skip it
			if(grid.is_marked(v)) continue;

		//	mark the vertex as processed
			grid.mark(v);

		//	number vertex
			aaVrtIndex[v] = n++;

		//	get position of vertex and write position to stream
			write_item_to_file(File, aaPos[v]);
		}
	}
}

/**
 * This function writes the vertices of a piece (the specified subsets) of
 * the grid to a vtk file.
 */
template <int TDim>
template <typename T>
void VTKOutput<TDim>::
write_points(VTKFileWriter& File,
             Grid::VertexAttachmentAccessor<Attachment<int> >& aaVrtIndex,
             const Grid::VertexAttachmentAccessor<Attachment<MathVector<TDim> > >& aaPos,
             Grid& grid, const T& iterContainer, const SubsetGroup& ssGrp, const int dim,
             const int numVert)
{
	if(!m_bWriteGrid){
		return;
	}

//	write starting xml tag for points
	File << VTKFileWriter::normal;
	File << "      <Points>\n";
	File << "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format="
		 <<	(m_bBinary ? "\"binary\"" : "\"ascii\"") << ">\n";
	int n = 3*sizeof(float) * numVert;
	if(m_bBinary)
		File << VTKFileWriter::base64_binary << n;

//	reset counter for vertices
	n = 0;

//	start marking of vertices
	grid.begin_marking();

//	switch dimension
	if(numVert > 0)
	for(size_t i = 0; i < ssGrp.size(); i++){
		switch(dim){
			case 0: write_points_elementwise<Vertex,T>(File, aaVrtIndex, aaPos, grid, iterContainer, ssGrp[i], n); break;
			case 1: write_points_elementwise<RegularEdge,T>(File, aaVrtIndex, aaPos, grid, iterContainer, ssGrp[i], n);
					write_points_elementwise<ConstrainingEdge,T>(File, aaVrtIndex, aaPos, grid, iterContainer, ssGrp[i], n); break;
			case 2: write_points_elementwise<Triangle,T>(File, aaVrtIndex, aaPos, grid, iterContainer, ssGrp[i], n);
					write_points_elementwise<Quadrilateral,T>(File, aaVrtIndex, aaPos, grid, iterContainer, ssGrp[i], n);
					write_points_elementwise<ConstrainingTriangle,T>(File, aaVrtIndex, aaPos, grid, iterContainer, ssGrp[i], n);
					write_points_elementwise<ConstrainingQuadrilateral,T>(File, aaVrtIndex, aaPos, grid, iterContainer, ssGrp[i], n); break;
			case 3:	write_points_elementwise<Tetrahedron,T>(File, aaVrtIndex, aaPos, grid, iterContainer, ssGrp[i], n);
					write_points_elementwise<Pyramid,T>(File, aaVrtIndex, aaPos, grid, iterContainer, ssGrp[i], n);
					write_points_elementwise<Prism,T>(File, aaVrtIndex, aaPos, grid, iterContainer, ssGrp[i], n);
					write_points_elementwise<Octahedron,T>(File, aaVrtIndex, aaPos, grid, iterContainer, ssGrp[i], n);
					write_points_elementwise<Hexahedron,T>(File, aaVrtIndex, aaPos, grid, iterContainer, ssGrp[i], n); break;
			default: UG_THROW("VTK::write_points: Dimension " << dim <<
			                        " is not supported.");
		}
	}

//	signal end of marking the grid
	grid.end_marking();

//	write closing tags
	File << VTKFileWriter::normal;
	File << "\n        </DataArray>\n";
	File << "      </Points>\n";
}

////////////////////////////////////////////////////////////////////////////////
// Cells
////////////////////////////////////////////////////////////////////////////////

template <int TDim>
void VTKOutput<TDim>::
write_cell_subset_names(VTKFileWriter& File, MGSubsetHandler &sh)
{
	File << VTKFileWriter::normal;
	File << "      <RegionInfo Name=\"regions\">\n";
	for(int i = 0; i < sh.num_subsets(); ++i){
		File << "        <Region Name=\"" << sh.get_subset_name(i) << "\"></Region>\n";
	}
	File << "      </RegionInfo>\n";
}

/**
 * This function writes the elements that are part of the specified subsets.
 */
template <int TDim>
template <typename T>
void VTKOutput<TDim>::
write_cells(VTKFileWriter& File,
            Grid::VertexAttachmentAccessor<Attachment<int> >& aaVrtIndex,
            Grid& grid, const T& iterContainer, const SubsetGroup& ssGrp, const int dim,
            const int numElem, const int numConn)
{
	File << VTKFileWriter::normal;

//	write opening tag to indicate that elements will be written
	File << "      <Cells>\n";

//	write connectivities of elements
	write_cell_connectivity(File, aaVrtIndex, grid, iterContainer, ssGrp, dim, numConn);

//	write offsets for elements (i.e. number of nodes counted up)
	write_cell_offsets(File, iterContainer, ssGrp, dim, numElem);

//	write a defined type for each cell
	write_cell_types(File, iterContainer, ssGrp, dim, numElem);

//	write closing tag
	File << VTKFileWriter::normal;
	File << "      </Cells>\n";
}

////////////////////////////////////////////////////////////////////////////////
// Connectivity
////////////////////////////////////////////////////////////////////////////////


template <int TDim>
template <class TElem, typename T>
void VTKOutput<TDim>::
write_cell_connectivity(VTKFileWriter& File,
                        Grid::VertexAttachmentAccessor<Attachment<int> >& aaVrtIndex,
                        Grid& grid, const T& iterContainer, const int si)
{
//	get reference element type
	typedef typename reference_element_traits<TElem>::reference_element_type
																ref_elem_type;

//	get reference element id
	static const ReferenceObjectID refID = ref_elem_type::REFERENCE_OBJECT_ID;

//	get iterators
	typedef typename IteratorProvider<T>::template traits<TElem>::const_iterator const_iterator;
	const_iterator iterBegin = IteratorProvider<T>::template begin<TElem>(iterContainer, si);
	const_iterator iterEnd = IteratorProvider<T>::template end<TElem>(iterContainer, si);

	if(m_bBinary)
		File << VTKFileWriter::base64_binary;
	else
		File << VTKFileWriter::normal;

//	loop all elements
	for( ; iterBegin != iterEnd; iterBegin++)
	{
	//	get element
		TElem* elem = *iterBegin;

	//	write ids of the element
		if(refID != ROID_PRISM && refID != ROID_OCTAHEDRON)
		{
			for(size_t i=0; i< (size_t) ref_elem_type::numCorners; i++)
			{
				Vertex* vert = elem->vertex(i);
				int id = aaVrtIndex[vert];
				File << id;
				if(!m_bBinary)
					File << ' ';
			}
		}
		else if(refID == ROID_PRISM)
		{
			int id = aaVrtIndex[elem->vertex(0)]; File << id;
			if(!m_bBinary) File << ' ';
			id = aaVrtIndex[elem->vertex(2)]; File << id;
			if(!m_bBinary) File << ' ';
			id = aaVrtIndex[elem->vertex(1)]; File << id;
			if(!m_bBinary) File << ' ';
			id = aaVrtIndex[elem->vertex(3)]; File << id;
			if(!m_bBinary) File << ' ';
			id = aaVrtIndex[elem->vertex(5)]; File << id;
			if(!m_bBinary) File << ' ';
			id = aaVrtIndex[elem->vertex(4)]; File << id;
			if(!m_bBinary) File << ' ';
		}
		//	case == ROID_OCTAHEDRON (handle by a splitting into a top and bottom pyramid)
		else
		{
		//	top pyramid
			int id = aaVrtIndex[elem->vertex(1)]; File << id;
			if(!m_bBinary) File << ' ';
			id = aaVrtIndex[elem->vertex(2)]; File << id;
			if(!m_bBinary) File << ' ';
			id = aaVrtIndex[elem->vertex(3)]; File << id;
			if(!m_bBinary) File << ' ';
			id = aaVrtIndex[elem->vertex(4)]; File << id;
			if(!m_bBinary) File << ' ';
			id = aaVrtIndex[elem->vertex(5)]; File << id;
			if(!m_bBinary) File << ' ';

		//	bottom pyramid
			id = aaVrtIndex[elem->vertex(1)]; File << id;
			if(!m_bBinary) File << ' ';
			id = aaVrtIndex[elem->vertex(2)]; File << id;
			if(!m_bBinary) File << ' ';
			id = aaVrtIndex[elem->vertex(3)]; File << id;
			if(!m_bBinary) File << ' ';
			id = aaVrtIndex[elem->vertex(4)]; File << id;
			if(!m_bBinary) File << ' ';
			id = aaVrtIndex[elem->vertex(0)]; File << id;
			if(!m_bBinary) File << ' ';
		}
	}
}

/**
 * This method writes the 'connectivity' for each element of a subset. The
 * connectivity are the indices of all vertex the element is formed of.
 */
template <int TDim>
template <typename T>
void VTKOutput<TDim>::
write_cell_connectivity(VTKFileWriter& File,
                        Grid::VertexAttachmentAccessor<Attachment<int> >& aaVrtIndex,
                        Grid& grid, const T& iterContainer, const SubsetGroup& ssGrp,
                        const int dim, const int numConn)
{
	File << VTKFileWriter::normal;
//	write opening tag to indicate that connections will be written
	File << "        <DataArray type=\"Int32\" Name=\"connectivity\" format="
		 <<	(m_bBinary ? "\"binary\"": "\"ascii\"") << ">\n";
	int n = sizeof(int) * numConn;

	if(m_bBinary)
		File << VTKFileWriter::base64_binary << n;
//	switch dimension
	if(numConn > 0)
	for(size_t i = 0; i < ssGrp.size(); i++){
		switch(dim)
		{
			case 0: break; // no elements -> nothing to do
			case 1: write_cell_connectivity<RegularEdge,T>(File, aaVrtIndex, grid, iterContainer, ssGrp[i]);
					write_cell_connectivity<ConstrainingEdge,T>(File, aaVrtIndex, grid, iterContainer, ssGrp[i]); break;
			case 2: write_cell_connectivity<Triangle,T>(File, aaVrtIndex, grid, iterContainer, ssGrp[i]);
					write_cell_connectivity<Quadrilateral,T>(File, aaVrtIndex, grid, iterContainer, ssGrp[i]);
					write_cell_connectivity<ConstrainingTriangle,T>(File, aaVrtIndex, grid, iterContainer, ssGrp[i]);
					write_cell_connectivity<ConstrainingQuadrilateral,T>(File, aaVrtIndex, grid, iterContainer, ssGrp[i]); break;
			case 3: write_cell_connectivity<Tetrahedron,T>(File, aaVrtIndex, grid, iterContainer, ssGrp[i]);
					write_cell_connectivity<Pyramid,T>(File, aaVrtIndex, grid, iterContainer, ssGrp[i]);
					write_cell_connectivity<Prism,T>(File, aaVrtIndex, grid, iterContainer, ssGrp[i]);
					write_cell_connectivity<Octahedron,T>(File, aaVrtIndex, grid, iterContainer, ssGrp[i]);
					write_cell_connectivity<Hexahedron,T>(File, aaVrtIndex, grid, iterContainer, ssGrp[i]); break;
			default: UG_THROW("VTK::write_cell_connectivity: Dimension " << dim <<
			                        " is not supported.");
		}
	}

//	write closing tag
	File << VTKFileWriter::normal;
	File << "\n        </DataArray>\n";
}

////////////////////////////////////////////////////////////////////////////////
// Offset
////////////////////////////////////////////////////////////////////////////////


template <int TDim>
template <class TElem, typename T>
void VTKOutput<TDim>::
write_cell_offsets(VTKFileWriter& File, const T& iterContainer, const int si, int& n)
{
//	get reference element
	typedef typename reference_element_traits<TElem>::reference_element_type
																ref_elem_type;

//	get iterators
	typedef typename IteratorProvider<T>::template traits<TElem>::const_iterator const_iterator;
	const_iterator iterBegin = IteratorProvider<T>::template begin<TElem>(iterContainer, si);
	const_iterator iterEnd = IteratorProvider<T>::template end<TElem>(iterContainer, si);

	if(m_bBinary)
		File << VTKFileWriter::base64_binary;
	else
		File << VTKFileWriter::normal;

//	loop all elements
	for( ; iterBegin != iterEnd; ++iterBegin)
	{
	//	handle octahedron separately by splitting into a top and bottom pyramid
		if(ref_elem_type::REFERENCE_OBJECT_ID != ROID_OCTAHEDRON)
		{
		//	increase counter of vertices
			n += ref_elem_type::numCorners;

		//	write offset
			File << n;
			if(!m_bBinary)
				File << ' ';
		}
		else
		{
		//	increase counter for top pyramid vertices
			n += 5;

		//	write offset for top pyramid
			File << n;
			if(!m_bBinary)
				File << ' ';

		//	increase counter for bottom pyramid vertices
			n += 5;

		//	write offset for bottom pyramid
			File << n;
			if(!m_bBinary)
				File << ' ';
		}
	}
}

/**
 * This method writes the 'offset' for each element of a subset.
 */
template <int TDim>
template <typename T>
void VTKOutput<TDim>::
write_cell_offsets(VTKFileWriter& File, const T& iterContainer, const SubsetGroup& ssGrp,
                   const int dim, const int numElem)
{
	File << VTKFileWriter::normal;
//	write opening tag indicating that offsets are going to be written
	File << "        <DataArray type=\"Int32\" Name=\"offsets\" format="
		 <<	(m_bBinary ? "\"binary\"": "\"ascii\"") << ">\n";
	int n = sizeof(int) * numElem;
	if(m_bBinary)
		File << VTKFileWriter::base64_binary << n;

	n = 0;
//	switch dimension
	if(numElem > 0)
	for(size_t i = 0; i < ssGrp.size(); i++){
		switch(dim)
		{
			case 0: break; // no elements -> nothing to do
			case 1: write_cell_offsets<RegularEdge,T>(File, iterContainer, ssGrp[i], n);
					write_cell_offsets<ConstrainingEdge,T>(File, iterContainer, ssGrp[i], n); break;
			case 2: write_cell_offsets<Triangle,T>(File, iterContainer, ssGrp[i], n);
					write_cell_offsets<Quadrilateral,T>(File, iterContainer, ssGrp[i], n);
					write_cell_offsets<ConstrainingTriangle,T>(File, iterContainer, ssGrp[i], n);
					write_cell_offsets<ConstrainingQuadrilateral,T>(File, iterContainer, ssGrp[i], n); break;
			case 3: write_cell_offsets<Tetrahedron,T>(File, iterContainer, ssGrp[i], n);
					write_cell_offsets<Pyramid,T>(File, iterContainer, ssGrp[i], n);
					write_cell_offsets<Prism,T>(File, iterContainer, ssGrp[i], n);
					write_cell_offsets<Octahedron,T>(File, iterContainer, ssGrp[i], n);
					write_cell_offsets<Hexahedron,T>(File, iterContainer, ssGrp[i], n); break;
			default: UG_THROW("VTK::write_cell_offsets: Dimension " << dim <<
			                        " is not supported.");
		}
	}

//	closing tag
	File << VTKFileWriter::normal;
	File << "\n        </DataArray>\n";
}

////////////////////////////////////////////////////////////////////////////////
// Types
////////////////////////////////////////////////////////////////////////////////


template <int TDim>
template <class TElem, typename T>
void VTKOutput<TDim>::
write_cell_types(VTKFileWriter& File, const T& iterContainer, const int si)
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
		case ROID_OCTAHEDRON: type = (char) 14; break; // use type id == 14 (PYRAMID) as we handle an octahedron by splitting into a top and bottom pyramid
		case ROID_HEXAHEDRON: type = (char) 12; break;
		default: UG_THROW("Element Type not known.");
	}

//	get iterators
	typedef typename IteratorProvider<T>::template traits<TElem>::const_iterator const_iterator;
	const_iterator iterBegin = IteratorProvider<T>::template begin<TElem>(iterContainer, si);
	const_iterator iterEnd = IteratorProvider<T>::template end<TElem>(iterContainer, si);

	if(m_bBinary)
		File << VTKFileWriter::base64_binary;
	else
		File << VTKFileWriter::normal;
//	loop all elements, write type for each element to stream
	for( ; iterBegin != iterEnd; ++iterBegin)
	{
	//	handle octahedral case separately by splitting into a top and bottom pyramid
		if(ref_elem_type::REFERENCE_OBJECT_ID != ROID_OCTAHEDRON)
		{
			if(m_bBinary)
				File << type;
			 else
				File << (int) type << ' ';
		}
		else
		{
			if(m_bBinary)
			{
				File << type; // top pyramid
				File << type; // bottom pyramid
			}
			else
			{
				File << (int) type << ' '; // top pyramid
				File << (int) type << ' '; // bottom pyramid
			}
		}
	}
}

/**
 * This method writes the 'type' for each element of a subset. The type is
 * a vtk-format specified number for a given element type.
 */
template <int TDim>
template <typename T>
void VTKOutput<TDim>::
write_cell_types(VTKFileWriter& File, const T& iterContainer, const SubsetGroup& ssGrp,
                 const int dim, int numElem)
{
	File << VTKFileWriter::normal;
//	write opening tag to indicate that types will be written
	File << "        <DataArray type=\"Int8\" Name=\"types\" format="
		 <<	(m_bBinary ? "\"binary\"": "\"ascii\"") << ">\n";
	if(m_bBinary)
		File << VTKFileWriter::base64_binary << numElem;

//	switch dimension
	if(numElem > 0)
	for(size_t i = 0; i < ssGrp.size(); i++)
	{
		switch(dim)
		{
			case 0: break; // no elements -> nothing to do
			case 1: write_cell_types<RegularEdge>(File, iterContainer, ssGrp[i]);
					write_cell_types<ConstrainingEdge>(File, iterContainer, ssGrp[i]); break;
			case 2: write_cell_types<Triangle>(File, iterContainer, ssGrp[i]);
					write_cell_types<Quadrilateral>(File, iterContainer, ssGrp[i]);
					write_cell_types<ConstrainingTriangle>(File, iterContainer, ssGrp[i]);
					write_cell_types<ConstrainingQuadrilateral>(File, iterContainer, ssGrp[i]); break;
			case 3: write_cell_types<Tetrahedron>(File, iterContainer, ssGrp[i]);
					write_cell_types<Pyramid>(File, iterContainer, ssGrp[i]);
					write_cell_types<Prism>(File, iterContainer, ssGrp[i]);
					write_cell_types<Octahedron>(File, iterContainer, ssGrp[i]);
					write_cell_types<Hexahedron>(File, iterContainer, ssGrp[i]); break;
			default: UG_THROW("VTK::write_cell_types: Dimension " << dim <<
			                        " is not supported.");
		}
	}

//	write closing tag
	File << VTKFileWriter::normal;
	File << "\n        </DataArray>\n";
}

////////////////////////////////////////////////////////////////////////////////
// Subset Indices
////////////////////////////////////////////////////////////////////////////////


template <int TDim>
template <class TElem, typename T>
void VTKOutput<TDim>::
write_cell_subsets(VTKFileWriter& File, const T& iterContainer, const int si, MGSubsetHandler& sh)
{
//	subset
	int subset;

//	get iterators
	typedef typename IteratorProvider<T>::template traits<TElem>::const_iterator const_iterator;
	const_iterator iterBegin = IteratorProvider<T>::template begin<TElem>(iterContainer, si);
	const_iterator iterEnd = IteratorProvider<T>::template end<TElem>(iterContainer, si);

	if(m_bBinary)
		File << VTKFileWriter::base64_binary;
	else
		File << VTKFileWriter::normal;
//	loop all elements, write type for each element to stream
	for( ; iterBegin != iterEnd; ++iterBegin)
	{
		subset = (char)sh.get_subset_index(*iterBegin);
		//iterContainer.get_subset_name(2);

		if(m_bBinary){
			File << (char) subset << ' ';
		}
		else{
			File << subset << ' ';
		}
	}
}

/**
 * This method writes the 'region' for each element of a subset. The region is
 * the ug4 subset index for a given element.
 */
template <int TDim>
template <typename T>
void VTKOutput<TDim>::
write_cell_subsets(VTKFileWriter& File, const T& iterContainer, const SubsetGroup& ssGrp,
                 const int dim, const int numElem, MGSubsetHandler& sh)
{
	File << VTKFileWriter::normal;
//	write opening tag to indicate that types will be written
	File << "        <DataArray type=\"Int8\" Name=\"regions\" format="
		 <<	(m_bBinary ? "\"binary\"": "\"ascii\"") << ">\n";
	if(m_bBinary)
		File << VTKFileWriter::base64_binary << numElem;

//	switch dimension
	if(numElem > 0)
	for(size_t i = 0; i < ssGrp.size(); i++)
	{
		switch(dim)
		{
			case 0: break; // no elements -> nothing to do
			case 1: write_cell_subsets<RegularEdge>(File, iterContainer, ssGrp[i], sh);
					write_cell_subsets<ConstrainingEdge>(File, iterContainer, ssGrp[i], sh); break;
			case 2: write_cell_subsets<Triangle>(File, iterContainer, ssGrp[i], sh);
					write_cell_subsets<Quadrilateral>(File, iterContainer, ssGrp[i], sh);
					write_cell_subsets<ConstrainingTriangle>(File, iterContainer, ssGrp[i], sh);
					write_cell_subsets<ConstrainingQuadrilateral>(File, iterContainer, ssGrp[i], sh); break;
			case 3: write_cell_subsets<Tetrahedron>(File, iterContainer, ssGrp[i], sh);
					write_cell_subsets<Pyramid>(File, iterContainer, ssGrp[i], sh);
					write_cell_subsets<Prism>(File, iterContainer, ssGrp[i], sh);
					write_cell_subsets<Octahedron>(File, iterContainer, ssGrp[i], sh);
					write_cell_subsets<Hexahedron>(File, iterContainer, ssGrp[i], sh); break;
			default: UG_THROW("VTK::write_cell_subsets: Dimension " << dim <<
			                        " is not supported.");
		}
	}

//	write closing tag
	File << VTKFileWriter::normal;
	File << "\n        </DataArray>\n";
}

////////////////////////////////////////////////////////////////////////////////
// Subset Indices
////////////////////////////////////////////////////////////////////////////////


template <int TDim>
template <class TElem, typename T>
void VTKOutput<TDim>::
write_cell_proc_ranks(VTKFileWriter& File, const T& iterContainer, const int si, MGSubsetHandler& sh)
{
//	subset
	int rank = 0;

#ifdef UG_PARALLEL
	rank = pcl::ProcRank();
#endif

//	get iterators
	typedef typename IteratorProvider<T>::template traits<TElem>::const_iterator const_iterator;
	const_iterator iterBegin = IteratorProvider<T>::template begin<TElem>(iterContainer, si);
	const_iterator iterEnd = IteratorProvider<T>::template end<TElem>(iterContainer, si);

	if(m_bBinary)
		File << VTKFileWriter::base64_binary;
	else
		File << VTKFileWriter::normal;
//	loop all elements, write type for each element to stream
	for( ; iterBegin != iterEnd; ++iterBegin)
	{

		if(m_bBinary){
			File << (char) rank << ' ';
		}
		else{
			File << rank << ' ';
		}
	}
}

/**
 * This method writes the proc ranks for each cell.
 */
template <int TDim>
template <typename T>
void VTKOutput<TDim>::
write_cell_proc_ranks(VTKFileWriter& File, const T& iterContainer, const SubsetGroup& ssGrp,
                 const int dim, const int numElem, MGSubsetHandler& sh)
{
	File << VTKFileWriter::normal;
//	write opening tag to indicate that types will be written
	File << "        <DataArray type=\"Int8\" Name=\"proc_ranks\" format="
		 <<	(m_bBinary ? "\"binary\"": "\"ascii\"") << ">\n";
	if(m_bBinary)
		File << VTKFileWriter::base64_binary << numElem;

//	switch dimension
	if(numElem > 0)
	for(size_t i = 0; i < ssGrp.size(); i++)
	{
		switch(dim)
		{
			case 0: break; // no elements -> nothing to do
			case 1: write_cell_proc_ranks<RegularEdge>(File, iterContainer, ssGrp[i], sh);
					write_cell_proc_ranks<ConstrainingEdge>(File, iterContainer, ssGrp[i], sh); break;
			case 2: write_cell_proc_ranks<Triangle>(File, iterContainer, ssGrp[i], sh);
					write_cell_proc_ranks<Quadrilateral>(File, iterContainer, ssGrp[i], sh);
					write_cell_proc_ranks<ConstrainingTriangle>(File, iterContainer, ssGrp[i], sh);
					write_cell_proc_ranks<ConstrainingQuadrilateral>(File, iterContainer, ssGrp[i], sh); break;
			case 3: write_cell_proc_ranks<Tetrahedron>(File, iterContainer, ssGrp[i], sh);
					write_cell_proc_ranks<Pyramid>(File, iterContainer, ssGrp[i], sh);
					write_cell_proc_ranks<Prism>(File, iterContainer, ssGrp[i], sh);
					write_cell_proc_ranks<Octahedron>(File, iterContainer, ssGrp[i], sh);
					write_cell_proc_ranks<Hexahedron>(File, iterContainer, ssGrp[i], sh); break;
			default: UG_THROW("VTK::write_cell_proc_ranks: Dimension " << dim <<
			                        " is not supported.");
		}
	}

//	write closing tag
	File << VTKFileWriter::normal;
	File << "\n        </DataArray>\n";
}

////////////////////////////////////////////////////////////////////////////////
// Nodal Value
////////////////////////////////////////////////////////////////////////////////

template <int TDim>
template <typename TElem, typename TFunction, typename TData>
void VTKOutput<TDim>::
write_nodal_data_elementwise(VTKFileWriter& File, TFunction& u, number time,
                             SmartPtr<UserData<TData, TDim> > spData,
                             Grid& grid, const int si)
{
//	get reference element
	typedef typename reference_element_traits<TElem>::reference_element_type
																ref_elem_type;
	static const ref_elem_type& refElem = Provider<ref_elem_type>::get();
	static const size_t numCo = ref_elem_type::numCorners;

	if(m_bBinary)
		File << VTKFileWriter::base64_binary;
	else
		File << VTKFileWriter::normal;

//	get iterators
	typedef typename IteratorProvider<TFunction>::template traits<TElem>::const_iterator const_iterator;
	const_iterator iterBegin = IteratorProvider<TFunction>::template begin<TElem>(u, si);
	const_iterator iterEnd = IteratorProvider<TFunction>::template end<TElem>(u, si);

	std::vector<MathVector<TDim> > vCorner(numCo);

	TData vValue[numCo];
	bool bNeedFct = spData->requires_grid_fct();

//	loop all elements
	for( ; iterBegin != iterEnd; ++iterBegin)
	{
	//	get element
		TElem *elem = *iterBegin;

	//	update the reference mapping for the corners
		CollectCornerCoordinates(vCorner, *elem, u.domain()->position_accessor(), true);

	//	get subset
		int theSI = si;
		if(si < 0) theSI = u.domain()->subset_handler()->get_subset_index(elem);

	//	get local solution if needed
		if(bNeedFct)
		{
		//	create storage
			LocalIndices ind;
			LocalVector locU;

		// 	get global indices
			u.indices(elem, ind);

		// 	adapt local algebra
			locU.resize(ind);

		// 	read local values of u
			GetLocalVector(locU, u);

		//	compute data
			try{
				(*spData)(vValue, &vCorner[0], time, theSI, elem,
							&vCorner[0], refElem.corners(), numCo, &locU, NULL);
			}
			UG_CATCH_THROW("VTK::write_nodal_data_elementwise: Cannot evaluate data.");
		}
		else
		{
		//	compute data
			try{
				(*spData)(vValue, &vCorner[0], time, theSI, numCo);
			}
			UG_CATCH_THROW("VTK::write_nodal_data_elementwise: Cannot evaluate data.");
		}

	//	loop vertices of element
		for(size_t co = 0; co < numCo; ++co)
		{
		//	get vertex of element
			Vertex* v = GetVertex(elem, co);

		//	if vertex has been handled before, skip
			if(grid.is_marked(v)) continue;

		//	mark as used
			grid.mark(v);

		//	loop all components
			write_item_to_file(File, vValue[co]);
		}
	}

}

/**
 * Writes the nodal data to the vtu file.
 */
template <int TDim>
template <typename TFunction, typename TData>
void VTKOutput<TDim>::
write_nodal_data(VTKFileWriter& File, TFunction& u, number time,
                 SmartPtr<UserData<TData, TDim> > spData,
                 const int numCmp,
                 const std::string& name,
                 Grid& grid, const SubsetGroup& ssGrp, const int dim, const int numVert)
{
	spData->set_function_pattern(u.function_pattern());

//	check that nodal data is possible
	if(!spData->continuous())
		UG_THROW("VTK: data with name '"<<name<<"' is scheduled for nodal output,"
				" but the data is not continuous. Cannot write it as nodal data.")

//	write opening tag
	File << VTKFileWriter::normal;
	File << "        <DataArray type=\"Float32\" Name=\""<<name<<"\" "
	"NumberOfComponents=\""<<numCmp<<"\" format="
		 <<	(m_bBinary ? "\"binary\"": "\"ascii\"") << ">\n";

	int n = sizeof(float) * numVert * numCmp;
	if(m_bBinary)
		File << VTKFileWriter::base64_binary << n;

//	start marking of grid
	grid.begin_marking();

//	switch dimension
	for(size_t i = 0; i < ssGrp.size(); i++)
	switch(dim)
	{
		case 1:	write_nodal_data_elementwise<RegularEdge,TFunction,TData>(File, u, time, spData, grid, ssGrp[i]);
				write_nodal_data_elementwise<ConstrainingEdge,TFunction,TData>(File, u, time, spData, grid, ssGrp[i]); break;
		case 2:	write_nodal_data_elementwise<Triangle,TFunction,TData>(File, u, time, spData, grid, ssGrp[i]);
				write_nodal_data_elementwise<Quadrilateral,TFunction,TData>(File, u, time, spData, grid, ssGrp[i]);
				write_nodal_data_elementwise<ConstrainingTriangle,TFunction,TData>(File, u, time, spData, grid, ssGrp[i]);
				write_nodal_data_elementwise<ConstrainingQuadrilateral,TFunction,TData>(File, u, time, spData, grid, ssGrp[i]); break;
		case 3:	write_nodal_data_elementwise<Tetrahedron,TFunction,TData>(File, u, time, spData, grid, ssGrp[i]);
				write_nodal_data_elementwise<Pyramid,TFunction,TData>(File, u, time, spData, grid, ssGrp[i]);
				write_nodal_data_elementwise<Prism,TFunction,TData>(File, u, time, spData, grid, ssGrp[i]);
				write_nodal_data_elementwise<Octahedron,TFunction,TData>(File, u, time, spData, grid, ssGrp[i]);
				write_nodal_data_elementwise<Hexahedron,TFunction,TData>(File, u, time, spData, grid, ssGrp[i]); break;
		default: UG_THROW("VTK::write_nodal_data: Dimension " << dim <<
		                        " is not supported.");
	}

//	end marking
	grid.end_marking();

//	write closing tag
	File << VTKFileWriter::normal;
	File << "\n        </DataArray>\n";
};


template <int TDim>
template <typename TElem, typename TFunction>
void VTKOutput<TDim>::
write_nodal_values_elementwise(VTKFileWriter& File, TFunction& u,
                               const std::vector<size_t>& vFct,
                               Grid& grid, const int si)
{
//	get reference element
	typedef typename reference_element_traits<TElem>::reference_element_type
																ref_elem_type;
	if(m_bBinary)
		File << VTKFileWriter::base64_binary;
	else
		File << VTKFileWriter::normal;

//	index vector
	std::vector<DoFIndex> vMultInd;

//	get iterators
	typedef typename IteratorProvider<TFunction>::template traits<TElem>::const_iterator const_iterator;
	const_iterator iterBegin = IteratorProvider<TFunction>::template begin<TElem>(u, si);
	const_iterator iterEnd = IteratorProvider<TFunction>::template end<TElem>(u, si);

//	loop all elements
	for( ; iterBegin != iterEnd; ++iterBegin)
	{
	//	get element
		TElem *elem = *iterBegin;

	//	loop vertices of element
		for(size_t co = 0; co < (size_t) ref_elem_type::numCorners; ++co)
		{
		//	get vertex of element
			Vertex* v = GetVertex(elem, co);

		//	if vertex has been handled before, skip
			if(grid.is_marked(v)) continue;

		//	mark as used
			grid.mark(v);

		//	loop all components
			for(size_t i = 0; i < vFct.size(); ++i)
			{
			//	get multi index of vertex for the function
				if(u.inner_dof_indices(v, vFct[i], vMultInd) != 1)
					UG_THROW("VTK:write_nodal_values_elementwise: "
							"The function component "<<vFct[i]<<" has "<<
							vMultInd.size()<<" DoFs in a vertex of an element of subset "
							<<si<<". To write a component to vtk, exactly "
							"one DoF must be given in any vertex.");

			//	flush stream
				write_item_to_file(File, DoFRef(u, vMultInd[0]));
			}



		//	fill with zeros up to 3d if vector type
			if(vFct.size() != 1) {
				for(size_t i = vFct.size(); i < 3; ++i) {
					write_item_to_file(File, 0.f);
				}
			}
		}
	}

}

template <int TDim>
template <typename TFunction>
void VTKOutput<TDim>::
write_nodal_values(VTKFileWriter& File, TFunction& u,
                   const std::vector<size_t>& vFct, const std::string& name,
                   Grid& grid, const SubsetGroup& ssGrp, const int dim, const int numVert)
{
	File << VTKFileWriter::normal;
//	write opening tag
	File << "        <DataArray type=\"Float32\" Name=\""<<name<<"\" "
	"NumberOfComponents=\""<<(vFct.size() == 1 ? 1 : 3)<<"\" format="
		 <<	(m_bBinary ? "\"binary\"": "\"ascii\"") << ">\n";

	int n = sizeof(float) * numVert * (vFct.size() == 1 ? 1 : 3);
	if(m_bBinary)
		File << VTKFileWriter::base64_binary << n;

//	start marking of grid
	grid.begin_marking();

//	switch dimension
	for(size_t i = 0; i < ssGrp.size(); i++)
	switch(dim)
	{
		case 0:	write_nodal_values_elementwise<Vertex>(File, u, vFct, grid, ssGrp[i]); break;
		case 1:	write_nodal_values_elementwise<RegularEdge>(File, u, vFct, grid, ssGrp[i]);
				write_nodal_values_elementwise<ConstrainingEdge>(File, u, vFct, grid, ssGrp[i]); break;
		case 2:	write_nodal_values_elementwise<Triangle>(File, u, vFct, grid, ssGrp[i]);
				write_nodal_values_elementwise<Quadrilateral>(File, u, vFct, grid, ssGrp[i]);
				write_nodal_values_elementwise<ConstrainingTriangle>(File, u, vFct, grid, ssGrp[i]);
				write_nodal_values_elementwise<ConstrainingQuadrilateral>(File, u, vFct, grid, ssGrp[i]); break;
		case 3:	write_nodal_values_elementwise<Tetrahedron>(File, u, vFct, grid, ssGrp[i]);
				write_nodal_values_elementwise<Pyramid>(File, u, vFct, grid, ssGrp[i]);
				write_nodal_values_elementwise<Prism>(File, u, vFct, grid, ssGrp[i]);
				write_nodal_values_elementwise<Octahedron>(File, u, vFct, grid, ssGrp[i]);
				write_nodal_values_elementwise<Hexahedron>(File, u, vFct, grid, ssGrp[i]); break;
		default: UG_THROW("VTK::write_nodal_values: Dimension " << dim <<
		                        " is not supported.");
	}

//	end marking
	grid.end_marking();

//	write closing tag
	File << VTKFileWriter::normal;
	File << "\n        </DataArray>\n";
};

template <int TDim>
template <typename TFunction>
void VTKOutput<TDim>::
write_nodal_values_piece(VTKFileWriter& File, TFunction& u, number time, Grid& grid,
                         const SubsetGroup& ssGrp, const int dim, const int numVert)
{
	if(!m_vSymbFct.empty()){
		for(std::map<std::string, std::vector<std::string> >::const_iterator iter =
				m_vSymbFct.begin(); iter != m_vSymbFct.end(); ++iter){
			const std::vector<std::string>& symbNames = (*iter).second;
			const std::string& vtkName = (*iter).first;

			bool bContinuous = true;
			for(size_t i = 0; i < symbNames.size(); ++i){
				size_t fct = u.fct_id_by_name(symbNames[i].c_str());
				if(!LocalFiniteElementProvider::continuous(u.local_finite_element_id(fct))){
					bContinuous = false; break;
				}
			}

			if(bContinuous){
				m_vSymbFctNodal[vtkName] = symbNames;
			}
		}
	}

//	check if something to do
	if(m_vSymbFctNodal.empty() && m_vScalarNodalData.empty() && m_vVectorNodalData.empty())
		return;

//	write opening tag to indicate point data
	File << VTKFileWriter::normal;
	File << "      <PointData>\n";

//	loop all selected symbolic names
	for(std::map<std::string, std::vector<std::string> >::const_iterator iter =
			m_vSymbFctNodal.begin(); iter != m_vSymbFctNodal.end(); ++iter)
	{
	//	get symb function
		const std::vector<std::string>& symbNames = (*iter).second;
		const std::string& vtkName = (*iter).first;

	//	create function group
		std::vector<size_t> fctGrp(symbNames.size());
		for(size_t i = 0; i < symbNames.size(); ++i)
			fctGrp[i] = u.fct_id_by_name(symbNames[i].c_str());

	//	check that all functions are contained in subset
		bool bContained = true;
		for(size_t i = 0; i < fctGrp.size(); ++i)
			for (size_t j = 0; j < ssGrp.size(); ++j)
				if(!u.is_def_in_subset(fctGrp[i], ssGrp[j]))
					bContained = false;

		if(!bContained) continue;

	//	write scalar value of this function
		try{
			write_nodal_values(File, u, fctGrp, vtkName, grid, ssGrp, dim, numVert);
		}
		UG_CATCH_THROW("VTK::write_piece: Failed write nodal Values.");
	}

//	loop all scalar data
	for(ScalarDataIterator iter = m_vScalarNodalData.begin();
			iter != m_vScalarNodalData.end(); ++iter)
	{
	//	get symb function
		SmartPtr<UserData<number,TDim> > spData = (*iter).second;
		const std::string& vtkName = (*iter).first;

	//	write scalar value of this data
		try{
			write_nodal_data<TFunction,number>
							(File, u, time, spData, 1, vtkName, grid, ssGrp, dim, numVert);
		}
		UG_CATCH_THROW("VTK::write_piece: Failed write nodal scalar Data.");
	}

//	loop all vector data
	for(VectorDataIterator	iter = m_vVectorNodalData.begin();
			iter != m_vVectorNodalData.end(); ++iter)
	{
	//	get symb function
		SmartPtr<UserData<MathVector<TDim>,TDim> > spData = (*iter).second;
		const std::string& vtkName = (*iter).first;

	//	write scalar value of this data
		try{
			write_nodal_data<TFunction,MathVector<TDim> >
							(File, u, time, spData, 3, vtkName, grid, ssGrp, dim, numVert);
		}
		UG_CATCH_THROW("VTK::write_piece: Failed write nodal vector Data.");
	}

//	loop all vector data
	for(MatrixDataIterator	iter = m_vMatrixNodalData.begin();
			iter != m_vMatrixNodalData.end(); ++iter)
	{
	//	get symb function
		SmartPtr<UserData<MathMatrix<TDim,TDim>,TDim> > spData = (*iter).second;
		const std::string& vtkName = (*iter).first;

	//	write scalar value of this data
		try{
			write_nodal_data<TFunction,MathMatrix<TDim,TDim> >
							(File, u, time, spData, 9, vtkName, grid, ssGrp, dim, numVert);
		}
		UG_CATCH_THROW("VTK::write_piece: Failed write nodal matrix Data.");
	}

//	write closing tag
	File << VTKFileWriter::normal;
	File << "      </PointData>\n";
}


////////////////////////////////////////////////////////////////////////////////
// Cell Value
////////////////////////////////////////////////////////////////////////////////

template <int TDim>
template <typename TElem, typename TFunction, typename TData>
void VTKOutput<TDim>::
write_cell_data_elementwise(VTKFileWriter& File, TFunction& u, number time,
                             SmartPtr<UserData<TData, TDim> > spData,
                             Grid& grid, const int si)
{
//	get reference element
	typedef typename reference_element_traits<TElem>::reference_element_type
																ref_elem_type;
	static const int refDim = reference_element_traits<TElem>::dim;
	static const ref_elem_type& refElem = Provider<ref_elem_type>::get();
	static const size_t numCo = ref_elem_type::numCorners;

	if(m_bBinary)
		File << VTKFileWriter::base64_binary;
	else
		File << VTKFileWriter::normal;

//	get iterators
	typedef typename IteratorProvider<TFunction>::template traits<TElem>::const_iterator const_iterator;
	const_iterator iterBegin = IteratorProvider<TFunction>::template begin<TElem>(u, si);
	const_iterator iterEnd = IteratorProvider<TFunction>::template end<TElem>(u, si);

	MathVector<refDim> localIP;
	MathVector<TDim> globIP;
	AveragePositions(localIP, refElem.corners(), numCo);

	std::vector<MathVector<TDim> > vCorner(numCo);

	TData value;
	bool bNeedFct = spData->requires_grid_fct();

//	loop all elements
	for( ; iterBegin != iterEnd; ++iterBegin)
	{
	//	get element
		TElem *elem = *iterBegin;

	//	update the reference mapping for the corners
		CollectCornerCoordinates(vCorner, *elem, u.domain()->position_accessor(), true);

	//	compute global integration points
		AveragePositions(globIP, &vCorner[0], numCo);

	//	get subset
		int theSI = si;
		if(si < 0) theSI = u.domain()->subset_handler()->get_subset_index(elem);

	//	get local solution if needed
		if(bNeedFct)
		{
		//	create storage
			LocalIndices ind;
			LocalVector locU;

		// 	get global indices
			u.indices(elem, ind);

		// 	adapt local algebra
			locU.resize(ind);

		// 	read local values of u
			GetLocalVector(locU, u);

		//	compute data
			try{
				(*spData)(value, globIP, time, theSI, elem,
							&vCorner[0], localIP, &locU);
			}
			UG_CATCH_THROW("VTK::write_cell_data_elementwise: Cannot evaluate data.");
		}
		else
		{
		//	compute data
			try{
				(*spData)(value, globIP, time, theSI);
			}
			UG_CATCH_THROW("VTK::write_cell_data_elementwise: Cannot evaluate data.");
		}

	//	handle octahedral case separately by splitting into top and bottom pyramid with the same values
		if(ref_elem_type::REFERENCE_OBJECT_ID != ROID_OCTAHEDRON)
		{
			write_item_to_file(File, value);
		}
		else
		{
			write_item_to_file(File, value); // top pyramid
			write_item_to_file(File, value); // bottom pyramid
		}
	}

}


template <int TDim>
template <typename TFunction, typename TData>
void VTKOutput<TDim>::
write_cell_data(VTKFileWriter& File, TFunction& u, number time,
                 SmartPtr<UserData<TData, TDim> > spData,
                 const int numCmp,
                 const std::string& name,
                 Grid& grid, const SubsetGroup& ssGrp, const int dim, const int numElem)
{
	spData->set_function_pattern(u.function_pattern());

//	write opening tag
	File << VTKFileWriter::normal;
	File << "        <DataArray type=\"Float32\" Name=\""<<name<<"\" "
	"NumberOfComponents=\""<<numCmp<<"\" format="
		 <<	(m_bBinary ? "\"binary\"": "\"ascii\"") << ">\n";

	int n = sizeof(float) * numElem * numCmp;
	if(m_bBinary)
		File << VTKFileWriter::base64_binary << n;

//	switch dimension
	for(size_t i = 0; i < ssGrp.size(); i++)
	switch(dim)
	{
		case 1:	write_cell_data_elementwise<RegularEdge,TFunction,TData>(File, u, time, spData, grid, ssGrp[i]);
				write_cell_data_elementwise<ConstrainingEdge,TFunction,TData>(File, u, time, spData, grid, ssGrp[i]); break;
		case 2:	write_cell_data_elementwise<Triangle,TFunction,TData>(File, u, time, spData, grid, ssGrp[i]);
				write_cell_data_elementwise<Quadrilateral,TFunction,TData>(File, u, time, spData, grid, ssGrp[i]);
				write_cell_data_elementwise<ConstrainingTriangle,TFunction,TData>(File, u, time, spData, grid, ssGrp[i]);
				write_cell_data_elementwise<ConstrainingQuadrilateral,TFunction,TData>(File, u, time, spData, grid, ssGrp[i]); break;
		case 3:	write_cell_data_elementwise<Tetrahedron,TFunction,TData>(File, u, time, spData, grid, ssGrp[i]);
				write_cell_data_elementwise<Pyramid,TFunction,TData>(File, u, time, spData, grid, ssGrp[i]);
				write_cell_data_elementwise<Prism,TFunction,TData>(File, u, time, spData, grid, ssGrp[i]);
				write_cell_data_elementwise<Octahedron,TFunction,TData>(File, u, time, spData, grid, ssGrp[i]);
				write_cell_data_elementwise<Hexahedron,TFunction,TData>(File, u, time, spData, grid, ssGrp[i]); break;
		default: UG_THROW("VTK::write_cell_data: Dimension " << dim <<
		                        " is not supported.");
	}

//	write closing tag
	File << VTKFileWriter::normal;
	File << "\n        </DataArray>\n";
};


template <int TDim>
template <typename TElem, typename TFunction>
void VTKOutput<TDim>::
write_cell_values_elementwise(VTKFileWriter& File, TFunction& u,
                               const std::vector<size_t>& vFct,
                               Grid& grid, const int si)
{
//	get reference element
	typedef typename reference_element_traits<TElem>::reference_element_type
																ref_elem_type;

	static const ref_elem_type& refElem = Provider<ref_elem_type>::get();
	static const ReferenceObjectID roid = ref_elem_type::REFERENCE_OBJECT_ID;
	static const int dim = ref_elem_type::dim;
	static const size_t numCo = ref_elem_type::numCorners;

//	index vector
	std::vector<DoFIndex> vMultInd;

//	get iterators
	typedef typename IteratorProvider<TFunction>::template traits<TElem>::const_iterator const_iterator;
	const_iterator iterBegin = IteratorProvider<TFunction>::template begin<TElem>(u, si);
	const_iterator iterEnd = IteratorProvider<TFunction>::template end<TElem>(u, si);

	if(m_bBinary)
		File << VTKFileWriter::base64_binary;
	else
		File << VTKFileWriter::normal;

	MathVector<dim> localIP;
	AveragePositions(localIP, refElem.corners(), numCo);

//	request for trial space
	try{

	// avoid passing this code, since LocalShapeFunctionSet might not exist for RefElem-Type
	if(iterBegin == iterEnd) return;

	std::vector<std::vector<number> > vvShape(vFct.size());
	std::vector<size_t> vNsh(vFct.size());

	for(size_t f = 0; f < vFct.size(); ++f)
	{
		const LFEID lfeID = u.local_finite_element_id(vFct[f]);
		const LocalShapeFunctionSet<dim>& lsfs
			 = LocalFiniteElementProvider::get<dim>(roid, lfeID);

		vNsh[f] = lsfs.num_sh();
		lsfs.shapes(vvShape[f], localIP);
	}

//	loop all elements
	for( ; iterBegin != iterEnd; ++iterBegin)
	{
	//	get element
		TElem *elem = *iterBegin;

	//	loop all components
		for(size_t f = 0; f < vFct.size(); ++f)
		{
		//	get multi index of vertex for the function
			if(u.dof_indices(elem, vFct[f], vMultInd) != vNsh[f])
				UG_THROW("VTK:write_cell_values_elementwise: "
						"Number of shape functions for component "<<vFct[f]<<
						" does not match number of DoFs");

			number ipVal = 0.0;
			for(size_t sh = 0; sh < vNsh[f]; ++sh)
			{
				ipVal += DoFRef(u, vMultInd[sh]) * vvShape[f][sh];
			}

		//	flush stream
			write_item_to_file(File, ipVal);

		//	flush stream a second time in case of an octahedron as we handle it by splitting into a top and bottom pyramid
			if(roid == ROID_OCTAHEDRON)
				write_item_to_file(File, ipVal); // bottom pyramid
		}

	//	fill with zeros up to 3d if vector type
		if(vFct.size() != 1){
			for(size_t i = vFct.size(); i < 3; ++i) {
				write_item_to_file(File, 0.f);

			//	flush stream a second time in case of an octahedron as we handle it by splitting into a top and bottom pyramid
				if(roid == ROID_OCTAHEDRON)
					write_item_to_file(File, 0.f); // bottom pyramid
			}
		}
	}
	}
	UG_CATCH_THROW("VTK: Could not find Shape function Set.");

}

template <int TDim>
template <typename TFunction>
void VTKOutput<TDim>::
write_cell_values(VTKFileWriter& File, TFunction& u,
                   const std::vector<size_t>& vFct, const std::string& name,
                   Grid& grid, const SubsetGroup& ssGrp, const int dim, const int numElem)
{
//	write opening tag
	File << VTKFileWriter::normal;
	File << "        <DataArray type=\"Float32\" Name=\""<<name<<"\" "
	"NumberOfComponents=\""<<(vFct.size() == 1 ? 1 : 3)<<"\" format="
		 <<	(m_bBinary ? "\"binary\"": "\"ascii\"") << ">\n";

	int n = sizeof(float) * numElem * (vFct.size() == 1 ? 1 : 3);
	if(m_bBinary)
		File << VTKFileWriter::base64_binary << n;

//	switch dimension
	for(size_t i = 0; i < ssGrp.size(); i++)
	switch(dim)
	{
		case 1:	write_cell_values_elementwise<RegularEdge>(File, u, vFct, grid, ssGrp[i]);
				write_cell_values_elementwise<ConstrainingEdge>(File, u, vFct, grid, ssGrp[i]); break;
		case 2:	write_cell_values_elementwise<Triangle>(File, u, vFct, grid, ssGrp[i]);
				write_cell_values_elementwise<Quadrilateral>(File, u, vFct, grid, ssGrp[i]);
				write_cell_values_elementwise<ConstrainingTriangle>(File, u, vFct, grid, ssGrp[i]);
				write_cell_values_elementwise<ConstrainingQuadrilateral>(File, u, vFct, grid, ssGrp[i]); break;
		case 3:	write_cell_values_elementwise<Tetrahedron>(File, u, vFct, grid, ssGrp[i]);
				write_cell_values_elementwise<Pyramid>(File, u, vFct, grid, ssGrp[i]);
				write_cell_values_elementwise<Prism>(File, u, vFct, grid, ssGrp[i]);
				write_cell_values_elementwise<Octahedron>(File, u, vFct, grid, ssGrp[i]);
				write_cell_values_elementwise<Hexahedron>(File, u, vFct, grid, ssGrp[i]); break;
		default: UG_THROW("VTK::write_cell_values: Dimension " << dim <<
		                        " is not supported.");
	}

//	write closing tag
	File << VTKFileWriter::normal;
	File << "\n        </DataArray>\n";
};

template <int TDim>
template <typename TFunction>
void VTKOutput<TDim>::
write_cell_values_piece(VTKFileWriter& File, TFunction& u, number time, Grid& grid,
                         const SubsetGroup& ssGrp, const int dim, const int numElem)
{
	if(!m_vSymbFct.empty()){
		for(std::map<std::string, std::vector<std::string> >::const_iterator iter =
				m_vSymbFct.begin(); iter != m_vSymbFct.end(); ++iter){
			const std::vector<std::string>& symbNames = (*iter).second;
			const std::string& vtkName = (*iter).first;

			bool bContinuous = true;
			for(size_t i = 0; i < symbNames.size(); ++i){
				size_t fct = u.fct_id_by_name(symbNames[i].c_str());
				if(!LocalFiniteElementProvider::continuous(u.local_finite_element_id(fct))){
					bContinuous = false; break;
				}
			}

			if(!bContinuous){
				m_vSymbFctElem[vtkName] = symbNames;
			}
		}
	}

//	check if something to do
	if(m_vSymbFctElem.empty() && m_vScalarElemData.empty() && m_vVectorElemData.empty())
		return;

//	write opening tag to indicate point data
	File << VTKFileWriter::normal;
	File << "      <CellData>\n";

//	loop all selected symbolic names
	for(ComponentsIterator iter = m_vSymbFctElem.begin();
			iter != m_vSymbFctElem.end(); ++iter)
	{
	//	get symb function
		const std::vector<std::string>& symbNames = (*iter).second;
		const std::string& vtkName = (*iter).first;

	//	create function group
		std::vector<size_t> fctGrp(symbNames.size());
		for(size_t i = 0; i < symbNames.size(); ++i)
			fctGrp[i] = u.fct_id_by_name(symbNames[i].c_str());

	//	check that all functions are contained in subset
		bool bContained = true;
		for(size_t i = 0; i < fctGrp.size(); ++i)
			for(size_t j = 0; j < ssGrp.size(); ++j)
				if(!u.is_def_in_subset(fctGrp[i], ssGrp[j]))
					bContained = false;

		if(!bContained) continue;

	//	write scalar value of this function
		try{
			write_cell_values(File, u, fctGrp, vtkName, grid, ssGrp, dim, numElem);
		}
		UG_CATCH_THROW("VTK::write_piece: Failed write cell Values.");
	}

//	loop all scalar data
	for(ScalarDataIterator iter = m_vScalarElemData.begin();
			iter != m_vScalarElemData.end(); ++iter)
	{
	//	get symb function
		SmartPtr<UserData<number,TDim> > spData = (*iter).second;
		const std::string& vtkName = (*iter).first;

	//	write scalar value of this data
		try{
			write_cell_data<TFunction,number>
							(File, u, time, spData, 1, vtkName, grid, ssGrp, dim, numElem);
		}
		UG_CATCH_THROW("VTK::write_piece: Failed to write cell scalar Data.");
	}

//	loop all vector data
	for(VectorDataIterator iter = m_vVectorElemData.begin();
			iter != m_vVectorElemData.end(); ++iter)
	{
	//	get symb function
		SmartPtr<UserData<MathVector<TDim>,TDim> > spData = (*iter).second;
		const std::string& vtkName =  (*iter).first;

	//	write scalar value of this data
		try{
			write_cell_data<TFunction,MathVector<TDim> >
							(File, u, time, spData, 3, vtkName, grid, ssGrp, dim, numElem);
		}
		UG_CATCH_THROW("VTK::write_piece: Failed to write cell vector Data.");
	}

//	loop all vector data
	for(MatrixDataIterator iter = m_vMatrixElemData.begin();
			iter != m_vMatrixElemData.end(); ++iter)
	{
	//	get symb function
		SmartPtr<UserData<MathMatrix<TDim,TDim>,TDim> > spData = (*iter).second;
		const std::string& vtkName =  (*iter).first;

	//	write scalar value of this data
		try{
			write_cell_data<TFunction,MathMatrix<TDim,TDim> >
							(File, u, time, spData, 9, vtkName, grid, ssGrp, dim, numElem);
		}
		UG_CATCH_THROW("VTK::write_piece: Failed to write cell matrix Data.");
	}

//	write closing tag
	File << VTKFileWriter::normal;
	File << "      </CellData>\n";
}

////////////////////////////////////////////////////////////////////////////////
// Grouping Files
////////////////////////////////////////////////////////////////////////////////

template <int TDim>
template <typename TFunction>
void VTKOutput<TDim>::
write_pvtu(TFunction& u, const std::string& filename,
           int si, int step, number time)
{
#ifdef UG_PARALLEL
//	File pointer
	FILE* file;

//	file name
	std::string name;

//	get and check number of procs (only for numProcs > 1 we write the pvtu)
	int numProcs = pcl::NumProcs();
	if(numProcs == 1) return;

//	check if this proc is output proc
	bool isOutputProc = GetLogAssistant().is_output_process();

//	max subset
	int maxSi = u.num_subsets() - 1;

//	only the master process writes this file
	if (isOutputProc)
	{
	//	get name for *.pvtu file
		pvtu_filename(name, filename, si, maxSi, step);

	//	open file
		file = fopen(name.c_str(), "w");
		if (file == NULL)
			UG_THROW("VTKOutput: Cannot print to file.");

	//	Write to file
		fprintf(file, "<?xml version=\"1.0\"?>\n");

	//	Write comment
		write_comment_printf(file);
	
		fprintf(file, "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\">\n");
		fprintf(file, "  <Time timestep=\"%.17g\"/>\n", time);
		fprintf(file, "  <PUnstructuredGrid GhostLevel=\"0\">\n");
		fprintf(file, "    <PPoints>\n");
		fprintf(file, "      <PDataArray type=\"Float32\" NumberOfComponents=\"3\"/>\n");
		fprintf(file, "    </PPoints>\n");

	// 	Node Data
		if(!m_vSymbFctNodal.empty() || !m_vScalarNodalData.empty() || !m_vVectorNodalData.empty())
		{
			fprintf(file, "    <PPointData>\n");
//			for(size_t sym = 0; sym < m_vSymbFctNodal.size(); ++sym)
			for(ComponentsIterator iter = m_vSymbFctNodal.begin();
					iter != m_vSymbFctNodal.end(); ++iter)
			{
			//	get symb function
				const std::vector<std::string>& symbNames = (*iter).second;
				const std::string& vtkName = (*iter).first;

			//	create function group
				std::vector<size_t> fctGrp(symbNames.size());
				for(size_t i = 0; i < symbNames.size(); ++i)
					fctGrp[i] = u.fct_id_by_name(symbNames[i].c_str());

			//	check that all functions are contained in subset
				bool bContained = true;
				for(size_t i = 0; i < fctGrp.size(); ++i)
					if(!u.is_def_in_subset(fctGrp[i], si))
						bContained = false;

				if(!bContained) continue;

				fprintf(file, "      <PDataArray type=\"Float32\" Name=\"%s\" "
							  "NumberOfComponents=\"%d\"/>\n",
							  vtkName.c_str(), (fctGrp.size() == 1 ? 1 : 3));
			}

		//	loop all scalar data
			for(ScalarDataIterator iter = m_vScalarNodalData.begin();
					iter != m_vScalarNodalData.end(); ++iter)
			{
			//	get symb function
				const std::string& vtkName = (*iter).first;

				fprintf(file, "      <PDataArray type=\"Float32\" Name=\"%s\" "
							  "NumberOfComponents=\"%d\"/>\n",
							  vtkName.c_str(), 1);
			}

		//	loop all vector data
			for(VectorDataIterator iter = m_vVectorNodalData.begin();
					iter != m_vVectorNodalData.end(); ++iter)
			{
			//	get symb function
				const std::string& vtkName = (*iter).first;

				fprintf(file, "      <PDataArray type=\"Float32\" Name=\"%s\" "
							  "NumberOfComponents=\"%d\"/>\n",
							  vtkName.c_str(), 3);
			}

		//	loop all matrix data
			for(MatrixDataIterator iter = m_vMatrixNodalData.begin();
					iter != m_vMatrixNodalData.end(); ++iter)
			{
			//	get symb function
				const std::string& vtkName = (*iter).first;

				fprintf(file, "      <PDataArray type=\"Float32\" Name=\"%s\" "
							  "NumberOfComponents=\"%d\"/>\n",
							  vtkName.c_str(), 9);
			}
			fprintf(file, "    </PPointData>\n");
		}

	// 	Elem Data
		if(!m_vSymbFctElem.empty() || !m_vScalarElemData.empty() || !m_vVectorElemData.empty() || m_bWriteSubsetIndices || m_bWriteProcRanks) //TODO: cleanup
		{
			fprintf(file, "    <PCellData>\n");
//			for(size_t sym = 0; sym < m_vSymbFctElem.size(); ++sym)
			for(ComponentsIterator iter = m_vSymbFctElem.begin();
					iter != m_vSymbFctElem.end(); ++iter)
			{
			//	get symb function
				const std::vector<std::string>& symbNames = (*iter).second;//m_vSymbFctElem[sym].first;
				const std::string& vtkName = (*iter).first;//m_vSymbFctElem[sym].second;

			//	create function group
				std::vector<size_t> fctGrp(symbNames.size());
				for(size_t i = 0; i < symbNames.size(); ++i)
					fctGrp[i] = u.fct_id_by_name(symbNames[i].c_str());

			//	check that all functions are contained in subset
				bool bContained = true;
				for(size_t i = 0; i < fctGrp.size(); ++i)
					if(!u.is_def_in_subset(fctGrp[i], si))
						bContained = false;

				if(!bContained) continue;

				fprintf(file, "      <PDataArray type=\"Float32\" Name=\"%s\" "
							  "NumberOfComponents=\"%d\"/>\n",
							  vtkName.c_str(), (fctGrp.size() == 1 ? 1 : 3));
			}

		//	loop all scalar data
			for(ScalarDataIterator iter = m_vScalarElemData.begin();
					iter != m_vScalarElemData.end(); ++iter)
			{
			//	get symb function
				const std::string& vtkName = (*iter).first;

				fprintf(file, "      <PDataArray type=\"Float32\" Name=\"%s\" "
							  "NumberOfComponents=\"%d\"/>\n",
							  vtkName.c_str(), 1);
			}

 //TODO: cleanup!!
			if(m_bWriteSubsetIndices){
				fprintf(file, "      <DataArray type=\"Int8\" Name=\"regions\" "
							  "NumberOfComponents=\"1\"/>\n");
			}
			if(m_bWriteProcRanks){
				fprintf(file, "      <DataArray type=\"Int8\" Name=\"proc_ranks\" "
							  "NumberOfComponents=\"1\"/>\n");
			}

		//	loop all vector data
			for(VectorDataIterator iter = m_vVectorElemData.begin();
					iter != m_vVectorElemData.end(); ++iter)
			{
			//	get symb function
				const std::string& vtkName = (*iter).first;

				fprintf(file, "      <PDataArray type=\"Float32\" Name=\"%s\" "
							  "NumberOfComponents=\"%d\"/>\n",
							  vtkName.c_str(), 3);
			}

		//	loop all vector data
			for(MatrixDataIterator iter = m_vMatrixElemData.begin();
					iter != m_vMatrixElemData.end(); ++iter)
			{
			//	get symb function
				const std::string& vtkName = (*iter).first;

				fprintf(file, "      <PDataArray type=\"Float32\" Name=\"%s\" "
							  "NumberOfComponents=\"%d\"/>\n",
							  vtkName.c_str(), 9);
			}
			fprintf(file, "    </PCellData>\n");
		}

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

	PCL_DEBUG_BARRIER_ALL();
#endif
}


template <int TDim>
template <typename TFunction>
void VTKOutput<TDim>::
write_time_pvd(const char* filename, TFunction& u)
{
//	File
	FILE* file;

//	filename
	std::string name;

// 	get some numbers
	bool isOutputProc = GetLogAssistant().is_output_process();
	int numProcs = 1;
#ifdef UG_PARALLEL
	numProcs = pcl::NumProcs();
#endif

//	get time steps
	std::vector<number>& vTimestep = m_mTimestep[filename];

//	change locale to ensure decimal . is really a .
	char* oldLocale = setlocale (LC_ALL, NULL);
	setlocale(LC_NUMERIC, "C");

	if (isOutputProc)
	{
	//	get file name
		pvd_filename(name, filename);

	//	open file
		file = fopen(name.c_str(), "w");
		if (file == NULL)
			UG_THROW("Cannot print to file.");

	// 	Write to file
		fprintf(file, "<?xml version=\"1.0\"?>\n");

	//	Write comment
		write_comment_printf(file);

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
				if(numProcs > 1) pvtu_filename(name, filename, -1, 0, step);

				name = FilenameWithoutPath(name);
				fprintf(file, "  <DataSet timestep=\"%.17g\" part=\"%d\" file=\"%s\"/>\n",
				        		vTimestep[step], 0, name.c_str());
			}
		}
		else
		{
			for(int step = 0; step < (int)vTimestep.size(); ++step)
				for(int si = 0; si < u.num_subsets(); ++si)
				{
					vtu_filename(name, filename, 0, si, u.num_subsets()-1, step);
					if(numProcs > 1) pvtu_filename(name, filename, si, u.num_subsets()-1, step);

					name = FilenameWithoutPath(name);
					fprintf(file, "  <DataSet timestep=\"%.17g\" part=\"%d\" file=\"%s\"/>\n",
					        	vTimestep[step], si, name.c_str());
				}
		}

	//	close file
		fprintf(file, "  </Collection>\n");
		fprintf(file, "</VTKFile>\n");
		fclose(file);
	}

// restore old locale
	setlocale(LC_NUMERIC, oldLocale);

	#ifdef UG_PARALLEL
		PCL_DEBUG_BARRIER_ALL();
	#endif
}


template <int TDim>
template <typename TFunction>
void VTKOutput<TDim>::
write_time_processwise_pvd(const char* filename, TFunction& u)
{
//	File
	FILE* file;

//	filename
	std::string name;

// 	get some numbers
	bool isOutputProc = GetLogAssistant().is_output_process();
	int numProcs = 1;
#ifdef UG_PARALLEL
	numProcs = pcl::NumProcs();
#endif

//	get time steps
	std::vector<number>& vTimestep = m_mTimestep[filename];

//	change locale to ensure decimal . is really a .
	char* oldLocale = setlocale (LC_ALL, NULL);
	setlocale(LC_NUMERIC, "C");

	if (isOutputProc && numProcs > 1)
	{
	//	adjust filename
		std::string procName = filename;
		procName.append("_processwise");
		pvd_filename(name, procName);

	//	open file
		file = fopen(name.c_str(), "w");
		if (file == NULL)
			UG_THROW("Cannot print to file.");

	// 	Write to file
		fprintf(file, "<?xml version=\"1.0\"?>\n");

	//	Write comment
		write_comment_printf(file);

		fprintf(file, "<VTKFile type=\"Collection\" version=\"0.1\">\n");
		fprintf(file, "  <Collection>\n");

	// 	include files from all procs
		for(int step = 0; step < (int)vTimestep.size(); ++step)
			for (int rank = 0; rank < numProcs; rank++)
				for(int si = 0; si < u.num_subsets(); ++si)
				{
					vtu_filename(name, filename, rank, si, u.num_subsets()-1, step);
					if(numProcs > 1) pvtu_filename(name, filename, si, u.num_subsets()-1, step);

					name = FilenameWithoutPath(name);
					fprintf(file, "  <DataSet timestep=\"%.17g\" part=\"%d\" file=\"%s\"/>\n",
					        vTimestep[step], rank, name.c_str());
				}

	//	end file
		fprintf(file, "  </Collection>\n");
		fprintf(file, "</VTKFile>\n");
		fclose(file);
	}

// restore old locale
	setlocale(LC_NUMERIC, oldLocale);

	#ifdef UG_PARALLEL
		PCL_DEBUG_BARRIER_ALL();
	#endif
}



template <int TDim>
template <typename TFunction>
void VTKOutput<TDim>::
write_time_pvd_subset(const char* filename, TFunction& u, int si)
{
//	File
	FILE* file;

//	filename
	std::string name;

// 	get some numbers
	bool isOutputProc = GetLogAssistant().is_output_process();
	int numProcs = 1;
#ifdef UG_PARALLEL
	numProcs = pcl::NumProcs();
#endif

//	get time steps
	std::vector<number>& vTimestep = m_mTimestep[filename];

//	change locale to ensure decimal . is really a .
	char* oldLocale = setlocale (LC_ALL, NULL);
	setlocale(LC_NUMERIC, "C");

	if (isOutputProc)
	{
	//	get file name
		pvd_filename(name, filename);

	//	open file
		file = fopen(name.c_str(), "w");
		if (file == NULL)
			UG_THROW("Cannot print to file.");

	// 	Write to file
		fprintf(file, "<?xml version=\"1.0\"?>\n");

	//	Write comment
		write_comment_printf(file);

		fprintf(file, "<VTKFile type=\"Collection\" version=\"0.1\">\n");
		fprintf(file, "  <Collection>\n");

	// 	include files from all procs
		for(int step = 0; step < (int)vTimestep.size(); ++step)
		{
			vtu_filename(name, filename, 0, si, u.num_subsets()-1, step);
			if(numProcs > 1) pvtu_filename(name, filename, si, u.num_subsets()-1, step);

			name = FilenameWithoutPath(name);
			fprintf(file, "  <DataSet timestep=\"%g\" part=\"%d\" file=\"%s\"/>\n",
						vTimestep[step], si, name.c_str());
		}

	//	close file
		fprintf(file, "  </Collection>\n");
		fprintf(file, "</VTKFile>\n");
		fclose(file);
	}

// restore old locale
	setlocale(LC_NUMERIC, oldLocale);

	#ifdef UG_PARALLEL
		PCL_DEBUG_BARRIER_ALL();
	#endif
}

}


