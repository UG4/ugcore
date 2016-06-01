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

#ifndef __H__UG__LIB_DISC__FUNCTION_SPACE__GRID_FUNCTION_UTIL__
#define __H__UG__LIB_DISC__FUNCTION_SPACE__GRID_FUNCTION_UTIL__

#include <boost/function.hpp>
#include "common/util/file_util.h"

#include "lib_algebra/cpu_algebra/sparsematrix_print.h"
#include "lib_algebra/operator/interface/matrix_operator.h"
#include "lib_algebra/operator/debug_writer.h"
#include "lib_algebra/operator/vector_writer.h"
#include "lib_disc/io/vtkoutput.h"
#include "lib_disc/spatial_disc/constraints/constraint_interface.h"
#include "lib_disc/dof_manager/dof_distribution.h"
#include <vector>
#include <string>
#include "lib_algebra/common/matrixio/matrix_io_mtx.h"
#include "lib_algebra/common/connection_viewer_output.h"
#include "lib_algebra/common/csv_gnuplot_output.h"
#include "lib_disc/spatial_disc/user_data/user_data.h"
#include "lib_disc/common/groups_util.h"
#include "lib_disc/common/geometry_util.h"
#include "lib_disc/function_spaces/integrate.h"
#include "lib_grid/tools/periodic_boundary_manager.h"

#include "grid_function.h"
#include "dof_position_util.h"

#ifdef UG_PARALLEL
#include "pcl/pcl.h"
#endif

namespace ug {

////////////////////////////////////////////////////////////////////////////////
// 	AverageComponent
////////////////////////////////////////////////////////////////////////////////

/**
 * Subtracts a given value from all DoFs of a given grid function associated with
 * a given element type.
 *
 * \tparam TGridFunction	grid function type
 * \tparam TBaseElem		type of the elements the DoFs are associated with
 */
template<typename TGridFunction, typename TBaseElem>
void SubtractValueFromComponent(SmartPtr<TGridFunction> spGF, size_t fct, number sub)
{
	typedef TGridFunction GF;
	typedef typename GF::template traits<TBaseElem>::const_iterator iter_type;

	iter_type iter = spGF->template begin<TBaseElem>();
	iter_type iterEnd = spGF->template end<TBaseElem>();

	PeriodicBoundaryManager* pbm = spGF->domain()->grid()->periodic_boundary_manager();

//  loop elems
	std::vector<DoFIndex> vMultInd;
	for(; iter != iterEnd; ++iter)
	{
	//	get element
		TBaseElem* elem = *iter;

	//	skip periodic ghosts
		if (pbm && pbm->is_slave(elem)) continue;

	//  get global indices
		spGF->inner_dof_indices(elem, fct, vMultInd);

	//	sum up value
		for(size_t i = 0; i < vMultInd.size(); ++i)
		{
			DoFRef(*spGF, vMultInd[i]) -= sub;
		}
	}
}

/**
 * Subtracts the average of a given grid function (computed by taking into account
 * the shape functions), shifted by the value of 'mean', from all the DoFs of this
 * grid function. Typically, it should result in a grid function with the average
 * equal to 'mean'.
 *
 * \remark Note that mean is not the desired average but the one multiplied by the area of the domain.
 *
 * \tparam TGridFunction	type of the grid function
 */
template<typename TGridFunction>
void AdjustMeanValue
(
	SmartPtr<TGridFunction> spGF, ///< the grid function
	const std::vector<std::string>& vCmp, ///< the component to adjust the average for
	number mean ///< the desired mean value
)
{
	typedef TGridFunction GF;
	PROFILE_FUNC_GROUP("gmg");

	if(vCmp.empty())
		return;

	if(spGF.invalid())
		UG_THROW("AverageComponent: expects a valid GridFunction.");

	ConstSmartPtr<DoFDistributionInfo> ddinfo =
								spGF->approx_space()->dof_distribution_info();

//	compute integral of components
	const number area = Integral(1.0, spGF);
	std::vector<number> vIntegral(vCmp.size(), 0.0);
	for(size_t f = 0; f < vCmp.size(); f++){
		const size_t fct = spGF->fct_id_by_name(vCmp[f].c_str());
		vIntegral[f] = Integral(spGF, vCmp[f].c_str(), NULL, ddinfo->lfeid(fct).order());
	}

//	subtract value
	for(size_t f = 0; f < vCmp.size(); f++)
	{
		const number sub = (vIntegral[f] - mean) / area;
		const size_t fct = spGF->fct_id_by_name(vCmp[f].c_str());

		if(ddinfo->max_fct_dofs(fct, VERTEX)) SubtractValueFromComponent<GF, Vertex>(spGF, fct, sub);
		if(ddinfo->max_fct_dofs(fct, EDGE)) SubtractValueFromComponent<GF, Edge>(spGF, fct, sub);
		if(ddinfo->max_fct_dofs(fct, FACE)) SubtractValueFromComponent<GF, Face>(spGF, fct, sub);
		if(ddinfo->max_fct_dofs(fct, VOLUME)) SubtractValueFromComponent<GF, Volume>(spGF, fct, sub);
	}
}

/**
 * Subtracts the average of a given grid function (computed by taking into account
 * the shape functions) from all the DoFs of this grid function. Typically, it
 * should result in a grid function with the zero average.
 *
 * \tparam TGridFunction	type of the grid function
 */
template<typename TGridFunction>
void AdjustMeanValue(SmartPtr<TGridFunction> spGF, const std::vector<std::string>& vCmp)
{
	AdjustMeanValue(spGF, vCmp, 0.0);
}

/**
 * Subtracts the average of a given grid function (computed by taking into account
 * the shape functions) from all the DoFs of this grid function. Typically, it
 * should result in a grid function with the zero average.
 *
 * \tparam TGridFunction	type of the grid function
 */
template<typename TGridFunction>
void AdjustMeanValue(SmartPtr<TGridFunction> spGF, const std::string& fcts)
{
	AdjustMeanValue(spGF, TokenizeTrimString(fcts), 0.0);
}

/**
 * Subtracts the average of a given grid function (computed by taking into account
 * the shape functions), shifted by the value of 'mean', from all the DoFs of this
 * grid function. Typically, it should result in a grid function with the average
 * equal to 'mean'.
 *
 * \remark Note that mean is not the desired average but the one multiplied by the area of the domain.
 *
 * \tparam TGridFunction	type of the grid function
 */
template<typename TGridFunction>
void AdjustMeanValue(SmartPtr<TGridFunction> spGF, const std::string& fcts, number mean)
{
	AdjustMeanValue(spGF, TokenizeTrimString(fcts), mean);
}

////////////////////////////////////////////////////////////////////////////////
// 	Summing up components
////////////////////////////////////////////////////////////////////////////////

/**
 * Sums values of a grid function at given elements (like vertices) in given subsets
 *
 * \tparam TGridFunction	the grid function type
 * \tparam TBaseElem		the given element type (like 'Vertex')
 */
template <typename TGridFunction, typename TBaseElem>
number SumGFValuesAt
(
	TGridFunction * u, ///< the grid function
	size_t fct, ///< index of the function
	SubsetGroup & ssGroup ///< the subsets
)
{
	typedef typename TGridFunction::template traits<TBaseElem>::const_iterator t_elem_iter;
	
	std::vector<DoFIndex> ind;
	number sum (1);
	
//	Loop the selected subsets
	sum = 0;
	for (size_t i = 0; i < ssGroup.size (); i++)
	{
		int ssi = ssGroup [i];
		
	//	Loop the vertices in the subset
		for (t_elem_iter vi = u->template begin<TBaseElem> (ssi);
					vi != u->template end<TBaseElem> (ssi); ++vi)
		{
		    TBaseElem * vert = *vi;
		    
		//	indices at this vertex
			u->inner_dof_indices (vert, fct, ind);
			
			if (ind.size () != 1)
				UG_THROW ("SumGridFuncValuesAt: The function must be scalar!");
		
		//	add the contribution
			sum += DoFRef (*u, ind [0]);
		}
	}
	
#ifdef UG_PARALLEL
	pcl::ProcessCommunicator procComm;
	sum = procComm.allreduce (sum, PCL_RO_SUM);
#endif
	
	return sum;
}

/**
 * Sums values of a grid function at vertices in given subsets.
 *
 * This version takes symbolic data.
 *
 * \tparam TGridFunction	the grid function type
 * \tparam TBaseElem		the given element type (like 'Vertex')
 */
template <typename TGridFunction, typename TBaseElem>
number SumGFValuesAt
(
	TGridFunction * u, ///< the grid function
	const char * fct_names, ///< index of the function
	const char * subset_names ///< the subsets
)
{
//	Get the function index
	std::vector<std::string> vfctNames;
	TokenizeString (fct_names, vfctNames);
	if (vfctNames.size () != 1)
		UG_THROW ("SumGridFuncValuesAt: Exactly one function name must be specified.");
	FunctionGroup fctGroup (u->function_pattern ());
	fctGroup.add (vfctNames [0]);

//	Get the subset indices
	std::vector<std::string> vssNames;
	TokenizeString (subset_names, vssNames);
	SubsetGroup ssGroup (u->domain()->subset_handler ());
	ssGroup.add (vssNames);
	
//	Compute the sum
	return SumGFValuesAt<TGridFunction, TBaseElem> (u, fctGroup[0], ssGroup);
}

////////////////////////////////////////////////////////////////////////////////
// 	Writing algebra
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
void WriteAlgebraIndices(std::string name, ConstSmartPtr<TDomain> domain,  ConstSmartPtr<DoFDistribution> dd)
{
	PROFILE_FUNC_GROUP("debug");

#ifdef UG_PARALLEL
	name = ConnectionViewer::GetParallelName(name, dd->layouts()->proc_comm());
#endif

	std::vector<size_t> fctIndex;
	std::vector<std::string> fctNames;

	ExtractAlgebraIndices<TDomain>(domain, dd, fctIndex);

	size_t numFct = dd->num_fct();
	if(numFct <= 1) return;

	fctNames.resize(numFct);
	for(size_t i=0; i<numFct; i++)
		fctNames[i] = dd->name(i);

	name.append(".indices");
	std::fstream file(name.c_str(), std::ios::out);

	//std::cout << "name is " << name << "\n";
	file << "NUMDOF " << fctNames.size() << "\n";
	for(size_t i=0; i<numFct; i++)
		file << fctNames[i] << "\n";

	for(size_t i=0; i<fctIndex.size(); i++)
		file << fctIndex[i] << "\n";
}

template<class TFunction>
void WriteMatrixToConnectionViewer(const char *filename,
		const typename TFunction::algebra_type::matrix_type &A,
		const TFunction &u) {

	PROFILE_FUNC_GROUP("debug");
//	check name
	if ( !FileTypeIs( filename, ".mat") ) {
		UG_THROW( "Only '.mat' format supported for matrices, but"
		          " filename is '" << filename << "'." );
	}

//	position array
	const static int dim = TFunction::domain_type::dim;
	std::vector<MathVector<dim> > vPos;
	ExtractPositions(u, vPos);

// 	write matrix
	if(vPos.empty())
		ConnectionViewer::WriteMatrixPar( filename, A, (MathVector<dim>*)NULL, dim );
	else
		ConnectionViewer::WriteMatrixPar( filename, A, &vPos[0], dim );

	WriteAlgebraIndices(filename, u.domain(),u.dof_distribution());

}

template<typename TGridFunction>
void SaveMatrixForConnectionViewer(
		TGridFunction& u,
		MatrixOperator<typename TGridFunction::algebra_type::matrix_type,
				typename TGridFunction::vector_type>& A,
		const char* filename) {
	PROFILE_FUNC_GROUP("debug");
//	forward
	WriteMatrixToConnectionViewer(filename, A.get_matrix(), u);
}

/**
 * \brief Save the assembled matrix of a matrix operator to MatrixMarket format
 * 
 * \param[in] filename name of the file; must end on '<tt>.mtx</tt>'
 * \param[in] A matrix operator of with \c CPUAlgebra matrix
 * \param[in] comment optional comment for the header of the MTX file
 * 
 * \note Until now only CPUAlgebra matrices are supported.
 */
// TODO extend MatrixIO to support other than CPUAlgebra
inline void SaveMatrixToMTX( const char *filename,
		MatrixOperator< CPUAlgebra::matrix_type, CPUAlgebra::vector_type >& A,
		std::string comment="%Generated with ug4." ) {
	PROFILE_FUNC_GROUP("debug");

	// check name
	if ( !FileTypeIs( filename, ".mtx" ) ) {
		UG_THROW( "Please use '.mtx' as file extension for MatrixMarket exchange files."
		          " (Filename is '" << filename << "')" );
	}
	
	// automatically detect the file mode to use
	MatrixIO::OpenMode openMode = (FileExists( filename )) ? MatrixIO::EXISTING : MatrixIO::NEW;
	
	// create MatrixIO object for handling the export
	MatrixIOMtx mtx( filename, openMode );
	mtx.write_from( A.get_matrix(), comment );
}

template<class TFunction>
void WriteVectorToConnectionViewer(const char *filename,
		const typename TFunction::algebra_type::vector_type &b,
		const TFunction &u,
		const typename TFunction::algebra_type::vector_type *pCompareVec = NULL) {
	PROFILE_FUNC_GROUP("debug");
//	check name
	if ( !FileTypeIs( filename, ".vec") ) {
		UG_THROW( "Only '.vec' format supported for vectors, but"
		          " filename is '" << filename << "'." );
	}

// 	get positions of vertices
	const static int dim = TFunction::domain_type::dim;
	std::vector<MathVector<dim> > vPos;
	ExtractPositions(u, vPos);

//	write vector
	ConnectionViewer::WriteVectorPar( filename, b, &vPos[0], dim, pCompareVec);
	WriteAlgebraIndices(filename, u.domain(),u.dof_distribution());
}

template<class TFunction>
void WriteVectorToConnectionViewer(
		const char *filename,
		const typename TFunction::algebra_type::matrix_type &A,
		const typename TFunction::algebra_type::vector_type &b,
		const TFunction &u,
		const typename TFunction::algebra_type::vector_type *pCompareVec = NULL) {
	PROFILE_FUNC_GROUP("debug");
//	get dimension
	const static int dim = TFunction::domain_type::dim;

//	check name
	if ( !FileTypeIs( filename, ".vec") ) {
		UG_THROW( "Only '.vec' format supported for vectors." );
	}

// 	get positions of vertices
	std::vector<MathVector<dim> > positions;
	ExtractPositions(u, positions);

//	write vector
	ConnectionViewer::WriteVectorPar( filename, A, b, &positions[0], dim, pCompareVec );
}

template<typename TGridFunction>
void SaveVectorForConnectionViewer(TGridFunction& b, const char* filename) {
	WriteVectorToConnectionViewer(filename, b, b);
}

template<typename TGridFunction>
void SaveVectorDiffForConnectionViewer(TGridFunction& b, TGridFunction& bCompare, const char* filename) {
	WriteVectorToConnectionViewer(filename, b, b, &bCompare);
}

template<typename TGridFunction>
void SaveVectorForConnectionViewer(
		TGridFunction& u,
		MatrixOperator<typename TGridFunction::algebra_type::matrix_type,
				typename TGridFunction::vector_type>& A,
		const char* filename) {
	WriteVectorToConnectionViewer(filename, A.get_matrix(), u, u);
}

template<typename TGridFunction>
void SaveVectorDiffForConnectionViewer(
		TGridFunction& u,
		TGridFunction& compareVec,
		MatrixOperator<typename TGridFunction::algebra_type::matrix_type,
				typename TGridFunction::vector_type>& A,
		const char* filename) {
//	forward
	WriteVectorToConnectionViewer(filename, A.get_matrix(), u, u, &compareVec);
}

// from connection_viewer_input.h
// with additional checks
template<typename vector_type>
bool ReadVector(std::string filename, vector_type &vec,int dim)
{
    Progress p;
	std::cout << " Reading std::vector from " <<  filename << "... ";
	std::fstream matfile(filename.c_str(), std::ios::in);
	if(matfile.is_open() == false) { std::cout << "failed.\n"; return false; }

	int version=-1, dimension=-1, gridsize;

	matfile >> version;
	matfile >> dimension;
	matfile >> gridsize;

	assert(version == 1);
	assert(dimension == dim);
	// todo check positions and not just size
	assert(gridsize == (int)vec.size());
	

	PROGRESS_START(prog, gridsize*2, "ReadVector " << dimension << "d from " << filename << " , " << gridsize << " x " << gridsize);
	for(int i=0; i<gridsize; i++)
	{
		if(i%100) { PROGRESS_UPDATE(prog, i); }
		if(matfile.eof())
		{
			std::cout << " failed.\n";
			assert(0);
			return false;
		}
		double x, y, z;
		matfile >> x >> y;
		if(dimension==3) matfile >> z;
	}

	int printStringsInWindow;
	matfile >> printStringsInWindow;

	// vec.resize(gridsize);
	bool bEOF = matfile.eof();
	while(!bEOF)
	{
		int from, to; double value;
		char c = matfile.peek();
		if(c == -1 || c == 'c' || c == 'v' || matfile.eof())
			break;

		matfile >> from >> to >> value;
		assert(from == to);
		vec[from] = value;
		if(from%100) { PROGRESS_UPDATE(prog, from); }
		bEOF = matfile.eof();
	}
	return true;
}

// load vector that has been saved in connection viewer format and write it
// into grid function
template<typename TGridFunction>
void LoadVector(TGridFunction& u,const char* filename){
	PROFILE_FUNC();
	typename TGridFunction::algebra_type::vector_type b;
	b.resize(u.num_indices());
	ReadVector(filename,b,TGridFunction::dim);
	u.assign(b);
}

// Same as before, but for comma separated value (CSV)
template<class TFunction>
void WriteVectorCSV(const char *filename,
		const typename TFunction::algebra_type::vector_type &b,
		const TFunction &u) {
	PROFILE_FUNC_GROUP("debug");
//	get dimension
	const static int dim = TFunction::domain_type::dim;

//	check name
	if ( !FileTypeIs( filename, ".csv") ) {
		UG_THROW( "Only '.csv' format supported for vectors, but"
		          " filename is '" << filename << "'." );
	}

//	extended filename
//	add p000X extension in parallel
#ifdef UG_PARALLEL
	std::string name(filename);
	size_t iExtPos = name.find_last_of(".");
	name.resize(iExtPos);
	int rank = pcl::ProcRank();
	char ext[20];
	sprintf(ext, "_p%05d.csv", rank);
	name.append(ext);
#endif

// 	get positions of vertices
	std::vector<MathVector<dim> > positions;
	ExtractPositions(u, positions);

//	write vector
	WriteVectorCSV( filename, b, &positions[0], dim );
}

template<typename TGridFunction>
void SaveVectorCSV(TGridFunction& b, const char* filename) {
	PROFILE_FUNC();
	WriteVectorCSV(filename, b, b);
}

/**
 * \brief Calculates the average of the pointwise difference of two functions on given subset
 * \details Iterates over all vertices of given subset, calculates difference
 *   of \c fct1 and \c fct2 and computes arithmetic mean of the differences at
 *   all vertices (grid points).
 * \note \c fct1 and \c fct2 must both be defined on \c subset !
 * \param[in] spGridFct GridFunction holding functions \c fct1 and \c fct2 defined on subset \c subset
 * \param[in] subset name of the subset to compare on
 * \param[in] fct1 name of the first function
 * \param[in] fct2 name of the second function
 * \return arithmetic average of the pointwise differences of \c fct1 and \c fct2 on
 *   all vertices of \c subset
 */
template<typename TDomain, typename TAlgebra>
number AverageFunctionDifference(
		SmartPtr< GridFunction<TDomain, TAlgebra> > spGridFct,
		std::string subset, std::string fct1, std::string fct2 )
{
	PROFILE_FUNC();
	// get subset index
	size_t subSetID = spGridFct->subset_id_by_name( subset.c_str() );

	// get function indices
	size_t fct1ID = spGridFct->fct_id_by_name( fct1.c_str() );
	size_t fct2ID = spGridFct->fct_id_by_name( fct2.c_str() );

	// create space for sum of difference
	number sum = 0.0;
	size_t numElements = 0;

	// loop over all vertices in given subset and compare values of fct1 and fct2
	typedef typename GridFunction<TDomain, TAlgebra>::template traits<Vertex>::const_iterator gridFctIterator;
	for( gridFctIterator iter = spGridFct->template begin<Vertex>((int)subSetID); 
	       iter != spGridFct->template end<Vertex>((int)subSetID); ++iter ) {
		// get dof_indices for the two functions on given subset
		std::vector< DoFIndex > indFct1, indFct2;
		spGridFct->template dof_indices<Vertex>( *iter, fct1ID, indFct1 );
		spGridFct->template dof_indices<Vertex>( *iter, fct2ID, indFct2 );

		// calculate the difference between the two functions at this grid point
		sum += DoFRef( *spGridFct, indFct1[0] ) - DoFRef( *spGridFct, indFct2[0] );
		numElements++;
	}
	
	// return overal arithmetic average of differences
	return sum / numElements;
}


////////////////////////////////////////////////////////////////////////////////
//	Grid Function Debug Writer
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain, typename TAlgebra>
class GridFunctionDebugWriter: public IDebugWriter<TAlgebra>
{
	///	dimension
	static const int dim = TDomain::dim;

public:
	///	type of matrix
	typedef TAlgebra algebra_type;

	///	type of vector
	typedef typename algebra_type::vector_type vector_type;

	///	type of matrix
	typedef typename algebra_type::matrix_type matrix_type;

	///	type of approximation space
	typedef ApproximationSpace<TDomain> approximation_space_type;

	typedef IDebugWriter<TAlgebra> super;
	using super::get_base_dir;

public:
	///	Constructor
	GridFunctionDebugWriter(
			SmartPtr<ApproximationSpace<TDomain> > spApproxSpace) :
			m_spApproxSpace(spApproxSpace), bConnViewerOut(
					true), bConnViewerIndices(false), bVTKOut(true), m_printConsistent(true)
	{
		IVectorDebugWriter<vector_type>::m_currentDim = dim;
		reset();
	}

	///	sets the grid level
	void set_grid_level(const GridLevel& gridLevel) {
		m_glFrom = gridLevel; m_glTo = gridLevel;
	}

	///	sets the grid level
	void set_grid_levels(const GridLevel& glFrom, const GridLevel& glTo) {
		m_glFrom = glFrom; m_glTo = glTo;
	}

	///	returns current grid level
	GridLevel grid_level() const {return m_glFrom;}

	///	sets to toplevel on surface
	void reset() {
		set_grid_level(GridLevel(GridLevel::TOP, GridLevel::SURFACE));
	}


	///	sets if writing to vtk
	void set_vtk_output(bool b) {bVTKOut = b;}

	///	sets if writing to conn viewer
	void set_conn_viewer_output(bool b) {bConnViewerOut = b;}

	/// sets if .indices file is written or conn viewer
	void set_conn_viewer_indices(bool b) { bConnViewerIndices = b; }

	///	sets if data shall be made consistent before printing
	void set_print_consistent(bool b) {m_printConsistent = b;}

	///	write vector
	virtual void write_vector(const vector_type& vec, const char* filename) {

		//	write to conn viewer
		if (bConnViewerOut)
			write_vector_to_conn_viewer(vec, filename);

		//	write to vtk
		if (bVTKOut)
			write_vector_to_vtk(vec, filename);
	}

	virtual void update_positions()
	{
		extract_positions(m_glFrom);
	}

	///	write matrix
	virtual void write_matrix(const matrix_type& mat, const char* filename) {
		PROFILE_FUNC_GROUP("debug");
		//	something to do ?
		if (!bConnViewerOut)
			return;

		update_positions();

		std::string name = get_base_dir() + "/" + filename;

		if (!FileExists(get_base_dir())) {
			UG_WARNING("GridFunctionDebugWriter::write_matrix: directory "
						<< get_base_dir() << "does not exist.");
			UG_WARNING("GridFunctionDebugWriter::write_matrix: using cwd "
						"as basedir.");
			name = "./"; name.append(filename);
		}

		if ( !FileTypeIs( filename, ".mat" ) ) {
			UG_THROW( "Only '.mat' format supported for matrices, but"
			          " filename is '" << filename << "'." );
		}

		//	write to file
		write_algebra_indices_CV(name);
		if(m_glFrom == m_glTo){
			if(mat.num_rows() != mat.num_cols())
				UG_THROW("DebugWriter: grid level the same, but non-square matrix.");

			const std::vector<MathVector<dim> >& vPos = this->template get_positions<dim>();
			if(vPos.empty())
				ConnectionViewer::WriteMatrixPar(name.c_str(), mat,(MathVector<dim>*)NULL, dim);
			else
				ConnectionViewer::WriteMatrixPar(name.c_str(), mat, &vPos[0], dim);
		}
		else{
			const std::vector<MathVector<dim> >& vFromPos = this->template get_positions<dim>();

			std::vector<MathVector<dim> > vToPos;
			ExtractPositions<TDomain>(m_spApproxSpace->domain(),
			                          m_spApproxSpace->dof_distribution(m_glTo),
			                          vToPos);

			if(mat.num_cols() == vFromPos.size() && mat.num_rows() == vToPos.size())
			{
				ConnectionViewer::WriteMatrixPar(name, mat, vFromPos, vToPos, dim);
			}
			else{
				UG_THROW("GridFunctionDebugWriter: Wrong size of matrix for writing"
						"matrix ("<<m_glTo<<" x "<<m_glFrom<<"), that would be a ("
						<<vToPos.size()<<" x "<<vFromPos.size()<<") matrix. But "
						"passed matrix of size ("<<mat.num_rows()<<" x "
						<<mat.num_cols()<<").");
			}
		}
	}

protected:

	void write_algebra_indices_CV(std::string name)
	{
		if(bConnViewerIndices)
			WriteAlgebraIndices<TDomain>(name, m_spApproxSpace->domain(), m_spApproxSpace->dof_distribution(m_glFrom));
	}

	///	reads the positions
	void extract_positions(const GridLevel& gridLevel) {
		//	extract positions for this grid function
		std::vector<MathVector<dim> >& vPos =
				this->template positions<dim>();

		vPos.clear();
		ExtractPositions<TDomain>(
				m_spApproxSpace->domain(),
				m_spApproxSpace->dof_distribution(gridLevel),
				vPos);
	}

	///	write vector
	virtual void write_vector_to_conn_viewer(const vector_type& vec,
			const char* filename)
	{
		PROFILE_FUNC_GROUP("debug");
		update_positions();

		std::string name = get_base_dir() + "/" + filename;

		if (!FileExists(get_base_dir())) {
			UG_WARNING("GridFunctionDebugWriter::write_vector_to_conn_viewer: directory "
						<< get_base_dir() << "does not exist.");
			UG_WARNING("GridFunctionDebugWriter::write_vector_to_conn_viewer: using cwd "
						"as basedir.");
			name = "./"; name.append(filename);
		}

		if ( !FileTypeIs( filename, ".vec" ) ) {
			UG_THROW( "Only '.vec' format supported for vectors, but"
			          " filename is '" << name << "'.");
		}

		//	write
		write_algebra_indices_CV(filename);
		const std::vector<MathVector<dim> >& vPos =
				this->template get_positions<dim>();
		if(vPos.empty())
			ConnectionViewer::WriteVectorPar(name.c_str(), vec, (MathVector<dim>*)NULL, dim);
		else
			ConnectionViewer::WriteVectorPar(name.c_str(), vec, &vPos[0], dim);
	}

	void write_vector_to_vtk(const vector_type& vec, const char* filename) {
		PROFILE_FUNC_GROUP("debug");
		//	check name
		std::string name = get_base_dir() + "/" + filename;

		if (!FileExists(get_base_dir())) {
			UG_WARNING("GridFunctionDebugWriter::write_vector_to_vtk: directory "
						<< get_base_dir() << "does not exist.");
			UG_WARNING("GridFunctionDebugWriter::write_vector_to_vtk: using cwd "
						"as basedir.");
			name = "./"; name.append(filename);
		}

		typedef GridFunction<TDomain, TAlgebra> TGridFunction;
		TGridFunction vtkFunc(
				m_spApproxSpace,
				m_spApproxSpace->dof_distribution(m_glFrom));
		vtkFunc.resize_values(vec.size());
		vtkFunc.assign(vec);
		VTKOutput<dim> out;
		out.print(filename, vtkFunc, m_printConsistent);
	}

protected:
	//	underlying approximation space
	SmartPtr<approximation_space_type> m_spApproxSpace;

	//	flag if write to conn viewer
	bool bConnViewerOut;
	bool bConnViewerIndices;

	//	flag if write to vtk
	bool bVTKOut;

	//	flag if data shall be made consistent before printing
	bool m_printConsistent;

	//	current grid level
	GridLevel m_glFrom, m_glTo;
};

////////////////////////////////////////////////////////////////////////////////
//	Grid Function Position Provider
////////////////////////////////////////////////////////////////////////////////

template<typename TGridFunction>
class GridFunctionPositionProvider: public IPositionProvider<
		TGridFunction::domain_type::dim> {
public:
	///	Constructor
	GridFunctionPositionProvider() :
			m_pGridFunc(NULL) {
	}
	GridFunctionPositionProvider(const TGridFunction& u) :
			m_pGridFunc(&u) {
	}
	void set_reference_grid_function(const TGridFunction& u) {
		m_pGridFunc = &u;
	}

	virtual bool get_positions(
			std::vector<MathVector<TGridFunction::domain_type::dim> >&vec) {
		UG_ASSERT(m_pGridFunc != NULL,
				"provide a grid function with set_reference_grid_function");
		ExtractPositions(*m_pGridFunc, vec);
		return true;
	}

protected:
	const TGridFunction *m_pGridFunc;
};

////////////////////////////////////////////////////////////////////////////////
//	Grid Function Vector Writer
////////////////////////////////////////////////////////////////////////////////

template<typename TGridFunction, typename TVector>
class GridFunctionVectorWriter: public IVectorWriter<TVector> {
public:
	typedef typename TVector::value_type value_type;
	typedef typename TGridFunction::domain_type domain_type;
	typedef TVector vector_type;

public:
	///	Constructor
	GridFunctionVectorWriter() :
			m_pGridFunc(NULL) {
	}

	void set_user_data(SmartPtr<UserData<number, domain_type::dim> > userData) {
		m_userData = userData;
	}

	void set_reference_grid_function(const TGridFunction& u) {
		m_pGridFunc = &u;
	}

	/*virtual double calculate(MathVector<3> &v, double time)
	 {
	 }*/

	virtual bool update(vector_type &vec) {
		PROFILE_FUNC_GROUP("debug");
		UG_ASSERT(m_pGridFunc != NULL,
				"provide a grid function with set_reference_grid_function");
		UG_ASSERT(m_userData.valid(), "provide user data with set_user_data");

		const TGridFunction &u = *m_pGridFunc;

		//	get position accessor

		const typename domain_type::position_accessor_type& aaPos =
				u.domain()->position_accessor();

		//	number of total dofs
		int nr = u.num_indices();

		//	resize positions
		vec.resize(nr);

		typedef typename TGridFunction::template traits<Vertex>::const_iterator const_iterator;

		//	loop all subsets
		for (int si = 0; si < u.num_subsets(); ++si) {
			//	loop all vertices
			const_iterator iter = u.template begin<Vertex>(si);
			const_iterator iterEnd = u.template end<Vertex>(si);

			for (; iter != iterEnd; ++iter) {
				//	get vertex
				Vertex* v = *iter;

				//	algebra indices vector
				std::vector<size_t> ind;

				//	load indices associated with vertex
				u.inner_algebra_indices(v, ind);

				number t = 0.0;

				number d;
				(*m_userData)(d, aaPos[v], t, si);

				//	write
				for (size_t i = 0; i < ind.size(); ++i) {
					const size_t index = ind[i];
					for (size_t alpha = 0; alpha < GetSize(vec[index]); ++alpha)
						BlockRef(vec[index], alpha) = d;
				}
			}
		}
		return true;
	}

protected:
	const TGridFunction *m_pGridFunc;
	SmartPtr<UserData<number, domain_type::dim> > m_userData;
};

////////////////////////////////////////////////////////////////////////////////
//	Grid Function Vector Writer Dirichlet 0
////////////////////////////////////////////////////////////////////////////////

template<typename TGridFunction>
class GridFunctionVectorWriterDirichlet0: public IVectorWriter<
		typename TGridFunction::algebra_type::vector_type> {
public:
	typedef typename TGridFunction::domain_type domain_type;

	typedef typename TGridFunction::approximation_space_type approximation_space_type;

	typedef typename TGridFunction::algebra_type algebra_type;
	typedef typename algebra_type::vector_type vector_type;
	typedef typename vector_type::value_type value_type;

public:
	///	Constructor
	GridFunctionVectorWriterDirichlet0() :
			m_pApproxSpace(NULL), m_spPostProcess(NULL), m_level(-1) {
	}

	void set_level(size_t level) {
		m_level = level;
	}

	void init(SmartPtr<IConstraint<algebra_type> > pp,
			approximation_space_type& approxSpace) {
		m_spPostProcess = pp;
		m_pApproxSpace = &approxSpace;
	}

	/*virtual double calculate(MathVector<3> &v, double time)
	 {
	 }*/

	virtual bool update(vector_type &vec) {
		PROFILE_FUNC_GROUP("debug");
		UG_ASSERT(m_spPostProcess.valid(), "provide a post process with init");
		UG_ASSERT(m_pApproxSpace != NULL, "provide approximation space init");

		size_t numDoFs = 0;
		GridLevel gl;
		if (m_level == (size_t) -1) {
			numDoFs = m_pApproxSpace->dof_distribution(GridLevel(GridLevel::TOP, GridLevel::SURFACE))->num_indices();
			gl = GridLevel();
		} else {
			numDoFs =
					m_pApproxSpace->dof_distribution(GridLevel(m_level, GridLevel::LEVEL, true))->num_indices();
			gl = GridLevel(m_level, GridLevel::LEVEL);
		}
		vec.resize(numDoFs);
		vec.set(1.0);

		m_spPostProcess->adjust_defect(vec, vec, m_pApproxSpace->dof_distribution(gl));

		return true;
	}

protected:
	approximation_space_type * m_pApproxSpace;
	SmartPtr<IConstraint<algebra_type> > m_spPostProcess;
	size_t m_level;
};

} // end namespace ug

#endif /* __H__UG__LIB_DISC__FUNCTION_SPACE__GRID_FUNCTION_UTIL__ */
