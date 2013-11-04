/*
 * grid_function_util.h
 *
 *  Created on: 17.08.2010
 *      Author: andreasvogel
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

#include "grid_function.h"
#include "dof_position_util.h"

#ifdef UG_PARALLEL
#include "pcl/pcl.h"
#endif

namespace ug {

template<typename TDomain>
void WriteAlgebraIndices(std::string name, ConstSmartPtr<TDomain> domain,  ConstSmartPtr<DoFDistribution> dd, const std::vector<int>* pvMapGlobalToPatch=NULL)
{
/*
	std::vector<size_t> fctIndex;
	std::vector<std::string> fctNames;

	ExtractAlgebraIndices<TDomain>(domain, dd, fctIndex, pvMapGlobalToPatch);

	size_t numFct = dd->num_fct();
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
		//*/
}

template<class TFunction>
void WriteMatrixToConnectionViewer(const char *filename,
		const typename TFunction::algebra_type::matrix_type &A,
		const TFunction &u) {

	PROFILE_FUNC();
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

	WriteAlgebraIndices(filename, u.domain(),u.dof_distribution(), NULL);

}

template<typename TGridFunction>
void SaveMatrixForConnectionViewer(
		TGridFunction& u,
		MatrixOperator<typename TGridFunction::algebra_type::matrix_type,
				typename TGridFunction::vector_type>& A,
		const char* filename) {
	PROFILE_FUNC();
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
	PROFILE_FUNC();

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
	PROFILE_FUNC();
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
}

template<class TFunction>
void WriteVectorToConnectionViewer(
		const char *filename,
		const typename TFunction::algebra_type::matrix_type &A,
		const typename TFunction::algebra_type::vector_type &b,
		const TFunction &u,
		const typename TFunction::algebra_type::vector_type *pCompareVec = NULL) {
	PROFILE_FUNC();
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
	PROFILE_FUNC();
	WriteVectorToConnectionViewer(filename, b, b);
}

template<typename TGridFunction>
void SaveVectorDiffForConnectionViewer(TGridFunction& b, TGridFunction& bCompare, const char* filename) {
	PROFILE_FUNC();
	WriteVectorToConnectionViewer(filename, b, b, &bCompare);
}

template<typename TGridFunction>
void SaveVectorForConnectionViewer(
		TGridFunction& u,
		MatrixOperator<typename TGridFunction::algebra_type::matrix_type,
				typename TGridFunction::vector_type>& A,
		const char* filename) {
	PROFILE_FUNC();
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
	PROFILE_FUNC();
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
	PROFILE_FUNC();
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
	int rank = pcl::GetProcRank();
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
 * \param[in] spGridFct GridFunction holding functions \c fct1 and \c fct2 defined on subset \subset
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
	// get subset index
	size_t subSetID = spGridFct->subset_id_by_name( subset.c_str() );

	// get function indices
	size_t fct1ID = spGridFct->fct_id_by_name( fct1.c_str() );
	size_t fct2ID = spGridFct->fct_id_by_name( fct2.c_str() );

	// create space for sum of difference
	number sum = 0.0;
	size_t numElements = 0;

	// loop over all vertices in given subset and compare values of fct1 and fct2
	typedef typename GridFunction<TDomain, TAlgebra>::template traits<VertexBase>::const_iterator gridFctIterator;
	for( gridFctIterator iter = spGridFct->template begin<VertexBase>((int)subSetID); 
	       iter != spGridFct->template end<VertexBase>((int)subSetID); ++iter ) {
		// get dof_indices for the two functions on given subset
		std::vector< DoFIndex > indFct1, indFct2;
		spGridFct->template dof_indices<VertexBase>( *iter, fct1ID, indFct1 );
		spGridFct->template dof_indices<VertexBase>( *iter, fct2ID, indFct2 );

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

	///  base directory for output
	std::string m_baseDir;

public:
	///	type of matrix
	typedef TAlgebra algebra_type;

	///	type of vector
	typedef typename algebra_type::vector_type vector_type;

	///	type of matrix
	typedef typename algebra_type::matrix_type matrix_type;

	///	type of approximation space
	typedef ApproximationSpace<TDomain> approximation_space_type;

public:
	///	Constructor
	GridFunctionDebugWriter(
			SmartPtr<ApproximationSpace<TDomain> > spApproxSpace) :
			m_baseDir("."), m_spApproxSpace(spApproxSpace), bConnViewerOut(
					true), bVTKOut(true), m_printConsistent(true), m_pvMapGlobalToPatch(NULL) {
		reset();
	}

	///	sets the grid level
	void set_grid_level(const GridLevel& gridLevel) {
		m_gridLevel = gridLevel; m_coarseGridLevel = gridLevel;
	}

	///	sets the grid level
	void set_grid_levels(const GridLevel& fineGridLevel, const GridLevel& coarseGridLevel) {
		m_gridLevel = fineGridLevel; m_coarseGridLevel = coarseGridLevel;
	}

	///	returns current grid level
	GridLevel grid_level() const {return m_gridLevel;}

	///	sets to toplevel on surface
	void reset() {
		set_grid_level(GridLevel(GridLevel::TOPLEVEL, GridLevel::SURFACE));
		m_pvMapGlobalToPatch = NULL;
	}

	void set_map_global_to_patch(const std::vector<int>* pvMapGlobalToPatch){
		m_pvMapGlobalToPatch = pvMapGlobalToPatch;
	}

	/// set the base directory for output files (.vec and .mat)
	inline void set_base_dir(const char* const baseDir) {m_baseDir = std::string(baseDir);}

	///	sets if writing to vtk
	void set_vtk_output(bool b) {bVTKOut = b;}

	///	sets if writing to conn viewer
	void set_conn_viewer_output(bool b) {bConnViewerOut = b;}

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

	///	write matrix
	virtual void write_matrix(const matrix_type& mat, const char* filename) {
		//	something to do ?
		if (!bConnViewerOut)
			return;

		std::string name(m_baseDir); name.append("/").append(filename);

		if (!FileExists(m_baseDir.c_str())) {
			UG_WARNING("GridFunctionDebugWriter::write_matrix: directory "
						<< m_baseDir << "does not exist.");
			UG_WARNING("GridFunctionDebugWriter::write_matrix: using cwd "
						"as basedir.");
			name = "./"; name.append(filename);
		}

		if ( !FileTypeIs( filename, ".mat" ) ) {
			UG_THROW( "Only '.mat' format supported for matrices, but"
			          " filename is '" << filename << "'." );
		}

		//	write to file
		extract_algebra_indices(name);
		if(m_gridLevel == m_coarseGridLevel){
			if(mat.num_rows() != mat.num_cols())
				UG_THROW("DebugWriter: grid level the same, but non-square matrix.");

			extract_positions(m_gridLevel);
			std::vector<MathVector<dim> >& vPos =
					this->template get_positions<dim>();
			if(vPos.empty())
				ConnectionViewer::WriteMatrixPar(name.c_str(), mat,(MathVector<dim>*)NULL, dim);
			else
				ConnectionViewer::WriteMatrixPar(name.c_str(), mat, &vPos[0], dim);
		}
		else{
			if(m_pvMapGlobalToPatch){
				UG_LOG("Cannot write on different levels using MapGlobalToPatch. Skipping.\n");
				return;
			}
			extract_positions(m_gridLevel);
			std::vector<MathVector<dim> >& vFinePos =
					this->template get_positions<dim>();
			std::vector<MathVector<dim> > vCoarsePos;

			ExtractPositions<TDomain>(m_spApproxSpace->domain(),
			                 m_spApproxSpace->dof_distribution(m_coarseGridLevel),
			                 vCoarsePos);

			if(mat.num_cols() == vFinePos.size())
				ConnectionViewer::WriteMatrixPar(name, mat, vFinePos, vCoarsePos, dim);
			else
				ConnectionViewer::WriteMatrixPar(name, mat, vCoarsePos, vFinePos, dim);
		}
	}

protected:

	void extract_algebra_indices(std::string name)
	{
		WriteAlgebraIndices<TDomain>(name, m_spApproxSpace->domain(), m_spApproxSpace->dof_distribution(m_gridLevel), m_pvMapGlobalToPatch);
	}

	///	reads the positions
	void extract_positions(const GridLevel& gridLevel) {
		//	extract positions for this grid function
		std::vector<MathVector<dim> >& vPos =
				this->template get_positions<dim>();

		vPos.clear();
		ExtractPositions<TDomain>(
				m_spApproxSpace->domain(),
				m_spApproxSpace->dof_distribution(gridLevel),
				vPos, m_pvMapGlobalToPatch);
	}

	///	write vector
	virtual void write_vector_to_conn_viewer(const vector_type& vec,
			const char* filename) {
		std::string name(m_baseDir); name.append("/").append(filename);

		if (!FileExists(m_baseDir.c_str())) {
			UG_WARNING("GridFunctionDebugWriter::write_vector_to_conn_viewer: directory "
						<< m_baseDir << "does not exist.");
			UG_WARNING("GridFunctionDebugWriter::write_vector_to_conn_viewer: using cwd "
						"as basedir.");
			name = "./"; name.append(filename);
		}

		if ( !FileTypeIs( filename, ".vec" ) ) {
			UG_THROW( "Only '.vec' format supported for vectors, but"
			          " filename is '" << name << "'.");
		}

		//	write
		extract_positions(m_gridLevel);
		std::vector<MathVector<dim> >& vPos =
				this->template get_positions<dim>();
		if(vPos.empty())
			ConnectionViewer::WriteVectorPar(name.c_str(), vec, (MathVector<dim>*)NULL, dim);
		else
			ConnectionViewer::WriteVectorPar(name.c_str(), vec, &vPos[0], dim);
	}

	void write_vector_to_vtk(const vector_type& vec, const char* filename) {
		//	check name
		std::string name(m_baseDir); name.append("/").append(filename);

		if(m_pvMapGlobalToPatch){
			UG_LOG("Cannot write to vtk using MapGlobalToPatch. Skipping.\n");
			return;
		}

		if (!FileExists(m_baseDir.c_str())) {
			UG_WARNING("GridFunctionDebugWriter::write_vector_to_vtk: directory "
						<< m_baseDir << "does not exist.");
			UG_WARNING("GridFunctionDebugWriter::write_vector_to_vtk: using cwd "
						"as basedir.");
			name = "./"; name.append(filename);
		}

		if (m_gridLevel.type() == GridLevel::LEVEL) {
			typedef GridFunction<TDomain, TAlgebra> TGridFunction;
			TGridFunction vtkFunc(
					m_spApproxSpace,
					m_spApproxSpace->dof_distribution(m_gridLevel));
			vtkFunc.resize_values(vec.size());
			vtkFunc.assign(vec);
			VTKOutput<dim> out;
			out.print(filename, vtkFunc, m_printConsistent);
		} else if (m_gridLevel.type() == GridLevel::SURFACE) {
			typedef GridFunction<TDomain, TAlgebra> TGridFunction;
			TGridFunction vtkFunc(
					m_spApproxSpace,
					m_spApproxSpace->dof_distribution(m_gridLevel));
			vtkFunc.resize_values(vec.size());
			vtkFunc.assign(vec);
			VTKOutput<dim> out;
			out.print(filename, vtkFunc, m_printConsistent);
		} else
			UG_THROW("Cannot find grid level.");
	}

protected:
	//	underlying approximation space
	SmartPtr<approximation_space_type> m_spApproxSpace;

	//	flag if write to conn viewer
	bool bConnViewerOut;

	//	flag if write to vtk
	bool bVTKOut;

	//	flag if data shall be made consistent before printing
	bool m_printConsistent;

	//	current grid level
	GridLevel m_gridLevel, m_coarseGridLevel;

	//	mapping to smaller index set (-1 indicates dropping)
	const std::vector<int>* m_pvMapGlobalToPatch;
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

		typedef typename TGridFunction::template traits<VertexBase>::const_iterator const_iterator;

		//	loop all subsets
		for (int si = 0; si < u.num_subsets(); ++si) {
			//	loop all vertices
			const_iterator iter = u.template begin<VertexBase>(si);
			const_iterator iterEnd = u.template end<VertexBase>(si);

			for (; iter != iterEnd; ++iter) {
				//	get vertex
				VertexBase* v = *iter;

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
		UG_ASSERT(m_spPostProcess.valid(), "provide a post process with init");
		UG_ASSERT(m_pApproxSpace != NULL, "provide approximation space init");

		size_t numDoFs = 0;
		GridLevel gl;
		if (m_level == (size_t) -1) {
			numDoFs = m_pApproxSpace->dof_distribution(GridLevel(GridLevel::TOPLEVEL, GridLevel::SURFACE))->num_indices();
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
