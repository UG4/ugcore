/*
 * hlibpro.h
 *
 *  Created on: 10.02.2011
 *      Author: iheppner, based on the coupling of ug 3.8 with HLibpro introduced
 *              by mrupp and blemke.
 */

/* Remarks of M. Rupp:
   * DONE: Integration von 'bFound' in 'get_element(size_t i, size_t j, bool &bFound)'
     and 'get_element(size_t i, size_t j)', siehe Mail von M. Rupp, 03032011 - nur
     was soll dann mit bFound in letzterer Funktion geschehen?

   * DONE: In 'build_crs_matrix()' nach M. Rupp 'm_bCRSSorted = matrix_type::rows_sorted'
     ergaenzt - wo geht das wie ein?
     Unter Umständen gibt es nämlich später Matrizen, die nicht die Spaltenindices sortieren.
*/

#ifndef __H__LIB_ALGEBRA__HLIB_PRO_OPERATOR__
#define __H__LIB_ALGEBRA__HLIB_PRO_OPERATOR__
#include <iostream>
#include <sstream>

#include "lib_algebra/operator/interface/operator_inverse.h"

#ifdef UG_PARALLEL
	#include "lib_algebra/parallelization/parallelization.h"
#endif

/* include all hlib headers - only recommended for small problems acc. to
   'hlibpro-user.pdf, p. 17! TODO: maybe include only the needed hlib headers */
#include <hlib.hh>

// to make functions of C interface available here:
#include "hlib-c.hh" // this is file 'hlibpro-<version>/src/include/hlib-c.hh'

#define CHECK_INFO { if ( info != HLIB_NO_ERROR ) \
{ char buf[1024]; hlib_error_desc( buf, 1024 ); \
printf( "\n%s\n\n", buf ); /*UG_ASSERT(0);*/ } }

#define HLIB_OUTPUT_PATH "."   

// additions for profiling
#include "common/profiler/profiler.h"
#define PROFILE_HLIB
#ifdef PROFILE_HLIB
	#define HLIB_PROFILE_FUNC()			PROFILE_FUNC_GROUP("algebra hlib")
	#define HLIB_PROFILE_BEGIN(name)	PROFILE_BEGIN_GROUP(name, "algebra hlib")
	#define HLIB_PROFILE_END()			PROFILE_END()
#else
	#define HLIB_PROFILE_FUNC()
	#define HLIB_PROFILE_BEGIN(name)
	#define HLIB_PROFILE_END()
#endif
// additions for profiling - end

namespace ug{

class CRSMatrix
{
	public:
		CRSMatrix() : m_bCRSSorted(true), m_nRows(0), m_nCols(0), m_nNonZeroes(0)
		{
		}

	//	method to display a crs matrix
		template<typename matrix_type>
		void build_crs_matrix(const matrix_type& A)
		{
			size_t A_num_rows = A.num_rows(); // number of rows of blocks
			size_t A_num_cols = A.num_cols(); // number of cols of blocks

			m_bCRSSorted = matrix_type::rows_sorted; // TODO: usage of this information?

			if(!block_traits<typename matrix_type::value_type>::is_static)
			{
                UG_LOG ("ERROR in 'HLIBSolver::init_lu': dynamic blocks not yet implemented!\n");
				m_nRows = m_nCols = 0; // used to detect error by caller
                return; // error handling in this case is done by caller!
			}

			const size_t block_num_rows = block_traits<typename matrix_type::value_type>::static_num_rows; // number of rows in a block (i.e. number of unknowns/equations)
			const size_t block_num_cols = block_traits<typename matrix_type::value_type>::static_num_cols; // number of cols in a block (i.e. number of unknowns/equations)
			UG_ASSERT(block_num_cols == block_num_rows, "only square blocks supported");

			m_nRows = A_num_rows*block_num_rows; // total number of rows
			m_nCols = A_num_cols*block_num_cols; // total number of cols

			UG_LOG("'build_crs_matrix()': matrix to convert consists of "
				   << A_num_rows << " x " << A_num_cols << " blocks with size "
				   << block_num_rows << " x " << block_num_cols << "\n");

		// count non zeroes
			m_nNonZeroes = 0;
			for(size_t row = 0; row < A_num_rows; ++row)
			{
				for(typename matrix_type::const_row_iterator conn = A.begin_row(row); conn != A.end_row(row); ++conn)
				{
					//	get corresponding block
					const typename matrix_type::value_type &val = conn.value();
					
					//	loop block components
					for(size_t rowCmp = 0; rowCmp < GetRows(val); ++rowCmp)
						for(size_t colCmp = 0; colCmp < GetCols(val); ++colCmp)
						{
							if(BlockRef(val, rowCmp, colCmp) != 0.0)
								m_nNonZeroes++;
						}
				}
			}
			
		// resize crs arrays
			m_row_ptr.resize(m_nRows+1);
			m_col_ind.resize(m_nNonZeroes);
			m_val.resize(m_nNonZeroes);

		// fill crs matrix arrays (assumption: static size of blocks (s.a.))
			m_nNonZeroes = 0;
			size_t rowIndex = 0;
			number a_ij;
			for(size_t row = 0; row < A_num_rows; ++row) // loop over rows of blocks ...
			{
				// .. and rows *in* each block (unfortunately we have to use this order of loops ...)
				for(size_t rowCmp = 0; rowCmp < block_num_rows; ++rowCmp)
				{
					m_row_ptr[rowIndex] = m_nNonZeroes;
					for(typename matrix_type::const_row_iterator conn = A.begin_row(row); conn != A.end_row(row); ++conn)
					{	
						size_t colIndex = conn.index()*block_num_cols; // offset for block column

						//	get corresponding block
						const typename matrix_type::value_type &val = conn.value();
					
						//	loop block components
						for(size_t colCmp = 0; colCmp < block_num_cols; ++colCmp)
						{
							a_ij = BlockRef(val, rowCmp, colCmp);
							if(a_ij != 0.0) {
								m_col_ind[m_nNonZeroes] = colIndex + colCmp;
								m_val[m_nNonZeroes] = a_ij;
								m_nNonZeroes++;
							}
						}
					}
					rowIndex++;
				}
			}

if(0)
{
			for(size_t row = 0; row < A_num_rows; ++row) // loop over rows of blocks
			{
				m_row_ptr[rowIndex] = m_nNonZeroes;
				for(typename matrix_type::const_row_iterator conn = A.begin_row(row); conn != A.end_row(row); ++conn)
				{	
					size_t colIndex = conn.index()*block_num_cols;

					//	get corresponding block
					const typename matrix_type::value_type &val = conn.value();
					
					//	loop block components
					for(size_t rowCmp = 0; rowCmp < GetRows(val); ++rowCmp)
						for(size_t colCmp = 0; colCmp < GetCols(val); ++colCmp)
						{
							a_ij = BlockRef(val, rowCmp, colCmp);
							if(a_ij != 0.0) {
								size_t ccc = colIndex + colCmp;
								m_col_ind[m_nNonZeroes] = ccc;
								m_val[m_nNonZeroes] = a_ij;
								m_nNonZeroes++;
							}
						}
				}
				rowIndex++;
			}
}
			m_row_ptr[m_nRows] = m_nNonZeroes; /* NOTE: Man darf hier nicht nach Konvention (see e.g. http://netlib.org/linalg/html_templates/node91.html) "m_nNonZeroes + 1" setzen,
												  sonst Absturz in 'TArray.hh' (07072010ih):
												  in "(TArray) operator[]" at "include/base/TArray.hh:182"
												  Error: out-of-bound error in array */

			/* some statistics: */
			//if(m_hlibVerb >= 2) {
				size_t sumOfCRSMem, fullMem;
				size_t bytesize;

				std::cout << "'build_crs_matrix()':" << std::endl;
				std::cout << "      matrix with " << rowIndex << " rows and " << m_nNonZeroes << " non zeroes in CRS format created (sparsity = "
						  << m_nNonZeroes/(number)(m_nRows*m_nCols) << " (= nNonZeroes/(m_nRows x m_nCols)))" << std::endl;

				fullMem = A_num_rows * A_num_cols * sizeof(number);
				fullMem *= block_num_rows * block_num_cols;
				std::cout << "      memory for dense matrix:                   \t"
						  << std::setprecision(3) << std::setw(10) << ((number) fullMem)/(1024.0 * 1024.0) << " Mbyte "
						  << "(" << fullMem << " byte = " << A_num_rows << "x" << A_num_cols << " * "
						  << block_num_rows << "x" << block_num_cols << " * " << sizeof(number)
						  << " byte (sizeof(number)))" << std::endl;

				bytesize = m_nNonZeroes * sizeof(number);
				sumOfCRSMem = bytesize;
				std::cout << "      memory for matrix in crs format:   val:    \t"
						  << std::setprecision(3) << std::setw(10) << ((number) bytesize)/(1024.0 * 1024.0) << " Mbyte "
						  << "(" << bytesize << " byte = " << m_nNonZeroes << " nNonZeroes * " << sizeof(number)
						  << " byte (sizeof(number)))\n";
		
				bytesize = m_nNonZeroes * sizeof(int);
				sumOfCRSMem += bytesize;
				std::cout << "                                         col_ind:\t"
						  << std::setprecision(3) << std::setw(10) << ((number) bytesize)/(1024.0 * 1024.0) << " Mbyte "
						  << "(" << bytesize << " byte = " << m_nNonZeroes << " nNonZeroes * " << sizeof(int)
						  << " byte (sizeof(int)))\n";
		
				bytesize = (m_nRows+1) * sizeof(int);
				sumOfCRSMem += bytesize;
				std::cout << "                                         row_ptr:\t"
						  << std::setprecision(3) << std::setw(10) << ((number) bytesize)/(1024.0 * 1024.0) << " Mbyte "
						  << "(" << bytesize << " byte = (" << m_nRows << " rows + 1) * " << sizeof(int)
						  << " byte (sizeof(int)))\n";
		
				std::cout << "                                         Sum:    \t"
						  << std::setprecision(3) << std::setw(10) << ((number) sumOfCRSMem)/(1024.0 * 1024.0) << " Mbyte "
						  << "(" << sumOfCRSMem << " byte)\n";

				std::cout << "      compression ratio = "
						  << std::setprecision(2) << std::setw(6) << 100.0 * ((number) sumOfCRSMem) / ((number) (fullMem))
						   << " percent: " << ((number) sumOfCRSMem) / (1024.0 * 1024.0) << " Mbyte compared to "
						  << ((number) (fullMem))   / (1024.0 * 1024.0) << " Mbyte "
						  << "(" << sumOfCRSMem << " byte compared to " << fullMem << " byte)\n" << std::endl;

				//} /* end verbose */

		} /* end 'build_crs_matrix()' */

	//	method to get an element of an crs matrix
		number get_element(size_t i, size_t j, bool &bFound)
		{
			/* Scan the column indices of entries in row 'i' */
			for(size_t k = m_row_ptr.at(i); k < (size_t)m_row_ptr.at(i+1); ++k) /* indices 'k' of row 'i' */
			{
				/* if j is the column index of a non-zero, then return it: */
				if((size_t)m_col_ind.at(k) == j)
				{
					bFound = true;
					return (m_val.at(k));
				}

				/* if we have passed j, then entry (i, j) must be a zero: */
				if(m_bCRSSorted && (size_t)m_col_ind.at(k) > j)
				{
					bFound = false;
					return (0.0);
				}
			}

			/* if we have reached the end of the non-zeroes without finding 'j',
			   it must be the index of a trailing 0: */
			return (0.0);
		}

		number get_element(size_t i, size_t j)
		{
			bool bFound;		/* TODO: what to do with 'bFound'??? */
			return get_element(i, j, bFound);
		}

	//	method to display a crs matrix (of course only reasonable for rather small matrices ...)
		void print_crs_matrix() //(CRSMatrix& C)
		{
			size_t i;

			std::cout << "'print_crs_matrix()': Displaying CRS matrix with '" << m_nNonZeroes
					  << "' non zero elements in '"
					  << m_nRows << "' rows:" << std::endl;

			//std::cout << "                      fromVec     = '" << fromVec
			//		  << "', toVec     = '" << toVec << std::endl;

			std::cout << "row_ptr:" << std::endl;
			for (i = 0; i < m_nRows+1; ++i)
				std::cout << std::setw(4)    << "   i: " << i << "\t" << m_row_ptr.at(i) << std::endl;

			std::cout << "col_ind:" << std::endl;
			for (i = 0; i < m_nNonZeroes; ++i)
				std::cout << std::setw(4)    << "   i: " << i <<  "\t" << m_col_ind.at(i) << std::endl;
	
			std::cout << "val:" << std::endl;
			for (i = 0; i < m_nNonZeroes; ++i)
				std::cout << std::setw(4)    << "   i: " << i <<  "\t"
						  << std::scientific << m_val.at(i) << std::endl;

			return;
		}

	public: // make it public ...
		bool m_bCRSSorted;
		size_t m_nRows, m_nCols, m_nNonZeroes;

		// data for compressed row storage
		std::vector<number> m_val;      /* pointer to non zero entries of original matrix (traversed row wise) */
		//std::vector<size_t> m_row_ptr;  /* list of value indexes where each row starts         */
		//std::vector<size_t> m_col_ind;  /* column indices corresp. to 'val'                    */
		// HLIB needs 'int' ...
		std::vector<int> m_row_ptr;  /* list of value indexes where each row starts         */
		std::vector<int> m_col_ind;  /* column indices corresp. to 'val'                    */

}; /* end class 'CRSMatrix' */

template <typename TAlgebra> //template <typename TAlgebra, int dim> // for dimension dependent handling of coordinates
class HLIBSolver
	: public IMatrixOperatorInverse<  typename TAlgebra::matrix_type,
	  	  	  	  	  	  	  	  	  typename TAlgebra::vector_type>,
	public DebugWritingObject<TAlgebra>
{
	public:
	// 	Algebra type
		typedef TAlgebra algebra_type;

	// 	Vector type
		typedef typename TAlgebra::vector_type vector_type;

	// 	Matrix type
		typedef typename TAlgebra::matrix_type matrix_type;

	///	Base type
		typedef IMatrixOperatorInverse<matrix_type,vector_type> base_type;

	protected:
		using base_type::convergence_check;
		using DebugWritingObject<TAlgebra>::write_debug;

	public:
		HLIBSolver() :
			m_pOperator(NULL),
			m_bIsExecutable(true),
			m_nMin(20),
			m_hlib_accuracy_H(1.0e-4), m_hlib_accuracy_LU(1.0e-4),
			m_hlibVerb(0),
			m_bClusterAlgebraically(true), m_bClusterWithNestedDissection(true)
			{
				sprintf(m_PSbaseName, "%s/%s", HLIB_OUTPUT_PATH, "hlib");
				//strcpy(m_PSbaseName, "hlib");
			};

		virtual const char* name() const {return "HLIBSolver";}

	//	sets nmin
		void set_hlib_nmin(size_t nmin)
		{
			m_nMin = nmin;
		}

	//	set verbositiy level of HLIB
		void set_hlib_accuracy_H ( const number acc )
		{
			m_hlib_accuracy_H = acc;
		}

	//	set verbositiy level of HLIB
		void set_hlib_accuracy_LU ( const number acc )
		{
			m_hlib_accuracy_LU = acc;
		}

	//	set verbositiy level of HLIB
		void set_hlib_verbosity ( const int verb )
		{
			m_hlibVerb = verb;
			hlib_set_verbosity(m_hlibVerb);
		}

	//	set clustering method
		void set_clustering_method(const char *ctype, const char *cmode)
		{

			// set clustering type, i.e., purely algebraic (default) or based on geometry
			if(strcmp("algebraic", ctype) == 0)
			{
				m_bClusterAlgebraically = true;
			} else if(strcmp("geometric", ctype) == 0)
			{
				m_bClusterAlgebraically = false;
			} else
			{
                UG_LOG ("ERROR in 'set_clustering_method()': clustering type '"
						<< ctype << "' not known! "
						"Use 'set_clustering_method(algebraic)' instead!\n");
				m_bIsExecutable = false;
                return;
			}

			// set clustering mode, i.e., with nested dissection (default) or without
			if(strcmp("nested dissection", cmode) == 0)
			{
				m_bClusterWithNestedDissection = true;
			} else
			{
				m_bClusterWithNestedDissection = false;
			}

			if (m_bClusterAlgebraically == false)
            {
				m_bIsExecutable = false;
                UG_LOG ("ERROR in 'set_clustering_method()': geometrical clustering not yet implemented! "
						"Use 'set_clustering_method(algebraic)' (default) instead!\n");
            }

			return;
		}

	//	set base name of HLIB PostScript files 
		void set_ps_basename(const char* nameIn)
		{

			sprintf(m_PSbaseName, "%s/%s", HLIB_OUTPUT_PATH, nameIn);
			//strcpy(m_PSbaseName,nameIn);

			return;
		}

	//protected:
		bool init_hlib(const matrix_type &A)
		{
	//	Convert matrix into CRS format
			HLIB_PROFILE_BEGIN(HLIBBuildCRSMatrix)
			m_CRSMatrix.build_crs_matrix(A);
			HLIB_PROFILE_END(); // end 'HLIB_PROFILE_BEGIN(HLIBBuildCRSMatrix)'
			if(m_CRSMatrix.m_nRows == 0)
			{
				UG_LOG("ERROR: In 'HLIBSolver::init_hlib': Cannot build CRS matrix!\n");
				return false;
			}

	//	Importing system matrix into HLIB, convert it to H-Matrix and prepare
	//	LU decomposition
			UG_LOG("'HLIBSolver::init_hlib()': Importing system matrix into HLIB ...\n");

			int                  info;
			hlib_cluster_t       ct;
			hlib_blockcluster_t  bct;
			int                  isSymmetric = 0;
			char psname[32];

	//	Importing CRS matrix into HLIB - NOTE: This doesn't convert 'm_hlib_A'
	//	into an H-Matrix! This will be done by 'hlib_matrix_build_sparse()' later!
			HLIB_PROFILE_BEGIN(HLIBImportCRSMatrix);
			m_hlib_A = hlib_matrix_import_crs(m_CRSMatrix.m_nRows,
											  m_CRSMatrix.m_nRows,
											  m_CRSMatrix.m_nNonZeroes,
											  &(m_CRSMatrix.m_row_ptr[0]),
											  &(m_CRSMatrix.m_col_ind[0]),
											  &(m_CRSMatrix.m_val[0]),
											  isSymmetric, &info); CHECK_INFO;
			HLIB_PROFILE_END(); // end 'HLIB_PROFILE_BEGIN(HLIBImportCRSMatrix)'

			if(m_hlibVerb >= 2) {
				sprintf(psname, "%s_A_pattern.ps", m_PSbaseName);
				UG_LOG("   printing pattern of A in '" << psname << "'\n");
				hlib_matrix_print_ps( m_hlib_A, psname, HLIB_MATIO_PATTERN, &info ); CHECK_INFO;
			}
			if(m_hlibVerb >= 3) {
				// print each entry of matrix - with its algebraic value, not only pattern!
				sprintf(psname, "%s_A_entry.ps", m_PSbaseName);
				UG_LOG("   printing entries of A in '" << psname << "'\n");
				hlib_matrix_print_ps( m_hlib_A, psname,   HLIB_MATIO_ENTRY, &info ); CHECK_INFO;
			}

	//	1. build cluster tree
			UG_LOG ("'init_hlib()': building cluster tree ...\n");

			HLIB_PROFILE_BEGIN(HLIBBuildClusterTree);

			if (m_bClusterAlgebraically) /* 1.1 purely algebraic clustering */
			{
				strcat(m_PSbaseName, "_alg");

				if (m_bClusterWithNestedDissection) {
					// nested dissection, suitable for LU decomposition
					strcat(m_PSbaseName, "_nd");
					ct = hlib_ct_build_alg_nd( m_hlib_A, &info ); CHECK_INFO;
				} else {
					// standard algebraic partitioning functions employing
					// bisection techniques, suitable for matrix inversion
					ct = hlib_ct_build_alg ( m_hlib_A, &info ); CHECK_INFO;
				}
			} else // 1.2 clustering based on coordinates
			{
				strcat(m_PSbaseName, "_geom");

				// TODO: not yet implemented! Import von Koordinaten evt. ueber eine vom Skript aus zu bestimmenden GridFunction machen!?
				/* 
				HLIB::TGeomBSPPartStrat part_strat( coord.get() );
				HLIB::TBSPCTBuilder ct_builder( & part_strat, m_nMin );
				ct = ct_builder.build( coord.get() ); // TODO: access to coordinates?
				*/

				hlib_coord_t co; // = hlib_coord_import ( m_CRSMatrix.m_nRows, dim, coords, NULL, &info );	

				if (m_bClusterWithNestedDissection) {
					strcat(m_PSbaseName, "_nd");
					ct = hlib_ct_build_bsp_nd( co, m_hlib_A, &info ); CHECK_INFO;
				} else { // use "binary space partitioning" to decompose the index set
					ct = hlib_ct_build_bsp   ( co,           &info ); CHECK_INFO;
				}
			}
			HLIB_PROFILE_END(); // end 'HLIB_PROFILE_BEGIN(HLIBBuildClusterTree)'

			/* postscript output for cluster tree: */
			if(m_hlibVerb >= 2) {
				sprintf(psname, "%s_ct.ps", m_PSbaseName);
				UG_LOG("   printing cluster tree in '"<< psname << "'\n");
				hlib_ct_print_ps( ct, psname, & info ); CHECK_INFO;
			}


	//	2. build block cluster tree using admissibility condition
			HLIB_PROFILE_BEGIN(HLIBBuildBlockClusterTree);
			UG_LOG ("'init_hlib()': building block cluster tree ...\n");
			bct = hlib_bct_build( ct, ct, &info ); CHECK_INFO;
			HLIB_PROFILE_END(); // end 'HLIB_PROFILE_BEGIN(HLIBBuildBlockClusterTree)'

			int sparsity_constant = hlib_bct_csp ( bct, &info ); CHECK_INFO;
			UG_LOG("   sparsity constant = '" << sparsity_constant << "'\n");

			/* postscript output for block cluster tree: */
			if(m_hlibVerb >= 2) {
				sprintf(psname, "%s_bct.ps", m_PSbaseName);
				UG_LOG("   printing block cluster tree in '"<< psname << "'\n");
				hlib_bct_print_ps( bct, psname, & info ); CHECK_INFO;
			}

	//	3. build low rank approximations for admissible blocks 
			UG_LOG ("'init_hlib()': building low rank approximations for admissible blocks ...\n");
			hlib_acc_t hlib_h_accuracy = hlib_acc_fixed( m_hlib_accuracy_H );

			HLIB_PROFILE_BEGIN(HLIBBuildLowRankApproximation);

			// convert 'm_hlib_A' into H-Matrix (save as 'm_hlib_LU' to avoid additional variable):
			m_hlib_LU  = hlib_matrix_build_sparse( bct, m_hlib_A, hlib_h_accuracy, &info ); CHECK_INFO;

			HLIB_PROFILE_END(); // end 'HLIB_PROFILE_BEGIN(HLIBBuildLowRankApproximation)'

			if(m_hlibVerb >= 2) {
				sprintf(psname, "%s_HA_pattern.ps", m_PSbaseName);
				UG_LOG("   printing pattern of system matrix as H-Matrix in '"<< psname << "'\n");
				hlib_matrix_print_ps( m_hlib_LU, psname, HLIB_MATIO_PATTERN, &info ); //CHECK_INFO;
			}
			if(m_hlibVerb >= 3) { // print each entry of matrix - algebraic values!
				sprintf(psname, "%s_HA_entry.ps", m_PSbaseName);
				UG_LOG("   printing entries of system matrix as H-Matrix in '"<< psname << "'\n");
				hlib_matrix_print_ps( m_hlib_LU, psname,   HLIB_MATIO_ENTRY,   &info ); //CHECK_INFO;
			}
	
	//	4. create LU decomposition (this overwrites the given matrix by a matrix
	//	   representing the inverse of the factors, e.g. (LU )^{-1}
			hlib_acc_t hlib_lu_accuracy = hlib_acc_fixed( m_hlib_accuracy_LU );

			HLIB_PROFILE_BEGIN(HLIBCreateLUDecomposition);

			hlib_matrix_inv_lu( m_hlib_LU, hlib_lu_accuracy, &info ); CHECK_INFO;

			HLIB_PROFILE_END(); // end 'HLIB_PROFILE_BEGIN(HLIBCreateLUDecomposition)'

			if(m_hlibVerb >= 2) {
				sprintf(psname, "%s_LU_svd.ps", m_PSbaseName);
				UG_LOG("   printing singular value decomposition of LU in '"<< psname << "'\n");
				hlib_matrix_print_ps( m_hlib_LU, psname, HLIB_MATIO_SVD, &info ); //CHECK_INFO;

				number lu_size = hlib_matrix_bytesize( m_hlib_LU, &info ) / (1024.0*1024.0);
				UG_LOG("   size of LU factor = "<< lu_size << " Mbyte\n");

				sprintf(psname, "%s_LU_pattern.ps", m_PSbaseName);
				UG_LOG("   printing pattern of LU decomposition in '"<< psname << "'\n");
				hlib_matrix_print_ps( m_hlib_LU, psname, HLIB_MATIO_PATTERN, &info ); //CHECK_INFO;
			}
			if(m_hlibVerb >= 3) { // print each entry of matrix - algebraic values!
				sprintf(psname, "%s_LU_entry.ps", m_PSbaseName);
				UG_LOG("   printing entries of LU decomposition in '"<< psname << "'\n");
				hlib_matrix_print_ps( m_hlib_LU, psname,   HLIB_MATIO_ENTRY,   &info ); //CHECK_INFO;
			}

			// check accuracy || I - A * (LU)^{-1} ||
			hlib_real_t hlib_lu_inv_accuracy_achieved = hlib_matrix_norm_inv_approx( m_hlib_A, m_hlib_LU, &info );
			UG_LOG("'init_hlib()': LU inversion error = " << std::scientific
				   << hlib_lu_inv_accuracy_achieved << "\n"); CHECK_INFO;
			if(hlib_lu_inv_accuracy_achieved > m_hlib_accuracy_LU)
			{
				//m_bIsExecutable = false;
				UG_LOG("WARNING: In 'HLIBSolver::init_hlib': "
					   "Cannot achieve prescribed LU inversion accuracy " <<
					   m_hlib_accuracy_LU <<
					   " (only " << hlib_lu_inv_accuracy_achieved << " reached)!\n");
				//return false;
			}

			return (true);
		} /* end 'init_hlib()' */

	public:
	//	set operator L, that will be inverted
		virtual bool init(MatrixOperator<matrix_type, vector_type>& Op)
		{

			if(m_bIsExecutable == false)
			{
				UG_LOG("ERROR: In 'HLIBSolver::apply': Solver is not executable "
					   "due to initialisation (i.e. wrong clustering method)!\n");
				return false;
			}

		// 	init hlib
			int argc = 0; char ** argv;
			argv = NULL;
			HLIB::INIT( &argc, &argv );
			
		// 	init hlib - by calling C interface function
			int info = 0;
			hlib_init( &argc, &argv, &info ); CHECK_INFO;
			hlib_set_verbosity(2);		/* verbosity level (see hlibpro-c.pdf, p. 13;
										   '9' seems to be the highest level; 22072010ih) */

		// 	remember operator
			m_pOperator = &Op;

		//	get matrix of Operator
			m_pMatrix = &m_pOperator->get_matrix();

		//	check that matrix exist
			if(m_pMatrix == NULL)
			{
				UG_LOG("ERROR in 'HLIBSolver::init': No Matrix given,\n"); return false;
			}

			UG_LOG("'HLIBSolver::init()': matrix to invert has "
				   << m_pMatrix->num_rows() << " x "
				   << m_pMatrix->num_cols() << " blocks\n")

		//	Debug output of matrix to convert
			if(m_pDebugWriter != NULL)
			{
                UG_LOG ("'HLIBSolver::init()': write matrix ...\n");
				m_pDebugWriter->write_matrix(*m_pMatrix, "HLIBMatrix");
			}

		//	init HLIB operator - convert matrix to crs, import to HLIB
			if(!init_hlib(*m_pMatrix))
			{
				UG_LOG("ERROR in 'HLIBSolver::init': Cannot init HLIB.\n"); return false;
			}

		//	we're done
			return true;
		}

	//	finalisation of HLIB
		virtual bool finalise()
		{
			int info = 0;
			hlib_done(&info); // calls 'HLIB::DONE()';
			CHECK_INFO;

		//	we're done
			return true;
		}

	// 	Compute u = L^{-1} * f
		virtual bool apply(vector_type& u, const vector_type& f)
		{
			convergence_check()->set_symbol('%');
			convergence_check()->set_name("HLIB Solver");

			size_t A_num_rows;    // number of rows of blocks
			int info = 0;

			UG_LOG("'HLIBSolver::apply()'\n");

			// The common checks ...
#ifdef UG_PARALLEL
			if(!f.has_storage_type(PST_ADDITIVE))
			{
				UG_LOG("ERROR: In 'HLIBSolver::apply': Inadequate storage format of Vector f.\n");
				return false;
			}
			if(!u.has_storage_type(PST_CONSISTENT))
			{
				UG_LOG("ERROR: In 'HLIBSolver::apply':Inadequate storage format of Vector u.\n");
				return false;
			}
#endif
			UG_ASSERT(f.size() == m_pMatrix->num_rows(),	"Vector and Row sizes have to match!");
			UG_ASSERT(u.size() == m_pMatrix->num_cols(), "Vector and Column sizes have to match!");
			UG_ASSERT(f.size() == u.size(),                         "Vector sizes have to match!");


			////////////////////////////////////////////////////////////////////////////////////
			// importing vectors of rhs and solution into HLIB
			////////////////////////////////////////////////////////////////////////////////////
			A_num_rows = m_pMatrix->num_rows(); // number of blocks per rows

			// allocate arrays - 'm_CRSMatrix.m_nRows' contains factor block size for systems of pdes!
			std::vector<number> f_array; f_array.resize(m_CRSMatrix.m_nRows);
			std::vector<number> u_array; u_array.resize(m_CRSMatrix.m_nRows);

			// import f from UG --> HLIB
			UG_LOG("'HLIBSolver::apply()': import rhs vec into HLIB\n");
			// copy elements
			//for(size_t rowIndex = 0; rowIndex < A_num_rows; ++rowIndex)
			for(size_t rowIndex = 0, k=0; rowIndex < f.size(); ++rowIndex)
			{
				for(size_t j=0; j<GetSize(f[rowIndex]); j++)
					f_array[k++] = BlockRef(f[rowIndex],j);
			}

			// import into HLIB
			hlib_vector_t hlib_f;
			hlib_f = hlib_vector_import_array(&(f_array[0]), m_CRSMatrix.m_nRows, &info); CHECK_INFO;
	

			// import u : UG --> HLIB
			UG_LOG("'HLIBSolver::apply()': import sol vec into HLIB\n");
			// copy elements
			//for(size_t rowIndex = 0; rowIndex < A_num_rows; ++rowIndex)
			for(size_t rowIndex = 0, k=0; rowIndex < u.size(); ++rowIndex)
			{
				for(size_t j=0; j<GetSize(u[rowIndex]); j++)
					u_array[k++] = BlockRef(u[rowIndex],j);
			}

			// import into HLIB
			hlib_vector_t hlib_u;
			hlib_u = hlib_vector_import_array(&(u_array[0]), m_CRSMatrix.m_nRows, &info); CHECK_INFO;

			////////////////////////////////////////////////////////////////////////////////////
			// Execute Solver
			////////////////////////////////////////////////////////////////////////////////////
			UG_LOG("'HLIBSolver::apply()': Executing solver ...\n");	
	
			hlib_solve_info_t solve_info;
			solve_info.converged = 0;
			solve_info.steps = -1;

			HLIB_PROFILE_BEGIN(HLIBSolvePrecond);
			//hlib_solve(m_hlib_A, hlib_u, hlib_f, &solve_info, &info); CHECK_INFO;
			hlib_solve_precond( m_hlib_A, m_hlib_LU, hlib_u, hlib_f, &solve_info, &info); CHECK_INFO;
			HLIB_PROFILE_END(); // end 'HLIB_PROFILE_BEGIN(HLIBSolvePrecond)'
	
			if ( solve_info.converged ) {
				// print results
				UG_LOG("...'hlib_solve()' converged in " << solve_info.steps
					   << " steps, res norm: " << solve_info.res_norm
					   << ", conv. rate: " << solve_info.conv_rate << "\n\n");
			} else {
				UG_LOG("ERROR in 'HLIBSolver::apply': Cannot apply HLIB solver.\n");
				UG_LOG("'hlib_solve()' not converged in " << solve_info.steps
					   << "  steps, res: " << solve_info.res_norm
					   << ", conv. rate: " << solve_info.conv_rate << "\n\n");
				return false;
			}

			////////////////////////////////////////////////////////////////////////////////////
			// Copy solution back: HLIB --> UG
			////////////////////////////////////////////////////////////////////////////////////
			UG_LOG("'HLIBSolver::apply()': copy HLIB solution back to UG solution vector\n");
			//for(size_t rowIndex = 0; rowIndex < A_num_rows; ++rowIndex)
			for(size_t rowIndex = 0, k=0; rowIndex < u.size(); ++rowIndex)
			{
				for(size_t j=0; j<GetSize(u[rowIndex]); j++)
					BlockRef(u[rowIndex],j) = u_array[k++];
			}

#ifdef UG_PARALLEL
			// todo: we set solution to consistent here, but that is only true for
			//			serial case. Handle parallel case.
			u.set_storage_type(PST_CONSISTENT);
#endif

			// clean up
			/////////////////////////////////////////////////////////////////////////////////////
			UG_LOG("'HLIBSolver::apply()': free hlib vectors\n");
			// free HLIB Vectors
			hlib_vector_free(hlib_u, &info); CHECK_INFO;
			hlib_vector_free(hlib_f, &info); CHECK_INFO;

			// TODO: It would be nice to also free the space allocated for CRS matrices, arrays ...
			// TODO: This has to be called in something like a "postprocess" method,
			// otherwise the licence will be released and isn't available anymore in repeated call of 'apply()'!
			// Or 'init()' has to be called before every call of 'apply()' ...
/*
			if(!finalise())
			{
				UG_LOG("ERROR in 'HLIBSolver::init': Cannot finalise HLIB.\n"); return false;
			}
*/			
		//	we're done
			return true;
		}

	// 	Compute u = L^{-1} * f AND return defect f := f - L*u
		virtual bool apply_return_defect(vector_type& u, vector_type& f)
		{
		//	solve u
			if(!apply(u, f)) return false;

		//	update defect
			if(!m_pMatrix->matmul_minus(f, u))
			{
				UG_LOG("ERROR in 'HLIBSolver::apply_return_defect': Cannot apply matmul_minus.\n");
				return false;
			}

		//	we're done
			return true;
		}

		void check_crs_matrix()
		{
			UG_LOG("'HLIBSolver::check_crs_matrix':\n");
			m_CRSMatrix.print_crs_matrix();
		}

	// 	Destructor
		virtual ~HLIBSolver() {};

	protected:
		// Operator to invert
		MatrixOperator<matrix_type, vector_type>* m_pOperator;

		// matrix to invert
		matrix_type* m_pMatrix;

		// CRS matrix to import into HLIB
		CRSMatrix m_CRSMatrix;

	//	execute this solver only if 'init()' has set this true
		bool m_bIsExecutable;

	// HLIB related
		size_t m_nMin;
		number m_hlib_accuracy_H;
		number m_hlib_accuracy_LU;
		unsigned int m_hlibVerb;

		bool m_bClusterAlgebraically;
		bool m_bClusterWithNestedDissection;

		char m_PSbaseName[128]; // basename for naming of postscript output

		hlib_matrix_t        m_hlib_A;
		hlib_matrix_t        m_hlib_LU;

		HLIB::TAlgCTBuilder* m_pCT_builder;

};

} // end namespace ug

#endif /* __H__LIB_ALGEBRA__HLIB_PRO_OPERATOR__ */
