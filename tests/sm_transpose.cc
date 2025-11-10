
#include "common/util/trace.h"
#include "lib_algebra/cpu_algebra/sparsematrix_impl.h"
#include "lib_algebra/graph_interface/sparsematrix_boost.h"

#include "common/log.cpp" // ?
#include "common/debug_id.cpp" // ?
#include "common/assert.cpp" // ?
#include "common/util/crc32.cpp" // ?
#include "common/util/ostream_buffer_splitter.cpp" // ?
#include "common/util/string_util.cpp" // ?

#include "boost/graph/graph_utility.hpp"
#include "common/stopwatch.h"

// sparse_matrix transpose test

template<class T>
void next_nonzero(T& it, T const& end)
{
	while(it != end){
		if(it.value()){
			break;
		}else{
			++it;
		}
	}
}

// operator==?
template<class T>
bool is_equal(ug::SparseMatrix<T> const& a, ug::SparseMatrix<T> const& b)
{
	typedef ug::SparseMatrix<T> M;
	typedef typename M::const_row_iterator const_row_iterator;
	if(&a == &b){ untested();
		return true;
	}else if(a.num_rows() != b.num_rows()){ untested();
		return false;
	}else if(a.num_cols() != b.num_cols()){ untested();
		return false;
	}

	for(int r=0; r<int(a.num_rows()); ++r){
		const_row_iterator ita = a.begin_row(r);
		const_row_iterator ea = a.end_row(r);
		const_row_iterator itb = b.begin_row(r);
		const_row_iterator eb = b.end_row(r);

		next_nonzero(ita, ea);
		next_nonzero(itb, eb);
		while(ita != ea && itb != eb){
			if(ita.value()!=itb.value()){ untested();
				return false;
			}else if(ita.index()!=itb.index()){ untested();
				return false;
			}else{
			}


			++ita; next_nonzero(ita, ea);
			++itb; next_nonzero(itb, eb);
		}

		if(ita!=ea){ untested();
			return false;
		}else if(itb!=eb){ untested();
			return false;
		}else{
		}
	}

	return true;
}

void test0(int N=7, int M=10)
{
	typedef ug::SparseMatrix<double> T;
	T A, B, C;

	A.resize_and_clear(N, M);

	for(int i=1; (i<N && i<M); ++i){
		A(i, i) = 2.;
		A(i, i-1) = 1.;
		A(i-1, i) = 1.;
	}

	if(N>6 && M>6){
		A(5, 5) = 0.;
		A(5, 0) = 0.;
		A(1, 6) = 1.;
	}else{
	}

	std::cout << A << "\n";
//	boost::print_graph(A);
	std::cout << "====\n";

	double t = -ug::get_clock_s();
	B.set_as_transpose_of(A);
	t += ug::get_clock_s();
	std::cout << B << " nnz " << B.total_num_connections() << " time " << t<< "\n";

	t = -ug::get_clock_s();
	C.set_as_transpose_of2(A);
	t += ug::get_clock_s();
	std::cout << C << " nnz " << C.total_num_connections() << " time " << t<< "\n";

//	assert(B==C);
	assert(is_equal(B, C));
	assert(is_equal(C, B));
}

void test3(int N=7)
{
	typedef ug::SparseMatrix<double> T;
	T A, B, C;

	A.resize_and_clear(N, N);

	for(int i=0; i<N; ++i){
		for(int j=0; j<N; ++j){
			A(i, j) = double(i+j);
		}
	}

	std::cout << A << "\n";
	// boost::print_graph(A);
	std::cout << "====\n";

	double t;

	t = -ug::get_clock_s();
	C.set_as_transpose_of2(A);
	t += ug::get_clock_s();
	std::cout << C << " nnz " << C.total_num_connections() << " time " << t << "\n";

	t = -ug::get_clock_s();
	B.set_as_transpose_of(A);
	t += ug::get_clock_s();
	std::cout << B << " nnz " << B.total_num_connections() << " time " << t << "\n";



//	assert(B==C);
	assert(is_equal(B, C));
	assert(is_equal(C, B));
}

int main()
{
	std::cout << "== test0\n";
	test0(2, 2);
	std::cout << "== test0 10000 10000\n";
	test0(10000, 10000);
	std::cout << "== test1\n";
	test0();
	std::cout << "== test2\n";
	test0(10, 7);
	std::cout << "== test3\n";
	test3(10);
	std::cout << "== test3 300\n";
	test3(300);
}
