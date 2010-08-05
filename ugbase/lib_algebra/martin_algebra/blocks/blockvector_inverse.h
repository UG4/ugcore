
///////////////////////////////////////////////////////////////////////////////////////
//!
//! smallInverse<int n>
//! A class to hold a inverse of a smallMatrix<n>
//! implemented with LAPACKs LU-Decomposition dgetrf
//! (uses double[n*n] for LU and interchange[n] for pivoting
//! functions:
//! setAsInverseOf(const blockDenseMatrix<n> &mat) : init as inverse of mat
//! blockVector<n> * smallInverse<n> = smallInverse<n> * blockVector<n>
//! = A^{-1} b
template<typename storage_type, int rows_, int cols_>
class smallInverse
{
private: // storage
	typedef typename storage_traits<storage_type, double, rows_, cols_>::array2_type array2_type;
	typedef typename storage_traits<storage_type, __CLPK_integer, rows_, 0>::array_type interchange_array_type;

	typedef blockVector<double, storage_type, rows_> vector_type;

	array2_type densemat;
	interchange_array_type interchange;

public:
	inline int num_cols() const
	{
		return densemat.num_cols();
	}

	inline int num_rows() const
	{
		return densemat.num_rows();
	}

///
public:
	//! initializes this object as inverse of mat
	void setAsInverseOf(const blockDenseMatrix<double, storage_type, rows_, cols_> &mat)
	{
		UG_ASSERT(mat.num_rows() == mat.num_cols(), "");
		__CLPK_integer rows = mat.num_rows();
		__CLPK_integer cols = mat.num_cols();

		densemat.resize(rows, cols);
		for(int r=0; r < rows; r++)
			for(int c=0; c < cols; c++)
				densemat[r + c*rows] = mat(r, c);

		interchange.resize(rows);

		__CLPK_integer info = 0;

		getrf(&rows, &cols, &densemat[0], &rows, &interchange[0], &info);
		UG_ASSERT(info == 0, "info is " << info << ( info > 0 ? ": SparseMatrix singular in U(i,i)" : ": i-th argument had had illegal value"));
	}

	//! calc dest = mat^{-1} * vec
	void apply(double *dest, const vector_type &vec) const
	{
		UG_ASSERT(num_rows() == num_cols() && num_cols() == vec.size(), "");
		for(int i=0; i<vec.size(); i++)
			dest[i] = vec(i);
		char trans ='N';
		__CLPK_integer nrhs = 1;
		__CLPK_integer dim = num_rows();
		__CLPK_integer info = 0;

		getrs(&trans,
				&dim,
				&nrhs,
				&(*const_cast<array2_type*> (&densemat))(0,0),
				&dim,
				&(*const_cast<interchange_array_type*> (&interchange))[0],
				dest,
				&dim,
				&info);
	}

	vector_type operator * (const vector_type &vec) const
	{
		vector_type erg(vec.size());
		apply(&erg(0), vec);
		return erg;
	}

// temporary prevention
	//! dest = this*vec . use this to prevent temporary variables
	void assign_mult(vector_type &dest, const vector_type &vec) const
	{
		apply(&dest(0), vec);
	}
	//! dest += this*vec . use this to prevent temporary variables
	void add_mult(vector_type &dest, const vector_type &vec) const
	{
		// we need one temporary variable
		// keep static so it gets reallocated only once or twice
		static std::vector<double> erg;
		erg.resize(vec.size());

		apply(&erg[0], vec);
		for(int i=0; i<vec.size(); i++)
			dest(i) += erg[i];
	}
	//! dest -= this*vec . use this to prevent temporary variables
	void sub_mult(vector_type &dest, const vector_type &vec) const
	{
		// we need one temporary variable
		// keep static so it gets reallocated only once or twice
		static std::vector<double> erg;
		erg.resize(vec.size());

		apply(&erg[0], vec);
		for(int i=0; i<vec.size(); i++)
			dest(i) -= erg[i];
	}
};

template<typename storage_type, int rows, int cols>
blockVector<double, storage_type, rows> operator * (const blockVector<double, storage_type, rows> &vec, const smallInverse<storage_type, rows, cols> &mat)
{
	return mat * vec;
}




template<typename value_type, typename storage_type, int rows, int cols>
void Invert(blockDenseMatrix<value_type, storage_type, rows, cols> &A)
{
	typename storage_traits<storage_type, int, rows_, cols_>::array_type
		interchange;
	interchange.set_size(max(A.num_cols(), A.num_rows()));
	int info = 0;
	int rows = A.num_rows();
	int cols = A.num_cols();

	double *ptr = &getAt(0,0);

	UG_ASSERT(info == 0, "info is " << info << ( info > 0 ? ": SparseMatrix singular in U(i,i)" : ": i-th argument had had illegal value"));

	// calc work size
	double worksize; int iWorksize = -1;
	getri(&rows, ptr, &rows, &interchange[0], &worksize, &iWorksize, &info);
	UG_ASSERT(info == 0, "");
	iWorksize = worksize;

	std::vector<double> work;
	work.resize(iWorksize);
	getri(&rows, ptr, &rows, &interchange[0], &work[0], &iWorksize, &info);
	UG_ASSERT(info == 0, "");
}

template<>
void Invert(const blockDenseMatrix<double, fixedStorage, 1, 1> &mat )
{
	mat(0,0) = 1/mat(0,0);
}

template<>
void Invert(const blockDenseMatrix<double, fixedStorage, 2, 2> &mat )
{
	double invD = 1.0/(mat(0, 0)*mat(1, 1) - mat(0, 1)*mat(1, 0));
	UG_ASSERT(invD != 0.0, "");
	swap(A(0,0), A(1,1));
	swap(A(1,0), A(0,1));
	A(0,0) *= invD;
	A(0,1) *= -invD;
	A(1,0) *= -invD;
	A(1,1) *= invD;
}

 template<typename value_type, typename storage_type, int rows_>
 inline void blockVector<value_type, storage_type, rows_>::operator /= (const blockDenseMatrix<value_type, storage_type, rows_, rows_> &mat )
 {
	 smallInverse<storage_type, rows_, rows_> inv;
	 inv.setAsInverseOf(mat);

	 inv.apply(&values[0], *this);
 }

// TODO: do this with template specialisation
// template<typename value_type, typename storage_type>
// inline void blockVector<value_type, storage_type, 1>::operator /= (const blockDenseMatrix<value_type, storage_type, 1, 1> &mat )

template<>
inline void blockVector<double, fixedStorage, 1>::operator /= (const blockDenseMatrix<double, fixedStorage, 1, 1> &mat )
{
	getAt(0) = getAt(0) / mat(0, 0);
}


template<>
inline void blockVector<double, fixedStorage, 2>::operator /= (const blockDenseMatrix<double, fixedStorage, 2, 2> &mat )
{
	double invD = 1.0/(mat(0, 0)*mat(1, 1) - mat(0, 1)*mat(1, 0));
	UG_ASSERT(invD != 0.0, "");
	double a = invD*(getAt(0) * mat(1,1) - getAt(1) * mat(0,1));
	double b = invD*(getAt(1) * mat(0,0) - getAt(0) * mat(1,0));
	getAt(0) = a;
	getAt(1) = b;
}
