/*
 *  SparseMatrix.hpp
 *  flexamg
 *
 *  Created by Martin Rupp on 04.11.09.
 *  Copyright 2009 . All rights reserved.
 *
 */

// construction etc.
//---------------------------------------------------

// constructor for empty SparseMatrix
template<typename entry_type>
SparseMatrix<entry_type>::SparseMatrix()
{
	name = "?";
	cols = rows = iTotalNrOfConnections = 0;
	iFragmentedMem = 0;
	cons = NULL;
	consmem = NULL; consmemsize = 0;
	iTotalNrOfConnections = 0;
	iNrOfConnections = NULL;
	fromlevel = tolevel = -1;
	bandwidth = 0;
	
	estimatedRowSize = 0;
	iMaxNrOfConnections = NULL;
	
	if(never_happens) 
	{
		print(); // force creation of this rountines for gdb.
		p();
		pr(0);
	}
}	

// destructor
template<typename entry_type>
SparseMatrix<entry_type>::~SparseMatrix()
{
	for(int i=0; i < rows; i++)
		safeSetConnections(i, NULL);
	
	delete [] iNrOfConnections;
	if(consmem)	delete [] consmem;
	delete [] cons;
}

// create: used to create the SparseMatrix
template<typename entry_type>
void SparseMatrix<entry_type>::create(int _rows, int _cols)
{
	ASSERT2(rows == 0 && cols == 0, *this << " not empty.");
	
	rows = _rows;
	cols = _cols;
	
	cons = new connection*[rows+1];
	memset(cons, 0, sizeof(connection*)*(rows+1));

	iNrOfConnections = new int[rows];
	memset(iNrOfConnections, 0, sizeof(int)*rows);
	
	iMaxNrOfConnections = new int[rows];
	memset(iMaxNrOfConnections, 0, sizeof(int)*rows);
	
	iTotalNrOfConnections = 0;
	bandwidth = 0;
}


/*void SparseMatrix<entry_type>::recreateWithMaxNrOfConnections(int newMax) const
 {
 // create new cons Memory
 connection *consmemNew = new connection[newMax];
 
 // adjust pointers
 int diff = consmemNew - consmem;
 for(int i=0; i<rows; i++)
 {
 if(cons[i] != NULL) 
 cons[i] += diff;
 }
 
 // copy, delete old and swap
 memcpy(consmemNew, consmem, sizeof(connection)*iTotalNrOfConnections);
 delete[] consmem;
 consmem = consmemNew;
 iMaxTotalNrOfConnections = newMax;
 }*/



template<typename entry_type>
inline const entry_type SparseMatrix<entry_type>::getDiag(int i) const
{
	ASSERT2(cons[i][0].iIndex == i, *this << " first entry has to be diagonal");
	// evtl anders, da nicht jede Matrix diageintrag
	return cons[i][0].dValue;
}

template<typename entry_type>
inline entry_type &SparseMatrix<entry_type>::getDiag(int i)
{
	ASSERT2(cons[i][0].iIndex == i, *this << " first entry has to be diagonal");
	return cons[i][0].dValue;
}

template<typename entry_type>
inline bool SparseMatrix<entry_type>::isUnconnected(int i) const
{
	ASSERT2(i < rows && i >= 0, *this << ": " << i << " out of bounds.");
	return iNrOfConnections[i] == 1;
}





////////////////////////////////////////////////////////////////////////////////////////////////////////////



template<typename entry_type>
void SparseMatrix<entry_type>::safeSetConnections(int row, connection *mem) const
{	
	if(cons[row] != NULL && cons[row] < consmem || cons[row] > consmem + consmemsize)
		delete[] cons[row];
	cons[row] = mem;
}

//!
//!
//! defrag
template<typename entry_type>
void SparseMatrix<entry_type>::defrag()
{
	ASSERT2(0, "this function is broken");
	iTotalNrOfConnections=0;
	for(int i=0; i<rows; i++)
		iTotalNrOfConnections+=iNrOfConnections[i];
	
	connection *consmemNew = new connection[iTotalNrOfConnections+3];
	connection *p= consmemNew;
	for(int i=0; i<rows; i++)
	{
		for(int k=0; k<iNrOfConnections[i]; k++)
			swap(p[k], cons[i][k]);
		safeSetConnections(i, p);
		p += iNrOfConnections[i];
		iMaxNrOfConnections[i] = iNrOfConnections[i];
	}
	delete[] consmem;
	consmem = consmemNew;
	iFragmentedMem = 0;
	consmemsize = iTotalNrOfConnections;
}

//!
//! getrow: returns a matrixrow object.
template<typename entry_type>
const matrixrow<entry_type> SparseMatrix<entry_type>::getrow(int i) const
{	
	return matrixrow<entry_type> (*this, i);
}

template<typename entry_type>
const matrixrow<entry_type> SparseMatrix<entry_type>::operator [] (int i) const
{		
	return matrixrow<entry_type> (*this, i);
}

// eliminateDirichletValues
//----------------------------
// eliminates Dirichlet Values by putting them on the rhs b.
template<typename entry_type>
template<typename Vector_type>
void SparseMatrix<entry_type>::eliminateDirichletValues(Vector_type &b)
{
	for(int i=0; i<rows; i++)
	{
		if(isUnconnected(i)) continue;		
		for(rowIterator conn = beginRow(i); !conn.isEnd(); ++conn)
		{
			int conindex = (*conn).iIndex;
			if(isUnconnected(conindex))
			{
				b[i] -= (*conn).dValue*b[conindex];
				(*conn).dValue = 0;
			}
		}
	}
}

// setDirichletRow
//----------------------------
// sets the row to i,i = 1.0.
template<typename entry_type>
void SparseMatrix<entry_type>::setDirichletRow(int row)
{
	ASSERT2(row >= 0 && row < rows, *this << ": row " << row << " out of bounds.");
	if(iNrOfConnections[row] > 0)
		cons[row][0] = 1.0;
	else
	{
		connection *c = new connection[1];
		safeSetConnections(row, c);
	}		
	iNrOfConnections[row] = 1;
	iMaxNrOfConnections[row] = 1;
}

template<typename entry_type>
void SparseMatrix<entry_type>::setDirichletRows(int *pRows, int nrows)
{
	for(int i=0; i<nrows; i++)
		setDirichletRow(pRows[i]);
}

// createAsTransposeOf
//-----------------------
// write in a empty SparseMatrix the transpose SparseMatrix of B.
template<typename entry_type>
void SparseMatrix<entry_type>::createAsTransposeOf(const SparseMatrix &B)
{
	ASSERT1(B.cols > 0 && B.rows > 0);
	create(B.cols, B.rows);
	fromlevel = B.tolevel;
	tolevel = B.fromlevel;
	
	// get length of each row
	for(int j=0; j < B.rows; j++)
		for(cRowIterator conn = B.beginRow(j); !conn.isEnd(); ++conn)
			if((*conn).dValue == 0) continue;
			else iNrOfConnections[(*conn).iIndex]++;
	
	int newTotal = 0;
	for(int i=0; i < rows; i++)
		newTotal += iNrOfConnections[i];
	
	
	
	int *nr = new int[rows];
	// init SparseMatrix data structure
	consmem = new connection[newTotal];
	consmemsize = newTotal;
	
	connection *p = consmem;
	for(int i=0; i < rows; i++)
	{
		nr[i] = 0;
		cons[i] = p;
		p += iNrOfConnections[i];
		iMaxNrOfConnections[i] = iNrOfConnections[i];
	}
	
	iTotalNrOfConnections = newTotal;
	iFragmentedMem = 0;
	
	bandwidth = 0;
	// write values
	for(int i=0; i < B.rows; i++)
		for(cRowIterator conn = B.beginRow(i); !conn.isEnd(); ++conn)
		{
			if((*conn).dValue == 0) continue;
			int ndx = (*conn).iIndex;
			ASSERT2(ndx >= 0 && ndx < rows, "connection " << (*conn) << " of SparseMatrix " << B << " out of range 0.." << rows);
			
			int k= nr[ndx];
			
			ASSERT2(k>=0 && k<iNrOfConnections[ndx], "k = " << k << ". precalculated nr of Connections " << iNrOfConnections[ndx] << " wrong?");
			//ASSERT(cons[ndx] + k < consmem+iMaxTotalNrOfConnections);
			//ASSERT(cons[ndx] + k >= consmem);
			
			cons[ndx][k].dValue = (*conn).dValue;
			cons[ndx][k].iIndex = i;
			if(bandwidth < i-ndx) bandwidth = i-ndx;
			nr[ndx]++;
		}
	
	delete[] nr;
}

// output functions
#define CONNECTION_VIEWER_VERSION 1
#define FLEXAMG_DIMENSIONS 3
// writeToFile: in somewhat SparseMatrix-market format:
// length \n then for each connection: from to value
template<typename entry_type>
void SparseMatrix<entry_type>::writeToFile(const char *filename) const
{
	fstream file(filename, ios::out);
	file << CONNECTION_VIEWER_VERSION << endl;
	file << FLEXAMG_DIMENSIONS << endl;
	
	writePosToStream(file);	
	file << 1 << endl;
	for(int i=0; i < rows; i++)
	{
		for(cRowIterator conn = beginRow(i); !conn.isEnd(); ++conn)
			if((*conn).dValue != 0.0)
				file << GetOriginalIndex(tolevel, i) << " " << GetOriginalIndex(fromlevel, (*conn).iIndex) << " " << ((*conn).dValue) <<		endl;
	}
}	

// writeToFile: in somewhat SparseMatrix-market format:
// length \n then for each connection: from to value
/*template<typename entry_type>
 void SparseMatrix<entry_type>::writeToFile2(const char *filename) const
 {
 fstream file(filename, ios::out);
 
 int level = min(fromlevel, tolevel);
 int m = max(rows, cols);
 cout << m << endl;
 for(int i=0; i < m; i++)
 {		
 int index = GetOriginalIndex(level, i)
 postype p = GetPosForIndex(index);
 file << p.x << " " << p.y << endl;
 }
 file << 1 << endl;
 for(int i=0; i < rows; i++)
 {
 for(cRowIterator conn = beginRow(i); !conn.isEnd(); ++conn)
 if((*conn).dValue != 0.0)
 file << GetOriginalIndex(tolevel, i) << " " << GetOriginalIndex(fromlevel, (*conn).iIndex) << " " << ((*conn).dValue) <<		endl;
 }
 }	*/

// print to console whole SparseMatrix
template<typename entry_type>
void SparseMatrix<entry_type>::print(const char * const text) const
{
	// cout << endl << "================ " << name;
	if(text) cout << " == " << text;
	cout << " == fromlevel: " << fromlevel << " tolevel: " << tolevel << ", " << rows << "x" << cols << " =================" << endl;
	for(int i=0; i < rows; i++)
		getrow(i).print();				
}


// print
//----------
template<typename entry_type>
void SparseMatrix<entry_type>::printrow(int row) const
{
	cout << row << " [" << GetOriginalIndex(tolevel, row) << "]: ";
	for(int i=0; i < iNrOfConnections[row]; i++)
	{
		connection &c = cons[row][i];
		if(c.dValue == 0.0) continue;
		cout << " ";
		cout << "(" << c.iIndex << "[" << GetOriginalIndex(fromlevel, c.iIndex) << "]-> " << c.dValue << ")";			
	}
	//cout << " SUM: " << sum() << endl;
	cout << endl;
}


template<typename entry_type>
void SparseMatrix<entry_type>::p() const
{
	print(NULL);
}

template<typename entry_type>
void SparseMatrix<entry_type>::pr(int row) const
{
	printrow(row);
}

template<typename entry_type>
void SparseMatrix<entry_type>::printtype() const
{ 
	cout << *this; 
}

template<typename entry_type>
void SparseMatrix<entry_type>::add(const submatrix<entry_type> &mat)
{
	connection *c = new connection[mat.getCols()];
	int nc;
	for(int i=0; i < mat.getRows(); i++)
	{
		nc = 0;
		for(int j=0; j < mat.getCols(); j++)
		{
			if(mat(i,j) != 0.0)
			{
				c[nc].iIndex = mat.getColIndex(j);
				c[nc].dValue = mat(i, j);
				nc++;
			}
		}
		if(nc > 0)
			addMatrixRow(mat.getRowIndex(i), c, nc);
	}
	delete[] c;
}

template<typename entry_type>
void SparseMatrix<entry_type>::set(const submatrix<entry_type> &mat)
{
	connection *c = new connection[mat.getCols()];
	int nc;
	for(int i=0; i < mat.getRows(); i++)
	{
		nc = 0;
		for(int j=0; j < mat.getCols(); j++)
		{
			if(mat(i,j) != 0.0)
			{
				c[nc].iIndex = mat.getColIndex(j);
				c[nc].dValue = mat(i, j);
				nc++;
			}
		}
		if(nc > 0)
			setMatrixRow(mat.getRowIndex(i), c, nc);
	}
	delete[] c;
}


template<typename entry_type>
void SparseMatrix<entry_type>::get(submatrix<entry_type> &mat) const
{
	int nr = mat.getNr();
	for(int i=0; i < mat.getRows(); i++)
	{
		int iindex = mat.getRowIndex(i), j = 0;
		
		rowIterator conn = beginRow(iindex);
		
		// diagonal
		mat(i,i) = (*conn).dValue;
		++conn;
		
		while(j < mat.getCols() && !conn.isEnd())
		{
			int jindex = mat.getColIndex(j);
			int cindex = (*conn).iIndex;
			
			if(cindex < jindex) 
				++conn;
			else if(cindex > jindex) 
			{
				mat(i, j) = 0.0;
				++j;
			}
			else
			{
				mat(i,j) = (*conn).dValue;
				++conn; ++j;
			}
		}
	}
}


template<typename entry_type>
int connection_compare(const void *a, const void *b)
{
	return ((typename SparseMatrix<entry_type>::connection*)a)->iIndex > ((typename SparseMatrix<entry_type>::connection*)b)->iIndex;
}

template<typename entry_type>
void sortConnections(typename SparseMatrix<entry_type>::connection *c, int nr, int row)
{
	// search diag
	if(c[0].iIndex != row)
		for(int i=1; i<nr; i++)
		{
			if(c[i].iIndex == row)
			{
				swap(c[i], c[0]);
				break;
			}			
		}
	if(nr-1 > 0)
		mergesort(c+1, nr-1, sizeof(typename SparseMatrix<entry_type>::connection), connection_compare<entry_type>);
}


template<typename entry_type>
void SparseMatrix<entry_type>::setMatrixRow(int row, connection *c, int nr)
{
	sortConnections<entry_type>(c, nr, row);
	connection *n;
	if(c[0].iIndex != row)
	{
		n = new connection[nr+1];		
		n[0].iIndex = row; n[0].dValue = 0.0;
		for(int i=0; i<nr; i++)
			n[i+1] = c[i];
		nr++;
	}
	else
	{
		n = new connection[nr];
		for(int i=0; i<nr; i++)
			n[i] = c[i];
	}
	for(int i=0; i<nr; i++)
	{
		bandwidth = max(bandwidth, abs(n[i].iIndex - row));
		ASSERT2(n[i].iIndex >= 0 && n[i].iIndex < getRows(), *this << " cannot have connection " << n[i] << ".");
	}
	iFragmentedMem += nr;
	
	safeSetConnections(row, n);
	iTotalNrOfConnections += nr - iNrOfConnections[row];
	iNrOfConnections[row] = nr;
	iMaxNrOfConnections[row] = nr;
}

template<typename entry_type>
void SparseMatrix<entry_type>::addMatrixRow(int row, connection *c, int nr)
{	
	connection *old = cons[row];	
	if(old == NULL)
	{
		
		setMatrixRow(row, c, nr);
		return;
	}
	
	//IFDEBUG( for(int i=0; i<nr; i++) ASSERT2(c[i].iIndex >= 0 && c[i].iIndex < getRows(), *A << " cannot have connection " << c[i] << "."); )
	
	int oldNrOfConnections = iNrOfConnections[row];
	// sort the connections
	sortConnections<entry_type>(c, nr, row);	
	
	int ic, skipped=0, iold=1;
	if(c[0].iIndex == row) 
	{
		old[0].dValue += c[0].dValue;
		ic = 1; 		
	}
	else
		ic = 0;
	
	while(ic < nr && iold < oldNrOfConnections)
	{
		if(c[ic].iIndex < old[iold].iIndex)
		{
			skipped++;
			ic++;
		}
		else if(c[ic].iIndex > old[iold].iIndex)
			iold++;
		else // "="
		{
			old[iold].dValue += c[ic].dValue;
			ic++;
			iold++;
		}
	}
	skipped += nr - ic;
	if(skipped == 0)  // everything already done
		return;
	
	
	// else realloc
	
	int iNewSize = oldNrOfConnections + skipped;
	connection *n = new connection[iNewSize];
	iFragmentedMem += iNewSize;
	n[0] = old[0];
	int i=1;
	iold=1;
	if(c[0].iIndex == row) 
		ic = 1; 		
	else
		ic = 0;
	
	while(ic < nr || iold < oldNrOfConnections)
	{
		if(iold >= oldNrOfConnections || c[ic].iIndex < old[iold].iIndex)
		{
			n[i++] = c[ic];
			ic++;
		}
		else if(ic >= nr || c[ic].iIndex > old[iold].iIndex)
		{
			n[i++] = old[iold];
			iold++;
		}
		else // "="
		{
			n[i++] = old[iold];
			ic++;
			iold++;
		}
	}
	
	safeSetConnections(row, n);
	iTotalNrOfConnections += iNewSize - oldNrOfConnections;
	iNrOfConnections[row] = iNewSize;
	iMaxNrOfConnections[row] = iNewSize;
}

template<typename entry_type>
void SparseMatrix<entry_type>::removezeros(int row)
{
	connection* con = cons[row];
	int nr = iNrOfConnections[row];
	// search from back for zero entries
	while(nr > 0 && con[nr-1].dValue == 0) nr--; // diagonaleintrag behalten
	for(int i=1; i < nr; i++)
		if(con[i].dValue == 0)
		{
			nr--;
			con[i] = con[nr];
			while(nr > 0 && con[nr-1].dValue == 0) nr--;
		}
	
	iNrOfConnections[row] = nr+1;
	iMaxNrOfConnections[row] = nr+1;
}


// res = A.T() * x
template<typename entry_type>
template<typename Vector_type>
void SparseMatrix<entry_type>::applyTransposed(Vector_type &res, const Vector_type &x) const
{
	res = 0.0;
	for(int i=0; i<rows; i++)
	{
		cRowIterator conn = beginRow(i);
		/*if((*conn).dValue != 0.0)
		{
			ASSERT1((*conn).iIndex  == i);
			res[(*conn).iIndex] += (*conn).dValue * x[i];
		}
		for(; !conn.isEnd(); ++conn)
			res[(*conn).iIndex] += (*conn).dValue * x[i];*/
		for(; !conn.isEnd(); ++conn)
		{
			if((*conn).dValue != 0.0)
				res[(*conn).iIndex] += (*conn).dValue * x[i];
		}
	}
}

// res = A*x
template<typename entry_type>
template<typename Vector_type>
void SparseMatrix<entry_type>::apply(Vector_type &res, const Vector_type &x) const
{
	for(int i=0; i<rows; i++)
	{
		cRowIterator conn = beginRow(i);
		/*if((*conn).dValue != 0.0)
		 {
		 ASSERT1((*conn).iIndex  == i);
		 res[(*conn).iIndex] += (*conn).dValue * x[i];
		 }
		 for(; !conn.isEnd(); ++conn)
		 res[(*conn).iIndex] += (*conn).dValue * x[i];*/
		res[i] = 0.0;
		for(; !conn.isEnd(); ++conn)
		{
			if((*conn).dValue != 0.0)
				res[i] += (*conn).dValue * x[(*conn).iIndex];
		}
	}
}

// res = res - A*x
template<typename entry_type>
template<typename Vector_type>
void SparseMatrix<entry_type>::matmul_minus(Vector_type &res, const Vector_type &x) const
{
	for(int i=0; i<rows; i++)
	{
		cRowIterator conn = beginRow(i);
		for(; !conn.isEnd(); ++conn)
		{
			if((*conn).dValue != 0.0)
				res[i] -= (*conn).dValue * x[(*conn).iIndex];
		}
	}
}





/// Template Expressions
//-------------------------------------------------------------------------------------
/*
 template<typename Operator> 
 struct Expression<double, Operator, SparseMatrix>
 { 
 const double& ld; 
 const SparseMatrix& r; 
 inline Expression(const double& ld_, const SparseMatrix &r_) : ld(ld_),r(r_) {} 
 
 inline double operator [] (int i) const
 {
 return Operator::apply( ld, r[i].getDiag() );
 } 
 
 inline int getLength() const	{	return r.getLength();	}
 void printtype() const	
 {
 cout << "(ID double[" << ld << "]";Operator::printtype();	r.printtype(); 	cout << ")";
 }
 }; */


template<typename entry_type, typename Vector_type>
Expression<SparseMatrix<entry_type>, Multiply_Operator<entry_type, typename Vector_type::entry_type>, Vector_type> 
	operator*(const SparseMatrix<entry_type> &l, const XD< Vector_type > &r)
{ 
	return Expression<SparseMatrix<entry_type>, Multiply_Operator<entry_type, typename Vector_type::entry_type>, Vector_type> (l, r.cast()); 
}


// todo: prevent x = A * x; mit feld forbiddenDestination
template<typename entry_type, typename Vector_type> 
class Expression<SparseMatrix<entry_type>, Multiply_Operator<entry_type, typename Vector_type::entry_type >, Vector_type > 
: public XD< Expression<SparseMatrix<entry_type>, Multiply_Operator<entry_type, typename Vector_type::entry_type >, Vector_type > >
{ 
public:	
	typedef typename Multiply_Operator<entry_type, typename Vector_type::entry_type>::ReturnType ReturnType;
	const SparseMatrix<entry_type>& l;
	const Vector_type & r; 
	inline Expression(const SparseMatrix<entry_type> &l_, const Vector_type &r_) : l(l_), r(r_) {} 
	
	inline ReturnType operator [] (int i) const
	{
		return l[i] * r;
	} 
	
	inline void copyTo(ReturnType &d, int i) const
	{
		l[i].copyToMult(d, r);
	}
	
	inline void addTo(ReturnType &d, int i) const
	{
		l[i].addToMult(d, r);
	}
	
	inline void substractFrom(ReturnType &d, int i) const
	{
		l[i].substractFromMult(d, r);
	}	
	
	inline int getLength() const	{	return l.getLength();	}
	
	// print routines
	friend ostream &operator<<(ostream &output, 
							   const Expression<SparseMatrix<entry_type>, Multiply_Operator<entry_type, typename Vector_type::entry_type >, Vector_type >  &ex)
	{
		output << "(" << ex.l	<< "*" << ex.r << ")"; 
		return output;
	}
	inline void printtype() const	{	cout << *this; }
}; 

/*
// todo: prevent x = A * x; mit feld forbiddenDestination
// x = r / A.Diag();
/*
template<> 
struct Expression<Vector, Divide_Operator, SparseMatrix<entry_type>::diagcomponent> 
{ 
	const Vector& l;
	const SparseMatrix<entry_type>::diagcomponent& r; 
	inline Expression(const Vector &l_, const SparseMatrix<entry_type>::diagcomponent &r_) : l(l_), r(r_) 
	{ ASSERT2(l.getLength() == r.getLength(), l << " has different length as " <<  r); } 
	
	inline double operator [] (int i) const
	{
		return l[i] / r[i];
	} 
	
	inline int getLength() const	{	return l.getLength();	}
	
	// print routines
	friend ostream &operator<<(ostream &output, const Expression<Vector, Divide_Operator, SparseMatrix<entry_type>::diagcomponent>  &ex)
	{
		output << "(" << ex.l << Divide_Operator::cTyp() << ex.r << ")"; 
		return output;
	}
	inline void printtype() const	{	cout << *this; }
}; 

Expression<Vector, Divide_Operator, SparseMatrix<entry_type>::diagcomponent> operator/(const Vector &l, const SparseMatrix<entry_type>::diagcomponent &r);
*/

		



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


template<typename entry_type, typename vec_type>
inline void multiplyCopyTo(vec_type &d, const matrixrow<entry_type> &r, const Vector<vec_type> &x)
{
	r.copyToMult(d, x);
}

template<typename entry_type, typename vec_type>
inline void multiplyAddTo(vec_type &d, const matrixrow<entry_type> &r, const Vector<vec_type> &x)
{
	r.addToMult(d, x);
}

template<typename entry_type, typename vec_type>
inline void multiplySubstractFrom(vec_type &d, const matrixrow<entry_type> &r, const Vector<vec_type> &x)
{
	r.substractFromMult(d, x);
}
