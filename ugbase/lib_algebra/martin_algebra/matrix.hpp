/*
 *  matrix.hpp
 *  flexamg
 *
 *  Created by Martin Rupp on 04.11.09.
 *  Copyright 2009 . All rights reserved.
 *
 */


template<typename mat_type>
inline const mat_type matrix<mat_type>::getDiag(int i) const
{
	// evtl anders, da nicht jede Matrix diageintrag
	return cons[i][0].dValue;
}

template<typename mat_type>
inline mat_type &matrix<mat_type>::getDiag(int i)
{
	ASSERT2(cons[i][0].iIndex == i, *this << " first entry has to be diagonal");
	return cons[i][0].dValue;
}

template<typename mat_type>
inline bool matrix<mat_type>::isUnconnected(int i) const
{
	return iNrOfConnections[i] == 1;
}

/*
 template<typename Operator> 
 struct Expression<double, Operator, matrix>
 { 
 const double& ld; 
 const matrix& r; 
 inline Expression(const double& ld_, const matrix &r_) : ld(ld_),r(r_) {} 
 
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


template<typename mat_type, typename R>
Expression<matrix<mat_type>, Multiply_Operator<mat_type, typename R::vec_type>, R> operator*(const matrix<mat_type> &l,const XD<R> &r)
{ 
	return Expression<matrix<mat_type>, Multiply_Operator<mat_type, typename R::vec_type>, R> (l, r.cast()); 
}


// todo: prevent x = A * x; mit feld forbiddenDestination
template<typename mat_type, typename vec_type> 
class Expression<matrix<mat_type>, Multiply_Operator<mat_type, vec_type >, Vector<vec_type> > 
: public XD< Expression<matrix<mat_type>, Multiply_Operator<mat_type, vec_type >, Vector<vec_type> > >
{ 
public:	
	typedef typename Multiply_Operator<mat_type, vec_type>::ReturnType ReturnType;
	const matrix<mat_type>& l;
	const Vector<vec_type> & r; 
	inline Expression(const matrix<mat_type> &l_, const Vector<vec_type> &r_) : l(l_), r(r_) {} 
	
	inline ReturnType operator [] (int i) const
	{
		return l[i] * r;
	} 
	
	inline int getLength() const	{	return l.getLength();	}
	
	// print routines
	friend ostream &operator<<(ostream &output, const Expression<matrix<mat_type>, Multiply_Operator<mat_type, vec_type>, Vector<vec_type> >  &ex)
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
struct Expression<Vector, Divide_Operator, matrix<mat_type>::diagcomponent> 
{ 
	const Vector& l;
	const matrix<mat_type>::diagcomponent& r; 
	inline Expression(const Vector &l_, const matrix<mat_type>::diagcomponent &r_) : l(l_), r(r_) 
	{ ASSERT2(l.getLength() == r.getLength(), l << " has different length as " <<  r); } 
	
	inline double operator [] (int i) const
	{
		return l[i] / r[i];
	} 
	
	inline int getLength() const	{	return l.getLength();	}
	
	// print routines
	friend ostream &operator<<(ostream &output, const Expression<Vector, Divide_Operator, matrix<mat_type>::diagcomponent>  &ex)
	{
		output << "(" << ex.l << Divide_Operator::cTyp() << ex.r << ")"; 
		return output;
	}
	inline void printtype() const	{	cout << *this; }
}; 

Expression<Vector, Divide_Operator, matrix<mat_type>::diagcomponent> operator/(const Vector &l, const matrix<mat_type>::diagcomponent &r);

#ifdef SPECIALIZE_EXPRESSION_TEMPLATES
inline void Vector::operator = (const Expression<matrix, Multiply_Operator, Vector> ex)
{
	IF_PRINTLEVEL(5) cout << *this << " = " << ex << endl;
	const matrix &A = ex.l;
	const Vector &v = ex.r;
	for(int i=0; i<length; i++)
		values[i] = A[i] * v;
}
#endif
*/


template<typename mat_type>
inline int matrixrow<mat_type>::getConNr(int index) const
{
	for(int i=0; i< getNrOfConnections(); i++)
	{
		if(A->cons[row][i].iIndex == index)
			return i;
	}
	return -1;
}

// the != 0.0 is bad but we need this for restriction, since A.cons[i][0].iIndex = i.
template<typename mat_type>
template<typename vec_type>
inline vec_type matrixrow<mat_type>::operator *(const Vector<vec_type> &x) const
{
	vec_type d;
	d = 0.0;
	
	citerator it(*this);
	if(!it.isEnd() && (*it).dValue == 0.0) 
		++it;
	
	for(; !it.isEnd(); ++it)
		d += (*it).dValue * x[(*it).iIndex];
	return d;
}

template<typename mat_type>
inline const typename matrixrow<mat_type>::connection &matrixrow<mat_type>::operator [] (int i) const
{
	ASSERT2(i < A->iNrOfConnections[row] && i >= 0, *this << " has no connection nr. " << i);
	/*ASSERT(A->cons[row]+i < A->consmem+A->iMaxTotalNrOfConnections
	 && A->cons[row]+i >= A->consmem);*/
	//ASSERT(A->consmem.isMemInChunk(A->cons[row][i]));		
	return A->cons[row][i];
}

template<typename mat_type>
inline typename matrixrow<mat_type>::connection &matrixrow<mat_type>::operator [] (int i)
{
	ASSERT2(i < A->iNrOfConnections[row] && i >= 0, *this << " has no connection nr. " << i);
	/*ASSERT(A->cons[row]+i < A->consmem+A->iMaxTotalNrOfConnections
	 && A->cons[row]+i >= A->consmem);*/
	//ASSERT(A->consmem.isMemInChunk(A->cons[row][i]));
	
	return A->cons[row][i];
}	


////////////////////////////////////////////////////////////////////////////////////////////////////////////


// constructor for empty matrix
template<typename mat_type>
matrix<mat_type>::matrix(const char *_name)
{
	length = iTotalNrOfConnections = 0;
	iFragmentedMem = 0;
	cons = NULL;
	consmem = NULL; consmemsize = 0;
	iTotalNrOfConnections = 0;
	iNrOfConnections = NULL;
	fromlevel = tolevel = -1;
	name = _name;
	bandwidth = 0;
}	

// constructor for matrix with length
template<typename mat_type>
matrix<mat_type>::matrix (int _length, const char *_name)
{
	iFragmentedMem = 0;
	fromlevel = 0;
	tolevel =0;
	length = 0;
	name = _name;
	fromlevel = tolevel = 0;
	bandwidth = 0;
	consmemsize = 0;
	create(_length);
}

// destructor
template<typename mat_type>
matrix<mat_type>::~matrix()
{
	for(int i=0; i<length; i++)
		saveSetConnections(i, NULL);
	
	delete [] iNrOfConnections;
	if(consmem)	delete [] consmem;
	delete [] cons;
}

// create: used to create the matrix
template<typename mat_type>
void matrix<mat_type>::create(int _length)
{
	ASSERT2(length == 0, *this << " not empty, has length " << length);
	
	length = _length;
	cons = new connection*[length+1];
	memset(cons, 0, sizeof(connection*)*(length+1));
	iNrOfConnections = new int[length];
	memset(iNrOfConnections, 0, sizeof(int)*length);
	
	iTotalNrOfConnections = 0;
	bandwidth = 0;
}


/*void matrix<mat_type>::recreateWithMaxNrOfConnections(int newMax) const
 {
 // create new cons Memory
 connection *consmemNew = new connection[newMax];
 
 // adjust pointers
 int diff = consmemNew - consmem;
 for(int i=0; i<length; i++)
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

template<typename mat_type>
void matrix<mat_type>::saveSetConnections(int row, connection *mem) const
{	
	if(cons[row] != NULL && cons[row] < consmem || cons[row] > consmem + consmemsize)
		delete[] cons[row];
	cons[row] = mem;
}

template<typename mat_type>
void matrix<mat_type>::defrag()
{
	iTotalNrOfConnections=0;
	for(int i=0; i<length; i++)
		iTotalNrOfConnections+=iNrOfConnections[i];
	
	connection *consmemNew = new connection[iTotalNrOfConnections+3];
	connection *p= consmemNew;
	for(int i=0; i<length; i++)
	{
		memcpy(p, cons[i], sizeof(connection)* iNrOfConnections[i]);
		saveSetConnections(i, p);
		p += iNrOfConnections[i];
	}
	delete[] consmem;
	consmem = consmemNew;
	iFragmentedMem = 0;
	consmemsize = iTotalNrOfConnections;
}

// getrow: returns a matrixrow object.
// can be used as A[i].print() or A[i]*x
template<typename mat_type>
const matrixrow<mat_type> matrix<mat_type>::getrow(int i) const
{	
	const matrixrow<mat_type> r(this, i);
	return r;
}

template<typename mat_type>
matrixrow<mat_type> matrix<mat_type>::getrow(int i)
{
	matrixrow<mat_type> r(this, i);
	
	return r;
}

template<typename mat_type>
matrixrow<mat_type> matrix<mat_type>::operator [] (int i) 
{		
	return getrow(i);
}

template<typename mat_type>
const matrixrow<mat_type> matrix<mat_type>::operator [] (int i) const
{		
	return getrow(i);
}

// eliminateDirichletValues
//----------------------------
// eliminates Dirichlet Values by putting them on the rhs b.
template<typename mat_type>
void matrix<mat_type>::eliminateDirichletValues(Vector<  matrix<mat_type>::vec_type> &b)
{
	for(int i=0; i<length; i++)
	{
		if(isUnconnected(i)) continue;
		matrixrow<mat_type> r(this, i);
		for(typename matrixrow<mat_type>::iterator conn(r); !conn.isEnd(); ++conn)
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
// createAsTransposeOf
//-----------------------
// write in a empty matrix the transpose matrix of B.
template<typename mat_type>
void matrix<mat_type>::createAsTransposeOf(const matrix &B, int length)
{
	create(length);
	
	// get length of each row
	for(int j=0; j<B.length; j++)
		for(typename matrixrow<mat_type>::citerator conn(B[j]); !conn.isEnd(); ++conn)
			if((*conn).dValue == 0) continue;
			else iNrOfConnections[(*conn).iIndex]++;
	
	int newTotal = 0;
	for(int i=0; i < length; i++)
		newTotal += iNrOfConnections[i];
	
	
	
	int *nr = new int[length];
	// init matrix data structure
	consmem = new connection[newTotal];
	consmemsize = newTotal;
	
	connection *p = consmem;
	for(int i=0; i < length; i++)
	{
		nr[i] = 0;
		cons[i] = p;
		p += iNrOfConnections[i];	
	}
	
	iTotalNrOfConnections = newTotal;
	iFragmentedMem = 0;
	
	bandwidth = 0;
	// write values
	for(int i=0; i<B.length; i++)
		for(typename matrixrow<mat_type>::citerator conn(B[i]); !conn.isEnd(); ++conn)
		{
			if((*conn).dValue == 0) continue;
			int ndx = (*conn).iIndex;
			ASSERT2(ndx >= 0 && ndx < length, "connection " << (*conn) << " of matrix " << B << " out of range 0.." << length);
			
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

// writeToFile: in somewhat matrix-market format:
// length \n then for each connection: from to value
template<typename mat_type>
void matrix<mat_type>::writeToFile(const char *filename) const
{
	fstream file(filename, ios::out);
	writePosToStream(file);
	file << 1 << endl;
	for(int i=0; i<length; i++)
	{
		for(typename matrixrow<mat_type>::citerator conn(getrow(i)); !conn.isEnd(); ++conn)
			if((*conn).dValue != 0.0)
				file << GetOriginalIndex(tolevel, i) << " " << GetOriginalIndex(fromlevel, (*conn).iIndex) << " " << ((*conn).dValue) <<		endl;
	}
}	

// print to console whole matrix
template<typename mat_type>
void matrix<mat_type>::print(const char * const text) const
{
	cout << endl << "================ " << name;
	if(text) cout << " == " << text;
	cout << " == fromlevel: " << fromlevel << " tolevel: " << tolevel << " length: " << length << " =================" << endl;
	for(int i=0; i<length; i++)
		getrow(i).print();				
}
template<typename mat_type>
void matrix<mat_type>::p()
{
	print(NULL);
}

template<typename mat_type>
void matrix<mat_type>::printtype() const
{ 
	cout << *this; 
}

template<typename mat_type>
void matrix<mat_type>::addSubmatrix	(submatrix<mat_type> &mat)
{
	int nr = mat.getNr();
	connection *c = new connection[nr];
	int nc;
	for(int i=0; i < nr; i++)
	{
		nc = 0;
		for(int j=0; j < nr; j++)
		{
			if(mat(i,j) != 0.0)
			{
				c[nc].iIndex = mat.getIndex(j);
				c[nc].dValue = mat(i, j);
				nc++;
			}
		}
		if(nc > 0)
			getrow( mat.getIndex(i) ).addMatrixRow(c, nc);
	}
	delete[] c;
}

template<typename mat_type>
int connection_compare(const void *a, const void *b)
{
	return ((typename matrix<mat_type>::connection*)a)->iIndex > ((typename matrix<mat_type>::connection*)b)->iIndex;
}

template<typename mat_type>
void sortConnections(typename matrix<mat_type>::connection *c, int nr, int row)
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
		mergesort(c+1, nr-1, sizeof(typename matrix<mat_type>::connection), connection_compare<mat_type>);
}


//						MATRIXROW
////////////////////////////////////////////////////////////////

// creation functions
/*
 void matrixrow<mat_type>::removezeros()
 {
 ASSERT2(A->cons[row+1] == NULL, *A << ": next row has to be uninitialized"); // nur von vorne nach hinten hinzufügen
 
 connection* con = A->cons[row];
 int nr = A->iNrOfConnections[row]-1;
 while(nr > 0 && con[nr].dValue == 0) nr--; // diagonaleintrag behalten
 for(int i=1; i<nr; i++)
 if(con[i].dValue == 0)
 {
 con[i] = con[nr--];
 while(nr > 0 && con[nr].dValue == 0) nr--;
 }
 
 A->iNrOfConnections[row] = nr+1;				
 }*/

template<typename mat_type>
void matrixrow<mat_type>::setMatrixRow(connection *c, int nr)
{
	sortConnections<mat_type>(c, nr, row);
	connection *n;
	if(c[0].iIndex != row)
	{
		n = new connection[nr+1];		
		n[0].iIndex = row; n[0].dValue = 0.0;
		memcpy(n+1, c, nr*sizeof(connection));
		nr++;
	}
	else
	{
		n = new connection[nr];
		memcpy(n, c, nr*sizeof(connection));
	}
	for(int i=0; i<nr; i++)
		A->bandwidth = max(A->bandwidth, abs(c[i].iIndex - row));
	A->iFragmentedMem += nr;
	
	A->saveSetConnections(row, n);
	A->iTotalNrOfConnections += nr - A->iNrOfConnections[row];
	A->iNrOfConnections[row] = nr;
}

template<typename mat_type>
void matrixrow<mat_type>::addMatrixRow(connection *c, int nr)
{	
	connection *old = A->cons[row];	
	if(old == NULL)
	{
		setMatrixRow(c, nr);
		return;
	}

	int iNrOfConnections = A->iNrOfConnections[row];
	// sort the connections
	sortConnections<mat_type>(c, nr, row);	

	int ic, skipped=0, iold=1;
	if(c[0].iIndex == row) 
	{
		old[0].dValue += c[0].dValue;
		ic = 1; 		
	}
	else
		ic = 0;

	while(ic < nr && iold < iNrOfConnections)
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

	int iNewSize = iNrOfConnections + skipped;
	connection *n = new connection[iNewSize];
	A->iFragmentedMem += iNewSize;
	n[0] = old[0];
	int i=1;
	iold=1;
	if(c[0].iIndex == row) 
	ic = 1; 		
	else
	ic = 0;

	while(ic < nr || iold < iNrOfConnections)
	{
		if(iold >= iNrOfConnections || c[ic].iIndex < old[iold].iIndex)
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
	 
	A->saveSetConnections(row, n);
	A->iTotalNrOfConnections += iNewSize - A->iNrOfConnections[row];
	A->iNrOfConnections[row] = iNewSize;

}


// print
//----------
template<typename mat_type>
void matrixrow<mat_type>::print() const
{
	cout << row << " [" << GetOriginalIndex(A->tolevel, row) << "]: ";
	for(int i=0; i < getNrOfConnections(); i++)
	{
		connection &c = A->cons[row][i];
		if(c.dValue == 0.0) continue;
		cout << " ";
		cout << "(" << c.iIndex << "[" << GetOriginalIndex(A->fromlevel, c.iIndex) << "]-> " << c.dValue << ")";			
	}
	//cout << " SUM: " << sum() << endl;
	cout << endl;
}

template<typename mat_type>
void matrixrow<mat_type>::p()
{
	print();
}

template<typename mat_type>
void matrixrow<mat_type>::printtype() const 
{
	cout << *this;
}