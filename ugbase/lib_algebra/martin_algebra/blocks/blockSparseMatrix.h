/*
 *  boundedSparse.h
 *  flexamg
 *
 *  Created by Martin Rupp on 16.12.09.
 *  Copyright 2009 G-CSC, University of Frankfurt. All rights reserved.
 *
 */

// THIS FILE IS NOT USED

// sparse
//		sparseMat< variableArray, variableArray >
// sparse with fixed nnz:
//		sparseMat< variableArray, fixedArray<5> >
// sparse with a fixed nr of rows:
//		sparseMat< fixedArray<5>, variableArray<5> >
// sparse with fixed nnz and fixed nr of rows:
//		sparseMat< fixedArray<5>, fixedArray<15> >


// todo: mat+mat, mat*vec, mat-mat, inverse_type << here we need cols. How to store ? perhaps in rowstart[rows+1]
///////////////////////////////////////////////////////////////////////////////////////
// connections are sorted < to < to

#pragma once
#if 0
template<typename row_storage_type, int rows_, typename nnz_storage_type, int nnzs_>
class sparseMat
{	
public:
	struct connection
	{
		int to;
		double val;
		friend ostream &operator << (ostream &out, const connection &c)
		{
			out << "[" << c.to << " -> " << c.val << "]" << endl;
			return out;
		}
	};
	
	// storage specific
	//--------------------
	typedef typename storage_traits<nnz_storage_type, connection, nnzs_, 0>::array_type connection_type;
	typedef typename storage_traits<row_storage_type, int, rows_+1, 0>::array_type rowstart_type;
	
	typedef sparseMat<row_storage_type, rows_, nnz_storage_type, nnzs_> matrix_type;
		
private:
	rowstart_type rowstart;
	connection_type conn;	
	
public:
	void setSize(int rows, int cols)
	{
		rowstart.setSize(rows+1);		
	}
	
	void setNNZ(int nnz)
	{
		conn.setSize(nnz);		
	}	
	
	
public:
	sparseMat() : rowstart(), conn()
	{
		for(int i=0; i < rowstart.size(); i++) rowstart[i] = 0;
	}
	
	int getRows() const
	{
		return rowstart.size()-1;
	}
	
	void addConnectionAt(int i, int row)
	{
		int nnzs = rowstart[getRows()];
		conn.setSize(nnzs+1, false);
		
		for(int j = nnzs-1; j > i; j--)
			conn[j+1] = conn[j];
		
		for(int j=r+1; j < rowstart.size(); j++)
			rowstart[j]++;		
	}

	int findEntry(int r, int c, bool bCreate=false)
	{
		int i;
		for(i=rowstart[r]; i < rowstart[r+1] && conn[i].to < c; i++) ;
		
		if(i < rowstart[r+1] && conn[i].to == c)
			return i;
		

		if(bCreate)
		{
			addConnecionAt(i, r);
			conn[i].to = c;
			conn[i].val = 0.0;
			return i;
		}
		else
			return -1;
	}
	
	
	
	void setAt(int r, int c, double a)
	{
		int i= findEntry(r, c, true);
		ASSERT2(i != -1, "entry not created ???");
		conn[i].val = a;
	}
	double getAt(int r, int c) const
	{
		int i= findEntry(r, c);
		if(i != -1) return conn[i].val;
		else return 0.0;
	}
	
	void operator += (const matrix_type &other)
	{
		for(int r=0; r< mat.getRows(); r++)
		{
			int c1 = rowstart[r]; c2 = other.rowstart[r];
			while(c1 < rowstart[r+1] && c2 < other.rowstart[r+1])
			{
				if(conn[c1].to == other.conn[c2].to)
				{
					conn[c1].val += other.conn[c2].val;
					c1++; c2++;
				}
				else if(conn[c1].to < other.conn[c2].to)
					c1++;
				else //if(conn[c1].to > other.conn[c2].to)
				{
					addConnectionAt(c1, r);
					conn[c1] = other.conn[c2];
					c1++; c2++;
				}
			}
		}		
	}
	
	matrix_type operator + (const matrix_type &other)
	{
		matrix_type mat(*this);
		mat += other;
		return mat;
	}
	
	void operator -= (const matrix_type &other)
	{
		for(int r=0; r< mat.getRows(); r++)
		{
			int c1 = rowstart[r]; c2 = other.rowstart[r];
			while(c1 < rowstart[r+1] && c2 < other.rowstart[r+1])
			{
				if(conn[c1].to == other.conn[c2].to)
				{
					conn[c1].val -= other.conn[c2].val;
					c1++; c2++;
				}
				else if(conn[c1].to < other.conn[c2].to)
					c1++;
				else //if(conn[c1].to > other.conn[c2].to)
				{
					addConnectionAt(c1, r);
					conn[c1].to = other.conn[c2].to;
					conn[c1].val = -other.conn[c2].val;
					c1++; c2++;
				}
			}
		}		
	}
	
	matrix_type operator - (const matrix_type &other)
	{
		matrix_type mat(*this);
		mat -= other;
		return mat;
	}
	
	friend ostream &operator << (ostream &out, const matrix_type &mat)
	{
		out << "sparse Mat, " << row_storage_type::getType() << " row storage, " << nnz_storage_type::getType() << " nnz storage." << endl;
		for(int r=0; r< mat.getRows(); r++)
		{
			cout << r << ": ";
			for(int i = mat.rowstart[r]; i < mat.rowstart[r+1]; i++)
				out << mat.conn[i] << " ";
		}
		return out;
	}
};
#endif