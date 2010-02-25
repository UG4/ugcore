/*
 *  graph.h
 *  flexamg
 *
 *  Created by Martin Rupp on 25.11.09.
 *  Copyright 2009 . All rights reserved.
 *
 */

#pragma once

//!
//! cgraph graph class
//! with fixed max nr of connections.
class cgraph
{
public:	
	//!
	//! construcotr
	//! @param n nr of nodes
	//! @param nnz max. nr of connections
	cgraph()
	{
		size = 0;
		conn = NULL; iNrOfConnections = NULL;
		consmem = NULL;
		iTotalNrOfConnections = 0;
		iMaxTotalNrOfConnections = 0;
	}
	
	cgraph(int n, int nnz)
	{
		size = n;
		conn = new int *[size+1];
		iNrOfConnections = new int[size];
		
		memset(conn, 0, sizeof(int*)*(size+1));
		memset(iNrOfConnections, 0, sizeof(int)*(size));
		
		iMaxTotalNrOfConnections = nnz;
		iTotalNrOfConnections = 0;
		
		consmem = new int[nnz];
		
		if(never_happens) { print(); pr(0);} 
	}
	
	//!
	//! destructor
	~cgraph()
	{
		if(iNrOfConnections) delete[] iNrOfConnections;
		if(consmem) delete[] consmem;
		if(conn) delete[] conn;
	}
	
	//!
	//! @param i init thi
	void init(int i)
	{
		ASSERT2(conn[i] == NULL, "already inited");
		if(i == 0)
			conn[i] = consmem;
		else 
		{
			ASSERT2(conn[i-1] != NULL, "prev not inited!");
			conn[i] = conn[i-1]+iNrOfConnections[i-1];
		}
	}
	
	void pr(int i)
	{
		cout << "graph row " << i << ":" << endl;
		for(int j=0; j<iNrOfConnections[i]; j++)
			cout << conn[i][j] << " "; cout << endl;
		cout.flush();
	}
	
	void print()
	{
		cout << "============= graph ================ " << endl;
		for(int i=0; i < size; i++)
		{
			cout << i << ": ";
			for(int j=0; j<iNrOfConnections[i]; j++)
				cout << conn[i][j] << " ";
			cout << endl;
		}
	}
	
	void writeToFile(const char *filename, int level) const
	{
		fstream file(filename, ios::out);
		file << CONNECTION_VIEWER_VERSION << endl;
		file << flexamg_dimensions << endl;
		
		writePosToStream(file);	
		file << 1 << endl;
		for(int i=0; i < size; i++)
		{
			for(cRowIterator conn = beginRow(i); !conn.isEnd(); ++conn)
				file << GetOriginalIndex(level, i) << " " << GetOriginalIndex(level, conn()) << " " << 1 << endl;
		}
	}	
	
	void transpose()
	{
		cgraph G;
		G.createAsTransposeOf(*this);
		swap(iNrOfConnections, G.iNrOfConnections);
		swap(iTotalNrOfConnections, G.iTotalNrOfConnections);
		swap(conn, G.conn);
		swap(consmem, G.consmem);
		swap(size, G.size);
		swap(iTotalNrOfConnections, G.iTotalNrOfConnections);
	}
	
	void createAsTransposeOf(const cgraph &other)
	{
		//ASSERT(conn == NULL && consmem == NULL);
		size = other.size;
		iNrOfConnections = new int[size+1];
		memset(iNrOfConnections, 0, sizeof(int)*size);			
		iTotalNrOfConnections = 0;
		for(int i=0; i<size; i++)
		{
			for(int j=0; j<other.iNrOfConnections[i]; j++)
				iNrOfConnections[other.conn[i][j]]++;
			iTotalNrOfConnections += other.iNrOfConnections[i];
		}
		
		conn = new int*[size+1];
		consmem = new int[iTotalNrOfConnections];
		conn[0] = consmem;
		for(int i=1; i < size+1; i++)
		{
			conn[i] = conn[i-1] + iNrOfConnections[i-1];
			iNrOfConnections[i-1] = 0;
		}
		
		for(int i=0; i < size; i++)
			for(int j=0; j<other.iNrOfConnections[i]; j++)
			{
				int to = other.conn[i][j];
				conn[ to ][ iNrOfConnections[to]++ ] = i;
			}
	}
	
	
	class cRowIterator 
	{
	public:
		const cgraph &C;
		int m_position;
		int row;
	public:
		inline cRowIterator(const cgraph &C_, int row_) : C(C_)
		{
			row = row_; 
			rewind(); 			
		}
		inline cRowIterator(const cRowIterator &other) : C(other.C) { row = other.row; m_position = other.m_position; }
		
		inline int operator () ()const {return C.conn[row][m_position];}
		
		inline void operator ++() {	m_position++; }
		
		inline void rewind() { m_position = 0;}
		inline int getPos() const{	return m_position;}
		
		inline bool isEnd() const { return m_position >=C.getNrOfConnections(row); }
	};
	
	cRowIterator beginRow(int row) const
	{
		return cRowIterator(*this, row);
	}
	
	int getNrOfConnections(int row) const
	{
		return iNrOfConnections[row];
	}
	
	
	void setConnection(int from, int to)
	{
		ASSERT2(conn[from] != NULL, "forgot an init?");
		ASSERT2(conn[from+1] == NULL, "only from back to front!")
		conn[from][iNrOfConnections[from]++] = to;
		iTotalNrOfConnections++;
		ASSERT2(iTotalNrOfConnections < iMaxTotalNrOfConnections, "too many connections, increase nnz!");
	}	
	
	int getLength() { return size; }
	
//private:
	int size;
	int **conn;
	int *consmem;

	int *iNrOfConnections;
	int iTotalNrOfConnections;
	int iMaxTotalNrOfConnections;
};