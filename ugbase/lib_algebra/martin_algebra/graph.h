/*
 *  graph.h
 *  flexamg
 *
 *  Created by Martin Rupp on 25.11.09.
 *  Copyright 2009 . All rights reserved.
 *
 */

class cgraph
	{
	public:
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
			
			if(never_happens) print();
		}
		
		~cgraph()
		{
			delete[] iNrOfConnections;
			delete[] consmem;
			delete[] conn;
		}
		
		void init(int from)
		{
			ASSERT2(conn[from] == NULL, "already inited");
			if(from == 0)
				conn[from] = consmem;
			else 
			{
				ASSERT2(conn[from-1] != NULL, "prev not inited!");
				conn[from] = conn[from-1]+iNrOfConnections[from-1];
			}
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
		
		void transpose()
		{
			int *newNrOfConnections = new int[size+1];
			memset(newNrOfConnections, 0, sizeof(int)*size);			
			for(int i=0; i<size; i++)
				for(int j=0; j<iNrOfConnections[i]; j++)
					newNrOfConnections[conn[i][j]]++;
			
			int **newConn = new int*[size];
			int *newconsmem = new int[iTotalNrOfConnections];
			newConn[0] = newconsmem;
			for(int i=1; i < size+1; i++)
			{
				newConn[i] = newConn[i-1] + newNrOfConnections[i-1];
				newNrOfConnections[i-1] = 0;
			}
			
			for(int i=0; i<size; i++)
				for(int j=0; j<iNrOfConnections[i]; j++)
				{
					int to = conn[i][j];
					newConn[ to ][ newNrOfConnections[to]++ ] = i;
				}
			
			delete[] iNrOfConnections;
			delete[] conn;
			delete[] consmem;
			
			iNrOfConnections = newNrOfConnections;
			conn = newConn;
			consmem = newconsmem;			
		}
		
		void setConnection(int from, int to)
		{
			ASSERT2(conn[from] != NULL, "forgot an init?");
			ASSERT2(conn[from+1] == NULL, "only from back to front!")
			conn[from][iNrOfConnections[from]++] = to;
			iTotalNrOfConnections++;
			ASSERT2(iTotalNrOfConnections < iMaxTotalNrOfConnections, "too many connections, increase nnz!");
		}	
		int size;
		int **conn;
		int *consmem;
		int *iNrOfConnections;
		int iTotalNrOfConnections;
		int iMaxTotalNrOfConnections;
	};