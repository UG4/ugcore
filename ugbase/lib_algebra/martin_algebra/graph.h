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