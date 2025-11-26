/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
 * Author: Daniel Jungblut
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#include <cmath>
#include "../ugmath_types.h"
#include "eigenvalues.h"

namespace ug
{

static void rot(number A[3][3], const number s, const number tau,
				const int i, const int j, const int k, const int l)
{
	
  number g = A[i][j];
  number h = A[k][l];
	
  A[i][j] = g - s * (h + g * tau);
  A[k][l] = h + s * (g - h * tau);												
}

// For symmetrische Matrizen:
bool CalculateEigenvalues(const ug::matrix33& mat, number& lambdaMinOut,
						number& lambdaMedOut, number& lambdaMaxOut,
						ug::vector3& evMinOut, ug::vector3& evMedOut,
						ug::vector3& evMaxOut)
{	
  
  number A[3][3];
  
  int i, j, ip, iq;
  number tresh, theta, tau, t, sm, s, h, g, c;
	
  int n = 3;
	
  number b[3];
  number d[3];
  number z[3];
	
  number V[3][3];

  // v wird als Einheitsmatrix initialisiert:
  // mat wird in A wird kopiert
  for(ip = 0; ip < 3; ip++) {
    for(iq = 0; iq < 3; iq++) {
	  A[ip][iq] = mat[ip][iq];
      V[ip][iq] = 0.0f;
    }
    V[ip][ip] = 1.0f;
  }
	
  // b und d werden mit der Diagonale von A initialisiert:
  // z wird mit 0 initialisiert:
  for(ip = 0; ip < 3; ip++) {
    b[ip] = d[ip] = A[ip][ip];
    z[ip] = 0;
  }
	
  //int nrot = 0;
	
  for(i = 1; i <= 50; i++) {
    sm = 0.0f;
    
    for(ip = 0; ip < 2; ip++) 
      for(iq = ip+1; iq < 3; iq ++)
        sm += fabs(A[ip][iq]);
		
    if(fabs(sm) <= 0.00001f) {
			
			// Berechnung fertig: 
			
      // Code zum sortieren der Eigenwerte einfügen:
			int max_index = 0;
			int med_index = 0;
			int min_index = 0;
			number max = 0.0f;
			
			for(int ii = 0; ii < 3; ii++)
				if(fabs(d[ii]) > max) {
					max_index = ii;
					max = fabs(d[ii]);
				}
			
			switch(max_index) {
					
			  case 0: if(fabs(d[1]) > fabs(d[2])) {med_index = 1; min_index = 2;} else {med_index = 2; min_index = 1;} break; 		
			  case 1: if(fabs(d[0]) > fabs(d[2])) {med_index = 0; min_index = 2;} else {med_index = 2; min_index = 0;} break; 
			  case 2: if(fabs(d[0]) > fabs(d[1])) {med_index = 0; min_index = 1;} else {med_index = 1; min_index = 0;} break; 
			  default: break;	
					
			}
			
			// Größter Eigenwert:
			lambdaMaxOut = d[max_index];
			evMaxOut[0] = V[0][max_index];
			evMaxOut[1] = V[1][max_index];
			evMaxOut[2] = V[2][max_index];
			
			// Mittlerer Eigenwert:
			lambdaMedOut = d[med_index];
			evMedOut[0] = V[0][med_index];
			evMedOut[1] = V[1][med_index];
			evMedOut[2] = V[2][med_index];
			
			// Kleinster Eigenwert:
			lambdaMinOut = d[min_index];
			evMinOut[0] = V[0][min_index];
			evMinOut[1] = V[1][min_index];
			evMinOut[2] = V[2][min_index];
			
			return true;   // Fertig!!!
			
    }
		
		
    if(i < 4)
      tresh = 0.2f * sm / 9.0f;
		
    else
      tresh = 0.0f;
		
    for(ip = 0; ip < (n-1); ip++) {
      for(iq = ip+1; iq < 3; iq++)  {
        g = 100.0f * fabs(A[ip][iq]);
				
        if((i > 4) && ((fabs(d[ip])+g) == fabs(d[ip])) && ((fabs(d[iq])+g) == fabs(d[iq])))
					A[ip][iq] = 0.0f;
				
        else if(fabs(A[ip][iq]) > tresh) {
          h = d[iq] - d[ip];
					
          if((fabs(h)+g) == fabs(h))
            t = A[ip][iq] / h;
          else {
            theta = 0.5f * h / A[ip][iq];
            t = 1.0f / (fabs(theta) + sqrt(1.0f + theta * theta));
            if(theta < 0.0f) t *= -1.0f;
          }
					
          c = 1.0f / sqrt(1.0f + t * t);
          s = t * c;
          tau = s / (1.0f + c);
          h = t * A[ip][iq];
          z[ip] -= h;
          z[iq] += h;
          d[ip] -= h;
          d[iq] += h;
          A[ip][iq] = 0.0f;
					
          for(j = 0; j < ip; j++)
            rot(A, s, tau, j, ip, j, iq);
					
          for(j = ip + 1; j < iq; j++)
            rot(A, s, tau, ip, j, j, iq);
					
          for(j = iq + 1; j < 3; j++)
            rot(A, s, tau, ip, j, iq, j);
					
          for(j = 0; j < 3; j++)
            rot(V, s, tau, j, ip, j, iq);
					
          //nrot ++;	
					
        }
      }
    }
		
    for(ip = 0; ip < 3; ip++) {
      b[ip] += z[ip];
      d[ip] = b[ip];
      z[ip] = 0.0f;
    }
		
  }
  
  return false;
}

}//	end of namespace

