/* Triangle/triangle intersection test routine,
 * by Tomas Moller, 1997.
 * See article "A Fast Triangle-Triangle Intersection Test",
 * Journal of Graphics Tools, 2(2), 1997
 *
 * int tri_tri_intersect(number V0[3],number V1[3],number V2[3],
 *                         number U0[3],number U1[3],number U2[3])
 *
 * parameters: vertices of triangle 1: V0,V1,V2
 *             vertices of triangle 2: U0,U1,U2
 * result    : returns 1 if the triangles intersect, otherwise 0
 *
 */

#include <cmath>
#include <cstring>
#include <cstdlib>
#include "math_util.h"

using namespace ug;

/* if USE_EPSILON_TEST is true then we do a check: 
         if |dv|<EPSILON then dv=0.0;
   else no check is done (which is less robust)
*/
//  change by sreiter: Note that a new variable snapThreshold was added to
//  tri_tri_intersect. EPSILON is now only used for the coplanarity check...
#define USE_EPSILON_TEST TRUE  
#define EPSILON 1.e-20


/* some macros */
#define CROSS(dest,v1,v2)                      \
              dest[0]=v1[1]*v2[2]-v1[2]*v2[1]; \
              dest[1]=v1[2]*v2[0]-v1[0]*v2[2]; \
              dest[2]=v1[0]*v2[1]-v1[1]*v2[0];

#define DOT(v1,v2) (v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])

#define SUB(dest,v1,v2)          \
            dest[0]=v1[0]-v2[0]; \
            dest[1]=v1[1]-v2[1]; \
            dest[2]=v1[2]-v2[2]; 

/* sort so that a<=b */
#define SORT(a,b)       \
             if(a>b)    \
             {          \
               number c; \
               c=a;     \
               a=b;     \
               b=c;     \
             }

#define ISECT(VV0,VV1,VV2,D0,D1,D2,isect0,isect1) \
              isect0=VV0+(VV1-VV0)*D0/(D0-D1);    \
              isect1=VV0+(VV2-VV0)*D0/(D0-D2);


#define COMPUTE_INTERVALS(VV0,VV1,VV2,D0,D1,D2,D0D1,D0D2,isect0,isect1) \
  if(D0D1>0.0f)                                         \
  {                                                     \
    /* here we know that D0D2<=0.0 */                   \
    /* that is D0, D1 are on the same side, D2 on the other or on the plane */ \
    ISECT(VV2,VV0,VV1,D2,D0,D1,isect0,isect1);          \
  }                                                     \
  else if(D0D2>0.0f)                                    \
  {                                                     \
    /* here we know that d0d1<=0.0 */                   \
    ISECT(VV1,VV0,VV2,D1,D0,D2,isect0,isect1);          \
  }                                                     \
  else if(D1*D2>0.0f || D0!=0.0f)                       \
  {                                                     \
    /* here we know that d0d1<=0.0 or that D0!=0.0 */   \
    ISECT(VV0,VV1,VV2,D0,D1,D2,isect0,isect1);          \
  }                                                     \
  else if(D1!=0.0f)                                     \
  {                                                     \
    ISECT(VV1,VV0,VV2,D1,D0,D2,isect0,isect1);          \
  }                                                     \
  else if(D2!=0.0f)                                     \
  {                                                     \
    ISECT(VV2,VV0,VV1,D2,D0,D1,isect0,isect1);          \
  }                                                     \
  else                                                  \
  {                                                     \
    /* triangles are coplanar */                        \
    return coplanar_tri_tri(N1,V0,V1,V2,U0,U1,U2);      \
  }

//	VV0, VV1, VV2, isect0 and isect1 have to be arrays of 3 numbers.
#define ISECT3(VV0,VV1,VV2,D0,D1,D2,isect0,isect1) \
              isect0[0]=VV0[0]+(VV1[0]-VV0[0])*D0/(D0-D1);    \
              isect0[1]=VV0[1]+(VV1[1]-VV0[1])*D0/(D0-D1);    \
              isect0[2]=VV0[2]+(VV1[2]-VV0[2])*D0/(D0-D1);    \
              isect1[0]=VV0[0]+(VV2[0]-VV0[0])*D0/(D0-D2);	  \
              isect1[1]=VV0[1]+(VV2[1]-VV0[1])*D0/(D0-D2);	  \
              isect1[2]=VV0[2]+(VV2[2]-VV0[2])*D0/(D0-D2);

//	VV0, VV1, VV2, isect0 and isect1 have to be arrays of 3 numbers.
#define COMPUTE_INTERVAL_POINTS(VV0,VV1,VV2,D0,D1,D2,D0D1,D0D2,isect0,isect1) \
  if(D0D1>0.0f)                                         \
  {                                                     \
    /* here we know that D0D2<=0.0 */                   \
    /* that is D0, D1 are on the same side, D2 on the other or on the plane */ \
    ISECT3(VV2,VV0,VV1,D2,D0,D1,isect0,isect1);          \
  }                                                     \
  else if(D0D2>0.0f)                                    \
  {                                                     \
    /* here we know that d0d1<=0.0 */                   \
    ISECT3(VV1,VV0,VV2,D1,D0,D2,isect0,isect1);          \
  }                                                     \
  else if(D1*D2>0.0f || D0!=0.0f)                       \
  {                                                     \
    /* here we know that d0d1<=0.0 or that D0!=0.0 */   \
    ISECT3(VV0,VV1,VV2,D0,D1,D2,isect0,isect1);          \
  }                                                     \
  else if(D1!=0.0f)                                     \
  {                                                     \
    ISECT3(VV1,VV0,VV2,D1,D0,D2,isect0,isect1);          \
  }                                                     \
  else if(D2!=0.0f)                                     \
  {                                                     \
    ISECT3(VV2,VV0,VV1,D2,D0,D1,isect0,isect1);          \
  }                                                     \
  else                                                  \
  {                                                     \
    /* triangles are coplanar */                        \
	/* we ignore this case since it we already returned.*/ \
  }



/* this edge to edge test is based on Franlin Antonio's gem:
   "Faster Line Segment Intersection", in Graphics Gems III,
   pp. 199-202 */ 
#define EDGE_EDGE_TEST(V0,U0,U1)                      \
  Bx=U0[i0]-U1[i0];                                   \
  By=U0[i1]-U1[i1];                                   \
  Cx=V0[i0]-U0[i0];                                   \
  Cy=V0[i1]-U0[i1];                                   \
  f=Ay*Bx-Ax*By;                                      \
  d=By*Cx-Bx*Cy;                                      \
  if((f>0 && d>=0 && d<=f) || (f<0 && d<=0 && d>=f))  \
  {                                                   \
    e=Ax*Cy-Ay*Cx;                                    \
    if(f>0)                                           \
    {                                                 \
      if(e>=0 && e<=f) return 2;                      \
    }                                                 \
    else                                              \
    {                                                 \
      if(e<=0 && e>=f) return 2;                      \
    }                                                 \
  }                                

#define EDGE_AGAINST_TRI_EDGES(V0,V1,U0,U1,U2) \
{                                              \
  number Ax,Ay,Bx,By,Cx,Cy,e,d,f;               \
  Ax=V1[i0]-V0[i0];                            \
  Ay=V1[i1]-V0[i1];                            \
  /* test edge U0,U1 against V0,V1 */          \
  EDGE_EDGE_TEST(V0,U0,U1);                    \
  /* test edge U1,U2 against V0,V1 */          \
  EDGE_EDGE_TEST(V0,U1,U2);                    \
  /* test edge U2,U1 against V0,V1 */          \
  EDGE_EDGE_TEST(V0,U2,U0);                    \
}

#define POINT_IN_TRI(V0,U0,U1,U2)           \
{                                           \
  number a,b,c,d0,d1,d2;                     \
  /* is T1 completly inside T2? */          \
  /* check if V0 is inside tri(U0,U1,U2) */ \
  a=U1[i1]-U0[i1];                          \
  b=-(U1[i0]-U0[i0]);                       \
  c=-a*U0[i0]-b*U0[i1];                     \
  d0=a*V0[i0]+b*V0[i1]+c;                   \
                                            \
  a=U2[i1]-U1[i1];                          \
  b=-(U2[i0]-U1[i0]);                       \
  c=-a*U1[i0]-b*U1[i1];                     \
  d1=a*V0[i0]+b*V0[i1]+c;                   \
                                            \
  a=U0[i1]-U2[i1];                          \
  b=-(U0[i0]-U2[i0]);                       \
  c=-a*U2[i0]-b*U2[i1];                     \
  d2=a*V0[i0]+b*V0[i1]+c;                   \
  if(d0*d1>0.0)                             \
  {                                         \
    if(d0*d2>0.0) return 2;                 \
  }                                         \
}

int coplanar_tri_tri(number N[3],number V0[3],number V1[3],number V2[3],
                     number U0[3],number U1[3],number U2[3])
{
   number A[3];
   short i0,i1;
   /* first project onto an axis-aligned plane, that maximizes the area */
   /* of the triangles, compute indices: i0,i1. */
//#pragma warning( disable : 4244 )
   A[0]=fabs(N[0]);
   A[1]=fabs(N[1]);
   A[2]=fabs(N[2]);
//#pragma warning( default : 4244 )
   if(A[0]>A[1])
   {
      if(A[0]>A[2])  
      {
          i0=1;      /* A[0] is greatest */
          i1=2;
      }
      else
      {
          i0=0;      /* A[2] is greatest */
          i1=1;
      }
   }
   else   /* A[0]<=A[1] */
   {
      if(A[2]>A[1])
      {
          i0=0;      /* A[2] is greatest */
          i1=1;                                           
      }
      else
      {
          i0=0;      /* A[1] is greatest */
          i1=2;
      }
    }               
                
    /* test all edges of triangle 1 against the edges of triangle 2 */
    EDGE_AGAINST_TRI_EDGES(V0,V1,U0,U1,U2);
    EDGE_AGAINST_TRI_EDGES(V1,V2,U0,U1,U2);
    EDGE_AGAINST_TRI_EDGES(V2,V0,U0,U1,U2);
                
    /* finally, test if tri1 is totally contained in tri2 or vice versa */
    POINT_IN_TRI(V0,U0,U1,U2);
    POINT_IN_TRI(U0,V0,V1,V2);

    return 0;
}


int tri_tri_intersect(number V0[3],number V1[3],number V2[3],
                      number U0[3],number U1[3],number U2[3],
					            number* ip1Out, number* ip2Out, const number snapThreshold)
{
  number E1[3],E2[3];
  number N1[3],N2[3],d1,d2;
  number du0,du1,du2,dv0,dv1,dv2;
  number D[3];
  number isect1[2], isect2[2];
  number du0du1,du0du2,dv0dv1,dv0dv2;
  short index;
  number vp0,vp1,vp2;
  number up0,up1,up2;
  number b,c,max;

  /* compute plane equation of triangle(V0,V1,V2) */
  SUB(E1,V1,V0);
  SUB(E2,V2,V0);
  CROSS(N1,E1,E2);
  d1=-DOT(N1,V0);
  /* plane equation 1: N1.X+d1=0 */

  /* put U0,U1,U2 into plane equation 1 to compute signed distances to the plane*/
  du0=DOT(N1,U0)+d1;
  du1=DOT(N1,U1)+d1;
  du2=DOT(N1,U2)+d1;

  /* coplanarity robustness check */
#if USE_EPSILON_TEST==TRUE
  if(fabs(du0)<EPSILON) du0=0.0;
  if(fabs(du1)<EPSILON) du1=0.0;
  if(fabs(du2)<EPSILON) du2=0.0;
#endif
  du0du1=du0*du1;
  du0du2=du0*du2;
  if(du0du1>0.0f && du0du2>0.0f) /* same sign on all of them + not equal 0 ? */
    return 0;                    /* no intersection occurs */

  /* compute plane of triangle (U0,U1,U2) */
  SUB(E1,U1,U0);
  SUB(E2,U2,U0);
  CROSS(N2,E1,E2);
  d2=-DOT(N2,U0);
  /* plane equation 2: N2.X+d2=0 */

  /* put V0,V1,V2 into plane equation 2 */
  dv0=DOT(N2,V0)+d2;
  dv1=DOT(N2,V1)+d2;
  dv2=DOT(N2,V2)+d2;

//#if USE_EPSILON_TEST==TRUE
  if(fabs(dv0)<snapThreshold) dv0=0.0;
  if(fabs(dv1)<snapThreshold) dv1=0.0;
  if(fabs(dv2)<snapThreshold) dv2=0.0;
//#endif

  dv0dv1=dv0*dv1;
  dv0dv2=dv0*dv2;
  if(dv0dv1>0.0f && dv0dv2>0.0f) /* same sign on all of them + not equal 0 ? */
    return 0;                    /* no intersection occurs */

  /* compute direction of intersection line */
  CROSS(D,N1,N2);

  /* compute and index to the largest component of D */
//#pragma warning( disable : 4244 )
  max=fabs(D[0]);
  index=0;
  b=fabs(D[1]);
  c=fabs(D[2]);
//#pragma warning( default : 4244 )
  if(b>max) max=b,index=1;
  if(c>max) max=c,index=2;

        /* this is the simplified projection onto L*/
        vp0=V0[index];
        vp1=V1[index];
        vp2=V2[index];

        up0=U0[index];
        up1=U1[index];
        up2=U2[index];

  /* compute interval for triangle 1 */
  COMPUTE_INTERVALS(vp0,vp1,vp2,dv0,dv1,dv2,dv0dv1,dv0dv2,isect1[0],isect1[1]);

  /* compute interval for triangle 2 */
  COMPUTE_INTERVALS(up0,up1,up2,du0,du1,du2,du0du1,du0du2,isect2[0],isect2[1]);

// those indices are used to find the maximum and minimum later on
  int ind1 = 1;
  int ind2 = 1;
  if(isect1[0] > isect1[1]){
  	std::swap(isect1[0], isect1[1]);
	ind1 = 0;
  }

  if(isect2[0] > isect2[1]){
  	std::swap(isect2[0], isect2[1]);
	ind2 = 0;
  }  
 // SORT(isect1[0],isect1[1]);
 // SORT(isect2[0],isect2[1]);
  if(isect1[1]<isect2[0] || isect2[1]<isect1[0]) return 0;

  if(ip1Out && ip2Out){
  	//	calculate the endpoints of the line segment that resembles the intersection.
	number tpa1[3]={0, 0, 0}, tpa2[3]={0, 0, 0}, tpb1[3]={0, 0, 0},
			tpb2[3]={0, 0, 0};

	/* compute interval-points for triangle 1 */
	COMPUTE_INTERVAL_POINTS(V0,V1,V2,dv0,dv1,dv2,dv0dv1,dv0dv2,tpa1,tpa2);

	/* compute interval-points for triangle 2 */
	COMPUTE_INTERVAL_POINTS(U0,U1,U2,du0,du1,du2,du0du1,du0du2,tpb1, tpb2);
	
	//choose the right ones and return
	if(isect1[0] > isect2[0])
		if(ind1 == 0)
			memcpy((char*)ip1Out, tpa2, sizeof(number)*3);
		else
			memcpy((char*)ip1Out, tpa1, sizeof(number)*3);
	else
		if(ind2 == 0)
			memcpy((char*)ip1Out, tpb2, sizeof(number)*3);
		else
			memcpy((char*)ip1Out, tpb1, sizeof(number)*3);
	
	if(isect1[1] < isect2[1])
		if(ind1 == 0)
			memcpy((char*)ip2Out, tpa1, sizeof(number)*3);
		else
			memcpy((char*)ip2Out, tpa2, sizeof(number)*3);
	else
		if(ind2 == 0)
			memcpy((char*)ip2Out, tpb1, sizeof(number)*3);
		else
			memcpy((char*)ip2Out, tpb2, sizeof(number)*3);
  }
  return 1;
}

namespace ug{
bool TriangleTriangleIntersection(const MathVector<3>& p0, const MathVector<3>& p1,
								  const MathVector<3>& p2, const MathVector<3>& q0,
								  const MathVector<3>& q1, const MathVector<3>& q2,
								  MathVector<3>* ip1Out, MathVector<3>* ip2Out,
                  number snapThreshold)
{
	return tri_tri_intersect((number*)&p0, (number*)&p1, (number*)&p2,
							 (number*)&q0, (number*)&q1, (number*)&q2,
							 (number*)ip1Out, (number*)ip2Out, snapThreshold) == 1;
}

}//	end of namespace
