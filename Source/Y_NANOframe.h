/* Copyright (C) 2000, Dr. Antonio Munjiza
 *
 * This code is provided as part of the book entitled "The Combined
 * Finite Discrete Element Method". It is distributed WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
 * PARTICULAR PURPOSE. Inclusion of a part or whole of this code into any other
 * commercial or research or other purpose code is not granted without author's
 * written explicit permission.
 * When results using whole or any part of this code
 * are published, Y code must be mentioned and acknowledgement to Dr Munjiza must be made.
 * Should you modify this source code, the Copyright (C) on the mdified code
 * as a whole belongs to Dr. A. Munjiza regardless of the extent or nature
 * of modifications.
 * Copyright (C) to whole of any code containing any part of this code
 * also belongs to Dr. A.Munjiza. 
 * Any code comprising any part of this source code
 * must be called Y program.
 * If you do not agree with this, you are not allowed to do 
 * any modifications to any part of this source code or included 
 * any part of it in any other program.
*/
#ifndef Y_NANOFRAME
#define Y_NANOFRAME

#define NeedFunctionPrototypes 1

#include <stdio.h>
#include <string.h>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#define CHR char
#define INT int
#define DBL double
#ifdef WIN32
    #define LONG __int64
#else
    #define LONG long long
#endif

#define YES 1
#define NO  0

#define D0          0.0
#define D1          1.0
#define D2          2.0
#define D3          3.0
#define D4          4.0
#define D5          5.0
#define D6          6.0
#define D7          7.0
#define D8          8.0
#define D9          9.0
#define D10         10.0
#define D24         24.0
#define D35         35.0
#define D45         45.0
#define D48         48.0
#define D50         50.0
#define D96         96.0
#define D100        100.0
#define D168        168.0
#define D624        624.0
#define D2000       2000.0
#define D4000       4000.0

#define MLAYER      4

#define DP1      0.1
#define DP2      0.2
#define DP5      0.5
#define DP7      0.7
#define DP15     1.5
#define DP45     4.5

#define DPP5     0.05

#define DPPP5    0.005

#define MAXIM(x,y) (((x)<(y))?(y):(x))
#define MINIM(x,y) (((x)>(y))?(y):(x))
#define DABS(x) (((x)<D0)?-(x):(x))
#define ABS(x) (((x)<0)?-(x):(x))
#define EXP exp
#define SQRT sqrt
#define SIGN(x) (((x)<D0)?-1:1)

#define MYPI 3.1415926535897932384626
#define DGP2 0.577350269189626
#define EPSILON 1.0e-10
#define BEPSILON 1.0e+15

#define MCHPO 30000     /* M_CHar_per_Object */
#define MBASE   80      /* M_BASE_characters */
#define MHASH  550      /* M_HASH_array      */
#define MPARA  300      /* M_PARAmeters      */

#define FILENULL ((FILE*)NULL)
#define CHR2NULL ((CHR **)NULL)
#define DBL1NULL ((DBL*)NULL)
#define DBL2NULL ((DBL**)NULL)
#define DBL3NULL ((DBL***)NULL)
#define INT1NULL ((INT*)NULL)
#define INT2NULL ((INT**)NULL)
#define INT3NULL ((INT***)NULL)

#define FREE(x) if ((x)!=NULL) free(x)
#define MALLOC(x) malloc(x)

/* =========STREAM======================================== 
****************************************************************************
*  This module contains macros to perform operations on STREAMS
*****************************************************************************/
#define FILEND(fptr) (feof(fptr))
extern CHR *DBL_S[20];
#define DBL_SR "%le" 
extern CHR *INT_S[20];
#define INT_SR   "%d"
 
#define CHR_S "%s"
#define CHRRETURN "\n"
#define CHRSPACE " "
#define CHRTERMINATE 0
#define INTr(file,x){fscanf((file),INT_SR,((x)));}
#define DBLr(file,x){fscanf((file),DBL_SR,((x)));}
#define CHRr(file,x){fscanf((file),CHR_S,((x)));}
#define INTw(file,x,ndigit)fprintf((file),(INT_S[(ndigit)]),((x))); 
#define SINTw(file,x,ndigit)sprintf((file),(INT_S[(ndigit)]),((x))); 
#define DBLw(file,x,ndigit) fprintf((file),(DBL_S[((ndigit))]),(x));
 
#define SCHRw(file,x){sprintf((file),CHR_S,(x));}
#define CHRw(file,x){fprintf((file),CHR_S,(x));}
#define CHRend(c,i){c[i]='\0';}
#define CHRwcr(file){fprintf((file),CHR_S,(CHRRETURN));}
#define CHRwsp(file){fprintf((file),CHR_S,(CHRSPACE));}
#define CHRcmp(c1,c2,n)(strncmp((c1),(c2),(n)))
#define CHRcpy(c1,c2)(strcpy((c1),(c2)))
#define CHRcat(c1,c2)(strcat((c1),(c2)))

#define TDBLNormalise(k,ddiv) ((DBL)((k))/(ddiv)); 
#define TINTDenormalise(dk,ddiv) ((INT)(dk*ddiv));
#define TDBLNormalDist(x) ((D3/DSQRT(D2*MYPI)*EXP(-DP45*(x)*(x))))
#define TINTCellnormalise(dcoor,dcmin,dcele) ((INT)((dcoor-dcmin)/dcele));

/* =========Small MATRIX======================================== 
****************************************************************************
*  This module contains macros to perform operations on  small matrices
*****************************************************************************/
   /* INVERSE A SMALL MATRIX - also return determinant  */
#define YMATINV3(m,minv,det)\
      {  det=m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1])-\
             m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0])+\
             m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0]);\
         minv[0][0]=(m[1][1]*m[2][2]-m[1][2]*m[2][1])/det;\
         minv[1][0]=(m[1][2]*m[2][0]-m[1][0]*m[2][2])/det;\
         minv[2][0]=(m[1][0]*m[2][1]-m[1][1]*m[2][0])/det;\
         minv[0][1]=(m[0][2]*m[2][1]-m[0][1]*m[2][2])/det;\
         minv[1][1]=(m[0][0]*m[2][2]-m[0][2]*m[2][0])/det;\
         minv[2][1]=(m[0][1]*m[2][0]-m[0][0]*m[2][1])/det;\
         minv[0][2]=(m[0][1]*m[1][2]-m[0][2]*m[1][1])/det;\
         minv[1][2]=(m[0][2]*m[1][0]-m[0][0]*m[1][2])/det;\
         minv[2][2]=(m[0][0]*m[1][1]-m[0][1]*m[1][0])/det;\
    }
#define YMATINV2(m,minv,det)\
      {  det=m[0][0]*m[1][1]-m[1][0]*m[0][1]; \
         minv[0][0]= m[1][1]/det; minv[1][0]=-m[1][0]/det;\
         minv[0][1]=-m[0][1]/det; minv[1][1]= m[0][0]/det;\
    }
/* =========vect3D======================================== 
****************************************************************************
*  This module contains macros to perform operations on 3D vectors.
*****************************************************************************/

/* Misc. Vector Math Macro Definitions */

extern DBL VTemp;
/* Assign; Null */
#define V3DAss(x1,y1,z1,x2,y2,z2) {(x1)=(x2);(y1)=(y2);z1=(z2);}
#define V3DNul(x1,y1,z1) {(x1)=R0;(y1)=R0;(z1)=R0;}
/* Vector Add; Subtract; Scale; Length square; Length; 
   Dot product; Cross product; Volume of triad */
#define V3DAdd(x1,y1,z1,x2,y2,z2,x3,y3,z3) \
              {(x1)=(x2)+(x3);(y1)=(y2)+(y3);(z1)=(z2)+(z3);}
#define V3DSub(x1,y1,z1,x2,y2,z2,x3,y3,z3) \
              {(x1)=(x2)-(x3);(y1)=(y2)-(y3);(z1)=(z2)-(z3);}
#define V3DSca(x1,y1,z1,s){(x1)=(x1)*(s);(y1)=(y1)*(s);(z1)=(z1)*(s);}
#define V3DDiv(x1,y1,z1,s){(x1)=(x1)/(s);(y1)=(y1)/(s);(z1)=(z1)/(s);}
#define V3DLe2(s,x1,y1,z1){(s)=(x1)*(x1)+(y1)*(y1)+(z1)*(z1);}
#define V3DLen(s,x1,y1,z1){(s)=SQRT((x1)*(x1)+(y1)*(y1)+(z1)*(z1));}
#define V3DDot(s,x1,y1,z1,x2,y2,z2) {(s)=((x1)*(x2))+((y1)*(y2))+((z1)*(z2));}
#define V3DCro(x1,y1,z1,x2,y2,z2,x3,y3,z3)\
              {(x1)=(y2)*(z3)-(z2)*(y3); \
               (y1)=(z2)*(x3)-(x2)*(z3); \
               (z1)=(x2)*(y3)-(y2)*(x3); }
#define V3DOmegaRot(x1,y1,z1,x2,y2,z2,x3,y3,z3)\
              {(x1)=(x3)+(y2)*(z3)-(z2)*(y3); \
               (y1)=(y3)+(z2)*(x3)-(x2)*(z3); \
               (z1)=(z3)+(x2)*(y3)-(y2)*(x3); }
#define V3DTranToGl(x1,y1,z1,a,b,c,x2,y2,z2,x3,y3,z3)\
              {(x1)=(a)*(x2)+(b)*(x3)+(c)*((y2)*(z3)-(z2)*(y3)); \
               (y1)=(a)*(y2)+(b)*(y3)+(c)*((z2)*(x3)-(x2)*(z3)); \
               (z1)=(a)*(z2)+(b)*(z3)+(c)*((x2)*(y3)-(y2)*(x3)); }
#define V3DTranToLoc(x1,y1,z1,a,b,c,x2,y2,z2,x3,y3,z3)\
              {(x1)=(a)*(x2)+(b)*(y2)+(c)*(z2);  \
               (y1)=(a)*(x3)+(b)*(y3)+(c)*(z3); \
(z1)=(a)*((y2)*(z3)-(z2)*(y3))+(b)*((x3)*(z2)-(x2)*(z3))+(c)*((x2)*(y3)-(y2)*(x3)); }
#define V3DVol(v,x2,y2,z2,x3,y3,z3,x1,y1,z1) \
              {(v)=((y2)*(z3)-(z2)*(y3))*(x1)+ \
                   ((z2)*(x3)-(x2)*(z3))*(y1)+ \
                   ((x2)*(y3)-(y2)*(x3))*(z1); }
  /*  Normalize a Vector; Halfway Vectors; */
#define V3DNor(s,x1,y1,z1) \
              {(s)=SQRT((x1)*(x1)+(y1)*(y1)+(z1)*(z1));   \
              if((s)>EPSILON)(x1)=(x1)/(s);  \
              if((s)>EPSILON)(y1)=(y1)/(s);  \
              if((s)>EPSILON)(z1)=(z1)/(s);  } 
#define V3DMid(x1,y1,z1,x2,y2,z2,x3,y3,z3) \                                   \
              {(x1)=(RP5)*((x2)+(x3)); \
               (y1)=(RP5)*((y2)+(y3)); \
               (z1)=(RP5)*((z2)+(z3)); }
/* Vector in Between two given vectors  */
#define V3DBet(x1,y1,z1,x2,y2,z2,x3,y3,z3,s2,s3)\
              {(x1)=(x2)*(s2)+(x3)*(s3);     \
               (y1)=(y2)*(s2)+(y3)*(s3);     \
               (z1)=(z2)*(s2)+(z3)*(s3);     }
/* =========vect2D======================================== 
****************************************************************************
*  This module contains macros to perform operations on 2D vectors.
*****************************************************************************/
/* Assign; Null */
#define V2DAss(x1,y1,x2,y2) {(x1)=(x2);(y1)=(y2)}
#define V2DNul(x1,y1) {(x1)=R0;(y1)=R0;} 
/* Vector Add; Subtract; Scale; Length square; Length; Dot product; Cross product;  */
#define V2DAdd(x1,y1,x2,y2,x3,y3) {(x1)=(x2)+(x3);(y1)=(y2)+(y3);}
#define V2DSub(x1,y1,x2,y2,x3,y3) {(x1)=(x2)-(x3);(y1)=(y2)-(y3);}
#define V2DSca(x1,y1,s) {(x1)=(x1)*s;(y1)=(y1)*s;}
#define V2DLe2(s,x1,y1) {(s)=(x1)*(x1)+(y1)*(y1);}
#define V2DLen(s,x1,y1) {(s)=SQRT((x1)*(x1)+(y1)*(y1);}
#define V2DDot(s,x1,y1,x2,y2) {s=(x1*x2)+(y1*y2);}
#define V2DCro(s,x2,y2,x3,y3){s=(x1*y2)-(y1*x2);} 
  
/*  Normalize a Vector; Compute a Vector Halfway Between Two Given Vectors; */ 
#define V2DNor(x1,y1,x2,y2){VTemp=SQRT((x1)*(x1)+(y1)*(y1));   \
              (x1)=(x2)/VTemp;(y1)=(y2)/VTemp;                 } 
#define V2DMid(x1,y1,x2,y2,x3,y3){(x1)=0.5*((x2)+(x3));(y1)=0.5*((y2)+(y3));}
                                               
/* Vector in Between two given vectors  */
#define V2DBet(x1,y1,x2,y2,x3,y3,s2,s3){(x1)=(x2)*(s2)+(x3)*(s3);     \
              (y1)=(y2)*(s2)+(y3)*(s3);     }


inline DBL DSQRT(DBL x)		/* Returns the Square Root of a DBL number */
{ return sqrt(x);
}

inline DBL DFABS(DBL x)		/* Returns the Square Root of a DBL number */
{ return fabs(x);
}

/* ALLOCATION SUBROUTINES		*/

inline INT *TalINT1(INT m1)
{ INT isize=m1*sizeof(INT);
  if(isize==0)return INT1NULL;
  return (INT*)MALLOC(isize);   
} 

inline CHR *TalCHR1(INT m1)
{ INT isize=m1*sizeof(CHR);
  if(isize==0)return NULL;
  return (CHR*)MALLOC(isize);   
} 

inline DBL *TalDBL1(INT m1)
{ INT isize=m1*sizeof(DBL);
  if(isize==0)return DBL1NULL;
  return (DBL*)MALLOC(isize);
}

inline INT **TalINT2(INT m2,INT m1)
{ INT isize,i2;
  INT     *p1;
  INT    **p2;
  void    *v1;

  isize=sizeof(INT*)*(m2+3)+
        sizeof(INT )*(m2*m1+3);
  if(isize==0)return INT2NULL;
  v1=MALLOC(isize);
  p2=(INT**)v1;
  p1=(INT*)v1;
  p1=p1+((m2+1)*sizeof(INT**))/sizeof(INT)+2;
  for(i2=0;i2<m2;i2++)
  { p2[i2]=p1+i2*m1; 
  }
  return p2; 
}

/* =========RandomNumberGenerator======================================== 
****************************************************************************
*  This module contains functions to get random numbers.
*****************************************************************************/

DBL TranNumUniDistGen(INT *i0dum);      /* UniformDistribution  */
DBL TranNumGaussDisGen(INT *i0dum);     /* NormalDistribution   */


/*************SPACESAVING by 90 FORMAT ************/

void codeCHRtoINT(  /* write array as coded */
#if NeedFunctionPrototypes 
  CHR *c1code,         /* coded array [0]=ndigit; [1]=nnum     */
  INT *i1num           /* INT array [0]=ndigit; [1]=nnum
                          the rest between 0 and base^ndigit    */
#endif
);

void codeINTtoCHR(  /* write array as coded */
#if NeedFunctionPrototypes 
  CHR *c1code,         /* coded array [0]=ndigit; [1]=nnum     */
  INT *i1num           /* INT array [0]=ndigit; [1]=nnum
                          the rest between 0 and base^ndigit    */
#endif
);

void ZCHRcat(       /* Appends a four digit number to a string  */
#if NeedFunctionPrototypes 
  CHR *c1,              /* String                               */
  INT inum              /* IntegerNumber                        */
#endif
);

/* =========The greatest common divisor Generator======================================== 
****************************************************************************
*  This module contains functions to get the greatest common divisor of two intergers.
*****************************************************************************/
INT Gcd(
#if NeedFunctionPrototypes 
  INT m,              
  INT n              
#endif
);

#endif
