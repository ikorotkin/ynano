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
#include "Y_NANOSetEquilibrium.h"
#include "Y_NANOtoolbox.h"
#include "Y_NANOcovepoint.h"
#include "Y_NANOconvert.h"

static void Equilibrate(Y_NANOContaiBinList *v0cbl)
{ INT i;                /* LoopsVariable                        */
  INT j;                /* LoopsVariable                        */
  INT natoms;           /* TotalNumberOfAtoms                   */
  DBL dmomen[3];        /* TotalMomentum                        */
  DBL A[3][3];          /* MatrixOfCoefficients                 */
  DBL Ainv[3][3];       /* InverseMatrix                        */
  DBL det;              /* DeterminantOfMatrixA                 */
  DBL b[3];             /* VectorOf....                         */
  DBL x1;               /* CrossProductXComponent               */
  DBL y1;               /* CrossProductYComponent               */
  DBL z1;               /* CrossProductZComponent               */
  DBL domega[3];        /* AngularVelocity                      */
  DBL dcgx;             /* CentreGravityXComponent              */
  DBL dcgy;             /* CentreGravityYComponent              */
  DBL dcgz;             /* CentreGravityZComponent              */
  DBL dvx;              /* DeltaVXDueToRotation                 */
  DBL dvy;              /* DeltaVYDueToRotation                 */
  DBL dvz;              /* DeltaVZDueToRotation                 */
  DBL dx;               /* X-CGX                                */
  DBL dy;               /* Y-CGY                                */
  DBL dz;               /* Z-CGZ                                */
  DBL dtmass;           /* TotalMass                            */
  Y_NANOobj *v0tmp;        /* TMPobject                            */
  Y_NANOcovepoint vptmp;   /* CoordinateVelocityPntTMP             */

  for(i=0;i<3;i++)
  { dmomen[i]=D0;
    b[i]=D0;
    domega[i]=D0;
    dtmass=D0;
    dcgx=D0;
    dcgy=D0;
    dcgz=D0;
    for(j=0;j<3;j++) A[i][j]=D0;
  }

  natoms=0;
  v0tmp=v0cbl->v0next;
  while(v0tmp!=v0cbl)
  { vtoolbox.Y_NANOtool(&convert,v0tmp,&vptmp);
    dtmass=dtmass+vptmp.dma;
    dcgx=dcgx+vptmp.dma*vptmp.dcx;
    dcgy=dcgy+vptmp.dma*vptmp.dcy;
    dcgz=dcgz+vptmp.dma*vptmp.dcz;
    dmomen[0]=dmomen[0]-vptmp.dma*vptmp.dvx;
    dmomen[1]=dmomen[1]-vptmp.dma*vptmp.dvy;
    dmomen[2]=dmomen[2]-vptmp.dma*vptmp.dvz;
    natoms=natoms+1;
    v0tmp=v0tmp->v0next;
  }

  dcgx=dcgx/dtmass;
  dcgy=dcgy/dtmass;
  dcgz=dcgz/dtmass;

  v0tmp=v0cbl->v0next;
  while(v0tmp!=v0cbl)
  { vtoolbox.Y_NANOtool(&convert,v0tmp,&vptmp);
    vptmp.dvx=vptmp.dvx+dmomen[0]/(natoms*vptmp.dma);
    vptmp.dvy=vptmp.dvy+dmomen[1]/(natoms*vptmp.dma);
    vptmp.dvz=vptmp.dvz+dmomen[2]/(natoms*vptmp.dma);
    vtoolbox.Y_NANOtool(&convert,&vptmp,v0tmp);
    dx=vptmp.dcx-dcgx;
    dy=vptmp.dcy-dcgy;
    dz=vptmp.dcz-dcgz;
    A[0][0]=A[0][0]+vptmp.dma*(dy*dy+dz*dz);
    A[1][1]=A[1][1]+vptmp.dma*(dx*dx+dz*dz);
    A[2][2]=A[2][2]+vptmp.dma*(dx*dx+dy*dy);
    A[0][1]=A[0][1]-vptmp.dma*dx*dy;
    A[0][2]=A[0][2]-vptmp.dma*dx*dz;
    A[1][2]=A[1][2]-vptmp.dma*dy*dz;
    V3DCro(x1,y1,z1,dx,dy,dz,vptmp.dvx,vptmp.dvy,vptmp.dvz);
    b[0]=b[0]+vptmp.dma*x1;
    b[1]=b[1]+vptmp.dma*y1;
    b[2]=b[2]+vptmp.dma*z1;
    v0tmp=v0tmp->v0next;
  }

  A[1][0]=A[0][1];
  A[2][0]=A[0][2];
  A[2][1]=A[1][2];
  YMATINV3(A,Ainv,det)
  for(i=0;i<3;i++)
  { for(j=0;j<3;j++)
    { domega[i]=domega[i]+Ainv[i][j]*b[j];
  } }

  v0tmp=v0cbl->v0next;
  while(v0tmp!=v0cbl)
  { vtoolbox.Y_NANOtool(&convert,v0tmp,&vptmp);
    dx=vptmp.dcx-dcgx;
    dy=vptmp.dcy-dcgy;
    dz=vptmp.dcz-dcgz;
    V3DCro(dvx,dvy,dvz,domega[0],domega[1],domega[2],dx,dy,dz);
    vptmp.dvx=vptmp.dvx-dvx;
    vptmp.dvy=vptmp.dvy-dvy;
    vptmp.dvz=vptmp.dvz-dvz;
    vtoolbox.Y_NANOtool(&convert,&vptmp,v0tmp);
    v0tmp=v0tmp->v0next;
  }
}

void SetEquilibrium(Y_NANOContaiBinList *v0cbl)
{ Equilibrate(v0cbl);
}

