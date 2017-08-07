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
#include "Y_NANOcontactBnCubAr.h"
#include "TMPPotEnerg.h"
#include "Y_NANOtoolbox.h"
#include "Y_NANOcontact.h"
#include "Y_NANObndCubical.h"
#include "Y_NANOatAr.h"
#include "Y_NANOsystem.h"


static void ContactWithCube(Y_NANObndCubical *v0bnd,Y_NANOatAr *v0at)
{ DBL dadb;             /* AtomDistanceBoundary                             */
  DBL dforc;            /* FORCe                                            */
  DBL dforcx;           /* FORCeX                                           */
  DBL dforcy;           /* FORCeY                                           */
  DBL dforcz;           /* FORCeZ                                           */
  DBL daccx;            /* AtomCenterCoordinateX                            */
  DBL daccy;            /* AtomCenterCoordinateY                            */
  DBL daccz;            /* AtomCenterCoordinateZ                            */
  DBL dcc1x;            /* CubePoint1CoordinateX                            */
  DBL dcc1y;            /* CubePoint1CoordinateY                            */
  DBL dcc1z;            /* CubePoint1CoordinateZ                            */
  DBL dcc2x;            /* CubePoint2CoordinateX                            */
  DBL dcc2y;            /* CubePoint2CoordinateY                            */
  DBL dcc2z;            /* CubePoint2CoordinateZ                            */
  DBL drefx;            /* ReferencePointX                                  */
  DBL drefy;            /* ReferencePointY                                  */
  DBL drefz;            /* ReferencePointZ                                  */
  DBL dpot;             /* AtomPotentialEnergy                              */
  DBL dtmp;             /* TMPVariable                                      */
  DBL dheatcurcont;     /* HeatCurrentContactPartModulus                    */
  INT ifcontact;        /* ContactFlag                                      */

/* Extract Atoms Variables                                                  */
  daccx=v0at->daccx;
  daccy=v0at->daccy;
  daccz=v0at->daccz;

/* Extract Cube Variables                                                   */
  dcc1x=v0bnd->dcc1x;
  dcc1y=v0bnd->dcc1y;
  dcc1z=v0bnd->dcc1z;
  dcc2x=v0bnd->dcc2x;
  dcc2y=v0bnd->dcc2y;
  dcc2z=v0bnd->dcc2z;

/****************************************************************************/
/* Contact with FACE A                                                      */
/****************************************************************************/

/* Reference Point */
  drefx=dcc1x;

  dadb=DABS(daccx-drefx);
  ifcontact=v0sys->pt2BNDArContactFunction(dadb,&dpot,&dforc);
  if(ifcontact==YES)
  { v0sys->dEPb=v0sys->dEPb+dpot;
    v0sys->iconbnd=v0sys->iconbnd+1;

/* Calculate The Total Pressure in The System                               */
    v0sys->dpress=v0sys->dpress+dforc;

    dforcx=dforc;

    v0at->davcx=v0at->davcx+dforcx/v0sys->v0ArDLJ->dama*v0sys->deltat;
    v0at->daepo=v0at->daepo+DP5*dpot;
    dtmp=v0at->davcxo;
    dheatcurcont=DP5*dadb*dforc*dtmp;
    v0sys->dheatcurcontx=v0sys->dheatcurcontx+dheatcurcont;
    v0sys->dTxx=v0sys->dTxx+DP5*dadb*dforc;
  }


/****************************************************************************/
/* Contact with FACE B                                                      */
/****************************************************************************/

/* Reference Point */
  drefz=dcc1z;

  dadb=DABS(daccz-drefz);
  ifcontact=v0sys->pt2BNDArContactFunction(dadb,&dpot,&dforc);
  if(ifcontact==YES)
  { v0sys->dEPb=v0sys->dEPb+dpot;
    v0sys->iconbnd=v0sys->iconbnd+1;

/* Calculate The Total Pressure in The System                               */
    v0sys->dpress=v0sys->dpress+dforc;

    dforcz=dforc;

    v0at->davcz=v0at->davcz+dforcz/v0sys->v0ArDLJ->dama*v0sys->deltat;
    v0at->daepo=v0at->daepo+DP5*dpot;
    dtmp=v0at->davczo;
    dheatcurcont=DP5*dadb*dforc*dtmp;
    v0sys->dheatcurcontz=v0sys->dheatcurcontz+dheatcurcont;
    v0sys->dTzz=v0sys->dTzz+DP5*dadb*dforc;
  }

/****************************************************************************/
/* Contact with FACE C                                                      */
/****************************************************************************/

/* Reference Point */
  drefx=dcc2x;

  dadb=DABS(drefx-daccx);
  ifcontact=v0sys->pt2BNDArContactFunction(dadb,&dpot,&dforc);
  if(ifcontact==YES)
  { v0sys->dEPb=v0sys->dEPb+dpot;
    v0sys->iconbnd=v0sys->iconbnd+1;

/* Calculate The Total Pressure in The System                               */
    v0sys->dpress=v0sys->dpress+dforc;

    dforcx=-dforc;

    v0at->davcx=v0at->davcx+dforcx/v0sys->v0ArDLJ->dama*v0sys->deltat;
    v0at->daepo=v0at->daepo+DP5*dpot;
    dtmp=-v0at->davcxo;
    dheatcurcont=DP5*dadb*dforc*dtmp;
    v0sys->dheatcurcontx=v0sys->dheatcurcontx+dheatcurcont;
    v0sys->dTxx=v0sys->dTxx-DP5*dadb*dforc;
  }

/****************************************************************************/
/* Contact with FACE D                                                      */
/****************************************************************************/

/* Reference Point */
  drefz=dcc2z;

  dadb=DABS(daccz-drefz);
  ifcontact=v0sys->pt2BNDArContactFunction(dadb,&dpot,&dforc);
  if(ifcontact==YES)
  { v0sys->dEPb=v0sys->dEPb+dpot;
    v0sys->iconbnd=v0sys->iconbnd+1;

/* Calculate The Total Pressure in The System                               */
    v0sys->dpress=v0sys->dpress+dforc;

    dforcz=-dforc;

    v0at->davcz=v0at->davcz+dforcz/v0sys->v0ArDLJ->dama*v0sys->deltat;
    v0at->daepo=v0at->daepo+DP5*dpot;
    dtmp=-v0at->davczo;
    dheatcurcont=DP5*dadb*dforc*dtmp;
    v0sys->dheatcurcontz=v0sys->dheatcurcontz+dheatcurcont;
    v0sys->dTzz=v0sys->dTzz-DP5*dadb*dforc;
  }


/****************************************************************************/
/* Contact with FACE E                                                      */
/****************************************************************************/

/* Reference Point */
  drefy=dcc1y;

  dadb=DABS(daccy-drefy);
  ifcontact=v0sys->pt2BNDArContactFunction(dadb,&dpot,&dforc);
  if(ifcontact==YES)
  { v0sys->dEPb=v0sys->dEPb+dpot;
    v0sys->iconbnd=v0sys->iconbnd+1;

/* Calculate The Total Pressure in The System                               */
    v0sys->dpress=v0sys->dpress+dforc;

    dforcy=dforc;

    v0at->davcy=v0at->davcy+dforcy/v0sys->v0ArDLJ->dama*v0sys->deltat;
    v0at->daepo=v0at->daepo+DP5*dpot;
    dtmp=v0at->davcyo;
    dheatcurcont=DP5*dadb*dforc*dtmp;
    v0sys->dheatcurconty=v0sys->dheatcurconty+dheatcurcont;
    v0sys->dTyy=v0sys->dTyy+DP5*dadb*dforc;
  }

/****************************************************************************/
/* Contact with FACE F                                                      */
/****************************************************************************/

/* Reference Point */
  drefy=dcc2y;

  dadb=DABS(daccy-drefy);
  ifcontact=v0sys->pt2BNDArContactFunction(dadb,&dpot,&dforc);
  if(ifcontact==YES)
  { v0sys->dEPb=v0sys->dEPb+dpot;
    v0sys->iconbnd=v0sys->iconbnd+1;

/* Calculate The Total Pressure in The System                               */
    v0sys->dpress=v0sys->dpress+dforc;

    dforcy=-dforc;

    v0at->davcy=v0at->davcy+dforcy/v0sys->v0ArDLJ->dama*v0sys->deltat;
    v0at->daepo=v0at->daepo+DP5*dpot;
    dtmp=-v0at->davcyo;
    dheatcurcont=DP5*dadb*dforc*dtmp;
    v0sys->dheatcurconty=v0sys->dheatcurconty+dheatcurcont;
    v0sys->dTyy=v0sys->dTyy-DP5*dadb*dforc;
  }
}

#define CDP5MASS 3.316762969 /* Constant 0.5*(atom mass) */

void ContactBnCubAr(Y_NANOobj *v0obj1, Y_NANOobj *v0obj2)
{ DBL dvelo2;           /* AtomVelocity^2                                   */
  DBL daEK;             /* AtomKineticEnergy                                */
  DBL daccx;            /* AtomCenterCoordinateX                            */
  DBL daccy;            /* AtomCenterCoordinateY                            */
  DBL daccz;            /* AtomCenterCoordinateZ                            */
                                                                            
  Y_NANObndCubical *v0bnd;/* Y_NANOpointerBouNDarycubical                         */
  Y_NANOatAr *v0at;        /* Y_NANOpointerATomargon                              */

  v0bnd=(Y_NANObndCubical *)v0obj1;
  v0at=(Y_NANOatAr *)v0obj2;

/* Extract Variables                                                        */
  daccx=v0at->daccx;
  daccy=v0at->daccy;
  daccz=v0at->daccz;

/* Calculate Kinetic Energy Of The Atom                                     */
  V3DLe2(dvelo2,v0at->davcx,v0at->davcy,v0at->davcz);
  daEK=CDP5MASS*dvelo2;
  v0at->daeki=daEK;
  v0sys->dEK=v0sys->dEK+daEK;

/* Resolve Contact between atom and boundary                                */
  ContactWithCube(v0bnd,v0at);
};


Y_NANOcontactBnCubAr::Y_NANOcontactBnCubAr(Y_NANOobj *v0obj1,Y_NANOobj *v0obj2,Y_NANOobj *v0obj3)
{ vtoolbox.Y_NANOregistertool(&contact,v0obj1,v0obj2,&ContactBnCubAr);
}

