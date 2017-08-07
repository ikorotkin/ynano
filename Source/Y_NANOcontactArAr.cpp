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
#include "Y_NANOcontactArAr.h"
#include "TMPPotEnerg.h"
#include "Y_NANOtoolbox.h"
#include "Y_NANOcontact.h"
#include "Y_NANOatAr.h"
#include "Y_NANOsystem.h"
#include "ContactBUF14_7ForAll.h"


INT ContactArArLJ12_6(DBL dad2,DBL *dpote,DBL *dforce)
{ INT ifcontact;        /* ContactFlag                                      */
  DBL dpot;             /* PotentialEnergy                                  */
  DBL dforc;            /* Force                                            */
  DBL dad;              /* DoubleAtomDistance                               */
  DBL dad6;             /* DoubleAtomDistance^6                             */
  DBL dad7;             /* DoubleAtomDistance^7                             */
  DBL dad12;            /* DoubleAtomDistance^12                            */
  DBL dad13;            /* DoubleAtomDistance^13                            */

  ifcontact=NO;
  dpot=D0;
  dforc=D0;

  if(dad2<v0sys->v0ArDLJ->dacr2)
  { dad=DSQRT(dad2);

    dad6=dad2*dad2*dad2;
    dad12=dad6*dad6;
    dad7=dad6*dad;
    dad13=dad12*dad;

/* Calculate Potential Energy Between The Atoms                             */
    dpot=v0sys->v0ArDLJ->dconstA/dad12-
         v0sys->v0ArDLJ->dconstB/dad6-v0sys->v0ArDLJ->ducoff-
         dad*v0sys->v0ArDLJ->dfcoff+v0sys->v0ArDLJ->dconstC;

    dforc=v0sys->v0ArDLJ->dconstD/dad13-v0sys->v0ArDLJ->dconstE/dad7-
          v0sys->v0ArDLJ->dfcoff;

    ifcontact=YES;
  }

  *dpote=dpot;
  *dforce=dforc;
  return ifcontact;
}


INT ContactArArBUF14_7(DBL dad2,DBL *dpote,DBL *dforce)
{ INT ifcontact;        /* ContactFlag                                      */

  ifcontact=ContactBUF14_7ForAll(v0sys->v0ArDBUF14_7->deps,
              v0sys->v0ArDBUF14_7->dconstA,v0sys->v0ArDBUF14_7->dconstB,
              v0sys->v0ArDBUF14_7->dconstC,v0sys->v0ArDBUF14_7->dconstD,
              dad2,v0sys->v0ArDBUF14_7->dacr,
              v0sys->v0ArDBUF14_7->dacr2,v0sys->v0ArDBUF14_7->ducoff,
              v0sys->v0ArDBUF14_7->dfcoff,dpote,dforce);

  return ifcontact;
}

void ContactArAr(Y_NANOobj *v0obj1, Y_NANOobj *v0obj2)
{ INT ifcontact;        /* ContactFlag                                      */
  DBL dad;              /* DoubleAtomDistance                               */
  DBL dad2;             /* DoubleAtomDistance^2                             */
  DBL dadx;             /* DoubleAtomDistanceX                              */
  DBL dady;             /* DoubleAtomDistanceY                              */
  DBL dadz;             /* DoubleAtomDistanceZ                              */
  DBL dad0x;            /* DoubleUnitVectorX                                */
  DBL dad0y;            /* DoubleUnitVectorY                                */
  DBL dad0z;            /* DoubleUnitVectorZ                                */
  DBL dad0xdt;          /* UnitVectorX*dt                                   */
  DBL dad0ydt;          /* UnitVectorY*dt                                   */
  DBL dad0zdt;          /* UnitVectorZ*dt                                   */
  DBL dforc;            /* DoubleFORCe                                      */
  DBL dpot;             /* PotentialEnergy                                  */
  DBL dscale;           /* ScalingFactor                                    */
  DBL dinvmass;         /* InverseOfMassAtom1                               */
  DBL dtmp;             /* TMPVariable                                      */
  DBL dheatcurcont;     /* HeatCurrentContactPartModulus                    */
  Y_NANOatAr *v0at1;       /* Y_NANOpointerATomargon1                             */
  Y_NANOatAr *v0at2;       /* Y_NANOpointerATomargon2                             */

  v0at1=(Y_NANOatAr *)v0obj1;
  v0at2=(Y_NANOatAr *)v0obj2;

  dadx=v0at2->daccx-v0at1->daccx;
  dady=v0at2->daccy-v0at1->daccy;
  dadz=v0at2->daccz-v0at1->daccz;
  V3DLe2(dad2,dadx,dady,dadz);

  ifcontact=v0sys->pt2ArArContactFunction(dad2,&dpot,&dforc);
  if(ifcontact==YES)
  { dad=DSQRT(dad2);

/* Calculate Versor                                                         */
    dscale=D1/dad;
    dad0x=dadx*dscale;
    dad0y=dady*dscale;
    dad0z=dadz*dscale;

    v0sys->iconato=v0sys->iconato+1;
    v0sys->dEPa=v0sys->dEPa+dpot;

/* Update Velocities                                                        */
    dad0xdt=dad0x*v0sys->deltat;
    dad0ydt=dad0y*v0sys->deltat;
    dad0zdt=dad0z*v0sys->deltat;

    dinvmass=D1/v0sys->v0ArDLJ->dama;
    v0at1->davcx=v0at1->davcx-dforc*dinvmass*dad0xdt;
    v0at1->davcy=v0at1->davcy-dforc*dinvmass*dad0ydt;
    v0at1->davcz=v0at1->davcz-dforc*dinvmass*dad0zdt;
    v0at1->daepo=v0at1->daepo+DP5*dpot;

    v0at2->davcx=v0at2->davcx+dforc*dinvmass*dad0xdt;
    v0at2->davcy=v0at2->davcy+dforc*dinvmass*dad0ydt;
    v0at2->davcz=v0at2->davcz+dforc*dinvmass*dad0zdt;
    v0at2->daepo=v0at2->daepo+DP5*dpot;

    V3DDot(dtmp,dad0x,dad0y,dad0z,
                (v0at1->davcxo+v0at2->davcxo),
                (v0at1->davcyo+v0at2->davcyo),
                (v0at1->davczo+v0at2->davczo));
    dheatcurcont=DP5*dad*dforc*dtmp;

    v0sys->dheatcurcontx=v0sys->dheatcurcontx+dheatcurcont*dad0x;
    v0sys->dheatcurconty=v0sys->dheatcurconty+dheatcurcont*dad0y;
    v0sys->dheatcurcontz=v0sys->dheatcurcontz+dheatcurcont*dad0z;

    v0sys->dTxx=v0sys->dTxx+DP5*dadx*dforc*dad0x;
    v0sys->dTyy=v0sys->dTyy+DP5*dady*dforc*dad0y;
    v0sys->dTzz=v0sys->dTzz+DP5*dadz*dforc*dad0z;
    v0sys->dTxy=v0sys->dTxy+DP5*dady*dforc*dad0x;
    v0sys->dTxz=v0sys->dTxz+DP5*dadz*dforc*dad0x;
    v0sys->dTyz=v0sys->dTyz+DP5*dadz*dforc*dad0y;
    v0sys->dTyx=v0sys->dTyx+DP5*dadx*dforc*dad0y;
    v0sys->dTzx=v0sys->dTzx+DP5*dadx*dforc*dad0z;
    v0sys->dTzy=v0sys->dTzy+DP5*dady*dforc*dad0z;
} };

Y_NANOcontactArAr::Y_NANOcontactArAr(Y_NANOobj *v0obj1,Y_NANOobj *v0obj2,Y_NANOobj *v0obj3)
{ vtoolbox.Y_NANOregistertool(&contact,v0obj1,v0obj2,&ContactArAr);
}
