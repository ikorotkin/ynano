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
#include "ContactBUF14_7ForAll.h"
#include "Y_NANOframe.h"

INT ContactBUF14_7ForAll(DBL deps,DBL dConstA,DBL dConstB,DBL dConstC,DBL dConstD,
                         DBL dR2,DBL dRco,DBL dRco2,DBL dUco,DBL dFco,
                         DBL *dpote,DBL *dforce)
{ INT ifcontact;    /* ContactFlag                                          */
  DBL dR;           /* InterAtomicDistance                                  */
  DBL dR6;          /* InterAtomicDistance^6                                */
  DBL dR7;          /* InterAtomicDistance^7                                */
  DBL dtmp;         /* TMPVariable                                          */
  DBL dtmp2;        /* TMPVariable                                          */
  DBL dtmp7;        /* TMPVariable                                          */
  DBL dA;           /* AuxiliaryVariable                                    */
  DBL dB;           /* AuxiliaryVariable                                    */
  DBL dC;           /* AuxiliaryVariable                                    */
  DBL dpot;         /* PotentialEnergy                                      */
  DBL dforc;        /* Force                                                */

  ifcontact=NO;
  dpot=D0;
  dforc=D0;

  if(dR2<dRco2)
  { dR=DSQRT(dR2);
    dR6=dR2*dR2*dR2;
    dR7=dR6*dR;
    dtmp=dR+dConstB;
    dtmp2=dtmp*dtmp;
    dtmp7=dtmp2*dtmp2*dtmp2*dtmp;

    dA=dConstA/dtmp7;
    dC=dConstC/(dR7+dConstD);
    dB=dC-D2;

    dpot=deps*dA*dB-dUco+dR*dFco-dRco*dFco;
    dforc=D7*deps*dA*(dR6*dC/(dR7+dConstD)+dB/(dR+dConstB))-dFco;

    ifcontact=YES;
  }
  *dpote=dpot;
  *dforce=dforc;

  return ifcontact;
}

