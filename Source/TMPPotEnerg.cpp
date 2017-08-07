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

#include "TMPPotEnerg.h"
#include "Y_NANOsystem.h"
#include "Y_NANOobj.h"
#include "Y_NANOatAr.h"

/****************************************************************/
/* For Full Lenard-Jones Potential 3-9                          */
/****************************************************************/
DBL PotBAr(DBL dad,DBL dad3,DBL dad9)
{ DBL dpote;            /* PotentialEnergy                      */

  dpote=v0sys->v0BnD->dBnC1/dad9-v0sys->v0BnD->dBnC2/dad3-
        v0sys->v0BnD->ducoff+
        (dad-v0sys->v0BnD->dacr)*v0sys->v0BnD->dfcoff;

  return dpote;
}

DBL PotBArWCA(DBL dad6,DBL dad12)
{ DBL dpote;            /* PotentialEnergy                      */

  dpote=v0sys->v0ArDLJ->dconstA/dad12-v0sys->v0ArDLJ->dconstB/dad6+
        v0sys->v0ArDLJ->deps;

  return dpote;
}

DBL PotArAr(DBL dad,DBL dad6,DBL dad12)
{ DBL dpote;            /* PotentialEnergy                      */

  dpote=v0sys->v0ArDLJ->dconstA/dad12-
        v0sys->v0ArDLJ->dconstB/dad6-v0sys->v0ArDLJ->ducoff+
        dad*v0sys->v0ArDLJ->dfcoff-v0sys->v0ArDLJ->dconstC;

  return dpote;
}

INT ContactBNDAr39(DBL dadb,DBL *dpote,DBL *dforce)
{ INT ifcontact;        /* IfContactFlag                                    */
  DBL dadb2;            /* AtomDistanceBoundary^2                           */
  DBL dadb3;            /* AtomDistanceBoundary^3                           */
  DBL dadb4;            /* AtomDistanceBoundary^4                           */
  DBL dadb8;            /* AtomDistanceBoundary^8                           */
  DBL dadb9;            /* AtomDistanceBoundary^9                           */
  DBL dadb10;           /* AtomDistanceBoundary^10                          */
  DBL dforc;            /* FORCe                                            */
  DBL dpot;             /* PotentialEnergy                                  */

  ifcontact=NO;
  dpot=D0;
  dforc=D0;

  if(dadb<v0sys->v0BnD->dacr)
  { dadb2=dadb*dadb;
    dadb3=dadb2*dadb;
    dadb4=dadb2*dadb2;
    dadb8=dadb4*dadb4;
    dadb9=dadb8*dadb;
    dadb10=dadb8*dadb2;

    dpot=PotBAr(dadb,dadb3,dadb9);

    dforc=v0sys->v0BnD->dBnC3/dadb10-v0sys->v0BnD->dBnC4/dadb4-
          v0sys->v0BnD->dfcoff;

    ifcontact=YES;
  }

  *dpote=dpot;
  *dforce=dforc;

  return ifcontact;
}


INT ContactBNDArWCA(DBL dadb,DBL *dpote,DBL *dforce)
{ INT ifcontact;        /* IfContactFlag                                    */
  DBL dadb2;            /* AtomDistanceBoundary^2                           */
  DBL dadb6;            /* AtomDistanceBoundary^6                           */
  DBL dadb7;            /* AtomDistanceBoundary^7                           */
  DBL dadb12;           /* AtomDistanceBoundary^12                          */
  DBL dadb13;           /* AtomDistanceBoundary^13                          */
  DBL dforc;            /* FORCe                                            */
  DBL dpot;             /* PotentialEnergy                                  */

  ifcontact=NO;
  dpot=D0;
  dforc=D0;

  if(dadb<v0sys->v0BnD->dacrWCA)
  { dadb2=dadb*dadb;
    dadb6=dadb2*dadb2*dadb2;
    dadb7=dadb6*dadb;
    dadb12=dadb6*dadb6;
    dadb13=dadb12*dadb;

    dpot=PotBArWCA(dadb6,dadb12);

    dforc=v0sys->v0ArDLJ->dconstD/dadb13-v0sys->v0ArDLJ->dconstE/dadb7;

    ifcontact=YES;
  }

  *dpote=dpot;
  *dforce=dforc;

  return ifcontact;
}

