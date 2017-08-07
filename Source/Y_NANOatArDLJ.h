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
#ifndef Y_NANOATARDLJ
#define Y_NANOATARDLJ

#include "Y_NANOframe.h"
#include "Y_NANOobj.h"

class Y_NANOatArDLJ:public Y_NANOobj
{ public:
  Y_NANOatArDLJ()
  { dama= 6.633525938;
    //deps=16.5679;
    deps  =16.25067; // oplsaa_097 (GROMACS)
    //dara= 3.4050;
    dara  = 3.401;   // oplsaa_097 (GROMACS)
    dara2=dara*dara;
    dara6=dara2*dara2*dara2;
    dara12=dara6*dara6;
    //dacr=2.5*dara;
    dacr  =3.0*dara; // 10.203 A
    dacr2=dacr*dacr;
    dacr6=dacr2*dacr2*dacr2;
    dacr7=dacr6*dacr;
    dacr12=dacr6*dacr6;
    dacr13=dacr6*dacr7;
    dfcoff=D48*deps*dara12/dacr13-D24*deps*dara6/dacr7;
    ducoff= D4*deps*dara12/dacr12- D4*deps*dara6/dacr6;
    dconstA=D4*deps*dara12;
    dconstB=D4*deps*dara6;
    dconstC=dacr*dfcoff;
    dconstD=D48*deps*dara12;
    dconstE=D24*deps*dara6;
  };

  DBL dama;             /* AtomMAss                             */
  DBL deps;             /* Epsilon                              */
  DBL dara;             /* AtomRAdius                           */
  DBL dara2;            /* AtomRAdius^2                         */
  DBL dara6;            /* AtomRAdius^6                         */
  DBL dara12;           /* AtomRAdius^12                        */
  DBL daep;             /* AtomEPsilon                          */
  DBL dacr;             /* AtomCutoffRadius                     */
  DBL dacr2;            /* AtomCutoffRadius^2                   */
  DBL dacr6;            /* AtomCutoffRadius^6                   */
  DBL dacr7;            /* AtomCutoffRadius^7                   */
  DBL dacr12;           /* AtomCutoffRadius^12                  */
  DBL dacr13;           /* AtomCutoffRadius^13                  */
  DBL dfcoff;           /* Force@CutOffDistance                 */
  DBL ducoff;           /* Potential@CutOffDist                 */
  DBL dconstA;          /* ConstantA                            */
  DBL dconstB;          /* ConstantB                            */
  DBL dconstC;          /* ConstantC                            */
  DBL dconstD;          /* ConstantD                            */
  DBL dconstE;          /* ConstantE                            */
};

#endif

