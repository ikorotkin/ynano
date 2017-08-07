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
#ifndef Y_NANOATARDBUF14_7
#define Y_NANOATARDBUF14_7

#include "Y_NANOframe.h"
#include "Y_NANOobj.h"

class Y_NANOatArDBuf14_7:public Y_NANOobj
{ public:
  Y_NANOatArDBuf14_7()
  { dama= 6.633525938;
    deps=19.78986329;
    dara= 3.4050;
    dara2=dara*dara;
    dara6=dara2*dara2*dara2;
    dara12=dara6*dara6;
    dacr12=dacr6*dacr6;
    dacr13=dacr6*dacr7;

    dRp=3.761;
    dRp2=dRp*dRp;
    dRp4=dRp2*dRp2;
    dRp7=dRp4*dRp2*dRp;

    dacr=2.5*dara;
    dacr2=dacr*dacr;
    dacr6=dacr2*dacr2*dacr2;
    dacr7=dacr6*dacr;

    dAco1=1.07*dRp/(dacr+0.07*dRp);
    dAco=dAco1*dAco1*dAco1*dAco1*dAco1*dAco1*dAco1;
    dCco=1.12*dRp7/(dacr7+0.12*dRp7);
    dBco=dCco-D2;

    ducoff=deps*dAco*dBco;
    dfcoff=D7*dAco*deps*(dacr6*dCco/(dacr7+0.12*dRp7)+dBco/(dacr+0.07*dRp));

    dconstA=pow((1.07*dRp),7);
    dconstB=0.07*dRp;
    dconstC=1.12*dRp7;
    dconstD=0.12*dRp7;
  };

  DBL dama;             /* AtomMAss                             */
  DBL deps;             /* Epsilon                              */
  DBL dara;             /* AtomRAdius                           */
  DBL dara2;            /* AtomRAdius^6                         */
  DBL dara6;            /* AtomRAdius^6                         */
  DBL dara12;           /* AtomRAdius^12                        */
  DBL daep;             /* AtomEPsilon                          */
  DBL dacr;             /* AtomCutoffRadius                     */
  DBL dacr2;            /* AtomCutoffRadius^2                   */
  DBL dacr6;            /* AtomCutoffRadius^6                   */
  DBL dacr7;            /* AtomCutoffRadius^7                   */
  DBL dacr12;           /* AtomCutoffRadius^12                  */
  DBL dacr13;           /* AtomCutoffRadius^13                  */
  DBL dRp;              /* RPrime                               */
  DBL dRp2;             /* RPrime^2                             */
  DBL dRp4;             /* RPrime^4                             */
  DBL dRp7;             /* RPrime^7                             */
  DBL dAco1;            /* CoefficientACutOff1                  */
  DBL dAco;             /* CoefficientACutOff                   */
  DBL dBco;             /* CoefficientBCutOff                   */
  DBL dCco;             /* CoefficientCCutOff                   */
  DBL dfcoff;           /* Force@CutOffDistance                 */
  DBL ducoff;           /* Potential@CutOffDist                 */
  DBL dconstA;          /* ConstantA                            */
  DBL dconstB;          /* ConstantB                            */
  DBL dconstC;          /* ConstantC                            */
  DBL dconstD;          /* ConstantD                            */
};

#endif

