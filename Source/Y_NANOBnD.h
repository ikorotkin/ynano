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
#ifndef Y_NANOBND
#define Y_NANOBND

#include "Y_NANOframe.h"
#include "Y_NANOobj.h"

#define DCOEF 0.99930625923385662755 /* 0.9718455053 -- Calculated as 4*(sigma/a)^3 where
                            sigma is the atom diameter and a is the
                            lattice constant for solid argon. In this
                            case a=5.4 10^-10 m.                            */

class Y_NANOBnD:public Y_NANOobj
{ public:
  Y_NANOBnD()
  {
    //deps=16.5679;
    deps  =16.25067; // oplsaa_097 (GROMACS)
    //dara= 3.4050;
    dara  = 3.401;   // oplsaa_097 (GROMACS)
    dara2=dara*dara;
    dara3=dara2*dara;
    dara9=dara3*dara3*dara3;
    //dacr=5.0*dara;
    dacr  =5.0*dara;
    dacrWCA=pow(D2,(D1/D6))*dara;
    dacr2=dacr*dacr;
    dacr3=dacr2*dacr;
    dacr4=dacr2*dacr2;
    dacr8=dacr4*dacr4;
    dacr9=dacr8*dacr;
    dacr10=dacr8*dacr2;
    dfcoff=DCOEF*((D4*deps*MYPI*dara9)/(D5*dacr10)-(D2*deps*MYPI*dara3)/dacr4);
    ducoff=DCOEF*((D4*MYPI*deps*dara9)/(D45*dacr9)-
                  (D2*MYPI*deps*dara3)/(D3*dacr3));
    dBnC1=(DCOEF*(D4*MYPI*deps*dara9))/D45;
    dBnC2=(DCOEF*(D2*MYPI*deps*dara3))/D3;
    dBnC3=(DCOEF*(D4*MYPI*deps*dara9))/D5;
    dBnC4=DCOEF*D2*MYPI*deps*dara3;
  };

  DBL deps;             /* Epsilon                              */
  DBL dara;             /* AtomRAdius                           */
  DBL dara2;            /* AtomRAdius^2                         */
  DBL dara3;            /* AtomRAdius^3                         */
  DBL dara9;            /* AtomRAdius^9                         */
  DBL dacrWCA;          /* AtomCutOffRadiusWCAPotential         */
  DBL dacr;             /* AtomCutoffRadius                     */
  DBL dacr2;            /* AtomCutoffRadius^2                   */
  DBL dacr3;            /* AtomCutoffRadius^3                   */
  DBL dacr4;            /* AtomCutoffRadius^4                   */
  DBL dacr8;            /* AtomCutoffRadius^8                   */
  DBL dacr9;            /* AtomCutoffRadius^9                   */
  DBL dacr10;           /* AtomCutoffRadius^10                  */
  DBL dfcoff;           /* Force@CutOffDistance                 */
  DBL ducoff;           /* Potential@CutOffDist                 */
  DBL dBnC1;            /* BoundaryConstant1                    */
  DBL dBnC2;            /* BoundaryConstant2                    */
  DBL dBnC3;            /* BoundaryConstant3                    */
  DBL dBnC4;            /* BoundaryConstant4                    */
};

#endif

