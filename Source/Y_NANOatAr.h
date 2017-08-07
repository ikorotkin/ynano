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
#ifndef Y_NANOATAR
#define Y_NANOATAR

#include "Y_NANOframe.h"
#include "Y_NANOobj.h"

class Y_NANOatAr:public Y_NANOobj
{ public:
   INT id;              /* AtomID                   */
   DBL daccx;           /* AtomCoordinateCurrentX   */
   DBL daccy;           /* AtomCoordinateCurrentY   */
   DBL daccz;           /* AtomCoordinateCurrentZ   */
   DBL davcx;           /* AtomVelocityCurrentX     */
   DBL davcy;           /* AtomVelocityCurrentY     */
   DBL davcz;           /* AtomVelocityCurrentZ     */
   DBL davcxo;          /* AtomVelocityCurrentXold  */
   DBL davcyo;          /* AtomVelocityCurrentYold  */
   DBL davczo;          /* AtomVelocityCurrentZold  */
   DBL daepo;           /* AtomPotentialEnergy      */
   DBL daeki;           /* AtomKineticEnergy        */
   DBL drad;            /* AtomRadius               */
   

  void Initia(Y_NANOobj *v0obj1, Y_NANOobj *v0obj2, Y_NANOobj *v0obj3);//initialize the atom Argon object binary list
};


#endif

