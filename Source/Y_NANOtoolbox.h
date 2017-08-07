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
#ifndef Y_NANOTOOLBOX
#define Y_NANOTOOLBOX

#include "Y_NANOframe.h"
#include "Y_NANOobj.h"
#include "Y_NANOsystem.h"

#define NFUN 20
#define NATOMS 20

class Y_NANOtoolbox
{ public:
   Y_NANOtoolbox();
   void Y_NANOtool(Y_NANOobj *vobj1);
   void Y_NANOtool(Y_NANOobj *vobj1,Y_NANOobj *vobj2);
   void Y_NANOtool(Y_NANOobj *vobj1,Y_NANOobj *vobj2,Y_NANOobj *vobj3);
   void Y_NANOregistertool(Y_NANOobj *vobj1,Y_NANOobj *vobj2,Y_NANOobj *vobj3,void (*pfx)(Y_NANOobj *,Y_NANOobj *));

  private:
   void (*TableofFunPointers[NFUN][NATOMS][NATOMS])(Y_NANOobj *,Y_NANOobj *);
};

void Y_NANOtoolerror(Y_NANOobj *vobj1, Y_NANOobj *vobj2);


#endif

