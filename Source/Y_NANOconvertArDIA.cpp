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
#include "Y_NANOconvertArDIA.h"
#include "Y_NANOtoolbox.h"
#include "Y_NANOconvert.h"
#include "Y_NANOatAr.h"
#include "Y_NANOdiapoint.h"

void ConvertArDIA(Y_NANOobj *v0obj1, Y_NANOobj *v0obj2)
{ Y_NANOatAr *v0at;        /* ATomargon                            */
  Y_NANOdiapoint *v0dp;    /* DiameterPoint                        */

  v0at=(Y_NANOatAr *)v0obj1;
  v0dp=(Y_NANOdiapoint *)v0obj2;
  
  v0dp->dadia=D2*v0sys->v0ArDLJ->dara;
  v0dp->dacdia=D2*v0sys->v0ArDLJ->dacr;
};

Y_NANOconvertArDIA::Y_NANOconvertArDIA(Y_NANOobj *v0obj1,Y_NANOobj *v0obj2,Y_NANOobj *v0obj3)
{ vtoolbox.Y_NANOregistertool(&convert,v0obj1,v0obj2,&ConvertArDIA);
}
