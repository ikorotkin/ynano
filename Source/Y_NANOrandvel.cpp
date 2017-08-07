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
#include "Y_NANOrandvel.h"
#include "Y_NANOatAr.h"
#include "Y_NANOsystem.h"

Y_NANOrandvel::Y_NANOrandvel()
{ idum=-1;
}

void Y_NANOrandvel::getrandvel(Y_NANOobj *v0obj)
{ Y_NANOatAr *v0atar;      /* ATomARgon                            */
  DBL dsig3;            /* ValueOf3Sigmas                       */

  dsig3=D3*SQRT(v0sys->dkbolt*v0sys->dtemp/v0sys->v0ArDLJ->dama);
  v0atar=(Y_NANOatAr *)v0obj;
  v0atar->davcx=dsig3*TranNumGaussDisGen(&idum);//random digits in the interval -1 to 1
  v0atar->davcy=dsig3*TranNumGaussDisGen(&idum);
  v0atar->davcz=dsig3*TranNumGaussDisGen(&idum);
  v0atar->davcxo=v0atar->davcx;
  v0atar->davcyo=v0atar->davcy;
  v0atar->davczo=v0atar->davcz;
} 
