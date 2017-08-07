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
#include "Y_NANOatAr.h"
#include "Y_NANOContaiBinList.h"
#include "Y_NANOprim.h"

void Y_NANOatAr::Initia(Y_NANOobj *v0obj1, Y_NANOobj *v0obj2, Y_NANOobj *v0obj3)
{ Y_NANOprim *v0pr;        /* PrimitiveObject                      */
  Y_NANOContaiBinList *v0cbl;/* ContainerBinaryList                */

  type=0;
  v0pr=(Y_NANOprim *)v0obj2;
  v0cbl=(Y_NANOContaiBinList *)v0obj1;
 
  if(v0cbl->v0next==v0cbl)
  { v0cbl->v0next=v0obj3;
    v0cbl->v0prev=v0obj3;
    v0prev=v0cbl;
    v0next=v0cbl;
  }
  else
  { v0prev=v0cbl->v0prev;
    v0cbl->v0prev->v0next=v0obj3;
    v0cbl->v0prev=v0obj3;
    v0next=v0cbl;
  }

  if(v0pr!=NULL)
  { id=v0pr->i1par[0];
    daccx=v0pr->d1par[0];
    daccy=v0pr->d1par[1];
    daccz=v0pr->d1par[2];
    davcx=v0pr->d1par[3];
    davcy=v0pr->d1par[4];
    davcz=v0pr->d1par[5];
    davcxo=v0pr->d1par[6];
    davcyo=v0pr->d1par[7];
    davczo=v0pr->d1par[8];
    daepo=v0pr->d1par[9];
    daeki=v0pr->d1par[10];
    //drad=1.7025;
    drad  =1.7005; // sigma/2
  }
}
