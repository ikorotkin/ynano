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
#include "Y_NANOArtoPrim.h"
#include "Y_NANOprim.h"
#include "Y_NANOatAr.h"

void ArtoPrim(Y_NANOobj *v0obj1, Y_NANOobj *v0obj2)
{ Y_NANOprim *v0pr;        /* Y_NANOPointertoPRimitive                */
  Y_NANOatAr *v0at;        /* Y_NANOPointertoATom                     */

  v0pr=(Y_NANOprim *)v0obj1;
  v0at=(Y_NANOatAr *)v0obj2;

  v0pr->type=v0at->type;
  v0pr->nipar=2;
  v0pr->i1par[0]=v0at->id;
  v0pr->ndpar=6;
  v0pr->d1par[0]=v0at->daccx;
  v0pr->d1par[1]=v0at->daccy;
  v0pr->d1par[2]=v0at->daccz;
  v0pr->d1par[3]=v0at->davcx;
  v0pr->d1par[4]=v0at->davcy;
  v0pr->d1par[5]=v0at->davcz;
  v0pr->d1par[6]=v0at->davcxo;
  v0pr->d1par[7]=v0at->davcyo;
  v0pr->d1par[8]=v0at->davczo;
  v0pr->d1par[9]=v0at->daepo;
  v0pr->d1par[10]=v0at->daeki;
};

Y_NANOArtoPrim::Y_NANOArtoPrim(Y_NANOobj *v0obj1,Y_NANOobj *v0obj2,Y_NANOobj *v0obj3)
{ vtoolbox.Y_NANOregistertool(&converttoprim,v0obj1,v0obj2,&ArtoPrim);
}

