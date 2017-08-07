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
#ifndef TMPPOTENERG
#define TMPPOTENERG

#include "Y_NANOframe.h"
#include "Y_NANOContaiBinList.h"
#include "Y_NANOsystem.h"

DBL PotArAr(DBL dad,DBL dad6,DBL dad12);

DBL PotBArWCA(DBL dad6,DBL dad12);
DBL PotBAr(DBL dad,DBL dad3,DBL dad9);
INT ContactBNDAr39(DBL dadb,DBL *dpote,DBL *dforce);
INT ContactBNDArWCA(DBL dadb,DBL *dpote,DBL *dforce);

#define FORCEBOUND(dforce,dr10,dr4)\
                {dforce=v0sys->v0BnD->dBnC3/dr10-\
                        v0sys->v0BnD->dBnC4/dr4-\
                        v0sys->v0BnD->dfcoff;}

extern Y_NANOsystem *v0sys;        /* Y_NANOsystemInstance        */

#endif

