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
#include "TMPPotEnerg.h"
#include "Y_NANOsystem.h"
#include "Y_NANOobj.h"
#include "Y_NANOatAr.h"

void TQuadratic(Y_NANOContaiBinList *v0cbinlist)
{ Y_NANOobj *v0con;        /* Contactor                            */
  Y_NANOobj *v0tar;        /* Target                               */

  v0con=v0cbinlist->v0next;
  while(v0con!=v0cbinlist)
  { v0tar=v0con->v0next;
    while(v0tar!=v0cbinlist)
    { vtoolbox.Y_NANOtool(&contact,v0con,v0tar);

      v0tar=v0tar->v0next;
    }
    vtoolbox.Y_NANOtool(&move_,v0con);
    v0con=v0con->v0next;
  }


}

