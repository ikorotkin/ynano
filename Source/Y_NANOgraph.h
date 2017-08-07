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
#ifndef Y_NANOGRAPH
#define Y_NANOGRAPH

#include "Y_NANOframe.h"
#include "Y_NANOobj.h"
#define MMPARA 3000		/* M_PARAmeters      */

class Y_NANOgraph:public Y_NANOobj
{ public:
   Y_NANOgraph()
   { nipar=0;
     ndpar=0;
     type=0;
   }
   INT nipar;	/* Number_of_Int_parameters - generalized object*/
   INT ndpar;	/* Number_of_Dbl_parameters - generalized object*/
   INT i1par[MMPARA];/* I1_PARameters_of_generalised_object			*/
   DBL d1par[MMPARA];/* D1_PARameters_of_generalised_object			*/
};

#endif

