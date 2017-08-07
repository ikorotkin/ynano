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
#ifndef Y_NANOIN
#define Y_NANOIN

#include "Y_NANOframe.h"
#include "Y_NANOobj.h"
#include "Y_NANOprim.h"
#include "Y_NANOsystem.h"

#include <string>
#include <fstream>
using namespace std;

class Y_NANOin
{ private:
   INT i1hat[MHASH]; 	/* I1_HAsh_Table				*/
   INT nprec;	        /* N_PRECision					*/
   CHR c1dig[MBASE];	/* C1_base_DIGit_characters			*/
   CHR c1word[MCHPO];	/* C1_WORD_characters				*/
   INT CHRtoTYPE(CHR *c1word, INT *i0);	/* translate CHR to TYPE number	*/
   INT CHRtoINT(CHR *c1word, INT *i0);/* translate CHR to INT number	*/
   INT IsTag(CHR *c1word, INT *i0);/* checks for TAG characters		*/
   DBL CHRtoDBL(CHR *c1word, INT *i0);/* translate CHR to DBL number	*/
   INT ReadFil(CHR *c1word);	/* Reads a file filtering rubbish	*/
   ifstream finpu;	/* FINput file					*/

  public:
   INT nipar;	/* Number_of_Int_parameters - generalized object	*/
   INT ndpar;	/* Number_of_Dbl_parameters - generalized object	*/
   INT ieoff;	/* IntegerEndOFFileflag					*/

   Y_NANOin(CHR *f0name);	/* Default constructor				*/
   void ReadObj(Y_NANOprim *v0prim);/* Reads objects from a file		*/
};

#endif

