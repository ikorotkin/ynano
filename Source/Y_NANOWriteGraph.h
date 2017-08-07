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
#ifndef Y_NANOWRGRPH
#define Y_NANOWRGRPH

#include "Y_NANOframe.h"
#include "Y_NANOgraph.h"
#include <string>
#include <fstream>
using namespace std;

class Y_NANOWriteGraph
{ private:
   INT *i1hat;                                      /* *I1HAshTable         */
   INT nprec;                                       /* NPRECision           */
   CHR *c1dig;                                      /* *C1baseDIGitcharacter*/
   CHR c1word[MCHPO];                               /* C1WORDcharacters     */
   void INTtoCHR(CHR *c1word, INT inum, INT *i0);   /* translateINTtoCHR    */
   void DBLtoCHR(CHR *c1word, DBL dnum, INT *i0);   /* translateDBLtoCHR    */
   ofstream foutp;                                  /* FOUTPutfile          */

  public:
   Y_NANOWriteGraph(CHR *f0name);                      /* Defaultconstructor   */
   ~Y_NANOWriteGraph();                                /* Defaultdestructor    */
   void WriteGrph(Y_NANOgraph *v0graph);
};

#endif


