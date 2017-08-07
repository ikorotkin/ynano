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
#ifndef Y_NANOCONTAIBINLIST
#define Y_NANOCONTAIBINLIST

#include "Y_NANOframe.h"
#include "Y_NANOobj.h"
#include "Y_NANOtoolbox.h"
#include "Y_NANOconvert.h"
#include "Y_NANOdblpoint.h"
#include "Y_NANOdiapoint.h"
#include "Y_NANOintpoint.h"
#include "Y_NANOintegerize.h"
#include "Y_NANOcontact.h"
#include "Y_NANOmove.h"

class Y_NANOContaiBinList:public Y_NANOobj
{ public:
   Y_NANOContaiBinList()
   { v0next=this;
     v0prev=this;
     dcell=-D1;
   }

   INT ncelx;			 /* CellsInXDirection		*/
   INT ncely;			 /* CellsInYDirection		*/
   INT ncelz;            /* CellsInZDirection		*/
   DBL dminx;            /* Xminimum                */
   DBL dminy;            /* Yminimum                */
   DBL dminz;            /* Zminimum                */
   DBL dmaxx;            /* Xmaximum                */
   DBL dmaxy;            /* Ymaximum                */
   DBL dmaxz;            /* Zmaximum                */
   DBL dmaxD;            /* maximumDiameter         */
   DBL dcell;            /* cellsize                */

   void TLContact();     /* LinearSearch            */
   void Quadratic();     /* QuadraticSearch         */

  private:
   void SetSpace();	     /* SetSpaceBoundaries		*/ 
   void TQSort();        /* QuickSort               */
   void TLSort();        /* LinearSort              */
   void TLSearch();      /* LinearSearch            */
   Y_NANOobj ****TalY_NANOobj4(INT m3,INT m2,INT m1) ;  /* AllocateMatrixofPointers */

};

extern Y_NANOtoolbox vtoolbox;         /* Y_NANOtoolbox               */
extern Y_NANOconvert convert;          /* Y_NANOconvertInstance       */
extern Y_NANOintegerize integerize;	/* Y_NANOintegerizeInstance	*/
extern Y_NANOcontact contact;          /* Y_NANOcontactInstance       */
extern Y_NANOmove move_;				/* Y_NANOmoveInstance			*/

#endif

