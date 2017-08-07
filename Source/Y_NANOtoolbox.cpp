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
#include "Y_NANOtoolbox.h"

void Y_NANOtoolerror(Y_NANOobj *vobj1, Y_NANOobj *vobj2)
{ cout<<"Function not registered!!"<<endl;
  exit(0);
}

Y_NANOtoolbox::Y_NANOtoolbox()
{ INT i1;
  INT i2;
  INT i3;
   
  for(i1=0;i1<NFUN;i1++)
  { for(i2=0;i2<NATOMS;i2++)
    { for(i3=0;i3<NATOMS;i3++)
      { TableofFunPointers[i1][i2][i3]=&Y_NANOtoolerror;
  } } } 
}


void Y_NANOtoolbox::Y_NANOtool(Y_NANOobj *vobj1)
{ INT itype;
     
  itype=vobj1->type;

/*  TableofFunPointers[itype][NULL][NULL](NULL,NULL);*/
  TableofFunPointers[itype][0][0](NULL,NULL);
}

void Y_NANOtoolbox::Y_NANOtool(Y_NANOobj *vobj1,Y_NANOobj *vobj2)
{ INT itype;
  INT jtype;
     
  itype=vobj1->type;
  jtype=vobj2->type;

/*  TableofFunPointers[itype][jtype][NULL](vobj2,NULL);*/
  TableofFunPointers[itype][jtype][0](vobj2,NULL);
}

void Y_NANOtoolbox::Y_NANOtool(Y_NANOobj *vobj1,Y_NANOobj *vobj2,Y_NANOobj *vobj3)
{ INT itype;
  INT jtype;
  INT ktype;

  itype=vobj1->type;
  jtype=vobj2->type;
  ktype=vobj3->type;

  TableofFunPointers[itype][jtype][ktype](vobj2,vobj3);
}

void Y_NANOtoolbox::Y_NANOregistertool(Y_NANOobj *vobj1,Y_NANOobj *vobj2,Y_NANOobj *vobj3,
                                   void (*pfx)(Y_NANOobj *,Y_NANOobj *))
{ INT itype;
  INT jtype;
  INT ktype;

  itype=0;
  jtype=0;
  ktype=0;

  if(vobj1!=NULL)
  { itype=vobj1->type;
  }
  if(vobj2!=NULL)
  { jtype=vobj2->type;
  }
  if(vobj3!=NULL)
  { ktype=vobj3->type;
  }

  TableofFunPointers[itype][jtype][ktype]=pfx;
}

