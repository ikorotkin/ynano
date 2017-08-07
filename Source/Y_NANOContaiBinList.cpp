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
#include "Y_NANOContaiBinList.h"
#include "Y_NANOtoolbox.h"
#include "Y_NANOconvert.h"
#include "Y_NANOdblpoint.h"
#include "Y_NANOdiapoint.h"
#include "Y_NANOintpoint.h"
#include "Y_NANOintegerize.h"
#include "Y_NANOcontact.h"
#include "Y_NANOmove.h"

#define ISAGTB(v0obj1,v0obj2)\
        (((v0obj1->iccz>v0obj2->iccz)||((v0obj1->iccz==v0obj2->iccz)&&\
         ((v0obj1->iccy>v0obj2->iccy)||((v0obj1->iccy==v0obj2->iccy)&&\
          (v0obj1->iccx>v0obj2->iccx))))))

#define ISALEB(v0obj1,v0obj2)\
        (((v0obj1->iccz<v0obj2->iccz)||((v0obj1->iccz==v0obj2->iccz)&&\
         ((v0obj1->iccy<v0obj2->iccy)||((v0obj1->iccy==v0obj2->iccy)&&\
          (v0obj1->iccx<=v0obj2->iccx))))))

#define ISALTB(v0obj1,v0obj2)\
        (((v0obj1->iccz<v0obj2->iccz)||((v0obj1->iccz==v0obj2->iccz)&&\
         ((v0obj1->iccy<v0obj2->iccy)||((v0obj1->iccy==v0obj2->iccy)&&\
          (v0obj1->iccx<v0obj2->iccx))))))

#define ISAGTCELL(v0obj1,ii,jj,kk)\
        (((v0obj1->iccz>(kk))||((v0obj1->iccz==(kk))&&\
         ((v0obj1->iccy>(jj))||((v0obj1->iccy==(jj))&&\
          (v0obj1->iccx>(ii)))))))

#define ISALECELL(v0obj1,ii,jj,kk)\
        (((v0obj1->iccz<(kk))||((v0obj1->iccz==(kk))&&\
         ((v0obj1->iccy<(jj))||((v0obj1->iccy==(jj))&&\
          (v0obj1->iccx<=(ii)))))))

#define ISALTCELL(v0obj1,ii,jj,kk)\
        (((v0obj1->iccz<(kk))||((v0obj1->iccz==(kk))&&\
         ((v0obj1->iccy<(jj))||((v0obj1->iccy==(jj))&&\
          (v0obj1->iccx<(ii)))))))

#define XXDSPSC     10  /* SpaceScale                           */

void Y_NANOContaiBinList::SetSpace()
{ DBL deltx;            /* SystemSizeXDirection                 */
  DBL delty;            /* SystemSizeYDirection                 */
  DBL deltz;            /* SystemSizeZDirection                 */
  Y_NANOobj *v0tmp;        /* TMPObject                            */
  Y_NANOdblpoint vdptmp;   /* DoublePointTMP                       */

  dmaxD=D0;
  dmaxx=EPSILON;
  dmaxy=EPSILON;
  dmaxz=EPSILON;
  dminx=BEPSILON;
  dminy=BEPSILON;
  dminz=BEPSILON;
  
  v0tmp=v0next;
  while(v0tmp!=this)
  { vtoolbox.Y_NANOtool(&convert,v0tmp,&vdptmp);//void ConvertArDbP(Y_NANOobj *v0obj1, Y_NANOobj *v0obj2)
    dmaxD=MAXIM(vdptmp.ddia,dmaxD);
    dmaxx=MAXIM(vdptmp.dx,dmaxx);
    dmaxy=MAXIM(vdptmp.dy,dmaxy);
    dmaxz=MAXIM(vdptmp.dz,dmaxz);
    dminx=MINIM(vdptmp.dx,dminx);
    dminy=MINIM(vdptmp.dy,dminy);
    dminz=MINIM(vdptmp.dz,dminz);
    v0tmp=v0tmp->v0next; 
  }
  dcell=dmaxD;
  deltx=dmaxx-dminx;
  delty=dmaxy-dminy;
  deltz=dmaxz-dminz;

  dminx=dminx-XXDSPSC*deltx;
  dminy=dminy-XXDSPSC*delty;
  dminz=dminz-XXDSPSC*deltz;
  dmaxx=dmaxx+XXDSPSC*deltx;
  dmaxy=dmaxy+XXDSPSC*delty;
  dmaxz=dmaxz+XXDSPSC*deltz;

  ncelx=((INT)((dmaxx-dminx)/dcell))+1;
  ncely=((INT)((dmaxy-dminy)/dcell))+1;
  ncelz=((INT)((dmaxz-dminz)/dcell))+1;

  v0tmp=v0next;
  while(v0tmp!=this)
  { vtoolbox.Y_NANOtool(&integerize,v0tmp,this);
    v0tmp=v0tmp->v0next; 
  }
}

#define Y_NANOobj4NULL ((Y_NANOobj****)NULL)

#define SUMPOINT(imed,jmed,kmed,v0obj1,v0obj2)\
      { imed=(v0obj1->iccx+v0obj2->iccx);\
        jmed=(v0obj1->iccy+v0obj2->iccy);\
        kmed=(v0obj1->iccz+v0obj2->iccz);}

#define HALFPOINT(imed,jmed,kmed,nx,ny)\
      { kmed=kmed+(jmed+imed/nx)/ny;\
        jmed=(jmed+imed/nx)%ny;\
        imed=imed%nx;\
        imed=(imed+((jmed+(kmed%2)*ny)%2)*nx)/2;\
        jmed=(jmed+(kmed%2)*ny)/2;\
        kmed=kmed/2;}

void Y_NANOContaiBinList::TQSort()
{ INT istuck;
  INT imed;             /* mediumintegerizedXcoord              */
  INT jmed;             /* mediumintegerizedYcoord              */
  INT kmed;             /* mediumintegerizedZcoord              */

  Y_NANOobj *firstobj[50]; /* arrayoffirstobjects                  */
  Y_NANOobj *lastobj[50];  /* arrayoflastobjects                   */
  Y_NANOobj *v0A;          /* parserpointerfromfirst               */
  Y_NANOobj *v0B;          /* parserpointerfromlast                */
  Y_NANOobj *v0tmp;        /* TMPobject                            */
  Y_NANOobj *v0beg;        /* RealBeginningOfArray                 */
  Y_NANOobj *v0end;        /* RealEndOfArray                       */
  Y_NANOobj *v0max;        /* MaximumObject                        */
  Y_NANOobj *v0min;        /* MinimumObject                        */


  SetSpace();           /* SetSystemSpace                       */
  firstobj[0]=v0next;
  lastobj[0]=v0prev;
  istuck=0;
  while(istuck>=0)
  { v0A=firstobj[istuck];/* GetArrayFromStuck                   */
    v0B=lastobj[istuck];
    v0beg=v0A;
    v0end=v0B;
    v0max=v0B;
    v0min=v0max;
    v0tmp=v0A;
    while(v0tmp!=v0B)
    { if(ISAGTB(v0tmp,v0max)) v0max=v0tmp;
      if(ISALEB(v0tmp,v0min)) v0min=v0tmp;
      v0tmp=v0tmp->v0next;
    }
    SUMPOINT(imed,jmed,kmed,v0min,v0max);
    HALFPOINT(imed,jmed,kmed,ncelx,ncely);
    while(v0A!=v0B)
    { while((v0A!=v0B) && ISALECELL(v0A,imed,jmed,kmed))
      { v0A=v0A->v0next;
      }
      while((v0A!=v0B) && ISAGTCELL(v0B,imed,jmed,kmed))
      { v0B=v0B->v0prev;
      }
      if(v0A!=v0B)
      { if(v0A==v0beg) v0beg=v0B;
        if(v0B==v0end) v0end=v0A;

        v0A->v0prev->v0next=v0B;
        v0B->v0next->v0prev=v0A;
        if(v0A->v0next!=v0B)
        { v0A->v0next->v0prev=v0B;
          v0B->v0prev->v0next=v0A;
          v0tmp=v0B->v0next;
          v0B->v0next=v0A->v0next;
          v0A->v0next=v0tmp;
          v0tmp=v0B->v0prev;
          v0B->v0prev=v0A->v0prev;
          v0A->v0prev=v0tmp;
        }
        else
        { v0A->v0next=v0B->v0next;
          v0B->v0next=v0A;
          v0B->v0prev=v0A->v0prev;
          v0A->v0prev=v0B;
        }
        v0tmp=v0A;
        v0A=v0B;
        v0B=v0tmp;
    } }
    if((v0beg!=v0A->v0prev) && (v0end!=v0A))
    { lastobj[istuck]=v0end;
      firstobj[istuck+1]=v0beg;
      lastobj[istuck+1]=v0A->v0prev;
      firstobj[istuck]=v0A;
      istuck=istuck+1;
    }
    else if((v0beg==v0A->v0prev) && (v0end!=v0A))
    { firstobj[istuck]=v0A;
      lastobj[istuck]=v0end;
    }
    else if((v0beg!=v0A->v0prev) && (v0end==v0A))
    { firstobj[istuck]=v0beg;
      lastobj[istuck]=v0A->v0prev;
    }
    else
    { istuck=istuck-1;
    }
} }

Y_NANOobj **** Y_NANOContaiBinList::TalY_NANOobj4(INT m3, INT m2, INT m1) 
{ INT isize,i2,i3;
  Y_NANOobj     **p1b;
  Y_NANOobj    ***p2b,  ***p2e, ***p2;
  Y_NANOobj   ****p3b, ****p3e;
  void    *v1;

  isize=sizeof(Y_NANOobj***)*(m3+3)+
        sizeof(Y_NANOobj** )*(m3*m2+3)+
        sizeof(Y_NANOobj*  )*(m3*m2*m1+3);
  if(isize==0)return Y_NANOobj4NULL;
  v1=MALLOC(isize);
  p3b=(Y_NANOobj****)v1; p3e=p3b+m3+1;
  p2b=(Y_NANOobj***)p3e; p2e=p2b+m3*m2+1;
  p1b=(Y_NANOobj**)p2e;
  p2=p2b;    
  for(i3=0;i3<m3;i3++)
  { p2=p2b+i3*m2;
    p3b[i3]=p2;
    for(i2=0;i2<m2;i2++)
    { p2[i2]=p1b+i3*m2*m1+i2*m1;
  } }
  return p3b; 
}

void Y_NANOContaiBinList::TLSort()
{ INT i;                /* LoopVariable                         */
  INT j;                /* LoopVariable                         */
  INT k;                /* LoopVariable                         */
  INT ixold;            /* OldIntegerizedCoordinateX            */
  INT iyold;            /* OldIntegerizedCoordinateY            */
  INT izold;            /* OldIntegerizedCoordinateZ            */
  INT ixmask;           /* XCoordCellOrderMask                  */
  INT iymask;           /* YCoordCellOrderMask                  */
  INT izmask;           /* ZCoordCellOrderMask                  */
  Y_NANOobj ****v0p;       /* OrderingMask                         */
  Y_NANOobj *v0tmp;        /* TMPobject                            */
  Y_NANOobj *v0con;        /* ContactorObjectP1Cell                */
  Y_NANOobj *v0conp;       /* ContactorPreviosuObject              */

  v0p=TalY_NANOobj4(3,3,3);
  v0con=v0next;
  v0conp=v0con;

  for(k=0;k<3;k++)      /* InitializeOrderingMask               */
  { for(j=0;j<3;j++)
    { for(i=0;i<3;i++)
      { v0p[k][j][i]=this;
  } } }

  while(v0con!=this)
  { ixold=v0con->iccx;
    iyold=v0con->iccy;
    izold=v0con->iccz;
    vtoolbox.Y_NANOtool(&integerize,v0con,this);
    if(ISALTB(v0con,v0conp)==1)
    { ixmask=v0con->iccx-ixold+1;
      iymask=v0con->iccy-iyold+1;
      izmask=v0con->iccz-izold+1;
      
      v0tmp=v0p[izmask][iymask][ixmask];/* AdvancePointer mask[ix][iy][iz]*/
      while(ISALTB(v0tmp,v0con)==1)
      { v0p[izmask][iymask][ixmask]=v0tmp;
        v0tmp=v0tmp->v0next;
      }
      v0con->v0next->v0prev=v0con->v0prev;
      v0con->v0prev->v0next=v0con->v0next;
      v0con->v0prev=v0p[izmask][iymask][ixmask];
      v0con->v0next=v0p[izmask][iymask][ixmask]->v0next;
      v0con->v0next->v0prev=v0con;
      v0p[izmask][iymask][ixmask]->v0next=v0con;
      v0con=v0conp->v0next;
    }
    else
    { v0conp=v0con;
      v0con=v0con->v0next;
    }
  }
  FREE(v0p);
}

void Y_NANOContaiBinList::TLSearch()
{ Y_NANOobj *v0con;        /* ContactorObjectP1Cell                */
  Y_NANOobj *v0tar;        /* TargetObject                         */
  Y_NANOobj *v0beg1;       /* ObjectImmediatelyBeforeP2Cell        */
  Y_NANOobj *v0beg2;       /* ObjectImmediatelyBeforeP2Cell        */
  Y_NANOobj *v0beg3;       /* ObjectImmediatelyBeforeP3Cell        */
  Y_NANOobj *v0beg4;       /* ObjectImmediatelyBeforeP4Cell        */
  Y_NANOobj *v0beg5;       /* ObjectImmediatelyBeforeP5Cell        */
  Y_NANOobj *v0end2;       /* ContactMaskEndObjectFor v0big2       */
  Y_NANOobj *v0end3;       /* ContactMaskEndObjectFor v0big3       */
  Y_NANOobj *v0end4;       /* ContactMaskEndObjectFor v0big4       */
  Y_NANOobj *v0end5;       /* ContactMaskEndObjectFor v0big5       */
  Y_NANOobj *v0next1;      /* ObjectInNextCellFor v0con            */
  Y_NANOobj *v0tmp;        /* TMPObject                            */

  v0con=v0next;
  v0next1=v0con;
  v0beg1=v0con;
  v0beg2=v0con;
  v0beg3=v0con;
  v0beg4=v0con;
  v0beg5=v0con;
  v0end2=v0con;
  v0end3=v0con;
  v0end4=v0con;
  v0end5=v0con;

  while(v0con!=this)
  { if(v0con==v0next1)
    { v0tmp=v0next1;    /* Advance v0next1                      */
      while((v0tmp!=this) && ISALEB(v0tmp,v0con))
      { v0tmp=v0tmp->v0next;
      }
      v0next1=v0tmp;
        
      v0tmp=v0beg1;     /* Advance v0beg1                       */
      while(ISALTCELL(v0tmp,v0con->iccx-1,v0con->iccy,v0con->iccz))
      { v0tmp=v0tmp->v0next;
      }
      v0beg1=v0tmp;
        
      v0tmp=v0beg2;     /* Advance v0beg2                       */
      while(ISALTCELL(v0tmp,v0con->iccx-1,v0con->iccy-1,v0con->iccz))
      { v0tmp=v0tmp->v0next;
      }
      v0beg2=v0tmp;

      v0tmp=v0beg3;     /* Advance v0beg3                       */
      while(ISALTCELL(v0tmp,v0con->iccx-1,v0con->iccy-1,v0con->iccz-1))
      { v0tmp=v0tmp->v0next;
      }
      v0beg3=v0tmp;
        
      v0tmp=v0beg4;     /* Advance v0beg4                       */
      while(ISALTCELL(v0tmp,v0con->iccx-1,v0con->iccy,v0con->iccz-1))
      { v0tmp=v0tmp->v0next;
      }
      v0beg4=v0tmp;
        
      v0tmp=v0beg5;     /* Advance v0beg5                       */
      while(ISALTCELL(v0tmp,v0con->iccx-1,v0con->iccy+1,v0con->iccz-1))
      { v0tmp=v0tmp->v0next;
      }
      v0beg5=v0tmp;
        
      v0tmp=v0end2;     /* Advance v0end2                       */
      while(ISALECELL(v0tmp,v0con->iccx+1,v0con->iccy-1,v0con->iccz))
      { v0tmp=v0tmp->v0next;
      }
      v0end2=v0tmp;
        
      v0tmp=v0end3;     /* Advance v0end3                       */
      while(ISALECELL(v0tmp,v0con->iccx+1,v0con->iccy-1,v0con->iccz-1))
      { v0tmp=v0tmp->v0next;
      }
      v0end3=v0tmp;

      v0tmp=v0end4;     /* Advance v0end4                       */
      while(ISALECELL(v0tmp,v0con->iccx+1,v0con->iccy,v0con->iccz-1))
      { v0tmp=v0tmp->v0next;
      }
      v0end4=v0tmp;

      v0tmp=v0end5;     /* Advance v0end5                       */
      while(ISALECELL(v0tmp,v0con->iccx+1,v0con->iccy+1,v0con->iccz-1))
      { v0tmp=v0tmp->v0next;
      }
      v0end5=v0tmp;
    }

    v0tar=v0beg1;       /* SearchForContacts                    */
    while(v0tar!=v0con)
    { vtoolbox.Y_NANOtool(&contact,v0con,v0tar);
      v0tar=v0tar->v0next;
    }

    v0tar=v0beg2;
    while(v0tar!=v0end2)
    { vtoolbox.Y_NANOtool(&contact,v0con,v0tar);
      v0tar=v0tar->v0next;
    }

    v0tar=v0beg3;
    while(v0tar!=v0end3)
    { vtoolbox.Y_NANOtool(&contact,v0con,v0tar);
      v0tar=v0tar->v0next;
    }

    v0tar=v0beg4;
    while(v0tar!=v0end4)
    { vtoolbox.Y_NANOtool(&contact,v0con,v0tar);
      v0tar=v0tar->v0next;
    }

    v0tar=v0beg5;
    while(v0tar!=v0end5)
    { vtoolbox.Y_NANOtool(&contact,v0con,v0tar);
      v0tar=v0tar->v0next;
    }
    v0con=v0con->v0next;
  }

  v0con=v0next;
  while(v0con!=this)
  { vtoolbox.Y_NANOtool(&move_,v0con);
    v0con=v0con->v0next;
  }
}

void Y_NANOContaiBinList::TLContact()
{ if(dcell<D0) TQSort();
  else         TLSort();
  TLSearch();
}


void Y_NANOContaiBinList::Quadratic()
{ Y_NANOobj *v0con;        /* ContactorObjectP1Cell                */
  Y_NANOobj *v0tar;        /* TargetObject                         */
  
  if(dcell<D0) TQSort();
  else         TLSort();

  v0con=v0next;
  while(v0con!=this)
  { v0tar=v0con->v0next;
    while(v0tar!=this)
    { vtoolbox.Y_NANOtool(&contact,v0con,v0tar);
      v0tar=v0tar->v0next;
    }
    vtoolbox.Y_NANOtool(&move_,v0con);
    v0con=v0con->v0next;
  }
}


